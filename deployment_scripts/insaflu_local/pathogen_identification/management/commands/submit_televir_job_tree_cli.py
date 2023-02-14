import os
import sys
from datetime import date
from typing import List

from django.contrib.auth.models import User
from django.core.management.base import BaseCommand
from managing_files.models import ProcessControler
from pathogen_identification.models import PIProject_Sample, Projects
from pathogen_identification.utilities.insaflu_cli import Insaflu_Cli
from pathogen_identification.utilities.tree_deployment import Tree_Progress
from pathogen_identification.utilities.utilities_pipeline import (
    Utility_Pipeline_Manager,
    Utils_Manager,
)
from utils.process_SGE import ProcessSGE


class Command(BaseCommand):
    help = "deploy run"

    def add_arguments(self, parser):
        parser.add_argument(
            "--user_name",
            type=str,
            help="user deploying the run (string)",
        )

        parser.add_argument(
            "--project_name",
            type=str,
            help="project to be run (string)",
        )

        parser.add_argument(
            "--technology",
            type=str,
            default=None,
            help="Project technology. used to create project if it does not exist",
        )

        parser.add_argument(
            "--fofn",
            "-f",
            type=str,
            default=None,
            help="file of fastq files. one file path per line.",
        )

        parser.add_argument(
            "--metadata",
            "-m",
            type=str,
            default=None,
            help="file of metadata. one file path per line.",
        )

        parser.add_argument(
            "-o",
            "--outdir",
            type=str,
            help="output directory",
        )

    @staticmethod
    def add_samples(user: User, project: Projects, **options):
        insaflu_cli = Insaflu_Cli()
        sample_list = []

        technology = project.technology

        if options["fofn"] is not None:
            sample_fofn = options["fofn"]
            sample = insaflu_cli.create_sample_from_fofn(sample_fofn, user, technology)
            sample_list.append(sample)

        if options["metadata"] is not None:
            metadata_fofn = options["metadata"]
            samples = insaflu_cli.create_sample_from_metadata(metadata_fofn, user)
            sample_list.extend(samples)

        for sample in sample_list:
            project_sample = insaflu_cli.piproject_sample_from_sample(
                sample, project, user
            )

    def process_arguments(self, *args, **options):

        try:
            user = User.objects.get(username=options["user_name"])
        except User.DoesNotExist:
            print("User does not exist")
            sys.exit(1)

        try:
            project = Projects.objects.get(name=options["project_name"], owner=user)
        except Projects.DoesNotExist:
            print(f"Project {options['project_name']} does not exist.")

            if options["technology"] is None:
                print("Please provide a technology for the project")
                sys.exit(1)
            else:
                technology = options["technology"]
                print(
                    f"Creating project {options['project_name']} with technology {technology}"
                )

            insaflu_cli = Insaflu_Cli()

            project = insaflu_cli.create_televir_project_if_not_exists(
                options["project_name"], user, technology
            )

        return user, project

    def handle(self, *args, **options):
        ###
        #### SETUP

        user, project = self.process_arguments(*args, **options)
        self.add_samples(user, project, **options)

        technology = project.technology

        ### PROCESS CONTROLER
        process_controler = ProcessControler()
        process_SGE = ProcessSGE()

        process_SGE.set_process_controler(
            user,
            process_controler.get_name_televir_project(project_pk=project.pk),
            ProcessControler.FLAG_RUNNING,
        )

        process = ProcessControler.objects.filter(
            owner__id=user.pk,
            name=process_controler.get_name_televir_project(project_pk=project.pk),
        )

        if process.exists():
            process = process.first()

            print(f"Process {process.pk} has been submitted")

        else:
            print("Process does not exist")

        ### UTILITIES
        utils = Utils_Manager()
        samples = PIProject_Sample.objects.filter(project=project)
        local_tree = utils.generate_project_tree(technology, project, user)
        local_paths = local_tree.get_all_graph_paths_explicit()

        tree_makeup = local_tree.makeup

        pipeline_tree = utils.generate_software_tree(technology, tree_makeup)

        ### MANAGEMENT
        matched_paths = {
            leaf: utils.utility_manager.match_path_to_tree_safe(path, pipeline_tree)
            for leaf, path in local_paths.items()
        }
        ### SUBMISSION

        pipeline_utils = Utility_Pipeline_Manager()

        reduced_tree = utils.tree_subset(pipeline_tree, list(matched_paths.values()))

        module_tree = pipeline_utils.compress_software_tree(reduced_tree)

        try:

            for project_sample in samples:
                if not project_sample.is_deleted:
                    deployment_tree = Tree_Progress(
                        module_tree, project_sample, project
                    )

                    print("leaves: ", module_tree.compress_dag_dict)

                    current_module = deployment_tree.get_current_module()
                    while current_module != "end":
                        # for x in range(0):

                        deployment_tree.deploy_nodes()
                        # deployment_tree.update_nodes()

                        current_module = deployment_tree.get_current_module()

            process_SGE.set_process_controler(
                user,
                process_controler.get_name_televir_project(project_pk=project.pk),
                ProcessControler.FLAG_FINISHED,
            )

        except Exception as e:
            print(e)
            process_SGE.set_process_controler(
                user,
                process_controler.get_name_televir_project(project_pk=project.pk),
                ProcessControler.FLAG_ERROR,
            )
            raise e
