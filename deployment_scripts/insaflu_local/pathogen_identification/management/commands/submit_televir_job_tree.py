import os
from datetime import date
from typing import List

from django.contrib.auth.models import User
from django.core.management.base import BaseCommand
from managing_files.models import ProcessControler
from pathogen_identification.models import (
    PIProject_Sample,
    Projects,
    SoftwareTree,
    SoftwareTreeNode,
)
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
            "--user_id",
            type=int,
            help="user deploying the run (pk)",
        )

        parser.add_argument(
            "--project_id",
            type=int,
            help="project to be run (pk)",
        )

        parser.add_argument(
            "-o",
            "--outdir",
            type=str,
            help="output directory",
        )

    def handle(self, *args, **options):
        ###
        #### SETUP

        user = User.objects.get(pk=options["user_id"])
        project = Projects.objects.get(pk=options["project_id"])
        technology = project.technology

        ### PROCESS CONTROLER
        process_controler = ProcessControler()
        process_SGE = ProcessSGE()

        process_SGE.set_process_controler(
            user,
            process_controler.get_name_televir_project(project_pk=project.pk),
            ProcessControler.FLAG_RUNNING,
        )

        ### UTILITIES
        utils = Utils_Manager()
        samples = PIProject_Sample.objects.filter(project=project)
        local_tree = utils.generate_project_tree(technology, project, user)
        local_paths = local_tree.get_all_graph_paths_explicit()

        tree_makeup = local_tree.makeup

        pipeline_tree = utils.generate_software_tree(technology, tree_makeup)
        pipeline_tree_index = utils.get_software_tree_index(technology, tree_makeup)
        pipeline_tree_query = SoftwareTree.objects.get(pk=pipeline_tree_index)

        ### MANAGEMENT
        submission_dict = {sample: [] for sample in samples if not sample.is_deleted}
        matched_paths = {
            leaf: utils.utility_manager.match_path_to_tree_safe(path, pipeline_tree)
            for leaf, path in local_paths.items()
        }
        available_paths = {
            leaf: path for leaf, path in matched_paths.items() if path is not None
        }

        available_path_nodes = {
            leaf: SoftwareTreeNode.objects.get(
                software_tree__pk=pipeline_tree_index, index=path
            )
            for leaf, path in available_paths.items()
        }

        ### SUBMISSION

        pipeline_utils = Utility_Pipeline_Manager()

        reduced_tree = utils.tree_subset(pipeline_tree, list(matched_paths.values()))

        module_tree = pipeline_utils.compress_software_tree(reduced_tree)

        for project_sample in samples:
            if not project_sample.is_deleted:
                deployment_tree = Tree_Progress(module_tree, project_sample, project)

                print("FIRST")
                print(len(deployment_tree.current_nodes))
                print(deployment_tree.get_current_module())
                deployment_tree.run_current_nodes()
                deployment_tree.update_nodes()
                print("SECOND")
                print(len(deployment_tree.current_nodes))
                print(deployment_tree.get_current_module())
                deployment_tree.deploy_nodes()
                print("THIRD")
                print(len(deployment_tree.current_nodes))
                print(deployment_tree.get_current_module())
                deployment_tree.deploy_nodes()
                print("FOURTH")
                print(len(deployment_tree.current_nodes))
                print(deployment_tree.get_current_module())
                deployment_tree.deploy_nodes()
                print("FIFTH")
                print(len(deployment_tree.current_nodes))
                print(deployment_tree.get_current_module())
                deployment_tree.deploy_nodes()

                print("SIXTH")
                print(len(deployment_tree.current_nodes))
                print(deployment_tree.get_current_module())
                deployment_tree.deploy_nodes()

        print(len(samples))
