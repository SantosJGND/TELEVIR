import os
from datetime import date
from typing import List

from django.contrib.auth.models import User
from django.core.management.base import BaseCommand
from managing_files.models import ProcessControler
from pathogen_identification.constants_settings import ConstantsSettings
from pathogen_identification.deployment_main import Run_Main_from_Leaf
from pathogen_identification.models import (
    ParameterSet,
    PIProject_Sample,
    Projects,
    SoftwareTree,
    SoftwareTreeNode,
)
from pathogen_identification.utilities.utilities_pipeline import (
    Parameter_DB_Utility,
    PipelineTree,
    Utility_Pipeline_Manager,
    Utils_Manager,
)
from utils.process_SGE import ProcessSGE


class Sample_Staging:
    """
    Class to stage samples for a project
    """

    def __init__(self, sample: PIProject_Sample):

        self.sample = sample
        self.is_deleted = self.sample.is_deleted


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
        available_paths_explicit = {
            z: g for z, g in local_paths.items() if z in available_paths.keys()
        }

        reduced_dag, reduced_nodes = local_tree.reduced_tree(
            list(available_paths_explicit.keys())
        )
        reduced_tree = utils.pipe_tree_from_dag_dict(
            reduced_dag, reduced_nodes, technology, tree_makeup
        )

        pipeline_utils = Utility_Pipeline_Manager()
        module_tree = pipeline_utils.compress_software_tree(reduced_tree)

        print(module_tree.compress_dag_dict)
        # for project_sample in samples:
        #    if not project_sample.is_deleted:
        #        for leaf, path in available_paths_explicit.items():
