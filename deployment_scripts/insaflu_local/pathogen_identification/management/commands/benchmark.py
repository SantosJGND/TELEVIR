#####
### generate tree
#####
import copy
import logging
import os
import shutil
import sys
from datetime import date
from threading import Thread
from tkinter.tix import Tree
from typing import List

import pandas as pd
from constants.constants import Televir_Metadata_Constants as Televir_Metadata
from constants.constants import TypePath
from constants.meta_key_and_values import MetaKeyAndValue
from django.conf import settings
from django.contrib.auth.models import User
from django.core.management.base import BaseCommand
from extend_user.models import Profile
from pathogen_identification.constants_settings import ConstantsSettings as PIConstants
from pathogen_identification.models import (
    ParameterSet,
    PIProject_Sample,
    Projects,
    SoftwareTree,
    SoftwareTreeNode,
)
from pathogen_identification.modules.remap_class import Mapping_Instance
from pathogen_identification.modules.run_main import RunMain_class
from pathogen_identification.utilities.benchmark_graph_utils import pipeline_tree
from pathogen_identification.utilities.insaflu_cli import Insaflu_Cli
from pathogen_identification.utilities.tree_deployment import Tree_Progress
from pathogen_identification.utilities.utilities_pipeline import (
    Parameter_DB_Utility,
    PipelineTree,
    Utility_Pipeline_Manager,
)
from settings.constants_settings import ConstantsSettings
from settings.models import Sample
from utils.process_SGE import ProcessSGE
from utils.utils import Utils


class Command(BaseCommand):
    help = "Populates the DBS"

    def add_arguments(self, parser):
        parser.add_argument(
            "--fofn",
            "-f",
            type=str,
            help="file of fastq files. one file path per line.",
        )
        parser.add_argument(
            "--fdir",
            "-d",
            type=str,
            help="dictory containing fofn files.",
        )

        parser.add_argument(
            "-t",
            "--tech",
            type=str,
            default="nanopore",
            help="technology. detemines parameters used. [default=nanopore]",
        )

        parser.add_argument(
            "-p", "--project", type=str, default="", help="Output directory."
        )

        parser.add_argument(
            "-u", "--user", type=str, default="admin", help="User to assign to project."
        )

        parser.add_argument(
            "--clean",
            action="store_true",
            default=False,
            help="move output reports to final output directory, intermediate files and config files to run output directories",
        )

        parser.add_argument(
            "--fdel",
            action="store_true",
            default=False,
            help="clean output repositories, keep only report files and assembly file. Recommend for large benchmarking runs.",
        )

        parser.add_argument(
            "--ref",
            type=str,
            required=False,
            default=None,
            help="reference genome for host depletion",
        )

        parser.add_argument(
            "--force",
            action="store_true",
            required=False,
            default=False,
            help="reference genome for host depletion",
        )

    def sanitize_input(self, options):
        if options["fofn"] is None and options["fdir"] is None:
            print("Error: must provide either fofn or fdir")
            sys.exit(1)

        if options["fofn"] is not None and options["fdir"] is not None:
            print("Error: must provide either fofn or fdir")
            sys.exit(1)

        if options["fofn"] is not None:
            if not os.path.exists(options["fofn"]):
                print("Error: fofn file does not exist")
                sys.exit(1)

        if options["fdir"] is not None:
            if not os.path.exists(options["fdir"]):
                print("Error: fofn file does not exist")
                sys.exit(1)

        if options["ref"] is not None:
            if not os.path.exists(options["ref"]):
                print("Error: reference file does not exist")
                sys.exit(1)

        technology = options["tech"]
        sample_fofn = options["fofn"]
        sample_fdir = options["fdir"]
        project_name = options["project"]
        user_name = options["user"]

        try:  # get user
            user = User.objects.get(username=user_name)
        except User.DoesNotExist:
            print("User does not exist")
            sys.exit(1)

        if technology.lower() in ["nanopore", "ont"]:
            technology = ConstantsSettings.TECHNOLOGY_minion
        elif technology.lower() in ["illumina", "pacbio", "Illumina/IonTorrent"]:
            technology = ConstantsSettings.TECHNOLOGY_illumina
        else:
            print("Technology not supported")
            sys.exit(1)

        if project_name == "":  # get project name
            project_name = "benchmark_" + date.today().strftime("%Y%m%d")

        return technology, sample_fofn, sample_fdir, project_name, user

    def generate_hardcoded_pipeline_tree(self, technology):

        pipe_tree = pipeline_tree()
        pipe_tree.param_input(technology)
        pipe_tree.create_pipe_tree()

        software_tree = PipelineTree(
            technology=technology,
            node_index=pipe_tree.node_index.reset_index().to_numpy().tolist(),
            edges=[x for x in pipe_tree.edge_dict if x[0] != "root"],
            leaves=pipe_tree.leaves,
            makeup=-1,
        )

        return software_tree

    def update_software_tree(self, software_tree: PipelineTree):
        parameter_util = Parameter_DB_Utility()
        utility_manager = Utility_Pipeline_Manager()

        if parameter_util.check_default_software_tree_exists(
            software_tree.technology, global_index=software_tree.makeup
        ):
            existing_pipeline_tree = parameter_util.query_software_default_tree(
                software_tree.technology, global_index=software_tree.makeup
            )

            tree_differences = utility_manager.compare_software_trees_given(
                existing_pipeline_tree, software_tree
            )

            if not tree_differences:
                parameter_util.update_software_tree(software_tree)
        else:

            parameter_util.update_software_tree(software_tree)

        software_tree_pk = parameter_util.query_software_tree_pk(software_tree)

        if not software_tree_pk is None:
            software_tree.software_tree_pk = software_tree_pk

        return software_tree

    def compress_software_tree(self, software_tree: PipelineTree):
        utility_manager = Utility_Pipeline_Manager()

        software_tree = utility_manager.compress_software_tree(software_tree)
        return software_tree

    def generate_modular_software_tree(self, technology):

        software_tree = self.generate_hardcoded_pipeline_tree(technology)

        self.update_software_tree(software_tree)

        software_tree = self.compress_software_tree(software_tree)

        return software_tree

    def handle(self, *args, **options):
        ###
        #

        insaflu_cli = Insaflu_Cli()

        technology, sample_fofn, sample_fdir, project_name, user = self.sanitize_input(
            options
        )

        project = insaflu_cli.create_televir_project_if_not_exists(
            project_name, user, technology
        )

        sample = insaflu_cli.create_sample_from_fofn(sample_fofn, user, technology)
        project_sample = insaflu_cli.piproject_sample_from_sample(sample, project, user)

        software_tree = self.generate_modular_software_tree(technology)
        print(software_tree.compress_dag_dict)
        print(software_tree.nodes_compress)
        print(software_tree.node_index)
        software_tree.node_index.to_csv("node_index.csv")

        deployment_tree = Tree_Progress(software_tree, project_sample, project)

        print("leaves: ", software_tree.compress_dag_dict)

        current_module = deployment_tree.get_current_module()
        while current_module != "end":
            # for x in range(0):

            print("NEXT")
            print(len(deployment_tree.current_nodes))
            print(deployment_tree.get_current_module())
            print([x.node_index for x in deployment_tree.current_nodes])
            print([x.children for x in deployment_tree.current_nodes])
            deployment_tree.deploy_nodes()

            current_module = deployment_tree.get_current_module()
