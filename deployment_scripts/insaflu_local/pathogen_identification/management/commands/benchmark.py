#####
### generate tree
#####
import logging
import os
import sys
from datetime import date
from tkinter.tix import Tree

import pandas as pd
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
from pathogen_identification.utilities.benchmark_graph_utils import pipeline_tree
from pathogen_identification.utilities.utilities_pipeline import (
    Parameter_DB_Utility,
    PipelineTree,
    Utility_Pipeline_Manager,
)
from settings.constants_settings import ConstantsSettings
from settings.models import Sample
from utils.process_SGE import ProcessSGE
from utils.utils import Utils


class Tree_Node:
    module: str
    name: str
    node_index: int
    children: list
    parameters: pd.DataFrame
    software_tree_pk: int

    def __init__(self, pipe_tree, node_index, software_tree_pk: int):
        node_metadata = pipe_tree.node_index.loc[node_index].node

        self.module = node_metadata[0]
        self.node_index = node_index
        self.branch = pipe_tree.nodes_df.loc[node_index].branch
        self.children = pipe_tree.edge_df[
            pipe_tree.edge_df.parent == node_index
        ].child.tolist()
        self.parameters = self.determine_params(pipe_tree)
        self.software_tree_pk = software_tree_pk

    def _is_node_leaf(self):
        return len(self.children) == 0

    def generate_software_tree_node_entry(self, pipe_tree):

        if not self._is_node_leaf():
            return
        node_metadata = pipe_tree.node_index.loc[self.node_index].node
        software_tree = SoftwareTree.objects.get(name=self.software_tree_pk).pk

        try:
            node = SoftwareTreeNode.objects.get(
                software_tree=software_tree,
                index=self.node_index,
                name=node_metadata[0],
                value=node_metadata[1],
                node_type=node_metadata[2],
                node_place=1,
            )
        except SoftwareTreeNode.DoesNotExist:
            node = None

        return node

    def setup_parameterset(self, project, sample):

        try:
            parameter_set = ParameterSet.objects.get(
                project=project, sample=sample, status=ParameterSet.STATUS_FINISHED
            )

        except ParameterSet.DoesNotExist:
            parameter_set = ParameterSet()
            parameter_set.project = project
            parameter_set.sample = sample
            parameter_set.status = ParameterSet.STATUS_FINISHED
            parameter_set.save()

            parameter_set.save()

        return parameter_set

    def determine_params(self, pipe_tree):
        arguments_list = []
        for node in self.branch:
            node_metadata = pipe_tree.node_index.loc[node].node
            arguments_list.append(node_metadata)

        arguments_df = pd.DataFrame(
            arguments_list, columns=["module", "software", "parameter"]
        )
        return arguments_df


class Tree_Progress:
    tree: PipelineTree
    current_nodes: list
    current_module: str

    sample: PIProject_Sample
    project: Projects

    def __init__(
        self, pipe_tree: PipelineTree, sample: PIProject_Sample, project: Projects
    ):
        pipe_tree.nodes_df = pd.DataFrame(
            pipe_tree.nodes_compress, columns=["node", "branch"]
        ).set_index("node")
        pipe_tree.edge_df = pd.DataFrame(
            pipe_tree.edge_compress, columns=["parent", "child"]
        )

        self.tree = pipe_tree
        self.current_nodes = [
            Tree_Node(pipe_tree, 0, software_tree_pk=pipe_tree.software_tree_pk)
        ]
        self.current_module = ""
        self.determine_current_module()
        self.sample = sample
        self.project = project

        self.logger = logging.getLogger(__name__)

    def get_current_module(self):
        return self.current_module

    def determine_current_module(self):
        self.current_module = self.current_nodes[0].module

    def update_nodes(self):
        new_nodes = []
        for node in self.current_nodes:
            children = node.children
            for child in children:
                new_nodes.append(Tree_Node(self.tree, child, node.software_tree_pk))

        if len(new_nodes) == 0:
            self.current_module = "end"
        else:
            self.current_nodes = new_nodes
        self.determine_current_module()


class Insaflu_Cli:
    def __init__(self):
        self.logger = logging.getLogger(__name__)

    @staticmethod
    def create_televir_project_if_not_exists(project_name, user, technology):
        try:
            project = Projects.objects.get(name=project_name)
        except Projects.DoesNotExist:
            project = Projects()
            project.name = project_name
            project.owner = user
            project.technology = technology
            project.save()
        return project

    @staticmethod
    def read_fofn(fofn):

        with open(fofn) as f:
            fastq_paths = f.read().splitlines()

        r1 = fastq_paths[0]
        if len(fastq_paths) > 1:
            r2 = fastq_paths[1]
        else:
            r2 = None

        return r1, r2

    def create_sample_from_fofn(self, fofn, user: User, technology: str):

        r1, r2 = self.read_fofn(fofn)
        name = os.path.basename(fofn).split(".")[0]

        try:
            sample = Sample.objects.get(name=name, owner=user)

        except Sample.DoesNotExist:

            sample = self.sample_save(name, user, r1, r2, technology)
            self.move_sample(sample, user)
            self.sample_preprocess(sample, user)

        return sample

    def sample_save(self, name, user, r1, r2, technology):
        utils = Utils()
        print(technology)

        if not os.path.exists(r1):
            raise FileNotFoundError(f"File {r1} does not exist")

        sample = Sample.objects.create(
            name=name,
            owner=user,
            path_name_1=r1,
            path_name_2=r2,
            type_of_fastq=int(technology == ConstantsSettings.TECHNOLOGY_minion),
            date_of_onset=date.today(),
            date_of_collection=date.today(),
        )

        sample.is_deleted = False
        sample.is_obsolete = False
        sample.file_name_1 = utils.clean_name(os.path.basename(sample.path_name_1.name))
        sample.is_valid_1 = True
        if sample.exist_file_2():
            sample.file_name_2 = utils.clean_name(
                os.path.basename(sample.path_name_2.name)
            )
            sample.is_valid_2 = True
        else:
            sample.is_valid_2 = False
        sample.has_files = True

        sample.manual_upload()

        sample.save()

        return sample

    def move_sample(self, sample: Sample, user: User):
        utils = Utils()

        sz_file_to = os.path.join(
            getattr(settings, "MEDIA_ROOT", None),
            utils.get_path_to_fastq_file(user.id, sample.id),
            sample.file_name_1,
        )
        utils.move_file(
            os.path.join(
                getattr(settings, "MEDIA_ROOT", None), sample.path_name_1.name
            ),
            sz_file_to,
        )
        sample.path_name_1.name = os.path.join(
            utils.get_path_to_fastq_file(user.id, sample.id),
            sample.file_name_1,
        )

        if sample.exist_file_2():
            sz_file_to = os.path.join(
                getattr(settings, "MEDIA_ROOT", None),
                utils.get_path_to_fastq_file(user.id, sample.id),
                sample.file_name_2,
            )
            utils.move_file(
                os.path.join(
                    getattr(settings, "MEDIA_ROOT", None), sample.path_name_2.name
                ),
                sz_file_to,
            )
            sample.path_name_2.name = os.path.join(
                utils.get_path_to_fastq_file(user.id, sample.id),
                sample.file_name_2,
            )
        sample.save()

    def sample_preprocess(self, sample: Sample, user: User):
        from managing_files.manage_database import ManageDatabase

        process_SGE = ProcessSGE()
        job_name = None

        try:
            (job_name_wait, job_name) = user.profile.get_name_sge_seq(
                Profile.SGE_PROCESS_clean_sample, Profile.SGE_SAMPLE
            )
            if sample.is_type_fastq_gz_sequencing():  ### default is Illumina
                taskID = process_SGE.set_run_trimmomatic_species(sample, user, job_name)
            else:  ### Minion, codify with other
                taskID = process_SGE.set_run_clean_minion(sample, user, job_name)

            sample.is_ready_for_projects = True
            sample.save()

        except Exception as e:
            self.logger.error("Fail to run: ProcessSGE - " + str(e))
            return

        ## refresh sample list for this user
        if not job_name is None:
            process_SGE.set_create_sample_list_by_user(user, [job_name])
        ###

        manageDatabase = ManageDatabase()

        manageDatabase.set_sample_metakey(
            sample,
            user,
            MetaKeyAndValue.META_KEY_Queue_TaskID,
            MetaKeyAndValue.META_VALUE_Queue,
            taskID,
        )

    def piproject_sample_from_sample(
        self, sample: Sample, project: Projects, user: User
    ):

        try:
            project_sample = PIProject_Sample.objects.get(
                project=project, sample=sample
            )
        except PIProject_Sample.DoesNotExist:
            project_sample_input = sample.file_name_1
            if sample.is_valid_2:
                project_sample_input += ";" + sample.file_name_2

            project_sample = PIProject_Sample.objects.create(
                project=project,
                sample=sample,
                name=sample.name,
                input=project_sample_input,
                technology=sample.type_of_fastq,
                report="report",
            )
            project_sample.save()

        return project_sample

    def create_benchmark_softwaretree_if_not_exists(self):

        try:
            benchmark_softwaretree = SoftwareTree.objects.filter(name="default")
        except SoftwareTree.DoesNotExist:
            benchmark_softwaretree = SoftwareTree.objects.create(
                name="default",
                owner=User.objects.get(username="admin"),
            )
            benchmark_softwaretree.save()

        return benchmark_softwaretree


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

            if len(tree_differences) > 0:
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
        run_manager = Tree_Progress(software_tree, project_sample, project)

        print(len(run_manager.current_nodes))
        print(run_manager.get_current_module())
        run_manager.update_nodes()
        print(len(run_manager.current_nodes))
        print(run_manager.get_current_module())
        run_manager.update_nodes()
        print(len(run_manager.current_nodes))
        print(run_manager.get_current_module())
        run_manager.update_nodes()
        print(len(run_manager.current_nodes))
        print(run_manager.get_current_module())
        run_manager.update_nodes()
        print(len(run_manager.current_nodes))
        print(run_manager.get_current_module())
        run_manager.update_nodes()
        print(len(run_manager.current_nodes))
        print(run_manager.get_current_module())
        run_manager.update_nodes()
        print(len(run_manager.current_nodes))
        print(run_manager.get_current_module())
        run_manager.update_nodes()
        print(len(run_manager.current_nodes))
        print(run_manager.get_current_module())
