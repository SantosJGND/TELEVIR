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
from pathogen_identification.utilities.update_DBs import (
    Update_Assembly,
    Update_Classification,
    Update_Remap,
    Update_RunMain_Initial,
    Update_RunMain_Secondary,
)
from pathogen_identification.utilities.utilities_pipeline import (
    Parameter_DB_Utility,
    PipelineTree,
    Utility_Pipeline_Manager,
)
from settings.constants_settings import ConstantsSettings
from settings.models import Sample
from utils.process_SGE import ProcessSGE
from utils.utils import Utils


class PathogenIdentification_Deployment_Manager:

    project: str
    prefix: str
    rdir: str
    threads: int
    run_engine: RunMain_class
    params = dict
    run_params_db = pd.DataFrame()
    pk: int = 0
    username: str
    prepped: bool = False
    sent: bool = False

    STATUS_ZERO = 0
    STATUS_PREPPED = 1
    STATUS_RUNNING = 2
    STATUS_SENT = 3

    def __init__(
        self,
        sample: PIProject_Sample,  # sample name
        project: Projects,
        username: str = "admin",
        technology: str = "ONT",
        deployment_root_dir: str = "/tmp/insaflu/insaflu_something",
        dir_branch: str = "deployment",
        threads: int = 3,
    ) -> None:

        self.username = username
        self.project = project.name
        self.sample = sample.name

        self.deployment_root_dir = deployment_root_dir
        self.dir_branch = dir_branch
        self.dir = os.path.join(self.deployment_root_dir, dir_branch)

        self.technology = technology
        self.install_registry = Televir_Metadata

        self.threads = threads
        self.file_r1 = sample.sample.get_fastq_available(TypePath.MEDIA_ROOT, True)
        if sample.sample.exist_file_2():
            self.file_r2 = sample.sample.get_fastq_available(TypePath.MEDIA_ROOT, False)
        else:
            self.file_r2 = ""

    def input_read_project_path(self, filepath) -> str:
        """copy input reads to project directory and return new path"""

        if not os.path.isfile(filepath):
            return ""

        rname = os.path.basename(filepath)
        new_rpath = os.path.join(self.dir, "reads") + "/" + rname
        shutil.copy(filepath, new_rpath)
        return new_rpath

    def configure(self) -> bool:
        """generate config dictionary for run_main, and copy input reads to project directory."""
        self.get_constants()

        self.generate_config_file()
        self.prep_test_env()

        new_r1_path = self.input_read_project_path(self.file_r1)
        new_r2_path = self.input_read_project_path(self.file_r2)

        self.config["sample_name"] = self.sample
        self.config["r1"] = new_r1_path
        self.config["r2"] = new_r2_path
        self.config["type"] = ["SE", "PE"][int(os.path.isfile(self.config["r2"]))]

        return True

    def get_constants(self):
        """set constants for technology"""
        if self.technology == ConstantsSettings.TECHNOLOGY_illumina:
            self.constants = PIConstants.CONSTANTS_ILLUMINA
        if self.technology == ConstantsSettings.TECHNOLOGY_minion:
            self.constants = PIConstants.CONSTANTS_ONT

    def generate_config_file(self):

        self.config = {
            "project": self.project,
            "source": self.install_registry.SOURCE,
            "deployment_root_dir": self.deployment_root_dir,
            "sub_directory": self.dir_branch,
            "directories": {},
            "threads": self.threads,
            "prefix": self.prefix,
            "project_name": self.project,
            "metadata": {
                x: os.path.join(self.install_registry.METADATA["ROOT"], g)
                for x, g in self.install_registry.METADATA.items()
            },
            "technology": self.technology,
            "bin": self.install_registry.BINARIES,
            "actions": {},
        }

        for dr, g in PIConstants.DIRS.items():
            self.config["directories"][dr] = os.path.join(self.dir, g)

        for dr, g in PIConstants.ACTIONS.items():
            self.config["actions"][dr] = g

        self.config.update(self.constants)

    def update_config(self):
        self.config["prefix"] = self.prefix

    def prep_test_env(self):
        """
        from main directory bearing scripts, params.py and main.sh, create metagenome run directory

        :return:
        """
        os.makedirs(self.dir, exist_ok=True)
        os.makedirs(
            os.path.join(PIConstants.media_directory, self.dir_branch),
            exist_ok=True,
        )
        os.makedirs(
            os.path.join(PIConstants.static_directory, self.dir_branch),
            exist_ok=True,
        )

        for directory in self.config["directories"].values():
            os.makedirs(directory, exist_ok=True)

    def close(self):
        if os.path.exists(self.dir):
            shutil.rmtree(self.dir)

    def import_params(self, run_params_db: pd.DataFrame):
        self.run_params_db = run_params_db

    def run_main_prep(self):

        if self.prepped:
            return

        self.run_engine = RunMain_class(self.config, self.run_params_db, self.username)
        self.run_engine.Prep_deploy()
        self.prepped = True

    def run_main(self):
        self.run_engine.Run_QC()
        self.run_engine.Run_PreProcess()
        self.run_engine.Sanitize_reads()
        self.run_engine.Run_Assembly()
        self.run_engine.Run_Classification()
        self.run_engine.Run_Remapping()

    def update_engine(self):
        self.update_config()
        self.run_engine.Update(self.config, self.run_params_db)

    def update_merged_targets(self, merged_targets: pd.DataFrame):
        self.run_engine.update_merged_targets(merged_targets)


class Tree_Node:
    module: str
    name: str
    node_index: int
    children: list
    parameters: pd.DataFrame
    software_tree_pk: int
    run_manager: PathogenIdentification_Deployment_Manager
    tree_node: SoftwareTreeNode
    parameter_set: ParameterSet

    def __init__(self, pipe_tree: PipelineTree, node_index: int, software_tree_pk: int):

        node_metadata = pipe_tree.node_index.loc[node_index].node

        self.module = node_metadata[0]
        self.node_index = node_index

        self.branch = pipe_tree.nodes_df.loc[node_index].branch
        self.children = pipe_tree.edge_df[
            pipe_tree.edge_df.parent == node_index
        ].child.tolist()

        self.parameters = self.determine_params(pipe_tree)
        self.software_tree_pk = software_tree_pk
        self.leaves = pipe_tree.leaves_from_node(node_index)

    def receive_run_manager(
        self, run_manager: PathogenIdentification_Deployment_Manager
    ):
        run_manager.prefix = f"run_leaf_{self.node_index}"
        run_manager.configure()
        run_manager.import_params(self.parameters)

        if self.node_index == 0:
            run_manager.run_main_prep()

        self.run_manager = run_manager

    def _is_node_leaf(self):
        return len(self.children) == 0

    def generate_software_tree_node_entry(self, pipe_tree: PipelineTree):

        if not self._is_node_leaf():
            return

        node_metadata = pipe_tree.node_index.loc[self.node_index].node
        software_tree = SoftwareTree.objects.get(pk=self.software_tree_pk)

        tree_node = SoftwareTreeNode.objects.filter(
            software_tree=software_tree,
            index=self.node_index,
            name=node_metadata[0],
            value=node_metadata[1],
            node_type=node_metadata[2],
        )

        try:
            tree_node = SoftwareTreeNode.objects.get(
                software_tree=software_tree,
                index=self.node_index,
                name=node_metadata[0],
                value=node_metadata[1],
                node_type=node_metadata[2],
            )
        except SoftwareTreeNode.DoesNotExist:
            tree_node = None

        return tree_node

    def setup_parameterset(
        self, project: Projects, sample: PIProject_Sample, node: SoftwareTreeNode
    ):

        try:
            parameter_set = ParameterSet.objects.get(
                project=project, sample=sample, status=ParameterSet.STATUS_FINISHED
            )

        except ParameterSet.DoesNotExist:
            parameter_set = ParameterSet()
            parameter_set.project = project
            parameter_set.sample = sample
            parameter_set.status = ParameterSet.STATUS_FINISHED
            parameter_set.leaf = node
            parameter_set.save()

        return parameter_set

    def register(self, project: Projects, sample: PIProject_Sample, tree: PipelineTree):

        tree_node = self.generate_software_tree_node_entry(tree)
        if tree_node is None:
            return False

        parameter_set = self.setup_parameterset(project, sample, tree_node)

        if parameter_set is None:
            return False

        self.tree_node = tree_node
        self.parameter_set = parameter_set

        return True

    def determine_params(self, pipe_tree):

        arguments_list = []
        for node in self.branch:

            node_metadata = pipe_tree.node_index.loc[node].node
            arguments_list.append(node_metadata)

        arguments_df = pd.DataFrame(
            arguments_list, columns=["parameter", "value", "flag"]
        )

        module_df = arguments_df[arguments_df.flag == "module"]
        module = module_df.parameter.values[0]
        software = module_df.value.values[0]
        parameters_df = arguments_df[arguments_df.flag == "param"]
        parameters_df["software"] = software
        parameters_df["module"] = module

        return parameters_df


class PathogenDeployment_Iterator:
    def __init__(self):
        self.deplyment_list = []

    def add_deployment(self, deployment: PathogenIdentification_Deployment_Manager):
        self.deplyment_list.append(deployment)

    def __getitem__(self, index):
        return self.deplyment_list[index]

    def __len__(self):
        return len(self.deplyment_list)


class TreeNode_Iterator:
    def __init__(self):
        self.node_list = []

    def add_node(self, node: Tree_Node):
        self.node_list.append(node)

    def __getitem__(self, index) -> Tree_Node:
        return self.node_list[index]

    def __len__(self):
        return len(self.node_list)


class Tree_Progress:
    tree: PipelineTree
    current_nodes: List[Tree_Node]
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

        self.logger = logging.getLogger(__name__)

        self.tree = pipe_tree
        self.sample = sample
        self.project = project

        self.initialize_nodes()
        self.determine_current_module()

    def setup_deployment_manager(self):
        utils = Utils()
        temp_dir = utils.get_temp_dir()

        prefix = f"{self.sample.sample.pk}_{self.sample.sample.name}"

        deployment_directory_structure = os.path.join(
            PIConstants.televir_subdirectory,
            f"{self.project.owner.pk}",
            f"{self.project.pk}",
            f"{self.sample.sample.pk}",
            prefix,
        )

        deployment_manager = PathogenIdentification_Deployment_Manager(
            self.sample,
            self.project,
            self.project.owner.username,
            self.project.technology,
            temp_dir,
            deployment_directory_structure,
            PIConstants.DEPLOYMENT_THREADS,
        )

        return deployment_manager

    def register_node_leaves(self, node: Tree_Node):

        if len(node.leaves) == 0:
            self.submit_node_run(node)

        for leaf in node.leaves:
            leaf_node = self.spawn_node_child(node, leaf)
            self.submit_node_run(leaf_node)

    def initialize_nodes(self):
        origin_node = Tree_Node(
            self.tree, 0, software_tree_pk=self.tree.software_tree_pk
        )

        run_manager = self.setup_deployment_manager()

        origin_node.receive_run_manager(run_manager)

        self.register_node_leaves(origin_node)

        self.current_nodes = [origin_node]
        # self.current_nodes.add_node(origin_node)

    def get_current_module(self):
        return self.current_module

    def determine_current_module(self):
        self.current_module = self.current_nodes[0].module

    def spawn_node_child(self, node: Tree_Node, child: int):

        new_node = Tree_Node(self.tree, child, node.software_tree_pk)
        run_manager_copy = copy.deepcopy(node.run_manager)
        new_node.receive_run_manager(run_manager_copy)
        new_node.run_manager.update_engine()

        return new_node

    def register_node(self, node: Tree_Node):
        print("registering node")

        if node.run_manager.sent:
            return False

        registraction_success = node.register(self.project, self.sample, self.tree)

        return registraction_success

    def update_node_dbs(self, node: Tree_Node, step="initial"):

        try:
            db_updated = Update_RunMain_Initial(
                node.run_manager.run_engine, node.parameter_set
            )
            if not db_updated:
                return False

            if (
                node.run_manager.run_engine.enrichment_performed
                or node.run_manager.run_engine.depletion_performed
            ):
                db_updated = Update_RunMain_Secondary(
                    node.run_manager.run_engine, node.parameter_set
                )
                if not db_updated:
                    return False

            if node.run_manager.run_engine.assembly_performed:
                db_updated = Update_Assembly(
                    node.run_manager.run_engine, node.parameter_set
                )
                if not db_updated:
                    return False

            if (
                node.run_manager.run_engine.read_classification_performed
                or node.run_manager.run_engine.contig_classification_performed
            ):
                db_updated = Update_Classification(
                    node.run_manager.run_engine, node.parameter_set, tag=step
                )
                if not db_updated:
                    return False

            if node.run_manager.run_engine.remapping_performed:
                node.run_manager.run_engine.export_sequences()
                node.run_manager.run_engine.export_final_reports()
                node.run_manager.run_engine.Summarize()
                node.run_manager.run_engine.generate_output_data_classes()
                node.run_manager.run_engine.export_logdir()
                db_updated = Update_Remap(
                    node.run_manager.run_engine, node.parameter_set
                )
                if not db_updated:
                    return False

            return True
        except Exception as e:
            self.logger.error(e)
            return False

    def submit_node_run(self, node: Tree_Node):

        print("Submitting node run: " + str(node.node_index))

        registration_success = self.register_node(node)
        print("Registration success: " + str(registration_success))
        if not registration_success:
            return

        self.update_node_dbs(node, step=node.module)

    def merge_node_targets(self):
        node_merged_targets = [
            n.run_manager.run_engine.merged_targets for n in self.current_nodes
        ]
        node_merged_targets = pd.concat(node_merged_targets, axis=0)
        node_merged_targets = node_merged_targets.drop_duplicates(subset=["taxid"])

        return node_merged_targets

    @staticmethod
    def get_remap_plans(nodes: List[Tree_Node]):

        for n in nodes:

            n.run_manager.run_engine.plan_remap_prep()

        return nodes

    def merge_node_targets_list(self, targetdf_list: List[Tree_Node]):

        node_merged_targets = [
            n.run_manager.run_engine.merged_targets for n in targetdf_list
        ]

        node_merged_targets = pd.concat(node_merged_targets, axis=0).reset_index()

        return node_merged_targets

    def get_node_node_targets(self, nodes_list: List[Tree_Node]):

        node_merged_targets = self.merge_node_targets_list(nodes_list)

        node_merged_targets = self.process_mapping_managerdf(node_merged_targets)

        return node_merged_targets

    @staticmethod
    def process_mapping_managerdf(df: pd.DataFrame):
        if "accid" in df.columns:
            df = df.drop_duplicates(subset=["accid"])
        if "protein_id" in df.columns:
            df = df.drop_duplicates(subset=["protein_id"])
        if "accession_id" in df.columns:
            df = df.drop_duplicates(subset=["accession_id"])

        return df

    def group_nodes_by_sample_sources(self):
        nodes_by_sample_sources = {}
        for node in self.current_nodes:
            sample_source = node.run_manager.run_engine.sample.sources_list()
            if sample_source not in nodes_by_sample_sources.keys():
                nodes_by_sample_sources[sample_source] = []
            nodes_by_sample_sources[sample_source].append(node)
        return nodes_by_sample_sources

    def group_nodes_by_source_and_parameters(self) -> List[List[Tree_Node]]:

        source_paramaters_combinations = {}

        def check_node_in_combination_dict(node: Tree_Node):
            sample_source = node.run_manager.run_engine.sample.sources_list()
            parameters = node.parameters
            for source, register in source_paramaters_combinations.items():
                if sample_source == register["source"] and parameters.equals(
                    register["parameters"]
                ):
                    return source

            return None

        def combination_dict_new_entry(node, sample_source, parameters):
            new_entry_idx = len(source_paramaters_combinations)
            source_paramaters_combinations[new_entry_idx] = {
                "source": sample_source,
                "parameters": parameters,
                "nodes": [node],
            }

        def update_combination_dict(node: Tree_Node):
            sample_source = node.run_manager.run_engine.sample.sources_list()
            parameters = node.parameters
            source = check_node_in_combination_dict(node)

            if source is None:
                combination_dict_new_entry(node, sample_source, parameters)
            else:
                source_paramaters_combinations[source]["nodes"].append(node)

        for node in self.current_nodes:
            sample_source = node.run_manager.run_engine.sample.sources_list()
            print("node: " + str(node.node_index))
            print("Sample source: " + str(sample_source))
            update_combination_dict(node)

        print("########### GROUPED NODES ###########")
        print(source_paramaters_combinations)

        grouped_nodes = [
            register["nodes"]
            for source, register in source_paramaters_combinations.items()
        ]

        return grouped_nodes

    def group_nodes_by_module(self):
        nodes_by_module = {}
        for node in self.current_nodes:
            if node.module not in nodes_by_module:
                nodes_by_module[node.module] = []
            nodes_by_module[node.module].append(node)
        return nodes_by_module

    def process_subject(self, volonteer: Tree_Node, group_targets: pd.DataFrame):
        original_targets = copy.deepcopy(
            volonteer.run_manager.run_engine.merged_targets
        )

        volonteer.run_manager.update_merged_targets(group_targets)

        volonteer.run_manager.run_main()

        volonteer.run_manager.update_merged_targets(original_targets)

        mapped_instances_shared = (
            volonteer.run_manager.run_engine.remap_manager.mapped_instances
        )

        return mapped_instances_shared

    def stacked_deployement(self, nodes_by_sample_sources: List[List[Tree_Node]]):

        current_nodes = []

        for nodes in nodes_by_sample_sources:

            if len(nodes) == 0:
                continue

            nodes = self.get_remap_plans(nodes)

            group_targets = self.get_node_node_targets(nodes)
            volonteer = nodes[0]

            mapped_instances_shared = self.process_subject(volonteer, group_targets)

            nodes = self.update_mapped_instances(nodes, mapped_instances_shared)

            current_nodes.extend(nodes)

        if len(current_nodes) > 0:
            self.current_nodes = current_nodes

    @staticmethod
    def update_mapped_instances(
        nodes_to_update: List[Tree_Node],
        mapped_instances_shared: List[Mapping_Instance],
    ):
        new_nodes = []

        for node in nodes_to_update:
            node.run_manager.run_engine.update_mapped_instances(mapped_instances_shared)
            new_nodes.append(node)

        return new_nodes

    def run_simplified_mapping(self):
        nodes_by_sample_sources = self.group_nodes_by_source_and_parameters()

        self.stacked_deployement(nodes_by_sample_sources)

    def update_nodes(self):
        new_nodes = []
        for node in self.current_nodes:
            children = node.children
            for child in children:
                new_node = self.spawn_node_child(node, child)
                new_nodes.append(new_node)

        print("updating nodes")
        print(len(new_nodes))
        self.current_nodes = new_nodes

        if len(new_nodes) == 0:
            self.current_module = "end"
        else:

            self.determine_current_module()

    def run_current_nodes(self):

        for node in self.current_nodes:
            self.run_node(node)

    def run_node(self, node: Tree_Node):
        node.run_manager.run_main()

    def run_nodes_sequential(self):
        self.run_current_nodes()
        self.register_current_nodes()
        self.update_nodes()

    def run_nodes_simply(self):
        self.run_simplified_mapping()
        self.register_current_nodes()
        self.update_nodes()

    def register_current_nodes(self):
        for node in self.current_nodes:

            self.register_node_leaves(node)

    def deploy_nodes(self):
        if self.current_module == "end":
            return
        if self.current_module == ConstantsSettings.PIPELINE_NAME_read_quality_analysis:
            self.update_nodes()
        if self.current_module == ConstantsSettings.PIPELINE_NAME_remapping:
            self.run_nodes_simply()
        else:
            self.run_nodes_sequential()

    def run_current_nodes_batch_parallel(self, batch=2):

        import multiprocessing as mp

        node_batches = [
            self.current_nodes[i : i + batch]
            for i in range(0, len(self.current_nodes), batch)
        ]

        for node_batch in node_batches:
            nproc = len(node_batch)
            pool = mp.Pool(nproc)
            drones = [
                pool.apply_async(node.run_manager.run_main) for node in node_batch
            ]

            for drone in drones:
                drone.get()

            for node in node_batch:
                self.register_node_leaves(node)

            pool.close()
            pool.join()


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

    def temp_reads(self, fofn):
        utils = Utils()
        r1, r2 = self.read_fofn(fofn)
        temp_dir = utils.get_temp_dir()

        new_r1 = os.path.join(temp_dir, os.path.basename(r1))
        new_r2 = None

        shutil.copy(r1, new_r1)

        if r2:
            new_r2 = os.path.join(temp_dir, os.path.basename(r2))
            shutil.copy(r2, new_r2)

        return new_r1, new_r2

    def create_sample_from_fofn(self, fofn, user: User, technology: str):

        r1, r2 = self.read_fofn(fofn)
        name = os.path.basename(fofn).split(".")[0]

        try:
            sample = Sample.objects.get(name=name, owner=user)

        except Sample.DoesNotExist:
            tmp_r1, tmp_r2 = self.temp_reads(fofn)
            sample = self.sample_save(name, user, tmp_r1, tmp_r2, technology)
            self.move_sample(sample, user)
            self.sample_preprocess(sample, user)

        return sample

    def sample_save(self, name, user, r1, r2, technology):
        utils = Utils()
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