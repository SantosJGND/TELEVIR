import argparse
import itertools as it
import json
import os
import shutil

import numpy as np
import pandas as pd
from background_task import background
from pathogen_detection.constants_settings import TestConstants
from pathogen_detection.run_main import RunMain_class
from pathogen_detection.update_DBs import (
    Update_QC_report,
    Update_Sample,
    Update_Sample_Runs,
)

from product.constants_settings import ConstantsSettings
from product.models import Fastq_Input, Processed, Submitted
from product.update_DBs import Update_project


def simplify_name(name):
    return (
        name.replace("_", "_")
        .replace("-", "_")
        .replace(" ", "_")
        .replace(".", "_")
        .lower()
    )


class AttrDict(dict):
    def __init__(self, *args, **kwargs):
        super(AttrDict, self).__init__(*args, **kwargs)
        self.__dict__ = self


class Dummy_deployment:

    project: str
    prefix: str
    rdir: str
    threads: int = 3
    run_engine: RunMain_class
    params = dict
    run_params_db = pd.DataFrame()
    pk: int = 0
    username: str

    def __init__(
        self,
        project: str = "test",
        prefix: str = "main",
        username: str = "admin",
        pk: int = 0,
        static_dir: str = "static",
    ) -> None:

        self.username = username
        self.project = project
        self.prefix = prefix
        self.rdir = os.path.join(ConstantsSettings.product_directory, username, project)
        self.dir = os.path.join(self.rdir, prefix)

        os.makedirs(self.dir, exist_ok=True)

        self.prefix = prefix
        self.pk = pk
        self.static_dir = static_dir

    def sample_main(
        self,
        sample=1,
        cols=["PREPROCESS", "ENRICHMENT", "ASSEMBLY", "CONTIG_CLASSIFICATION"],
    ):
        """
        sample module / software combinations from dictionaries in params.py. random.
        :param sample: how many combinations to sample. corresponds to number of metclass run directories to be created.
        :param cols: keys to sample combinations of in SOFTWARE dict in params.py
        :return: data frame.
        """
        if len(cols) == 0:
            cols = list(self.params_dict.SOFTWARE.keys())

        venue = [self.params_dict.SOFTWARE[x] for x in cols]
        venues = list(it.product(*venue))

        if sample > 0 and sample < len(venues):

            vex = np.random.choice(list(range(len(venues))), sample, replace=False)
            venues = [venues[x] for x in vex]

        venues = pd.DataFrame(venues)
        venues.columns = cols
        #
        return venues

    def params_extract(self, show, modules=[], sample=1):
        """
        takes list of software, which might or not have entries in the argument dictionaries.
        """

        if len(modules) == 0:
            modules = list(self.params_dict.SOFTWARE.keys())

        relate = []
        new_features = []
        nvens = []
        #
        for ix, soft in enumerate(show):
            soft_module = modules[ix]
            if soft in self.modules_to_stores[soft_module].keys():
                for c, g in self.modules_to_stores[soft_module][soft].items():
                    relate.append([soft_module, soft])
                    new_features.append(c)
                    nvens.append(g)
        #
        nvens = list(it.product(*nvens))
        if sample:
            vex = np.random.choice(list(range(len(nvens))), sample, replace=False)
            nvens = [nvens[x] for x in vex]

        relate = pd.DataFrame(relate, columns=["module", "software"])
        nvens = pd.DataFrame(nvens).reset_index(drop=True)
        nvens.columns = new_features

        return nvens, relate

    def generate_combinations(self, ncomb: int = 0, modules: list = []):

        hdconf = self.sample_main(sample=ncomb, cols=modules)
        params2 = {}
        paramCombs = [
            self.params_extract(hdconf.iloc[idx], modules=modules, sample=0)
            for idx in range(hdconf.shape[0])
        ]
        linked_dbs = [x[1] for x in paramCombs]
        paramCombs = [x[0] for x in paramCombs]

        return hdconf, linked_dbs, paramCombs

    @staticmethod
    def extract_parameters(
        parameters_df_list: list, linked_db_list: list, common_index: tuple
    ):
        """
        extract parameters from list of data frames.
        :param parameters_df_list: list of data frames.
        :param linked_db_list: list of data frames.
        :param common_index: index of common index in data frames.

        :return: data frame.
        """

        params = (
            parameters_df_list[common_index[0]]
            .loc[[common_index[1]]]
            .reset_index(drop=True)
        )
        params = pd.DataFrame([params.columns, params.loc[0]]).T
        params = pd.concat(
            (linked_db_list[common_index[0]], params), axis=1
        ).reset_index(drop=True)
        params.columns = ["module", "software", "param", "value"]

        return params

    def parameters_generate_random(self):
        hdconf, linked_dbs, paramCombs = self.generate_combinations(ncomb=1)

        params = self.extract_parameters(paramCombs, linked_dbs, common_index=(0, 0))
        return params

    def get_params_ont(self):

        with open(ConstantsSettings.ont_params_json, "r") as f:
            params = json.load(f)
            params = AttrDict(params)

        self.params_dict = params

    def get_params_illumina(self):

        with open(ConstantsSettings.illumina_params_json, "r") as f:
            params = json.load(f)
            params = AttrDict(params)

        self.params_dict = params

    def generate_config_file(self):

        self.config = {
            "project": self.project,
            "source": self.params_dict.SOURCE,
            "directories": {
                "root": self.rdir,
            },
            "static_dir": self.static_dir,
            "threads": self.threads,
            "prefix": self.prefix,
            "project_name": self.project,
            "metadata": {
                x: os.path.join(self.params_dict.METADATA["ROOT"], g)
                for x, g in self.params_dict.METADATA.items()
            },
            "technology": self.params_dict.DATA_TYPE,
            "bin": self.params_dict.BINARIES,
            "actions": {},
        }

        for dr, g in self.params_dict.DIRS.items():
            self.config["directories"][dr] = os.path.join(self.dir, g)

        for dr, g in self.params_dict.ACTIONS.items():
            self.config["actions"][dr] = True
            self.config["actions"]["DEPLETE"] = False
            self.config["actions"]["CLEAN"] = False

        self.config.update(self.params_dict.CONSTANTS)

        self.modules_to_stores = {
            "PREPROCESS": self.params_dict.ARGS_QC,
            "ENRICHMENT": self.params_dict.ARGS_ENRICH,
            "ASSEMBLY": self.params_dict.ARGS_ASS,
            "CONTIG_CLASSIFICATION": self.params_dict.ARGS_CLASS,
            "READ_CLASSIFICATION": self.params_dict.ARGS_CLASS,
            "REMAPPING": self.params_dict.ARGS_REMAP,
        }

        self.run_params_db = self.parameters_generate_random()

    def configure_ont(self, sample_path: str) -> None:
        self.get_params_ont()
        self.generate_config_file()
        self.prep_test_env()

        sample_name = os.path.basename(sample_path)
        new_sample_path = os.path.join(self.dir, "reads") + "/" + sample_name

        shutil.copy(sample_path, new_sample_path)

        self.config["sample_name"] = simplify_name(sample_name)
        self.config["r1"] = new_sample_path
        self.config["r2"] = "none"
        self.config["type"] = "SE"
        self.config["technology"] = "nanopore"

    def configure_illumina(self, r1_path: str, r2_path: str = "") -> None:

        self.get_params_illumina()
        self.generate_config_file()
        self.prep_test_env()

        r1_name = os.path.basename(r1_path)
        new_r1_path = os.path.join(self.dir, "reads") + "/" + r1_name
        shutil.copy(r1_path, new_r1_path)

        if r2_path:
            r2_name = os.path.basename(r1_path)
            new_r2_path = os.path.join(self.dir, "reads") + "/" + r2_name
            shutil.copy(r2_path, new_r2_path)

        self.config["sample_name"] = simplify_name(r1_name)
        self.config["r1"] = new_r1_path
        self.config["r2"] = new_r2_path
        self.config["type"] = ["SE", "PE"][
            int(os.path.isfile(self.config["r2"] is not None))
        ]
        self.config["technology"] = "illumina"

    def prep_test_env(self):
        """
        from main directory bearing scripts, params.py and main.sh, create metagenome run directory

        :return:
        """
        #

        for dir in self.config["directories"].values():
            os.makedirs(dir, exist_ok=True)

    def close(self):
        if os.path.exists(self.dir):
            shutil.rmtree(self.dir)

    def run_main_prep(self):

        os.makedirs(self.static_dir, exist_ok=True)

        self.run_engine = RunMain_class(self.config, self.run_params_db, self.username)


class Run_Main_from_Fastq_Input:

    user: str = "test"
    file_r1: str
    file_r2: str
    technology: str
    email: str
    name: str
    project_name: str
    description: str
    date_created: str
    date_modified: str
    pk: int

    container: Dummy_deployment

    def __init__(self, pk):
        print(f"pk init: {pk}")
        input_data = Fastq_Input.objects.get(pk=pk)

        self.file_r1 = input_data.file_r1.path
        if input_data.file_r2:
            self.file_r2 = input_data.file_r2.path
        self.technology = input_data.technology
        self.user = input_data.name
        self.project_name = input_data.project_name
        self.date_created = input_data.date_created
        self.date_modified = input_data.date_modified
        self.pk = input_data.pk

        self.static_dir = os.path.join(
            ConstantsSettings.static_directory_product,
            self.user,
            self.project_name,
            simplify_name(os.path.basename(self.file_r1)),
        )

        self.static_dir = os.path.join(
            ConstantsSettings.static_directory, self.static_dir
        )

        self.container = Dummy_deployment(
            project=self.project_name,
            prefix=f"{self.project_name}_{self.pk}",
            username=self.user,
            static_dir=self.static_dir,
        )

    def configure(self):

        Update_project(self.container.dir, user=self.user, submit_index=self.pk)

        if self.technology == "illumina":
            self.container.configure_illumina(self.file_r1, self.file_r2)
        else:
            self.container.configure_ont(self.file_r1)

    def Deploy(self):
        self.container.run_main_prep()
        self.container.run_engine.Run()
        self.container.run_engine.Summarize()
        self.container.run_engine.generate_output_data_classes()

    def Update_dbs(self):

        Update_Sample(
            self.container.run_engine.sample,
        )

        Update_QC_report(self.container.run_engine.sample)
        Update_Sample_Runs(self.container.run_engine)

    def register_submission(self):
        try:
            submitted = Submitted.objects.get(pk=self.pk)
        except Submitted.DoesNotExist:
            fastq_input = Fastq_Input.objects.get(pk=self.pk)
            submitted = Submitted(
                fastq_input=fastq_input,
            )

            submitted.save()

    def register_completion(self):

        try:
            processed = Processed.objects.get(pk=self.pk)
        except Processed.DoesNotExist:
            fastq_input = Fastq_Input.objects.get(pk=self.pk)
            processed = Processed(
                fastq_input=fastq_input,
            )
            processed.save()

    def Submit(self):

        self.register_submission()
        self.configure()
        self.Deploy()
        self.Update_dbs()
        self.container.close()
        self.register_completion()
