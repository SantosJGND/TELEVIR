import itertools as it
import json
import os
import pickle
import random
import shutil

import numpy as np
import pandas as pd
from django.test import TestCase
from pathogen_detection.constants_settings import TestConstants
from pathogen_detection.object_classes import (
    Read_class,
    RunCMD,
    Sample_runClass,
    Software_detail,
)
from pathogen_detection.preprocess_class import Preprocess
from pathogen_detection.run_main import RunMain_class, get_bindir_from_binaries


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

    def __init__(self, project: str = "test", prefix: str = "main") -> None:
        self.project = project
        self.prefix = prefix
        self.rdir = os.path.join(TestConstants.Test_Temp_Directory, project)
        self.dir = os.path.join(self.rdir, prefix)

        os.makedirs(self.dir, exist_ok=True)
        self.static_dir = os.path.join(TestConstants.Test_Temp_Directory, "static")
        self.prefix = prefix

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

        with open(TestConstants.ont_params_json, "r") as f:
            params = json.load(f)
            params = AttrDict(params)

        self.params_dict = params

    def get_params_illumina(self):

        with open(TestConstants.illumina_params_json, "r") as f:
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
            self.config["actions"][dr] = False

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

        self.config["sample_name"] = sample_name
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

        self.config["sample_name"] = r1_name
        self.config["r1"] = new_r1_path
        self.config["r2"] = new_r2_path
        self.config["type"] = ["SE", "PE"][int(os.path.isfile(self.config["r2"]))]
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

        self.run_engine = RunMain_class(self.config, self.run_params_db)


class test_runcmd(TestCase):

    container = Dummy_deployment(prefix="runCMD")

    def setUp(self):

        self.container.configure_ont(TestConstants.ont_fastq_gz_file_path)

        self.cmd = RunCMD(
            get_bindir_from_binaries(self.container.config["bin"], "PREPROCESS"),
            self.container.config["directories"]["log_dir"],
            prefix="test",
            task="test",
        )

    def tearDown(self) -> None:
        self.container.close()
        return super().tearDown()

    def test_runcmd_config(self):
        self.assertEqual(self.cmd.logger.level, 50)
        self.assertEqual(self.cmd.logger.propagate, False)
        self.assertEqual(
            self.cmd.logfile,
            os.path.join(
                self.container.config["directories"]["log_dir"], "test_test.log"
            ),
        )
        self.assertEqual(
            self.cmd.logdir, self.container.config["directories"]["log_dir"]
        )
        self.assertEqual(self.cmd.prefix, "test")

    def test_flag_error(self):
        for flag in self.cmd.error_flags:
            self.assertEqual(self.cmd.flag_error(flag), True)


class test_read_class(TestCase):

    container: Dummy_deployment = Dummy_deployment(prefix="Read_class")

    def setUp(self):

        self.container.configure_ont(TestConstants.ont_fastq_gz_file_path)

        self.r1 = Read_class(
            self.container.config["r1"],
            self.container.config["directories"]["PREPROCESS"],
            self.container.config["directories"]["reads_enriched_dir"],
            self.container.config["directories"]["reads_depleted_dir"],
            bin=get_bindir_from_binaries(self.container.config["bin"], "PREPROCESS"),
        )
        self.cmd = RunCMD(
            get_bindir_from_binaries(self.container.config["bin"], "PREPROCESS")
        )

    def tearDown(self) -> None:

        self.container.close()

        return super().tearDown()

    def test_read_class_configuration(self):

        self.assertEqual(
            self.r1.prefix,
            os.path.splitext(os.path.basename(TestConstants.ont_fastq_gz_file_path))[0],
        )

        self.assertEqual(
            self.r1.clean,
            os.path.join(
                self.container.config["directories"]["PREPROCESS"],
                os.path.splitext(self.container.config["sample_name"])[0]
                + ".clean.fastq.gz",
            ),
        )

        self.assertEqual(
            self.r1.enriched,
            os.path.join(
                self.container.config["directories"]["reads_enriched_dir"],
                os.path.splitext(self.container.config["sample_name"])[0]
                + ".enriched.fastq.gz",
            ),
        )

        self.assertEqual(
            self.r1.depleted,
            os.path.join(
                self.container.config["directories"]["reads_depleted_dir"],
                os.path.splitext(self.container.config["sample_name"])[0]
                + ".depleted.fastq.gz",
            ),
        )

    def test_readclass_update(self):

        copy_child = lambda obj: pickle.loads(pickle.dumps(obj))

        r1_copy = copy_child(self.r1)

        r1_copy.update(
            self.container.config["directories"]["PREPROCESS"],
            self.container.config["directories"]["reads_enriched_dir"],
            self.container.config["directories"]["reads_depleted_dir"],
        )

        self.assertEqual(self.r1.current, r1_copy.current)
        self.assertEqual(self.r1.enriched, r1_copy.enriched)
        self.assertEqual(self.r1.depleted, r1_copy.depleted)

    def test_fake_quality_check(self):
        temp_path = os.path.join(self.container.dir, "temp.fastq.gz")
        self.r1.copy(temp_path)
        temp_r1 = Read_class(
            temp_path,
            self.container.config["directories"]["PREPROCESS"],
            self.container.config["directories"]["reads_enriched_dir"],
            self.container.config["directories"]["reads_depleted_dir"],
            bin=get_bindir_from_binaries(self.container.config["bin"], "PREPROCESS"),
        )
        temp_r1.fake_quality()

        cmd_qual_zgrep = [
            "zcat",
            temp_r1.current,
            "|",
            "awk",
            "'NR % 4 == 0'",
        ]

        qual_grep = self.cmd.run_bash_return(cmd_qual_zgrep).decode().split("\n")
        self.assertEqual(list(set(qual_grep[0])), ["3"])
        os.remove(temp_path)

    def test_read_filter_move(self):
        temp_path = os.path.join(self.container.dir, "temp.fastq.gz")
        r1_reads = self.r1.get_read_list()
        r1_reads_sample = random.sample(r1_reads, 10)

        temp_r1 = Read_class(
            temp_path,
            self.container.config["directories"]["PREPROCESS"],
            self.container.config["directories"]["reads_enriched_dir"],
            self.container.config["directories"]["reads_depleted_dir"],
            bin=get_bindir_from_binaries(self.container.config["bin"], "PREPROCESS"),
        )

        self.r1.read_filter_move(self.r1.current, r1_reads_sample, temp_r1.current)
        temp_r1.exists = True

        temp_reads = temp_r1.get_read_list()
        self.assertEqual(len(temp_reads), len(r1_reads_sample))

    def test_enrich(self):
        r1_reads = self.r1.get_read_list()
        r1_reads_sample = random.sample(r1_reads, 10)

        self.r1.enrich(r1_reads_sample)

        enriched_reads = self.r1.get_read_list()
        self.assertEqual(len(enriched_reads), self.r1.read_number_enriched)
        self.assertEqual(self.r1.current, self.r1.enriched)
        self.assertEqual(self.r1.current_status, "enriched")
        self.assertEqual(self.r1.read_number_enriched, len(r1_reads_sample))

    def test_deplete(self):
        r1_reads = self.r1.get_read_list()
        r1_reads_sample = random.sample(r1_reads, 10)

        self.r1.deplete(r1_reads_sample)

        depleted_reads = self.r1.get_read_list()

        self.assertEqual(len(depleted_reads), self.r1.read_number_depleted)
        self.assertEqual(self.r1.current, self.r1.depleted)
        self.assertEqual(self.r1.current_status, "depleted")
        self.assertEqual(self.r1.read_number_depleted, len(r1_reads_sample))

    def test_read_number_count(self):

        read_list = self.r1.get_read_list()
        read_number = self.r1.current_fastq_read_number()
        self.assertEqual(read_number, len(read_list))


class test_sample_class(TestCase):

    container: Dummy_deployment = Dummy_deployment(prefix="Sample_runClass")

    def setUp(self):

        self.container.configure_ont(TestConstants.ont_fastq_gz_file_path)

        self.r1 = Read_class(
            self.container.config["r1"],
            self.container.config["directories"]["PREPROCESS"],
            self.container.config["directories"]["reads_enriched_dir"],
            self.container.config["directories"]["reads_depleted_dir"],
            bin=get_bindir_from_binaries(self.container.config["bin"], "PREPROCESS"),
        )

        self.r2 = Read_class(
            self.container.config["r2"],
            self.container.config["directories"]["PREPROCESS"],
            self.container.config["directories"]["reads_enriched_dir"],
            self.container.config["directories"]["reads_depleted_dir"],
            bin=get_bindir_from_binaries(self.container.config["bin"], "PREPROCESS"),
        )

        self.sample = Sample_runClass(
            self.r1,
            self.r2,
            self.container.config["sample_name"],
            self.container.config["project_name"],
            self.container.config["technology"],
            self.container.config["type"],
            0,
            ",".join(
                [os.path.basename(self.r1.current), os.path.basename(self.r2.current)]
            ),
            bin=get_bindir_from_binaries(self.container.config["bin"], "PREPROCESS"),
            threads=self.container.config["threads"],
        )

        self.container.configure_ont(TestConstants.ont_fastq_gz_file_path)

        self.cmd = RunCMD(
            get_bindir_from_binaries(self.container.config["bin"], "PREPROCESS")
        )

    def tearDown(self) -> None:
        self.container.close()
        return super().tearDown()

    def test_sample_init(self):
        self.assertEqual(self.sample.r1.current, self.r1.current)
        self.assertEqual(self.sample.r2.current, self.r2.current)
        self.assertEqual(self.sample.r1.current_status, "raw")
        self.assertEqual(self.sample.r2.current_status, "raw")
        self.assertEqual(self.sample.sample_name, self.container.config["sample_name"])
        self.assertEqual(
            self.sample.project_name, self.container.config["project_name"]
        )
        self.assertEqual(self.sample.technology, self.container.config["technology"])
        self.assertEqual(self.sample.type, self.container.config["type"])
        self.assertEqual(self.sample.combinations, 0)
        self.assertEqual(self.sample.threads, self.container.config["threads"])
        self.assertEqual(
            self.sample.reads_before_processing,
            self.r1.read_number_raw + self.r2.read_number_raw,
        )

    def test_sample_clean_unique(self):
        unique_reads = self.r1.get_read_list() + self.r2.get_read_list()
        unique_reads = list(set(unique_reads))

        self.sample.clean_unique()

        post_process_reads = self.r1.get_read_list() + self.r2.get_read_list()
        self.assertEqual(len(post_process_reads), len(unique_reads))

    def test_trimmomatic_sort(self):
        unique_reads = self.r1.get_read_list() + self.r2.get_read_list()
        unique_reads = list(set(unique_reads))

        self.sample.trimmomatic_sort()

        post_process_reads = self.r1.get_read_list() + self.r2.get_read_list()
        self.assertEqual(len(post_process_reads), len(unique_reads))


class test_pathid_Preprocess(TestCase):

    container_ont = Dummy_deployment(prefix="Preprocess_ont")
    container_illumina = Dummy_deployment(prefix="Preprocess_illumina")
    preprocess_engine = Preprocess
    cmd: RunCMD

    def setUp(self):

        self.container_ont.configure_ont(TestConstants.ont_fastq_gz_file_path)

        self.container_ont.run_main_prep()

        self.preprocess_engine_ont = Preprocess(
            self.container_ont.run_engine.sample.r1.current,
            self.container_ont.run_engine.sample.r2.current,
            self.container_ont.run_engine.filtered_reads_dir,
            self.container_ont.run_engine.type,
            self.container_ont.run_engine.preprocess_method,
            self.container_ont.run_engine.sample.r1.clean,
            self.container_ont.run_engine.sample.r2.clean,
            self.container_ont.run_engine.threads,
            self.container_ont.run_engine.subsample,
        )
        ##########
        ##########
        self.container_illumina.configure_illumina(
            TestConstants.illumina_fastq_gz_file_r1_path,
            TestConstants.illumina_fastq_gz_file_r2_path,
        )

        self.container_illumina.run_main_prep()

        self.preprocess_engine_illumina = Preprocess(
            self.container_illumina.run_engine.sample.r1.current,
            self.container_illumina.run_engine.sample.r2.current,
            self.container_illumina.run_engine.filtered_reads_dir,
            self.container_illumina.run_engine.type,
            self.container_illumina.run_engine.preprocess_method,
            self.container_illumina.run_engine.sample.r1.clean,
            self.container_illumina.run_engine.sample.r2.clean,
            self.container_illumina.run_engine.threads,
            self.container_illumina.run_engine.subsample,
        )

        self.cmd = RunCMD(
            get_bindir_from_binaries(self.container_ont.config["bin"], "PREPROCESS")
        )

    def tearDown(self) -> None:
        self.container_ont.close()
        # self.container_illumina.close()
        return super().tearDown()

    def test_check_file_not_empty(self):
        result = self.preprocess_engine_ont.check_gz_file_not_empty(
            self.container_ont.run_engine.r1.current
        )
        self.assertEqual(result, True)

        neg_file = self.container_ont.dir + "/temp_test.txt"
        open(neg_file, "w").close()

        negative_result = self.preprocess_engine_ont.check_gz_file_not_empty(neg_file)

        self.assertFalse(negative_result)

    def test_fastqc_input_SE(self):
        self.preprocess_engine_ont.fastqc_input("input_data")
        self.assertTrue(os.path.exists(self.preprocess_engine_ont.input_qc_report))

    def test_fastqc_input_PE(self):
        self.preprocess_engine_illumina.fastqc_input("input_data")
        self.assertTrue(os.path.exists(self.preprocess_engine_illumina.input_qc_report))

    def test_preprocess_illumina_SE(self):
        self.preprocess_engine_illumina.trimmomatic_SE()

        self.assertTrue(
            os.path.exists(self.preprocess_engine_illumina.preprocess_name_fastq_gz)
        )

    def test_preprocess_illumina_PE(self):
        self.preprocess_engine_illumina.trimmomatic_PE()

        self.assertTrue(
            os.path.exists(self.preprocess_engine_illumina.preprocess_name_fastq_gz)
        )
        self.assertTrue(
            os.path.exists(self.preprocess_engine_illumina.preprocess_name_r2_fastq_gz)
        )

    def test_preprocess_ont(self):
        self.preprocess_engine_ont.preprocess_QC()

        self.assertTrue(
            os.path.exists(self.preprocess_engine_ont.preprocess_name_fastq_gz)
        )
