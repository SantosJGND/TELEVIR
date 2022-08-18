import json
import os
import random
import shutil

from django.test import TestCase
from pathogen_detection.constants_settings import TestConstants
from pathogen_detection.object_classes import Read_class, RunCMD
from pathogen_detection.run_main import get_bindir_from_binaries
from pyrsistent import T
from this import d


class AttrDict(dict):
    def __init__(self, *args, **kwargs):
        super(AttrDict, self).__init__(*args, **kwargs)
        self.__dict__ = self


class Dummy_deployment:

    rdir: str
    threads: int = 3

    def __init__(self, project: str = "test", prefix: str = "main") -> None:
        self.rdir = os.path.join(TestConstants.Test_Temp_Directory, project)
        self.prefix = prefix

        with open(TestConstants.ont_params_json, "r") as f:
            params = json.load(f)
            params = AttrDict(params)

        self.config = {
            "project": project,
            "source": params.SOURCE,
            "directories": {
                "root": self.rdir,
            },
            "static_dir": self.static_dir,
            "threads": self.threads,
            "prefix": self.prefix,
            "project_name": project,
            "metadata": {
                x: os.path.join(params.METADATA["ROOT"], g)
                for x, g in params.METADATA.items()
            },
            "technology": params.DATA_TYPE,
            "bin": params.BINARIES,
        }

        for dr, g in params.DIRS.items():
            self.config["directories"][dr] = self.dir + g
        for dr, g in self.actions.items():
            self.config["actions"][dr] = False
        self.config.update(params.CONSTANTS)

    def configure_ont(self, sample_path: str) -> None:
        sample_name = os.path.basename(sample_path)
        shutil.copy(sample_path, os.path.join(self.dir, "reads"))
        sample_path = os.path.join(self.dir, "reads", sample_name)

        self.config["sample_name"] = sample_name
        self.config["r1"] = sample_path
        self.config["r2"] = "none"
        self.config["type"] = "SE"
        self.config["technology"] = "nanopore"

    def configure_illumina(self, r1_path: str, r2_path: str = "") -> None:
        sample_name = os.path.basename(sample_path)
        shutil.copy(sample_path, os.path.join(self.dir, "reads"))
        sample_path = os.path.join(self.dir, "reads", sample_name)
        self.config["sample_name"] = sample_name
        self.config["r1"] = r1_path
        self.config["r2"] = r2_path
        self.config["type"] = ["SE", "PE"][int(os.path.isfile(self.config["r2"]))]
        self.config["technology"] = "illumina"

    def prep_test_env(self, rdir=""):
        """
        from main directory bearing scripts, params.py and main.sh, create metagenome run directory

        :return:
        """
        #
        if not self.rdir:
            rdir = os.getcwd()

        self.dir = rdir + "{}/".format(self.prefix)

        for dir in self.config["directories"].values():
            os.system("mkdir -p " + self.dir + dir)

class test_runcmd(TestCase):
    
    container= Dummy_deployment()

    def SetUp(self):
        self.cmd= RunCMD(
            get_bindir_from_binaries(self.container.config["bin"], "PREPROCESS"),
            self.container.config["directories"]["log_dir"],
            prefix= "test", 
            task= "test"
        )

    def test_runcmd_config(self):
        self.assertEqual(self.cmd.logger.level, "CRITICAL")
        self.assertEqual(self.cmd.logger.propagate, False)
        self.assertEqual(self.cmd.logfile, os.path.join(self.container.config["directories"]["log_dir"], "test_test.log"))
        self.assertEqual(self.cmd.logdir, self.container.config["directories"]["log_dir"])
        self.assertEqual(self.cmd.prefix, "test")
    
    def test_flag_error_success(self):



class test_read_class(TestCase):

    container: Dummy_deployment = Dummy_deployment()

    def SetUp(self):

        self.cmd = RunCMD(
            get_bindir_from_binaries(self.container.config["bin"], "PREPROCESS")
        )

        self.container.configure_ont(TestConstants.ont_fastq_gz_file_path)

        self.r1 = Read_class(
            self.container.config["r1"],
            self.container.config["directories"]["PREPROCESS"],
            self.container.config["directories"]["reads_enriched_dir"],
            self.container.config["directories"]["reads_depleted_dir"],
            bin=get_bindir_from_binaries(self.container.config["bin"], "PREPROCESS"),
        )

    def test_read_class_configuration(self):

        self.assertEqual(self.r1.filepath, TestConstants.ont_fastq_gz_file_path)
        self.assertEqual(self.r1.current, TestConstants.ont_fastq_gz_file_path)
        self.assertEqual(
            self.r1.prefix,
            os.path.splitext(os.path.basename(TestConstants.ont_fastq_gz_file_path)),
        )

        self.assertEqual(
            self.r1.clean,
            os.path.join(
                self.container.config["directories"]["PREPROCESS"],
                self.container.prefix + ".clean.fastq.gz",
            ),
        )

        self.assertEqual(
            self.r1.enriched,
            os.path.join(
                self.container.config["directories"]["PREPROCESS"],
                self.container.prefix + ".enriched.fastq.gz",
            ),
        )

        self.assertEqual(
            self.r1.depleted,
            os.path.join(
                self.container.config["directories"]["PREPROCESS"],
                self.container.prefix + ".depleted.fastq.gz",
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

        self.assertEqual(self.r1.current, r1_copy.clean)
        self.assertEqual(self.r1.enriched, r1_copy.enriched)
        self.assertEqual(self.r1.depleted, r1_copy.depleted)

    def test_fake_quality_check(self):
        temp_path = os.path.join(self.container.dir, "temp.fastq")
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
            temp_r1.current
            "|"
            "awk",
            "'NR % 4 == 0'",
            "|",
            "head",
        ]

        qual_grep = self.cmd.run_bash_return(cmd_qual_zgrep)
        self.assertEqual(list(set(qual_grep[0].strip()))[0], "3")
        os.remove(temp_path)

    def test_read_filter_move(self):
        temp_path = os.path.join(self.container.dir, "temp.fastq")
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

        temp_reads = temp_r1.get_read_list()
        self.assertEqual(len(temp_reads), len(r1_reads_sample))

    def test_enrich(self):
        r1_reads = self.r1.get_read_list()
        r1_reads_sample = random.sample(r1_reads, 10)

        self.r1.enrich(r1_reads_sample)

        enriched_reads = self.r1.get_read_list()
        self.assertEqual(len(enriched_reads), len(r1_reads_sample))
        self.assertEqual(self.r1.current, self.r1.enriched)
        self.assertEqual(self.r1.current_status, "enriched")
        self.assertEqual(self.r1.read_number_enriched, len(r1_reads_sample))

    def test_deplete(self):
        r1_reads = self.r1.get_read_list()
        r1_reads_sample = random.sample(r1_reads, 10)

        self.r1.deplete(r1_reads_sample)

        depleted_reads = self.r1.get_read_list()
        self.assertEqual(len(depleted_reads), len(r1_reads_sample))
        self.assertEqual(self.r1.current, self.r1.depleted)
        self.assertEqual(self.r1.current_status, "depleted")
        self.assertEqual(self.r1.read_number_depleted, len(r1_reads_sample))

    def test_read_number_count(self):

        read_list = self.r1.get_read_list()
        read_number = self.r1.current_fastq_read_number()
        self.assertEqual(read_number, len(read_list))


class test_pathid_deployment(TestCase):
    pass
