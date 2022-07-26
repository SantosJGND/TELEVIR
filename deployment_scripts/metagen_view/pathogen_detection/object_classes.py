import logging
import os
import subprocess
from dataclasses import dataclass, field
from random import randint
from typing import Type

import pandas as pd

from pathogen_detection.utilities import fastqc_parse


class RunCMD:
    """
    Run command line commands.
    """

    def __init__(self, bin):
        """
        Initialize.
        """
        if bin:
            if bin[-1] != "/":
                bin += "/"

        self.bin = bin
        self.logs = []

        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(logging.INFO)
        self.logger.addHandler(logging.StreamHandler())

    @staticmethod
    def flag_error(subprocess_errorlog):
        """
        Check if error in subprocess.
        """

        if subprocess_errorlog:
            try:
                subprocess_errorlog = subprocess_errorlog.decode("utf-8")
                if "[error]" in subprocess_errorlog.lower():
                    return True
            except Exception as e:
                return False

        return False

    def run(self, cmd):
        """
        Run software.
        """

        if isinstance(cmd, list):
            cmd = " ".join(cmd)

        print(f"running command: {self.bin}{cmd}")

        proc_prep = subprocess.Popen(
            f"{self.bin}{cmd}",
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        out, err = proc_prep.communicate()

        if self.flag_error(err):
            self.logger.error(f"errror in command: {self.bin}{cmd}")
            raise Exception(err.decode("utf-8"))

        self.logs.append(cmd)
        self.logs.append(out)

    def run_python(self, cmd):
        """
        Run software with python.
        """

        if isinstance(cmd, list):
            cmd = " ".join(cmd)

        # self.logger.info(f"running command: python {self.bin}{cmd}")
        proc_prep = subprocess.Popen(
            f"python {self.bin}{cmd}",
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        out, err = proc_prep.communicate()

        if self.flag_error(err):
            self.logger.error(f"errror in command: python {self.bin}{cmd}")
            raise Exception(err.decode("utf-8"))

        self.logs.append(cmd)
        self.logs.append(out)

    def run_java(self, cmd):
        """
        Run software with java.
        """

        if isinstance(cmd, list):
            cmd = " ".join(cmd)

        # print(f"running command: java -cp {self.bin} {cmd}")
        proc_prep = subprocess.Popen(
            f"java -cp {self.bin} {cmd}",
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        out, err = proc_prep.communicate()

        if self.flag_error(err):
            self.logger.error(f"errror in command: java -cp {cmd}")
            raise Exception(err.decode("utf-8"))

        self.logs.append(cmd)
        self.logs.append(out)

    def run_bash(self, cmd):
        """
        Run bash command.
        """

        if isinstance(cmd, list):
            cmd = " ".join(cmd)

        proc_prep = subprocess.Popen(
            f"{cmd}",
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
        )
        out, err = proc_prep.communicate()

        if self.flag_error(err):
            self.logger.error(f"errror in command: {cmd}")
            raise Exception(err.decode("utf-8"))

        self.logs.append(cmd)
        self.logs.append(out)

    def run_bash_return(self, cmd):
        """
        Run bash command and return stdout.
        """

        if isinstance(cmd, list):
            cmd = " ".join(cmd)

        proc_prep = subprocess.Popen(
            f"{cmd}", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )
        out, err = proc_prep.communicate()

        if self.flag_error(err):
            self.logger.error(f"errror in command: {cmd}")
            raise Exception(err.decode("utf-8"))

        self.logs.append(cmd)
        self.logs.append(out)
        return out


class Read_class:
    def __init__(
        self, filepath, clean_dir: str, enriched_dir: str, depleted_dir: str, bin: str
    ):
        """
        Initialize.

        Args:
            filepath: path to file.
            clean_dir: path to clean directory.
            enriched_dir: path to enriched directory.
            depleted_dir: path to depleted directory.
            bin: path to bin directory.

        """
        self.cmd = RunCMD(bin)

        self.exists = os.path.isfile(filepath)

        self.filepath = filepath
        self.current = filepath
        self.prefix = self.determine_read_name(filepath)
        self.clean = os.path.join(clean_dir, self.prefix + ".clean.fastq.gz")
        self.enriched = os.path.join(enriched_dir, self.prefix + ".enriched.fastq.gz")
        self.depleted = os.path.join(depleted_dir, self.prefix + ".depleted.fastq.gz")
        self.current_status = "raw"
        self.read_number_raw = self.get_current_fastq_read_number()
        self.read_number_clean = 0
        self.read_number_enriched = 0
        self.read_number_depleted = 0

    def determine_read_name(self, filepath):

        if not self.exists:
            return "none"

        if "gz" not in filepath:
            self.cmd.run("bgzip -f %s" % filepath)

        filename = os.path.basename(filepath)

        filename = filename.replace(".gz", "")

        filename = filename.replace(".fasta", "").replace(".fa", "")

        return filename

    def read_filter_move(self, input: str, read_list: list, output: str = ""):

        if not self.exists:
            return

        temp_reads_keep = os.path.join(
            os.path.dirname(input), f"temp_{randint(1,1999)}.fq.gz"
        )
        with open(temp_reads_keep, "w") as f:
            f.write("\n".join(read_list))

        self.read_filter(input, output, temp_reads_keep)

        os.remove(temp_reads_keep)

    def read_filter_inplace(self, input: str, read_list: str):

        if not self.exists:
            return

        output_dir = os.path.dirname(input)
        tempreads = os.path.join(output_dir, f"temp_{randint(1,1999)}.fq.gz")

        self.read_filter(input, tempreads, read_list)

        if os.path.isfile(tempreads) and os.path.getsize(tempreads):
            os.remove(input)
            os.rename(tempreads, input)

    def read_filter(self, input: str, output: str, read_list: str):
        """
        filter read file using exisiting lsit of reads.
        Args:

            input: path to input file.
            output: path to output file.
            read_list: path to file containing read list.
        """
        if not self.exists:
            return

        cmd = "seqtk subseq %s %s | gzip > %s" % (input, read_list, output)

        self.cmd.run(cmd)

    def enrich(self, read_list):
        """
        filter reads and set current status to enriched.
        """
        self.read_filter_move(self.filepath, read_list, self.enriched)

        if os.path.isfile(self.enriched) and os.path.getsize(self.enriched):
            self.is_enriched()

    def deplete(self, read_list):
        """
        filter reads and aset current status to depleted.
        """
        self.read_filter_move(self.filepath, read_list, self.depleted)

        if os.path.isfile(self.depleted) and os.path.getsize(self.depleted):
            self.is_depleted()

    def is_clean(self):
        """
        Set current status to clean.
        """
        self.current = self.clean
        self.read_number_clean = self.get_current_fastq_read_number()
        self.current_status = "clean"
        self.filepath = os.path.dirname(self.current)

    def is_enriched(self):
        """
        Set current status to enriched.
        """
        self.current = self.enriched
        self.read_number_enriched = self.get_current_fastq_read_number()
        self.current_status = "enriched"
        self.filepath = os.path.dirname(self.current)

    def is_depleted(self):
        """
        Set current status to depleted.
        """
        self.current = self.depleted
        self.read_number_depleted = self.get_current_fastq_read_number()
        self.current_status = "depleted"
        self.filepath = os.path.dirname(self.current)

    def get_current_fastq_read_number(self):
        """
        Get number of reads in current fastq file."""

        if not self.exists:
            return 0

        cmd = "zcat %s | wc -l" % self.current
        rnumber = self.cmd.run_bash_return(cmd).decode("utf-8")
        return int(rnumber) // 4

    def current_fastq_read_number(self):
        """
        Get cached number of reads in current fastq file."""

        if self.current_status == "raw":
            return self.read_number_raw
        elif self.current_status == "clean":
            return self.read_number_clean
        elif self.current_status == "enriched":
            return self.read_number_enriched
        elif self.current_status == "depleted":
            return self.read_number_depleted
        else:
            return 0

    def __str__(self):
        return self.filepath


class Sample_runClass:

    report: str
    qcdata: dict
    reads_before_processing: int = 0
    reads_after_processing: int = 0
    sample_dir: str = ""

    def __init__(
        self,
        r1: Type[Read_class],
        r2: Type[Read_class],
        sample_name: str,
        project_name: str,
        technology: str,
        type: str,
        combinations: str,
        input: str,
        qc_soft: str,
        input_fastqc_report: os.PathLike,
        processed_fastqc_report: os.PathLike,
    ) -> None:
        self.r1 = r1
        self.r2 = r2
        self.sample_name = sample_name
        self.project_name = project_name
        self.technology = technology
        self.type = type
        self.combinations = combinations
        self.input = input
        self.qc_soft = qc_soft
        self.input_fastqc_report = input_fastqc_report
        self.processed_fastqc_report = processed_fastqc_report

        QCdir = os.path.dirname(self.r1.clean)

        self.get_qc_data(QCdir)

        self.reads_before_processing = self.r1.read_number_raw + self.r2.read_number_raw
        self.reads_after_processing = (
            self.r1.get_current_fastq_read_number()
            + self.r2.get_current_fastq_read_number()
        )

    def QC_summary(self, QCdir: str):
        """get QC summary
        runs fastqc parse on zip files in QCdir

        :param QCdir:
        :return: qc_summary
        """
        qc_summary = {
            "input": fastqc_parse(os.path.join(QCdir, "input_data.zip")),
            "processed": fastqc_parse(os.path.join(QCdir, "processed_data.zip")),
        }

        return qc_summary

    def get_qc_data(self, reads_dir: str = "reads/clean/"):
        """get qc data for sample_class. Update sample_class.qc_data.

        :param sample_class:
        :return: None
        """

        self.qcdata = self.QC_summary(reads_dir)


class Software_detail:
    def __init__(self, module, args_df: pd.DataFrame, config: dict, prefix: str):
        """

        Args:
            module: name of module.
            args_df: dataframe containing module arguments.
            config: dictionary containing module configuration.
            prefix: prefix of module.
        """
        if module not in args_df.module.unique():
            print(f"module {module} not found")
            self.module = "None"
            self.name = "None"
            self.args = "None"
            self.db = "None"
            self.bin = "None"
            self.dir = "None"
            self.output_dir = "None"
        else:
            method_details = args_df[(args_df.module == module)]
            self.module = module
            self.name = method_details.software.values[0]

            try:
                self.args = method_details[
                    method_details.param.str.contains("ARGS")
                ].value.values[0]
            except IndexError:
                self.args = ""

            try:
                db_name = method_details[
                    method_details.param.str.contains("DB")
                ].value.values[0]
                if ".gz" in db_name:
                    self.db = os.path.join(config["source"]["REF_FASTA"], db_name)
                else:
                    self.db = os.path.join(config["source"]["DBDIR_MAIN"], db_name)

            except IndexError:
                self.db = ""

            try:
                self.bin = os.path.join(
                    config["bin"]["ROOT"], config["bin"]["software"][self.name], "bin"
                )
            except KeyError:
                try:
                    self.bin = os.path.join(
                        config["bin"]["ROOT"], config["bin"][module]["default"], "bin"
                    )
                except KeyError:
                    self.bin = ""

            try:
                self.dir = config["directories"][module]
            except KeyError:
                self.dir = ""

            self.output_dir = os.path.join(self.dir, f"{self.name}.{prefix}")

    def __str__(self):
        return f"{self.module}:{self.name}:{self.args}:{self.db}:{self.bin}"


@dataclass
class Remap_Target:
    accid: str
    acc_simple: str
    taxid: str
    file: str
    run_prefix: str
    description: str = ""
    accid_in_file: list = field(default_factory=lambda: [])

    def __post_init__(self):
        self.name = f"{self.run_prefix}_{self.acc_simple}_{self.taxid}_{os.path.splitext(os.path.basename(self.file))[0]}"


@dataclass(frozen=True)
class Run_detail_report:
    max_depth: int
    max_depthR: int
    max_gaps: int
    max_prop: int
    max_mapped: int
    input: str
    processed: str
    processed_percent: str
    sift_preproc: bool
    sift_remap: bool
    sift_removed_pprc: str
    processing_final: str
    processing_final_percent: str
    merged: bool
    merged_number: int
    merged_files: str


@dataclass(frozen=True)
class Contig_classification_results:
    performed: bool
    method: str
    args: str
    db: str
    classification_number: int
    classification_minhit: int
    success: bool


@dataclass(frozen=True)
class Read_classification_results:
    performed: bool
    method: str
    args: str
    db: str
    classification_number: int
    classification_minhit: int
    success: bool


@dataclass(frozen=True)
class Remap_main:
    performed: bool
    success: bool
    method: str
    found_total: int
    coverage_min: int
    coverage_max: int


@dataclass(frozen=True)
class Assembly_results:
    performed: bool
    assembly_soft: str
    assembly_args: str
    assembly_number: int
    assembly_min: int
    assembly_mean: int
    assembly_max: int
    assembly_trim: int
