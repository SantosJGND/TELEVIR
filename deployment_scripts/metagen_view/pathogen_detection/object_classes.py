import gzip
import logging
import os
import shutil
import subprocess
import time
from dataclasses import dataclass, field
from random import randint
from typing import Type

import pandas as pd
from numpy import ERR_CALL

from pathogen_detection.constants_settings import ConstantsSettings
from pathogen_detection.utilities import fastqc_parse


class RunCMD:
    """
    Run command line commands.
    """

    def __init__(self, bin, logdir: str = "", prefix: str = "run", task: str = "NONE"):
        """
        Initialize.
        """
        if bin:
            if bin[-1] != "/":
                bin += "/"

        self.bin = bin
        self.logs = []
        self.task = task

        self.logger = logging.getLogger(f"{prefix}_{task}")
        self.logger.setLevel(logging.ERROR)
        self.logger.addHandler(logging.StreamHandler())
        self.logger.propagate = False

        self.logfile = os.path.join(logdir, f"{prefix}_{task}.log")
        self.logdir = logdir
        self.prefix = prefix

    def temp_script_log(self):
        """
        Create temporary script.
        """

        seed = str(randint(1000000, 9999999))

        bash_script = f"temp_{self.task}_{seed}.sh"
        bash_log = f"temp_{self.task}_{seed}.log"
        bash_flag = f"temp_{self.task}_{seed}.flag"

        return bash_script, bash_log, bash_flag

    def flag_error(self, subprocess_errorlog, cmd: str = ""):
        """
        Check if error in subprocess.
        """

        if subprocess_errorlog:
            try:
                subprocess_errorlog = subprocess_errorlog.decode("utf-8")
                if "Killed" in subprocess_errorlog:
                    print(f"Killed {cmd}")
                if "[error]" in subprocess_errorlog.lower():
                    return True

            except Exception as e:
                return False

        return False

    @staticmethod
    def process_cmd_log(cmd_out):
        """
        Process command log.
        """

        if not cmd_out:
            return ""

        if isinstance(cmd_out, bytes):
            try:  # python 3
                cmd_out = cmd_out.decode()
            except Exception as e:
                print("log error:", e)
                print("out:", cmd_out)
                cmd_out = ""

        if len(cmd_out) > 300:
            cmd_out = cmd_out[:70] + "..." + cmd_out[-70:].strip()

        return cmd_out

    def output_disposal(self, cmd: str, err: str, out: str, exec_time: float, bin: str):
        if self.logdir:
            with open(os.path.join(self.logdir, self.logfile), "a") as f:
                software = cmd.split(" ")[0]
                f.write(f"exec\t{software}\t{exec_time}\n")
                f.write(f"bin\t{bin}\n")
                f.write(f"{cmd}\n")

                out = self.process_cmd_log(out)
                err = self.process_cmd_log(err)

                f.write(f"out\t{out}\n")
                f.write(f"err\t{err}\n")

        else:
            self.logs.append(cmd)
            self.logs.append(out)

    def run(self, cmd):
        """
        Run software.
        """

        if isinstance(cmd, list):
            cmd = " ".join(cmd)

        self.logger.info(f"running: {self.bin}{cmd}")

        start_time = time.perf_counter()
        proc_prep = subprocess.Popen(
            f"{self.bin}{cmd}",
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        out, err = proc_prep.communicate()

        exec_time = time.perf_counter() - start_time

        if self.flag_error(err, cmd):
            self.logger.error(f"errror in command: {self.bin}{cmd}")
            raise Exception(err.decode("utf-8"))

        self.output_disposal(cmd, err, out, exec_time, self.bin)

    def run_python(self, cmd):
        """
        Run software with python.
        """

        if isinstance(cmd, list):
            cmd = " ".join(cmd)

        self.logger.info(f"running command: python {self.bin}{cmd}")

        start_time = time.perf_counter()

        proc_prep = subprocess.Popen(
            f"python {self.bin}{cmd}",
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        out, err = proc_prep.communicate()

        exec_time = time.perf_counter() - start_time

        if self.flag_error(err):
            self.logger.error(f"errror in command: {self.bin}{cmd}")
            raise Exception(err.decode("utf-8"))

        self.output_disposal(cmd, err, out, exec_time, self.bin)

    def run_java(self, cmd):
        """
        Run software with java.
        """

        if isinstance(cmd, list):
            cmd = " ".join(cmd)

        # print(f"running command: java -cp {self.bin} {cmd}")

        start_time = time.perf_counter()

        proc_prep = subprocess.Popen(
            f"java -cp {self.bin} {cmd}",
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        out, err = proc_prep.communicate()

        exec_time = time.perf_counter() - start_time

        if self.flag_error(err):
            self.logger.error(f"errror in command: {self.bin}{cmd}")
            raise Exception(err.decode("utf-8"))

        self.output_disposal(cmd, err, out, exec_time, self.bin)

    def run_bash(self, cmd):
        """
        Run bash command.
        """

        if isinstance(cmd, list):
            cmd = " ".join(cmd)

        start_time = time.perf_counter()
        proc_prep = subprocess.Popen(
            f"{cmd}",
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
        )
        out, err = proc_prep.communicate()

        exec_time = time.perf_counter() - start_time

        if self.flag_error(err):
            self.logger.error(f"errror in command: {self.bin}{cmd}")
            raise Exception(err.decode("utf-8"))

        self.output_disposal(cmd, err, out, exec_time, "")

    def run_script_software(self, cmd):
        """
        Run bash script.
        """

        if isinstance(cmd, list):
            cmd = " ".join(cmd)

        start_time = time.perf_counter()

        bash_script, bash_log, bash_flag = self.temp_script_log()
        bash_script = os.path.join(self.logdir, bash_script)
        bash_log = os.path.join(self.logdir, bash_log)
        bash_flag = os.path.join(self.logdir, bash_flag)

        with open(bash_script, "w") as f:
            f.write("#!/bin/bash")
            f.write("\n")
            f.write(f"{self.bin}{cmd}")
            f.write("\n")
            f.write("touch " + bash_flag)

        proc_prep = subprocess.Popen(
            "bash " + bash_script + " &> " + bash_log,
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
        )

        out, err = proc_prep.communicate()

        found_flag = False

        while not found_flag:
            time.sleep(1)
            found_flag = os.path.exists(bash_flag)

        err = open(bash_log).read()

        exec_time = time.perf_counter() - start_time

        if self.flag_error(err):
            self.logger.error(f"errror in command: {self.bin}{cmd}")

        os.remove(bash_script)
        os.remove(bash_log)
        os.remove(bash_flag)

        self.output_disposal(cmd, err, out, exec_time, "")

    def run_script(self, cmd):
        """
        Run bash script.
        """

        if isinstance(cmd, list):
            cmd = " ".join(cmd)

        start_time = time.perf_counter()

        bash_script, bash_log, bash_flag = self.temp_script_log()
        bash_script = os.path.join(self.logdir, bash_script)
        bash_log = os.path.join(self.logdir, bash_log)
        bash_flag = os.path.join(self.logdir, bash_flag)

        with open(bash_script, "w") as f:
            f.write("#!/bin/bash")
            f.write("\n")
            f.write(cmd)
            f.write("\n")
            f.write("touch " + bash_flag)

        proc_prep = subprocess.Popen(
            "bash " + bash_script + " &> " + bash_log,
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
        )

        out, err = proc_prep.communicate()

        found_flag = False

        while not found_flag:
            time.sleep(1)
            found_flag = os.path.exists(bash_flag)

        err = open(bash_log).read()
        out = open(bash_log).read()

        exec_time = time.perf_counter() - start_time

        if self.flag_error(err):
            self.logger.error(f"errror in command: {self.bin}{cmd}")

        os.remove(bash_script)
        os.remove(bash_log)
        os.remove(bash_flag)

        self.output_disposal(cmd, err, out, exec_time, "")

    def run_script_return(self, cmd):
        """
        Run bash script.
        """

        if isinstance(cmd, list):
            cmd = " ".join(cmd)

        start_time = time.perf_counter()

        bash_script, bash_log, bash_flag = self.temp_script_log()
        bash_script = os.path.join(self.logdir, bash_script)
        bash_log = os.path.join(self.logdir, bash_log)
        bash_flag = os.path.join(self.logdir, bash_flag)

        with open(bash_script, "w") as f:
            f.write("#!/bin/bash")
            f.write("\n")
            f.write(cmd)
            f.write("\n")
            f.write("touch " + bash_flag)

        proc_prep = subprocess.Popen(
            "bash " + bash_script + " &> " + bash_log,
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
        )

        out, err = proc_prep.communicate()

        found_flag = False

        while not found_flag:
            time.sleep(1)
            found_flag = os.path.exists(bash_flag)

        err = open(bash_log).read()
        out = open(bash_log).read()

        exec_time = time.perf_counter() - start_time

        if self.flag_error(err):
            self.logger.error(f"errror in command: {self.bin}{cmd}")

        os.remove(bash_script)
        os.remove(bash_log)
        os.remove(bash_flag)

        self.output_disposal(cmd, err, out, exec_time, "")

        return out.strip()

    def run_bash_return(self, cmd):
        """
        Run bash command and return stdout.
        """

        if isinstance(cmd, list):
            cmd = " ".join(cmd)

        start_time = time.perf_counter()

        proc_prep = subprocess.Popen(
            f"{cmd}", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )
        out, err = proc_prep.communicate()

        if self.flag_error(err):
            self.logger.error(f"errror in command: {cmd}")
            raise Exception(err.decode("utf-8"))

        exec_time = time.perf_counter() - start_time

        if self.flag_error(err):
            self.logger.error(f"errror in command: {self.bin}{cmd}")
            raise Exception(err.decode("utf-8"))

        self.output_disposal(cmd, err, "", exec_time, "")

        return out


class Read_class:
    filepath: str
    current: str
    prefix: str
    clean: str
    enriched: str
    depleted: str
    current_status: str
    read_number_raw: int = 0
    read_number_clean: int = 0
    read_number_enriched: int = 0
    read_number_depleted: int = 0
    read_number_current: int = 0

    def __init__(
        self,
        filepath,
        prefix,
        clean_dir: str,
        enriched_dir: str,
        depleted_dir: str,
        bin: str,
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
        self.filepath = filepath
        self.current = filepath

        self.exists = os.path.isfile(filepath) and os.path.getsize(filepath) > 100

        if self.exists:
            run_filepath = os.path.join(
                os.path.dirname(clean_dir), prefix + ".input.fastq.gz"
            )
            shutil.copy(
                self.current,
                run_filepath,
            )
            self.current = run_filepath
            self.filepath = run_filepath

            self.read_number_raw = self.get_current_fastq_read_number()

            if self.read_number_raw == 0:
                self.exists = False

        self.prefix = self.determine_read_name(filepath)
        self.clean = os.path.join(clean_dir, self.prefix + ".clean.fastq.gz")
        self.enriched = os.path.join(enriched_dir, self.prefix + ".enriched.fastq.gz")
        self.depleted = os.path.join(depleted_dir, self.prefix + ".depleted.fastq.gz")
        self.current_status = "raw"

        self.read_number_clean = 0
        self.read_number_enriched = 0
        self.read_number_depleted = 0
        self.read_number_filtered = 0

    def update(self, prefix, clean_dir: str, enriched_dir: str, depleted_dir: str):

        self.clean = os.path.join(clean_dir, self.prefix + ".clean.fastq.gz")
        self.enriched = os.path.join(enriched_dir, self.prefix + ".enriched.fastq.gz")
        self.depleted = os.path.join(depleted_dir, self.prefix + ".depleted.fastq.gz")

        if self.exists:
            self.read_number_raw = self.get_current_fastq_read_number()
            run_filepath = os.path.join(
                os.path.dirname(clean_dir), self.prefix + ".input.fastq.gz"
            )

            if self.current != run_filepath:
                shutil.copy(
                    self.current,
                    run_filepath,
                )
                self.current = run_filepath
                self.filepath = run_filepath

            if self.current_status == "clean":
                self.clean = self.current
            if self.current_status == "enriched":
                self.enriched = self.current
            if self.current_status == "depleted":
                self.depleted = self.current

        else:
            print(self.current, self.current_status, prefix, "does not exist")
            # raise Exception("File does not exist")

    def copy(self, filepath):
        """
        Copy file.
        """
        shutil.copy(self.filepath, filepath)

    def determine_read_name(self, filepath):

        if not self.exists:
            return "none"

        if "gz" not in filepath:
            self.cmd.run("bgzip -f %s" % filepath)

        filename = os.path.basename(filepath)

        filename = os.path.splitext(filename)[0]
        filename = os.path.splitext(filename)[0]

        return filename

    def get_read_list(self) -> str:
        """
        Get list of reads.
        """

        if not self.exists:
            return []

        read_grep_cmd = [
            "zcat",
            self.current,
            "|",
            "sed",
            "-n",
            "1~4p",
            "|",
            "awk",
            "'{ print $1 }'",
            "|",
            "sed",
            "-e",
            "'s/^@'//g",
        ]

        read_list = self.cmd.run_bash_return(read_grep_cmd).decode().split("\n")
        read_list = [x for x in read_list if x]
        return read_list

    def fake_quality(self, quality: str = "30"):

        if not self.exists:
            return

        tempf = os.path.join(os.path.dirname(self.current), "temp.fq.gz")

        cmd = [
            "seqtk",
            "seq",
            "-F",
            quality,
            self.current,
            "|",
            "gzip",
            "-c",
            ">",
            tempf,
        ]
        self.cmd.run(cmd)

        if os.path.isfile(tempf):
            os.remove(self.current)

            shutil.copy(tempf, self.current)

    def filter_cigar_strings(self):

        if not self.exists:
            return

        tempf = os.path.join(os.path.dirname(self.current), "temp.fq.gz")

        with gzip.open(self.current, "rb") as f_in, gzip.open(tempf, "wb") as f_out:
            counter = 0
            lines = []
            for line in f_in:

                if len(lines) == 4:

                    if lines[-1][0].decode() != "@":
                        for ln in lines:
                            f_out.write(ln)

                    else:
                        print(lines)
                    counter = 0
                    lines = []

                lines.append(line)
                counter += 1

        if os.path.isfile(tempf):
            os.remove(self.current)
            os.rename(tempf, self.current)

    def read_filter_move(self, input: str, read_list: list, output: str = ""):

        if not self.exists:
            return

        temp_reads_keep = os.path.join(
            os.path.dirname(output), f"keep_temp_{randint(1,1999)}.lst"
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

        self.cmd.run_script(cmd)

    def enrich(self, read_list):
        """
        filter reads and set current status to enriched.
        """

        if len(read_list) > 0:
            self.read_filter_move(self.current, read_list, self.enriched)
            self.is_enriched()

    def deplete(self, read_list):
        """
        filter reads and aset current status to depleted.
        """

        if len(read_list) > 0:
            self.read_filter_move(self.current, read_list, self.depleted)
            self.is_depleted()

    def is_clean(self):
        """
        Set current status to clean.
        """
        if not self.exists:
            return
        shutil.copy(self.current, self.clean)
        self.current = self.clean
        self.read_number_clean = self.get_current_fastq_read_number()
        self.current_status = "clean"
        self.filepath = os.path.dirname(self.current)

    def is_enriched(self):
        """
        Set current status to enriched.
        """
        if not self.exists:
            return
        shutil.copy(self.current, self.enriched)
        self.current = self.enriched
        self.read_number_enriched = self.get_current_fastq_read_number()
        self.current_status = "enriched"
        self.filepath = os.path.dirname(self.current)
        self.read_number_filtered = self.read_number_enriched

    def is_depleted(self):
        """
        Set current status to depleted.
        """
        if not self.exists:
            return
        shutil.copy(self.current, self.depleted)
        self.current = self.depleted
        self.read_number_depleted = self.get_current_fastq_read_number()
        self.current_status = "depleted"
        self.filepath = os.path.dirname(self.current)
        self.read_number_filtered = self.read_number_depleted

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

    def move_to_static(self, main_static, sub_static):
        """
        Move current fastq file to static folder.
        """

        new_file = os.path.basename(self.current)
        new_file = os.path.join(
            ConstantsSettings.static_directory, main_static, sub_static, new_file
        )
        if self.exists and not os.path.isfile(self.current):
            print("File does not exist: %s" % self.current)
            print(self.clean, self.enriched, self.depleted)

        if not os.path.exists(new_file) and os.path.isfile(self.current):
            shutil.copy(self.current, new_file)

        self.current = new_file

    def __str__(self):
        return self.filepath


class Sample_runClass:

    r1: Type[Read_class]
    r2: Type[Read_class]
    report: str
    qcdata: dict
    reads_before_processing: int = 0
    reads_after_processing: int = 0
    sample_dir: str = ""
    qc_soft: str = ""
    input_fastqc_report: os.PathLike
    processed_fastqc_report: os.PathLike
    threads: int = 1
    user_name: str

    def __init__(
        self,
        r1: Type[Read_class],
        r2: Type[Read_class],
        sample_name: str,
        project_name: str,
        user_name: str,
        technology: str,
        type: str,
        combinations: str,
        input: str,
        bin: str,
        threads: int = 1,
    ) -> None:
        self.cmd = RunCMD(bin)
        self.r1 = r1
        self.r2 = r2
        self.sample_name = sample_name
        self.project_name = project_name
        self.user_name = user_name
        self.technology = technology
        self.type = type
        self.combinations = combinations
        self.input = input
        self.threads = threads

        self.QCdir = os.path.dirname(self.r1.clean)

        self.reads_before_processing = self.r1.read_number_raw + self.r2.read_number_raw

    def current_total_read_number(self):
        return (
            self.r1.get_current_fastq_read_number()
            + self.r2.get_current_fastq_read_number()
        )

    def QC_summary(self):
        """get QC summary
        runs fastqc parse on zip files in QCdir

        :param QCdir:
        :return: qc_summary
        """
        qc_summary = {
            "input": fastqc_parse(os.path.join(self.QCdir, "input_data.zip")),
            "processed": fastqc_parse(os.path.join(self.QCdir, "processed_data.zip")),
        }

        return qc_summary

    def get_qc_data(self, reads_dir: str = "reads/clean/"):
        """get qc data for sample_class. Update sample_class.qc_data.

        :param sample_class:
        :return: None
        """

        self.qcdata = self.QC_summary()

    def get_fake_qc_data(self):
        """get fake qc data for sample_class. Update sample_class.qc_data.

        :param sample_class:
        :return: None
        """

        fake_df = pd.DataFrame(
            [
                ["Total_sequences", 1000000],
                ["Sequence length", 5000],
                ["%GC", 50],
                ["Q20", 90],
            ],
            columns=["measure", "value"],
        ).set_index("measure")

        self.qcdata = {
            "input": fake_df,
            "processed": fake_df,
        }

    def fake_quality_strings(self):
        """
        Fake quality strings for fastqc.
        """
        self.r1.fake_quality()
        self.r2.fake_quality()

    def filter_cigar_strings(self):
        """
        Filter reads with cigar strings.
        """
        self.r1.filter_cigar_strings()
        self.r2.filter_cigar_strings()

    def clean_unique(self):
        if self.type == "SE":
            self.clean_unique_SE()
        else:
            self.clean_unique_PE()

    def clean_unique_SE(self):
        WHERETO = os.path.dirname(self.r1.current)
        unique_reads = os.path.join(WHERETO, "unique_reads.lst")

        cmd = [
            "zgrep",
            "'^@'",
            self.r1.current,
            "|",
            "awk",
            "'{print $1}'",
            "|",
            "sed",
            "'s/^@//g'",
            "|",
            "sort",
            "|",
            "uniq",
            ">",
            unique_reads,
        ]

        self.cmd.run_script(cmd)
        if not os.path.exists(unique_reads) or os.path.getsize(unique_reads) == 0:
            print(
                f"No unique reads found in {self.r1.current}, skipping unique read cleaning"
            )
            return

        self.r1.read_filter_inplace(self.r1.current, unique_reads)

    def clean_unique_PE(self):

        WHERETO = os.path.dirname(self.r1.current)
        common_reads = os.path.join(WHERETO, f"common_reads_{randint(0,10000)}.lst")

        cmd_find_common = [
            f"seqkit common -n -i {self.r1.current} {self.r1.current} | paste - - - - | cut -f1 | sort | uniq | sed 's/^@//g' > {common_reads}"
        ]

        self.cmd.run_script(cmd_find_common)

        if os.path.getsize(common_reads) == 0:
            return

        self.r1.read_filter_inplace(self.r1.current, common_reads)
        self.r2.read_filter_inplace(self.r2.current, common_reads)

    def trimmomatic_sort(self):
        if self.type == "SE":
            self.trimmomatic_sort_SE()
        else:
            self.trimmomatic_sort_PE()

    def trimmomatic_sort_SE(self):

        tempdir = os.path.dirname(self.r1.current)
        tempfq = os.path.join(tempdir, f"temp{randint(1,1999)}")

        cmd_trimsort = [
            "trimmomatic",
            "SE",
            "-threads",
            f"{self.threads}",
            self.r1.current,
            tempfq,
            "MINLEN:20",
        ]

        self.cmd.run(cmd_trimsort)

        if tempfq in os.listdir(tempdir):
            os.remove(self.r1.current)
            os.rename(tempfq, self.r1.current)

    def trimmomatic_sort_PE(self):
        if self.type == "SE":
            return

        tempdir = os.path.dirname(self.r1.current)
        tempfq = os.path.join(tempdir, f"temp{randint(1,1999)}")

        cmd_trimsort = [
            f"trimmomatic PE -phred33 -threads {self.threads} {self.r1.current} {self.r2.current} -baseout {tempfq}.fastq.gz MINLEN:20"
        ]

        self.cmd.run(cmd_trimsort)

        if tempfq + "_1P.fastq.gz" in os.listdir(tempdir):
            os.remove(self.r1.current)
            os.remove(self.r2.current)
            os.rename(tempfq + "_1P.fastq.gz", self.r1.current)
            os.rename(tempfq + "_2P.fastq.gz", self.r2.current)


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
            self.module = "None"
            self.name = "None"
            self.args = "None"
            self.db = "None"
            self.db_name = "None"
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
                    self.db_name = db_name
                else:
                    self.db = os.path.join(config["source"]["DBDIR_MAIN"], db_name)
                    self.db_name = db_name

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
    reads: bool = False
    contigs: bool = False

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
    success: int
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
