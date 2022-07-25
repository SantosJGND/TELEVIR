import logging
import os
import time
from random import randint
from typing import Type

import numpy as np
import pandas as pd

from pathogen_detection.assembly_class import Assembly_class
from pathogen_detection.classification_class import Classifier
from pathogen_detection.metadata_handler import Metadata_handler
from pathogen_detection.object_classes import (
    Assembly_results,
    Contig_classification_results,
    Read_class,
    Read_classification_results,
    Remap_main,
    Remap_Target,
    Run_detail_report,
    RunCMD,
    Sample_runClass,
    Software_detail,
)
from pathogen_detection.preprocess_class import Preprocess
from pathogen_detection.remap_class import Remapping


def get_bindir_from_binaries(binaries, key, value: str = ""):

    if value == "":
        try:
            return os.path.join(binaries["ROOT"], binaries[key]["default"], "bin")
        except KeyError:
            return ""
    else:
        try:
            return os.path.join(binaries["ROOT"], binaries[key][value], "bin")
        except KeyError:
            return ""


class RunDetail_main:
    suprun: str = None

    threads: int
    config: dict
    prefix: str
    ## input

    type: str
    r1_suffix: str
    r2_suffix: str

    r1: Type[Read_class]
    r2: Type[Read_class]

    sample: Type[Sample_runClass]
    ##  metadata
    metadata_tool: Type[Metadata_handler]
    sift_query: str
    max_remap: int
    taxid_limit: int

    ## actions
    quality_control: bool
    sift: bool
    depletion: bool
    enrichment: bool
    assembly: bool
    classification: bool
    remapping: bool
    house_cleaning: bool

    ## methods
    preprocess_method: Software_detail
    depletion_method: Software_detail
    enrichment_method: Software_detail
    assembly_method: Software_detail
    assembly_classification_method: Software_detail
    read_classification_method: Software_detail
    remapping_method: Software_detail

    ## directories.
    root: str

    input_reads_dir: str
    filtered_reads_dir: str
    depleted_reads_dir: str

    log_dir: str
    ## output content
    report: pd.DataFrame

    def __init__(self, config: dict, method_args: pd.DataFrame):

        self.logger_level_main = logging.INFO
        self.logger_level_detail = logging.CRITICAL
        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(self.logger_level_main)
        self.logger.addHandler(logging.StreamHandler())
        self.runtime = 0
        self.start_time = time.perf_counter()
        self.exec_time = 0

        self.method_args = method_args
        self.config = config
        self.cmd = RunCMD(get_bindir_from_binaries(config["bin"], "PREPROCESS"))
        self.threads = config["threads"]
        self.prefix = config["prefix"]
        self.type = config["type"]
        self.suprun = self.prefix

        self.r1 = Read_class(
            config["r1"],
            config["directories"]["PREPROCESS"],
            config["directories"]["reads_enriched_dir"],
            config["directories"]["reads_depleted_dir"],
            bin=get_bindir_from_binaries(config["bin"], "PREPROCESS"),
        )
        self.r2 = Read_class(
            config["r2"],
            config["directories"]["PREPROCESS"],
            config["directories"]["reads_enriched_dir"],
            config["directories"]["reads_depleted_dir"],
            bin=get_bindir_from_binaries(config["bin"], "PREPROCESS"),
        )

        self.min_scaffold_length = config["assembly_contig_min_length"]
        self.minimum_coverage = config["minimum_coverage_threshold"]
        self.maximum_coverage = 100000
        ### metadata
        self.metadata_tool = Metadata_handler(
            self.config["metadata"], sift_query=config["sift_query"]
        )

        self.max_remap = config["max_output_number"]
        self.taxid_limit = config["taxid_limit"]

        ### actions
        self.subsample = False
        self.quality_control = config["actions"]["QCONTROL"]
        self.sift = config["actions"]["SIFT"]
        self.depletion = config["actions"]["DEPLETE"]
        self.enrichment = config["actions"]["ENRICH"]
        self.assembly = config["actions"]["ASSEMBLY"]
        self.classification = config["actions"]["CLASSIFY"]
        self.remapping = config["actions"]["REMAPPING"]
        self.house_cleaning = config["actions"]["CLEANING"]

        ### methods
        self.preprocess_method = Software_detail(
            "PREPROCESS",
            self.method_args,
            config,
            self.prefix,
        )
        self.assembly_method = Software_detail(
            "ASSEMBLY",
            method_args,
            config,
            self.prefix,
        )
        self.depletion_method = Software_detail(
            "DEPLETION",
            method_args,
            config,
            self.prefix,
        )

        self.enrichment_method = Software_detail(
            "ENRICHMENT",
            method_args,
            config,
            self.prefix,
        )

        self.contig_classification_method = Software_detail(
            "CONTIG_CLASSIFICATION",
            method_args,
            config,
            self.prefix,
        )
        self.read_classification_method = Software_detail(
            "READ_CLASSIFICATION",
            method_args,
            config,
            self.prefix,
        )

        self.remapping_method = Software_detail(
            "REMAPPING",
            method_args,
            config,
            self.prefix,
        )

        ### drones
        print("empty drones")
        self.depletion_drone = Classifier(
            Software_detail("NONE", method_args, config, self.prefix),
            logging_level=self.logger_level_detail,
        )
        self.enrichment_drone = Classifier(
            Software_detail("NONE", method_args, config, self.prefix),
            logging_level=self.logger_level_detail,
        )

        ### directories.
        self.root = config["directories"]["root"]
        self.filtered_reads_dir = config["directories"]["PREPROCESS"]
        self.log_dir = config["directories"]["log_dir"]

        self.report = pd.DataFrame()
        self.rclass_summary = pd.DataFrame()
        self.aclass_summary = pd.DataFrame()
        self.merged_targets = pd.DataFrame()

        ###

    def Update(self, config: dict, method_args: pd.DataFrame):

        self.method_args = method_args
        # with open(config_json) as json_file:
        #    config = json.load(json_file)

        self.config = config
        self.prefix = config["prefix"]
        self.type = config["type"]
        print("prefix:", self.prefix)
        print("type:", self.type)
        self.start_time = time.perf_counter()

        ### actions
        self.subsample = False
        self.quality_control = config["actions"]["QCONTROL"]
        self.sift = config["actions"]["SIFT"]
        self.depletion = config["actions"]["DEPLETE"]
        self.enrichment = config["actions"]["ENRICH"]
        self.assembly = config["actions"]["ASSEMBLY"]
        self.classification = config["actions"]["CLASSIFY"]
        self.remapping = config["actions"]["REMAPPING"]
        self.house_cleaning = config["actions"]["CLEANING"]

        self.preprocess_method = Software_detail(
            "PREPROCESS",
            self.method_args,
            config,
            self.prefix,
        )
        self.assembly_method = Software_detail(
            "ASSEMBLY",
            method_args,
            config,
            self.prefix,
        )

        self.depletion_method = Software_detail(
            "DEPLETION",
            method_args,
            config,
            self.prefix,
        )

        self.enrichment_method = Software_detail(
            "ENRICHMENT",
            method_args,
            config,
            self.prefix,
        )

        self.contig_classification_method = Software_detail(
            "CONTIG_CLASSIFICATION",
            method_args,
            config,
            self.prefix,
        )
        self.read_classification_method = Software_detail(
            "READ_CLASSIFICATION",
            method_args,
            config,
            self.prefix,
        )

        self.remapping_method = Software_detail(
            "REMAPPING",
            method_args,
            config,
            self.prefix,
        )

    def Update_exec_time(self):
        """
        Update the execution time of the pipeline.
        """
        self.exec_time = self.exec_time + time.perf_counter() - self.start_time


class RunMethods(RunDetail_main):
    def __init__(self, config_json: os.PathLike, method_args: pd.DataFrame):
        super().__init__(config_json, method_args)
        self.mapped_instances = []

    def deploy_QC(self):
        self.preprocess_drone = Preprocess(
            self.r1.current,
            self.r2.current,
            self.filtered_reads_dir,
            self.type,
            self.preprocess_method,
            self.r1.clean,
            self.r2.clean,
            self.threads,
            self.subsample,
            logging_level=self.logger_level_detail,
        )

        print("r1 reads: ", self.r1.get_current_fastq_read_number())
        print("r2 reads: ", self.r2.get_current_fastq_read_number())

        self.preprocess_drone.run()

    def deploy_HD(self):
        self.depletion_drone = Classifier(
            self.depletion_method,
            self.r1.current,
            type=self.type,
            r2=self.r2.current,
            prefix=self.prefix,
            threads=self.threads,
            bin=get_bindir_from_binaries(self.config["bin"], "REMAPPING"),
            logging_level=self.logger_level_detail,
        )

        self.depletion_drone.run()

    def deploy_EN(self):
        self.enrichment_drone = Classifier(
            self.enrichment_method,
            self.r1.current,
            type=self.type,
            r2=self.r2.current,
            prefix=self.prefix,
            threads=self.threads,
            bin=get_bindir_from_binaries(self.config["bin"], "REMAPPING"),
            logging_level=self.logger_level_detail,
        )

        self.enrichment_drone.run()

    def deploy_ASSEMBLY(self):
        self.assembly_drone = Assembly_class(
            self.r1.current,
            self.assembly_method,
            self.type,
            min_scaffold_length=self.min_scaffold_length,
            r2=self.r2.current,
            prefix=self.prefix,
            threads=self.threads,
            bin=get_bindir_from_binaries(self.config["bin"], "REMAPPING"),
            logging_level=self.logger_level_detail,
        )

        self.assembly_drone.run()

    def deploy_CONTIG_CLASSIFICATION(self):

        self.contig_classification_drone = Classifier(
            self.contig_classification_method,
            self.assembly_drone.assembly_file_fasta_gz,
            r2="",
            prefix=self.prefix,
            threads=self.threads,
            bin=get_bindir_from_binaries(self.config["bin"], "REMAPPING"),
            logging_level=self.logger_level_detail,
        )

        self.contig_classification_drone.run()

    def deploy_READ_CLASSIFICATION(self):

        self.read_classification_drone = Classifier(
            self.read_classification_method,
            self.r1.current,
            type=self.type,
            r2=self.r2.current,
            prefix=self.prefix,
            threads=self.threads,
            bin=get_bindir_from_binaries(self.config["bin"], "REMAPPING"),
            logging_level=logging.INFO,  # self.logger_level_detail,  #
        )

        self.read_classification_drone.run()

    def reference_map(self, remap_target: Type[Remap_Target]):
        rdir = os.path.join(
            self.remapping_method.dir,
            remap_target.name,
            "reference",
        )

        target_remap_drone = Remapping(
            self.r1.current,
            remap_target,
            self.remapping_method,
            self.assembly_drone.assembly_file_fasta_gz,
            self.type,
            self.prefix,
            rdir,
            self.threads,
            r2=self.r2.current,
            minimum_coverage=self.minimum_coverage,
            bin=get_bindir_from_binaries(self.config["bin"], "REMAPPING"),
            logging_level=self.logger_level_detail,
            cleanup=self.house_cleaning,
        )

        target_remap_drone.run_remap()

        return target_remap_drone

    def assembly_map(self, reference_remap: Type[Remapping]):

        if len(reference_remap.mapped_contigs) == 0:
            return None

        output_directory = os.path.join(
            self.remapping_method.dir,
            reference_remap.target.name,
            "assembly",
        )

        assembly_target = Remap_Target(
            "none",
            "assembly",
            "none",
            self.assembly_drone.assembly_file_fasta_gz,
            self.prefix,
            "description",
            reference_remap.mapped_contigs,
        )

        assembly_remap_drone = Remapping(
            reference_remap.mapped_subset_r1,
            assembly_target,
            self.remapping_method,
            self.assembly_drone.assembly_file_fasta_gz,
            self.type,
            self.prefix,
            output_directory,
            self.threads,
            r2=reference_remap.mapped_subset_r2,
            minimum_coverage=self.minimum_coverage,
            bin=get_bindir_from_binaries(self.config["bin"], "REMAPPING"),
            logging_level=self.logger_level_detail,
            cleanup=self.house_cleaning,
        )

        assembly_remap_drone.run_remap()

        return assembly_remap_drone

    def deploy_REMAPPING(self):

        for remap_target in self.remap_targets:

            target_remap_drone = self.reference_map(remap_target)

            mapped_instance = {
                "reference": target_remap_drone,
                "assembly": self.assembly_map(target_remap_drone),
            }

            self.mapped_instances.append(mapped_instance)

    def select_targets(self):

        self.metadata_tool.get_metadata()

        print(self.read_classification_drone.classification_report)

        rdata = self.metadata_tool.results_process(
            self.read_classification_drone.classification_report
        )
        cdata = self.metadata_tool.results_process(
            self.contig_classification_drone.classification_report
        )

        merged_targets = self.metadata_tool.merge_reports(
            cdata,
            rdata,
            max_remap=self.max_remap,
        )

        self.aclass_summary = cdata
        self.rclass_summary = rdata
        self.merged_targets = merged_targets

        #######
        #######

        remap_targets, remap_absent = self.metadata_tool.generate_mapping_targets(
            merged_targets,
            prefix=self.prefix,
            taxid_limit=self.taxid_limit,
            fasta_main_dir=self.config["source"]["REF_FASTA"],
        )

        self.remap_targets = remap_targets
        self.remap_absent_taxid_list = remap_absent

    def clean_unique(self):
        if self.type == "SE":
            return

        WHERETO = os.path.dirname(self.r1.current)
        common_reads = os.path.join(WHERETO, "common_reads.lst")

        cmd_find_common = [
            f"seqkit common -n -i {self.r1.current} {self.r1.current} | paste - - - - | cut -f1 | uniq | sed 's/^@//g' > {common_reads}"
        ]

        self.cmd.run(cmd_find_common)
        if os.path.getsize(common_reads) == 0:
            self.logger.info("No common reads found")
            return

        self.r1.read_filter_inplace(self.r1.current, common_reads)
        self.r2.read_filter_inplace(self.r2.current, common_reads)

    def trimmomatic_sort(self):
        if self.type == "SE":
            return

        tempdir = os.path.dirname(self.r1.current)
        tempfq = os.path.join(tempdir, f"temp{randint(1,1999)}")

        cmd_trimsort = [
            f"trimmomatic PE -phred33 -threads {self.threads} {self.r1.current} {self.r2.current} -baseout {tempfq}.fastq.gz MINLEN:20"
        ]

        self.cmd.run(cmd_trimsort)

        if tempfq + ".fastq.gz" in os.listdir(tempdir):
            os.remove(self.r1.current)
            os.remove(self.r2.current)
            os.rename(tempfq + "_1P.fastq.gz", self.r1.current)
            os.rename(tempfq + "_2P.fastq.gz", self.r2.current)

    def Summarize(self):

        self.logger.info(f"prefix: {self.prefix}")
        with open(os.path.join(self.log_dir, self.prefix + "_latest.fofn"), "w") as f:
            f.write(self.r1.current + "\n")
            if self.type == "PE":
                f.write(self.r2.current + "\n")

        with open(os.path.join(self.log_dir, "reads_latest.stats"), "w") as f:
            f.write(f"CLEAN\t{self.r1.read_number_clean}\n")
            f.write(f"ENRICHED\t{self.r1.read_number_enriched}\n")

    def merge_mapping_reports(self):

        full_report = []
        ntax_cols = [
            "suffix",
            "taxid",
            "refseq",
            "description",
            "rclass",
            "aclass",
        ]

        for instance in self.mapped_instances:
            success = "none"
            apres = False
            if not instance["reference"].remapping_successful:
                self.logger.info(
                    f"No mapping output found for reference {instance['reference'].target.name}"
                )
                continue

            ###

            mapped = instance["reference"].number_of_reads_mapped

            if instance["assembly"]:
                apres = True

            if mapped and apres:
                success = "reads and contigs"
            elif mapped:
                success = "reads"
            elif apres:
                success = "contigs"

            ntax = [
                [
                    self.prefix,
                    instance["reference"].target.taxid,
                    instance["reference"].target.accid,
                    instance["reference"].target.description,
                    True,
                    apres,
                ]
            ]

            def simplify_taxid(x):
                return (
                    x.replace(";", "_")
                    .replace(":", "_")
                    .replace(".", "_")
                    .replace("|", "_")
                )

            ntax = pd.DataFrame(ntax, columns=ntax_cols)

            ntax = pd.concat((ntax, instance["reference"].report), axis=1)
            ntax["mapped"] = mapped
            ntax["mapped_prop"] = 100 * (mapped / self.sample.reads_after_processing)
            ntax["ref_prop"] = 100 * (mapped / self.sample.reads_before_processing)
            ntax["refdb"] = instance["reference"].target.file
            ntax["ID"] = instance["reference"].target.accid
            ntax["simple_id"] = ntax["ID"].apply(simplify_taxid)
            ntax["unique_id"] = ntax["ID"].apply(simplify_taxid)
            ntax["contig_length"] = instance["reference"].reference_fasta_length
            ntax["contig_string"] = instance["reference"].reference_fasta_string
            ntax["success"] = success

            ntax["refa_dotplot_exists"] = instance["reference"].dotplot_exists
            ntax["covplot_exists"] = instance["reference"].coverage_plot_exists
            ntax["refa_dotplot_path"] = instance["reference"].dotplot
            ntax["covplot_path"] = instance["reference"].coverage_plot
            ntax["bam_path"] = instance["reference"].read_map_sorted_bam
            ntax["bam_index_path"] = instance["reference"].read_map_sorted_bam_index
            ntax["reference_path"] = instance["reference"].reference_file
            ntax["reference_index_path"] = instance["reference"].reference_fasta_index
            ntax["reference_assembly_paf"] = instance["reference"].assembly_map_paf

            if apres:
                ntax["mapped_scaffolds_path"] = instance[
                    "assembly"
                ].reference_fasta_index
                ntax["mapped_scaffolds_index_path"] = instance[
                    "assembly"
                ].assembly_map_paf
            else:
                ntax["mapped_scaffolds_path"] = ""
                ntax["mapped_scaffolds_index_path"] = ""

            ntax = ntax.sort_values(["taxid", "Hdepth"])

            full_report.append(ntax)

        if len(full_report) > 0:

            self.report = pd.concat(full_report, axis=0)
            self.clean_final_report()
        else:
            self.report = pd.DataFrame()

    def clean_final_report(self):

        self.report.ngaps = self.report.ngaps.fillna(0)

    def get_assembly_stats(self):
        """Assembly contig summary_stats"""

        if self.assembly_drone.contig_summary.shape[0] > 0:
            assembly_min = (
                self.assembly_drone.contig_summary["contig_length"].fillna(0).min()
            )
            assembly_max = (
                self.assembly_drone.contig_summary["contig_length"].fillna(0).max()
            )
            assembly_mean = (
                self.assembly_drone.contig_summary["contig_length"].fillna(0).mean()
            )
        else:
            assembly_min = 0
            assembly_max = 0
            assembly_mean = 0

        return assembly_min, assembly_max, assembly_mean

    def get_classification_input(self):
        """
        Get the classification read number input for classification.
        if sift was performed, then final processing is grabbed from the output ofsift.
        otherwise, final_processed_reads == post_processed_reads.

        """
        post_processed_reads = self.sample.reads_after_processing

        final_processing_reads = (
            self.r1.current_fastq_read_number() + self.r2.current_fastq_read_number()
        )

        return post_processed_reads, final_processing_reads

    def generate_output_data_classes(self):
        ###
        processed_reads = (
            self.sample.r1.read_number_clean + self.sample.r2.read_number_clean
        )

        post_processed_reads, final_processing_reads = self.get_classification_input()

        post_percent = (int(post_processed_reads) / processed_reads) * 100
        final_processing_percent = (final_processing_reads / processed_reads) * 100

        assembly_number = self.assembly_drone.contig_summary.shape[0]
        assembly_min, assembly_max, assembly_mean = self.get_assembly_stats()

        ###
        minhit_assembly = self.aclass_summary["counts"].min()
        if not minhit_assembly or not self.aclass_summary.shape[0]:
            minhit_assembly = 0

        minhit_reads = self.rclass_summary["counts"].min()
        if np.isnan(minhit_reads):
            minhit_reads = 0

        ###
        max_gaps = self.report.ngaps.max()
        if np.isnan(max_gaps):
            max_gaps = 0

        max_prop = self.report.ref_prop.max()
        if np.isnan(max_prop):
            max_prop = 0

        max_mapped = self.report.mapped.max()
        if np.isnan(max_mapped):
            max_mapped = 0

        files = list(set([t["reference"].target.file for t in self.mapped_instances]))

        self.run_detail_report = Run_detail_report(
            self.report.Hdepth.max(),
            self.report.HdepthR.max(),
            max_gaps,
            max_prop,
            max_mapped,
            f"{processed_reads:,}",
            f"{post_processed_reads:,}",
            f"{post_percent:.2f}",
            False,
            self.sift,
            f"{self.metadata_tool.sift_report.loc[0]['removed']:,}",
            f"{final_processing_reads:,}",
            round(final_processing_percent, 3),
            self.remapping,
            self.merged_targets.taxid.nunique(),
            ", ".join(files),
        )

        self.contig_classification_results = Contig_classification_results(
            True,
            self.contig_classification_drone.classifier_method.name,
            self.aclass_summary.shape[0],
            minhit_assembly,
            self.aclass_summary.shape[0] > 0,
        )

        self.read_classification_results = Read_classification_results(
            True,
            self.read_classification_drone.classifier_method.name,
            self.rclass_summary.shape[0],
            minhit_reads,
            self.rclass_summary.shape[0] > 0,
        )

        self.assembly_report = Assembly_results(
            True,
            self.assembly_drone.assembly_method.name,
            assembly_number,
            f"{assembly_min:,}",
            f"{int(assembly_mean):,}",
            f"{assembly_max:,}",
            f"{int(self.min_scaffold_length):,}",
        )

        self.remap_main = Remap_main(
            True,
            self.report.shape[0] > 0,
            self.remapping_method.name,
            len(self.mapped_instances),
            self.minimum_coverage,
            self.maximum_coverage,
        )

    def export_reports(self):

        self.full_report = os.path.join(self.root, f"{self.prefix}_full_report.tsv")
        self.assembly_classification_summary = os.path.join(
            self.root, f"{self.prefix}_aclass_summary.tsv"
        )
        self.read_classification_summary = os.path.join(
            self.root, f"{self.prefix}_rclass_summary.tsv"
        )
        self.merged_classification_summary = os.path.join(
            self.root, f"{self.prefix}_mclass_summary.tsv"
        )

        ### main report
        self.report.to_csv(
            self.full_report,
            index=False,
            sep="\t",
            header=True,
        )

        ### contig classification report
        self.aclass_summary.to_csv(
            self.assembly_classification_summary,
            index=False,
            sep="\t",
            header=True,
        )

        ### read classification report
        self.rclass_summary.to_csv(
            self.read_classification_summary,
            index=False,
            sep="\t",
            header=True,
        )

        ### merged classification report

        self.merged_targets.to_csv(
            self.merged_classification_summary,
            index=False,
            sep="\t",
            header=True,
        )


class RunMain_class(RunMethods):
    def __init__(self, config_json: os.PathLike, method_args: pd.DataFrame):
        super().__init__(config_json, method_args)

    def Run(self):
        print("quality control: ", self.quality_control)
        print("enrichment: ", self.enrichment)
        print("depletion: ", self.depletion)
        print("remapping: ", self.remapping)
        print("assembly: ", self.assembly)
        print("classification: ", self.classification)
        print("sift: ", self.sift)

        if self.quality_control:
            self.deploy_QC()

            self.r1.is_clean()
            self.r2.is_clean()

        if self.enrichment:
            self.deploy_EN()

            self.r1.enrich(self.enrichment_drone.classified_reads_list)
            self.r2.enrich(self.enrichment_drone.classified_reads_list)

        if self.depletion:
            self.deploy_HD()

            self.r1.deplete(self.depletion_drone.classified_reads_list)
            self.r2.deplete(self.depletion_drone.classified_reads_list)

        if self.enrichment or self.depletion:
            self.clean_unique()
            self.trimmomatic_sort()

        if self.assembly:
            self.deploy_ASSEMBLY()

        if self.classification:
            self.deploy_READ_CLASSIFICATION()
            self.deploy_CONTIG_CLASSIFICATION()

        if self.remapping:
            self.select_targets()
            self.deploy_REMAPPING()

        self.Update_exec_time()
