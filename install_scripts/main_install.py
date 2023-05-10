#!/usr/bin/python3
import datetime
import logging
import os
import sys
from distutils.command.install import install
from pathlib import Path
from re import I
from typing import Iterable, List

import dnaio
import xopen
from fastq_filter import fastq_records_to_file, file_to_fastq_records


def fastq_records_to_files(records: Iterable[dnaio.Sequence], filepath_template: str = "test",
                           compression_level: int = 2, max_file_size: int = 1000000000, margin: int = 1000000) -> List[str]:
    """
    split fastq records into files of max_file_size"""

    files = []
    file_index = 0
    filepath = filepath_template + f"_{file_index}" + ".fastq.gz"

    output_h = xopen.xopen(filepath, mode='wb', threads=0,
                           compresslevel=compression_level)

    print(f"writing to {filepath}")

    for record in records:

        if os.path.getsize(filepath) > (max_file_size - margin):
            output_h.close()
            files.append(filepath)
            file_index += 1
            filepath = filepath_template + f"_{file_index}" + ".fastq.gz"
            output_h = xopen.xopen(filepath, mode='wb', threads=0,
                                   compresslevel=compression_level)
            print(f"writing to {filepath}")

        header = ">" + record.name + "\n"
        header = header.encode("ascii")
        output_h.write(header)
        sequence = record.sequence + "\n"
        output_h.write(sequence.encode("ascii"))

    output_h.close()
    files.append(filepath)
    return files


def bioinf_splitext(filepath: str):

    basename, ext = os.path.splitext(filepath)
    if ext == ".gz":
        basename, ext = os.path.splitext(basename)

    return basename, ext


def predict_files_split(filepath: str, max_file_size: int = 1000000000, file_template: str = "test"):
    """
    predict files that will be created by split_file"""
    file_size = os.path.getsize(filepath)
    n_files = file_size/max_file_size
    # round up
    n_files = int(n_files) + 1
    files_predict = [
        file_template + f"_{i}" + ".fastq.gz" for i in range(n_files)
    ]
    return files_predict


def check_proceed(filepath: str, file_template: str, max_file_size: int = 1000000000):
    """
    check if needed to continue"""
    files_predict = predict_files_split(
        filepath, max_file_size=max_file_size, file_template=file_template)

    files_exist = [(os.path.exists(file) and os.path.getsize(file) > 100)
                   for file in files_predict]
    if all(files_exist):
        return False
    else:
        return True


def split_file(filepath: str, max_file_size=1000000000) -> List[str]:
    """
    split file into smaller files"""

    if os.path.getsize(filepath) < max_file_size:
        return [filepath]

    else:
        records = file_to_fastq_records(filepath)
        basename, _ = bioinf_splitext(filepath)
        if check_proceed(filepath, basename, max_file_size=max_file_size):
            print(f"splitting {filepath}")
            return fastq_records_to_files(records, filepath_template=basename, max_file_size=max_file_size)

        else:
            return predict_files_split(filepath, max_file_size=max_file_size, file_template=basename)


def process_file_list(file_list: List[str], max_file_size: int = 1000000000) -> List[str]:
    """
    split files in file_list into smaller files"""

    output_files = []
    for filepath in file_list:
        output_files += split_file(filepath, max_file_size=max_file_size)

    return output_files


def process_nuc_fasta_dict(nuc_fasta_dict: dict, max_file_size: int = 1000000000) -> dict:
    """
    split files in file_list into smaller files"""

    output_files = {}
    for software, file_list in nuc_fasta_dict.items():
        output_files[software] = process_file_list(
            file_list, max_file_size=max_file_size)

    return output_files


def get_args_install():
    """
    get user defined arguments.
    """
    try:
        import argparse

        parser = argparse.ArgumentParser(description="parse arguments")
        parser.add_argument(
            "--envs", action="store_true", default=False, help="Install environments"
        ),
        parser.add_argument(
            "--seqdl",
            action="store_true",
            default=False,
            help="download sequence databases",
        )
        parser.add_argument(
            "--soft",
            action="store_true",
            default=False,
            help="Install software databases",
        )

        parser.add_argument(
            "--taxdump",
            default="metadata/taxdump.tar.gz",
            help=(
                "path to ncbi taxdump. Extract manually and provide due to \
                corruption when using wget or curl. \
                    Find at https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz. \
                        Suggest browser download."
            ),
        )

        parser.add_argument(
            "--nanopore",
            action="store_true",
            default=False,
            help="Install software specific to 3d generation sequencing technologies.",
        )

        parser.add_argument(
            "--test",
            action="store_true",
            default=False,
            help="test software installation",
        )

        parser.add_argument(
            "--organism",
            default="viral",
            choices=["viral", "bacterial", "archaeal", "fungal"],
            help=("organism to install. options: viral, bacteria, archaea, fungi"),
        )

        parser.add_argument(
            "--home",
            default="",
            help="path to install reference databases, sequences.",
        )

        args = parser.parse_args()

    except TypeError as e:
        print("check report args")
        print(e)
        sys.exit(1)

    return args


class software_item:
    def __init__(self, name, path, database, installed, env_path) -> None:
        self.name = name
        self.path = path
        self.database = database
        self.installed = installed
        self.env_path = env_path

    def __repr__(self) -> str:
        return f"({self.name}, {self.path}, {self.database}, {self.installed}, {self.env_path})"


class database_item:
    def __init__(self, name, path, installed) -> None:
        self.name = name
        self.path = path
        self.installed = installed

    def __repr__(self) -> str:
        return f"({self.name}, {self.path}, {self.installed})"


class main_setup:
    """
    prepare metagenomics run environments, databases.
    """

    # internal log. read this to know what is available.
    available_classifiers = [
        "kraken2",
        "centrifuge",
        # "metaphlan2",
        # "humann2",
        "diamond",
        # "vsearch",
        "blast",
        "blastp",
        "kaiju",
        "kuniq",
        "desamba",
        "clark",
        "fastviromeexplorer",
    ]

    available_assemblers = [
        "spades",
        "raven",
        "velvet",
        "flye",
    ]

    available_preprocessors = ["trimmomatic", "nanofilt"]

    available_remapers = [
        "bwa",
        "minimap2",
        "rematch",
        "snippy",
        "bowtie2",
    ]

    installed_software = []
    installed_databases = []

    def __init__(
        self,
        env_install,
        setup_dl,
        setup_install,
        repository,
        pdir="",
        ENVS_PARAMS="",
        INSTALL_PARAMS="",
        install_config="full",
        install_type="local",
    ) -> None:
        if not ENVS_PARAMS:
            from install_source import ENVS_PARAMS
        if not INSTALL_PARAMS:
            from install_source import INSTALL_PARAMS

        self.ENVS_PARAMS = ENVS_PARAMS
        self.INSTALL_PARAMS = INSTALL_PARAMS
        self.wdir = setup_dl(INSTALL_PARAMS)
        self.env_install_class = env_install(ENVS_PARAMS)
        self.setup_install_class = setup_install
        self.install_config = install_config

        if not pdir:
            pdir = os.getcwd()
        if pdir[-1] != "/":
            pdir += "/"
        self.pdir = pdir

        self.setup_config()
        self.utilities = repository(
            db_path=self.wdir.home, install_type=install_type)

    @ staticmethod
    def software_install_string(software_name: str):
        date = datetime.datetime.now().strftime("%Y-%m-%d")
        return f"{software_name} installed on {date}"

    @ staticmethod
    def database_install_string(database_name: str):
        date = datetime.datetime.now().strftime("%Y-%m-%d")
        return f"{database_name} installed on {date}"

    def setup_config(self):
        if self.install_config == "full":

            try:
                from install_scripts.config import Televir_Layout_full
            except ModuleNotFoundError as e:
                print("check install config.py exists")
                sys.exit()

            self.layout = Televir_Layout_full()

        elif self.install_config == "docker":

            try:
                from install_scripts.config import Televir_Layout_docker
            except ModuleNotFoundError as e:
                print("check install config.py exists")
                sys.exit()

            self.layout = Televir_Layout_docker()

        else:

            try:
                from install_scripts.config import Televir_Layout_minimal
            except ModuleNotFoundError as e:
                print("check install config.py exists")
                sys.exit()

            self.layout = Televir_Layout_minimal()

    def user_input(self):
        args = get_args_install()
        self.install_envs = args.envs
        self.seqdl = args.seqdl
        self.soft = args.soft
        self.nanopore = args.nanopore
        self.taxdump = args.taxdump
        self.setup_install_class.test = args.test
        self.organism = args.organism
        self.wdir.organism = self.organism

        if args.home:
            self.wdir.home = args.home
            self.wdir.dbdir = args.home + "ref_db/"
            self.wdir.seqdir = args.home + "ref_fasta/"
            self.wdir.metadir = args.home + "metadata/"
            self.ENVS_PARAMS["ENVSDIR"] = args.home + "envs/"

    def object_input(
        self,
        envs,
        seqdl,
        soft,
        nanopore=False,
        taxdump="",
        test=False,
        organism="viral",
    ):
        self.install_envs = envs
        self.seqdl = seqdl
        self.soft = soft
        self.nanopore = nanopore
        self.taxdump = taxdump
        self.setup_install_class.test = test
        self.organism = organism
        self.wdir.organism = organism

    def env_prepare_conda(self, ENVS_PARAMS):
        """install environments described in ENVS_PARAMS script.

        :param ENVS_PARAMS: dictionary of environment parameters
        :type ENVS_PARAMS: dict
        :return: None
        """
        envprep = self.env_install_class
        envprep.prep_dir()
        envprep.conda_install()

    def env_prepare(self, ENVS_PARAMS):
        """install environments described in ENVS_PARAMS script.

        :param ENVS_PARAMS: dictionary of environment parameters
        :type ENVS_PARAMS: dict
        :return: None
        """
        envprep = self.env_install_class
        envprep.prep_dir()
        if self.layout.install_flye:
            envprep.flye_install()

        if self.layout.install_clark:
            envprep.clark_install()

        if self.layout.install_fastviromeexplorer:
            envprep.fve_install()

        if self.layout.install_desamba:
            envprep.deSAMBA_install()

        os.system(
            f"cp {self.pdir}bin/krakenuniq-download* {ENVS_PARAMS['ENVSDIR']}hostDepletion/hostdep_env/bin/"
        )
        os.system(
            f"cp {self.pdir}bin/centrifuge-download {ENVS_PARAMS['ENVSDIR']}hostDepletion/hostdep_env/bin/"
        )
        os.system(
            f"cp {self.pdir}bin/rsync_from_ncbi.pl {ENVS_PARAMS['ENVSDIR']}kraken2/kraken_env/libexec/"
        )

    def setup_envs_conda(self):

        if self.install_envs:
            self.env_prepare_conda(self.ENVS_PARAMS)

    def setup_envs(self):

        if self.install_envs:
            self.env_prepare(self.ENVS_PARAMS)

    def setup_dir(self):
        """Download reference databases, store according to
        description in INSTALL_PARAMS
        """
        if self.seqdl or self.soft:

            self.wdir.mkdirs()
            logging.info("Environment directory: %s", self.wdir.envs["ROOT"])
            logging.info("Fasta directory: %s", self.wdir.seqdir)

    def prep_dl(self):
        """
        download prot sequences and get taxids.
        """

        if self.layout.install_request_sequences:
            request_success = self.wdir.install_requests()
            if request_success:
                self.installed_databases.append(
                    self.database_install_string("requests")
                )
                self.utilities.add_database(
                    self.utilities.database_item(
                        "requests",
                        self.wdir.fastas["nuc"]["requests"][0],
                        True,
                    )
                )

        if self.layout.install_refseq_prot:
            success_refprot = self.wdir.refseq_prot_dl()
            if success_refprot:
                self.installed_databases.append(
                    self.database_install_string("refseq_prot")
                )

                self.utilities.add_database(
                    self.utilities.database_item(
                        "refseq_prot",
                        self.wdir.fastas["prot"]["refseq_prot"],
                        True,
                    )
                )

        if self.layout.install_refseq_gen:
            success_refnuc = self.wdir.refseq_gen_dl()
            if success_refnuc:
                self.installed_databases.append(
                    self.database_install_string("refseq_gen")
                )

                self.utilities.add_database(
                    self.utilities.database_item(
                        "refseq_gen",
                        self.wdir.fastas["nuc"]["refseq"][0],
                        True,
                    )
                )

        if self.layout.install_swissprot:
            swissprot_dl = self.wdir.swissprot_dl()
            if swissprot_dl:
                self.installed_databases.append(
                    self.database_install_string("swissprot")
                )

                self.utilities.add_database(
                    self.utilities.database_item(
                        "swissprot",
                        self.wdir.fastas["prot"]["swissprot"],
                        True,
                    )
                )

        if self.layout.install_hg38:

            success_hg38 = self.wdir.download_hg38()
            if success_hg38:
                self.installed_databases.append(
                    self.database_install_string("hg38"))
                self.utilities.add_database(
                    self.utilities.database_item(
                        "hg38",
                        self.wdir.fastas["host"]["hg38"],
                        True,
                    )
                )

        if self.layout.install_grc38:

            success_hg38 = self.wdir.download_grc38()
            if success_hg38:
                self.installed_databases.append(
                    self.database_install_string("grc38"))
                self.utilities.add_database(
                    self.utilities.database_item(
                        "grc38",
                        self.wdir.fastas["host"]["grc38"],
                        True,
                    )
                )

        if self.layout.install_request_sequences:
            request_success = self.wdir.install_requests()
            if request_success:
                self.installed_databases.append(
                    self.database_install_string("requests")
                )
                self.utilities.add_database(
                    self.utilities.database_item(
                        "requests",
                        self.wdir.fastas["nuc"]["requests"][0],
                        True,
                    )
                )

        if self.organism == "viral":

            if self.layout.install_virosaurus:
                success_virosaurus = self.wdir.virosaurus_dl()
                if success_virosaurus:
                    self.installed_databases.append(
                        self.database_install_string("virosaurus")
                    )

                    self.utilities.add_database(
                        self.utilities.database_item(
                            "virosaurus",
                            self.wdir.fastas["nuc"]["virosaurus"][0],
                            True,
                        )
                    )

            if self.layout.install_rvdb:
                success_rvdb = self.wdir.RVDB_dl()
                if success_rvdb:
                    self.installed_databases.append(
                        self.database_install_string("rvdb")
                    )

                    self.utilities.add_database(
                        self.utilities.database_item(
                            "rvdb",
                            self.wdir.fastas["prot"]["rvdb"],
                            True,
                        )
                    )

        # self.wdir.index_fasta_files()

    def dl_metadata_prot(self):
        """
        generate metadata files for prot databases.
        """
        self.wdir.prot_metadata()
        self.wdir.generate_main_protacc_to_taxid()

    def dl_metadata_nuc(self):
        """
        generate metadata files for nuc databases."""
        self.wdir.nuc_metadata()

    def db_generate_intrinsic(
        self, INSTALL_PARAMS, prepdl, nanopore=False, taxdump="", test=False
    ):
        """
        Install databases for each software envisaged.
        Could be changed so that the software in question are also drawn from config file,
        if installation methods are stored in a relational dictionary themselves beforehand.


        :param INSTALL_PARAMS: dictionary of software installation parameters.
        :type INSTALL_PARAMS: dict
        :param prepdl: boolean, whether to prepare databases.
        :type prepdl: bool
        :param nanopore: boolean, whether to install nanopore software.
        :type nanopore: bool
        :param taxdump: path to taxdump file.
        :type taxdump: str
        :param test: boolean, whether to test installation.
        :type test: bool
        :return: None
        """

        # install databases using organism name only.
        sofprep = self.setup_install_class(
            INSTALL_PARAMS, taxdump=taxdump, test=test, organism=self.organism
        )
        logging.info("database directory %s", sofprep.dbdir)

        sofprep.install_prep()

        logging.info("install prepped")

        ######################### centrifuge ###############################

        if self.layout.install_centrifuge:
            install_success = sofprep.centrifuge_install(dbname=self.organism)
            centlib = f"refseq-{self.organism}.dust.fna.gz"
            if os.path.isfile(sofprep.dbs["centrifuge"]["fasta"]):
                os.system(
                    f"mv {sofprep.dbs['centrifuge']['fasta']} {prepdl.seqdir}{centlib}"
                )

            else:
                if not os.path.isfile(f"{prepdl.seqdir}{centlib}"):
                    logging.info("centrifuge database not found.")
            if os.path.isfile(f"{prepdl.seqdir}{centlib}"):
                prepdl.fastas["nuc"][
                    "centrifuge"
                ] = [f"{prepdl.seqdir}{centlib}"]  # add to fastas dict

                if install_success:
                    self.installed_software.append(
                        self.software_install_string("centrifuge")
                    )

                    self.utilities.add_software(
                        self.utilities.software_item(
                            "centrifuge",
                            sofprep.dbs["centrifuge"]["db"],
                            "default",
                            True,
                            sofprep.envs["ROOT"] + sofprep.envs["centrifuge"],
                        )
                    )

        if self.layout.install_centrifuge:
            install_success = sofprep.centrifuge_install(
                dbname="bacteria", threads="16")
            centlib = f"refseq-bacteria.dust.fna.gz"
            if os.path.isfile(sofprep.dbs["centrifuge"]["fasta"]):
                os.system(
                    f"mv {sofprep.dbs['centrifuge']['fasta']} {prepdl.seqdir}{centlib}"
                )

            else:
                if not os.path.isfile(f"{prepdl.seqdir}{centlib}"):
                    logging.info("centrifuge database not found.")
            if os.path.isfile(f"{prepdl.seqdir}{centlib}"):
                prepdl.fastas["nuc"][
                    "centrifuge_bacteria"
                ] = [f"{prepdl.seqdir}{centlib}"]  # add to fastas dict

                if install_success:
                    self.installed_software.append(
                        self.software_install_string("centrifuge")
                    )

                    self.utilities.add_software(
                        self.utilities.software_item(
                            "centrifuge",
                            sofprep.dbs["centrifuge"]["db"],
                            "bacteria",
                            True,
                            sofprep.envs["ROOT"] + sofprep.envs["centrifuge"],
                        )
                    )

        ########################## clark ##################################
        if self.layout.install_clark:
            success_install = sofprep.clark_install(dbname=self.organism)
            if success_install:
                self.installed_software.append(
                    self.software_install_string("clark"))

                self.utilities.add_software(
                    self.utilities.software_item(
                        "clark",
                        sofprep.dbs["clark"]["db"],
                        "default",
                        True,
                        sofprep.envs["ROOT"] + sofprep.envs["clark"],
                    )
                )

        ########################## kraken2 ###############################

        if self.layout.install_kraken2:
            success_install = sofprep.kraken2_two_strategies_install(
                dbname=self.organism)
            krlib = f"kraken2-{self.organism}-library.fna.gz"
            if os.path.isfile(sofprep.dbs["kraken2"]["fasta"]):
                os.system(
                    f"mv {sofprep.dbs['kraken2']['fasta']} {prepdl.seqdir}{krlib}"
                )
            else:
                if not os.path.isfile(f"{prepdl.seqdir}{krlib}"):
                    logging.info("kraken2 database fasta not found.")

            if os.path.isfile(f"{prepdl.seqdir}{krlib}"):
                prepdl.fastas["nuc"]["kraken2"] = [f"{prepdl.seqdir}{krlib}"]

            if success_install:

                self.installed_software.append(
                    self.software_install_string("kraken2"))

                self.utilities.add_software(
                    self.utilities.software_item(
                        "kraken2",
                        sofprep.dbs["kraken2"]["db"],
                        "default",
                        True,
                        sofprep.envs["ROOT"] + sofprep.envs["kraken2"],
                    )
                )

        if self.layout.install_kraken2:
            success_install = sofprep.kraken2_download_install(
                dbname="bacteria")
            krlib = f"kraken2-bacteria-library.fna.gz"
            if os.path.isfile(sofprep.dbs["kraken2"]["fasta"]):
                os.system(
                    f"mv {sofprep.dbs['kraken2']['fasta']} {prepdl.seqdir}{krlib}"
                )
            else:
                if not os.path.isfile(f"{prepdl.seqdir}{krlib}"):
                    logging.info("kraken2 database fasta not found.")

            if os.path.isfile(f"{prepdl.seqdir}{krlib}"):
                prepdl.fastas["nuc"]["kraken2-bacteria"] = [
                    f"{prepdl.seqdir}{krlib}"]

            if success_install:

                self.installed_software.append(
                    self.software_install_string("kraken2"))

                self.utilities.add_software(
                    self.utilities.software_item(
                        "kraken2",
                        sofprep.dbs["kraken2"]["db"],
                        "bacteria",
                        True,
                        sofprep.envs["ROOT"] + sofprep.envs["kraken2"],
                    )
                )

        if self.layout.install_kraken2:
            success_install = sofprep.kraken2_download_install(
                dbname="standard")
            krlib = f"kraken2-standard-library.fna.gz"
            if os.path.isfile(sofprep.dbs["kraken2"]["fasta"]):
                os.system(
                    f"mv {sofprep.dbs['kraken2']['fasta']} {prepdl.seqdir}{krlib}"
                )
            else:
                if not os.path.isfile(f"{prepdl.seqdir}{krlib}"):
                    logging.info("kraken2 database fasta not found.")

            if os.path.isfile(f"{prepdl.seqdir}{krlib}"):
                prepdl.fastas["nuc"]["kraken2-standard"] = [
                    f"{prepdl.seqdir}{krlib}"]

            if success_install:

                self.installed_software.append(
                    self.software_install_string("kraken2"))

                self.utilities.add_software(
                    self.utilities.software_item(
                        "kraken2",
                        sofprep.dbs["kraken2"]["db"],
                        "standard",
                        True,
                        sofprep.envs["ROOT"] + sofprep.envs["kraken2"],
                    )
                )

        if self.layout.install_kraken2:
            success_install = sofprep.kraken2_download_install(
                dbname="ribo16s")
            krlib = f"kraken2-ribo16s-library.fna.gz"
            if os.path.isfile(sofprep.dbs["kraken2"]["fasta"]):
                os.system(
                    f"mv {sofprep.dbs['kraken2']['fasta']} {prepdl.seqdir}{krlib}"
                )
            else:
                if not os.path.isfile(f"{prepdl.seqdir}{krlib}"):
                    logging.info("kraken2 database fasta not found.")

            if os.path.isfile(f"{prepdl.seqdir}{krlib}"):
                prepdl.fastas["nuc"]["kraken2-ribo16s"] = [
                    f"{prepdl.seqdir}{krlib}"]

            if success_install:

                self.installed_software.append(
                    self.software_install_string("kraken2"))

                self.utilities.add_software(
                    self.utilities.software_item(
                        "kraken2",
                        sofprep.dbs["kraken2"]["db"],
                        "ribo16s",
                        True,
                        sofprep.envs["ROOT"] + sofprep.envs["kraken2"],
                    )
                )

        ########################## krakenuniq ###############################
        if self.layout.install_krakenuniq:
            success_install = sofprep.kuniq_install(dbname=self.organism)
            if success_install:
                self.installed_software.append(
                    self.software_install_string("krakenuniq")
                )
                self.utilities.add_software(
                    self.utilities.software_item(
                        "krakenuniq",
                        sofprep.dbs["krakenuniq"]["db"],
                        "default",
                        True,
                        sofprep.envs["ROOT"] + sofprep.envs["krakenuniq"],
                    )
                )

            success_install = sofprep.kuniq_install(dbname="fungi")
            if success_install:
                self.installed_software.append(
                    self.software_install_string("krakenuniq")
                )
                self.utilities.add_software(
                    self.utilities.software_item(
                        "krakenuniq",
                        sofprep.dbs["krakenuniq"]["db"],
                        "fungi",
                        True,
                        sofprep.envs["ROOT"] + sofprep.envs["krakenuniq"],
                    )
                )

        # install viral specific databases
        if self.organism == "viral":

            if self.layout.install_kaiju:
                sofprep.kaiju_viral_install()
                self.installed_software.append(
                    self.software_install_string("kaiju"))

                self.utilities.add_software(
                    self.utilities.software_item(
                        "kaiju",
                        sofprep.dbs["kaiju"]["db"],
                        "default",
                        True,
                        sofprep.envs["ROOT"] + sofprep.envs["kaiju"],
                    )
                )

        # install host dbs.

        return sofprep

    def db_generate_external(
        self, prepdl, sofprep, nanopore=False, taxdump="", test=False
    ):

        # sofprep = self.setup_install_class(
        #    INSTALL_PARAMS, taxdump=taxdump, test=test, organism=self.organism
        # )

        logging.info("install prepped")

        for fname, fpath in prepdl.fastas["host"].items():
            bwa_install = sofprep.bwa_install(dbname=fname, reference=fpath)
            if bwa_install:
                self.installed_software.append(
                    self.software_install_string("bwa"))
                self.utilities.add_software(
                    self.utilities.software_item(
                        "bwa",
                        sofprep.dbs["bwa"]["fasta"],
                        fname,
                        True,
                        sofprep.envs["ROOT"] + sofprep.envs["bwa"],
                    )
                )

            if self.layout.install_bowtie2:
                bowtie2_install = sofprep.bowtie2_index(
                    dbname=fname, reference=fpath)
                if bowtie2_install:
                    self.installed_software.append(
                        self.software_install_string("bowtie2")
                    )
                    self.utilities.add_software(
                        self.utilities.software_item(
                            "bowtie2",
                            sofprep.dbs["bowtie2"]["db"],
                            fname,
                            True,
                            sofprep.envs["ROOT"] + sofprep.envs["bowtie2"],
                        )
                    )

            self.utilities.add_software(
                self.utilities.software_item(
                    "minimap2",
                    fpath,
                    fpath,
                    True,
                    sofprep.envs["ROOT"] + sofprep.envs["bwa"],
                )
            )

        # install prot databases using local files.
        for fname, fdb in prepdl.fastas["prot"].items():

            if self.layout.install_diamond:
                install_success = sofprep.diamond_install(dbname=fname, db=fdb)

                if install_success:
                    self.installed_software.append(
                        self.software_install_string("diamond")
                    )

                    self.utilities.add_software(
                        self.utilities.software_item(
                            "diamond",
                            sofprep.dbs["diamond"]["db"],
                            fname,
                            True,
                            sofprep.envs["ROOT"] + sofprep.envs["diamond"],
                        )
                    )

            if self.layout.install_blast:
                if fname == "refseq":

                    success_install = sofprep.blast_install(
                        reference=fdb,
                        dbname=f"refseq_{self.organism}_prot",
                        nuc=False,
                        taxid_map=sofprep.metadir + "acc2taxid.prot.map",
                        args="-parse_seqids",
                        title=f"refseq {self.organism} prot",
                    )

                    if success_install:
                        self.installed_software.append(
                            self.software_install_string("blastp")
                        )

                        self.utilities.add_software(
                            self.utilities.software_item(
                                "blastp",
                                sofprep.dbs["blast"]["db"],
                                f"refseq_{self.organism}_prot",
                                True,
                                sofprep.envs["ROOT"] + sofprep.envs["blast"],
                            )
                        )

        # install nuc databases using local files.
        for fname, fd_list in prepdl.fastas["nuc"].items():

            for fdb in fd_list:

                bwa_install = sofprep.bwa_install(dbname=fname, reference=fdb)
                if bwa_install:
                    self.installed_software.append(
                        self.software_install_string("bwa"))
                    self.utilities.add_software(
                        self.utilities.software_item(
                            "bwa",
                            sofprep.dbs["bwa"]["fasta"],
                            fname,
                            True,
                            sofprep.envs["ROOT"] + sofprep.envs["bwa"],
                        )
                    )

                self.utilities.add_software(
                    self.utilities.software_item(
                        "minimap2",
                        fdb,
                        fdb,
                        True,
                        sofprep.envs["ROOT"] + sofprep.envs["bwa"],
                    )
                )

                if self.layout.install_fastviromeexplorer:
                    install_success = sofprep.fve_install(
                        reference=fdb,
                        dbname=fname,
                        virus_list=sofprep.metadir + f"{fname}-list.txt",
                        list_create=True,
                    )
                    if install_success:
                        self.installed_software.append(
                            self.software_install_string("fastviromeexplorer")
                        )

                        self.utilities.add_software(
                            self.utilities.software_item(
                                "fastviromeexplorer",
                                sofprep.dbs["fastviromeexplorer"]["db"],
                                fname,
                                True,
                                sofprep.envs["ROOT"] +
                                sofprep.envs["fastviromeexplorer"],
                            )
                        )

                if self.layout.install_blast:

                    if fname == "refseq":
                        install_success = sofprep.blast_install(
                            reference=fdb,
                            dbname=f"refseq_{self.organism}_genome",
                            nuc=True,
                            taxid_map=sofprep.metadir + "acc2taxid.nuc.map",
                            args="-parse_seqids",
                            title=f"refseq {self.organism} genome",
                        )

                        if install_success:
                            self.installed_software.append(
                                self.software_install_string("blastn")
                            )

                            self.utilities.add_software(
                                self.utilities.software_item(
                                    "blastn",
                                    sofprep.dbs["blast"]["db"],
                                    f"refseq_{self.organism}_genome",
                                    True,
                                    sofprep.envs["ROOT"] +
                                    sofprep.envs["blast"],
                                )
                            )

                if nanopore:

                    if self.layout.install_desamba:
                        install_success = sofprep.deSAMBA_install(
                            reference=fdb,
                            dbname=fname,
                        )

                        if install_success:
                            self.installed_software.append(
                                self.software_install_string("desamba")
                            )

                            self.utilities.add_software(
                                self.utilities.software_item(
                                    "desamba",
                                    sofprep.dbs["desamba"]["db"],
                                    fname,
                                    True,
                                    sofprep.envs["ROOT"] +
                                    sofprep.envs["desamba"],
                                )
                            )

    def setup_soft(self):

        if self.seqdl or self.soft:

            self.utilities.reset_tables()

            # repdl = self.setup_dir(self.INSTALL_PARAMS)
            self.prep_dl()
            logging.info("Downloading databases and software")
            self.dl_metadata_prot()

            if self.soft:

                sofprep = self.db_generate_intrinsic(
                    self.INSTALL_PARAMS,
                    self.wdir,
                    nanopore=self.nanopore,
                    taxdump=self.taxdump,
                )

                self.db_generate_external(
                    self.wdir,
                    sofprep,
                    nanopore=self.nanopore,
                    taxdump=self.taxdump,
                )

                self.wdir.fastas["nuc"] = process_nuc_fasta_dict(
                    self.wdir.fastas["nuc"], max_file_size=4000000000
                )

                self.dl_metadata_nuc()

            if self.seqdl:
                self.utilities.dump_database(self.wdir.home)

            if self.soft:
                self.utilities.dump_software(self.wdir.home)

    def register_install_logs(self):

        with open(os.path.join(self.wdir.home, "install_log.txt"), "w") as install_log:
            if self.installed_software:
                install_log.write("Installed software:\n")
                for software in self.installed_software:
                    install_log.write(f"{software}\n")
            if self.installed_databases:
                install_log.write("Installed databases:\n")
                for database in self.installed_databases:
                    install_log.write(f"{database}\n")

    def setup_deploy(self):

        envprep = self.env_install_class
        envprep.prep_dir()

        envprep.install_deployment_software()


if __name__ == "__main__":

    from modules.db_install import setup_dl, setup_install
    from modules.env_install import env_install
    from modules.utility_manager import Utility_Repository

    logging.basicConfig(
        format="%(asctime)s - %(message)s",
        datefmt="%d-%b-%y %H:%M:%S",
        level=logging.INFO,
    )

    metagen_prep = main_setup(env_install, setup_dl,
                              setup_install, Utility_Repository)
    metagen_prep.user_input()

    metagen_prep.setup_envs_conda()
    metagen_prep.setup_envs()
    metagen_prep.setup_dir()
    metagen_prep.setup_soft()
