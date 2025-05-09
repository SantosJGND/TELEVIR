#!/usr/bin/python3
import datetime
import logging
import os
import sys

from install_scripts.modules.db_install import setup_dl, setup_install
from install_scripts.modules.env_install import env_install
from install_scripts.modules.utility_manager import Utility_Repository
from install_scripts.utils.parse_utils import (
    process_nuc_fasta_dict,
)
from install_scripts.config import TelevirLayout


class LayoutWithReport(TelevirLayout):
    """
    Layout with report.
    """

    def report_config(self):
        """
        Print the configuration.
        """
        print("Configuration:")

        attributes = [
            attr
            for attr in dir(self)
            if not callable(getattr(self, attr)) and not attr.startswith("__")
        ]

        for attr in attributes:
            print(f"{attr} : {getattr(self, attr)}")

        print("#" * 20)
        print("\n")


def get_args_install():
    """
    get user defined arguments.
    """
    try:
        import argparse

        parser = argparse.ArgumentParser(description="parse arguments")
        parser.add_argument(
            "--envs", action="store_true", default=False, help="Install environments"
        )
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


class main_setup:
    """
    prepare metagenomics run environments, databases.
    """

    installed_software = []
    installed_databases = []

    wdir: setup_dl
    utilities: Utility_Repository
    env_manager: env_install
    setup_install_class: setup_install

    def __init__(
        self,
        env_install_engine: env_install,
        dl_engine: setup_dl,
        install_engine: setup_install,
        repository: Utility_Repository,
        pdir="",
        ENVS_PARAMS="",
        INSTALL_PARAMS="",
    ) -> None:
        self.ENVS_PARAMS = ENVS_PARAMS
        self.INSTALL_PARAMS = INSTALL_PARAMS
        self.wdir = dl_engine
        self.env_manager = env_install_engine
        self.setup_install_class = install_engine

        if not pdir:
            pdir = os.getcwd()
        if pdir[-1] != "/":
            pdir += "/"
        self.pdir = pdir

        self.setup_config()
        self.utilities = repository

    @staticmethod
    def software_install_string(software_name: str):
        date = datetime.datetime.now().strftime("%Y-%m-%d")
        return f"{software_name} installed on {date}"

    @staticmethod
    def database_install_string(database_name: str):
        date = datetime.datetime.now().strftime("%Y-%m-%d")
        return f"{database_name} installed on {date}"

    def setup_config(self):

        try:
            self.layout = LayoutWithReport()
            self.layout.report_config()
        except ModuleNotFoundError as e:
            print("check install config.py exists")
            print(e)
            sys.exit()

    def user_input(self):
        """
        get user defined arguments."""
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
        envprep = self.env_manager
        envprep.prep_dir()
        envprep.conda_install()

    def env_prepare(self, ENVS_PARAMS):
        """install environments described in ENVS_PARAMS script.

        :param ENVS_PARAMS: dictionary of environment parameters
        :type ENVS_PARAMS: dict
        :return: None
        """
        envprep = self.env_manager
        envprep.prep_dir()
        # if self.layout.install_flye:
        #    envprep.flye_install()

        # if self.layout.install_clark:
        #    envprep.clark_install()

        if self.layout.install_fastviromeexplorer:
            envprep.fve_install()

        if self.layout.install_voyager_viral:
            envprep.voyager_install()

        # if self.layout.install_desamba:
        #    envprep.deSAMBA_install()

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

        if self.layout.install_refseq_16s:
            success_16s = self.wdir.refseq_16s_dl()
            if success_16s:
                self.installed_databases.append(
                    self.database_install_string("refseq_16s")
                )

                self.utilities.add_database(
                    self.utilities.database_item(
                        "refseq_16s",
                        self.wdir.fastas["filter"]["refseq_16s"],
                        True,
                    )
                )

        if self.layout.install_ribo16s:
            success_16s = self.wdir.silva_16s_dl(fname="arb-silva_ribo16s")
            if success_16s:
                self.installed_databases.append(
                    self.database_install_string("arb-silva_ribo16s")
                )

                self.utilities.add_database(
                    self.utilities.database_item(
                        "arb-silva_ribo16s",
                        self.wdir.fastas["filter"]["arb-silva_ribo16s"],
                        True,
                    )
                )

        for host_name in self.layout.HOSTS_TO_INSTALL:
            success_install = self.wdir.download_host(host_name)
            if success_install:
                self.installed_databases.append(self.database_install_string(host_name))

                self.utilities.add_database(
                    self.utilities.database_item(
                        host_name,
                        self.wdir.fastas["host"][host_name],
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

    def db_generate_intrinsic(self, prepdl: setup_dl):
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
        sofprep = self.setup_install_class

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
                prepdl.fastas["nuc"]["centrifuge"] = [
                    f"{prepdl.seqdir}{centlib}"
                ]  # add to fastas dict

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

        if self.layout.install_centrifuge_bacteria:
            install_success = sofprep.centrifuge_install(
                dbname="bacteria", threads="16"
            )
            centlib = "refseq-bacteria.dust.fna.gz"
            if os.path.isfile(sofprep.dbs["centrifuge"]["fasta"]):
                os.system(
                    f"mv {sofprep.dbs['centrifuge']['fasta']} {prepdl.seqdir}{centlib}"
                )

            else:
                if not os.path.isfile(f"{prepdl.seqdir}{centlib}"):
                    logging.info("centrifuge database not found.")
            if os.path.isfile(f"{prepdl.seqdir}{centlib}"):
                prepdl.fastas["nuc"]["centrifuge_bacteria"] = [
                    f"{prepdl.seqdir}{centlib}"
                ]  # add to fastas dict

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

        ########################### clark ##################################

        # if self.layout.install_clark:
        #    success_install = sofprep.clark_install(dbname=self.organism)
        #    if success_install:
        #        self.installed_software.append(self.software_install_string("clark"))
        #
        #        self.utilities.add_software(
        #            self.utilities.software_item(
        #                "clark",
        #                sofprep.dbs["clark"]["db"],
        #                "default",
        #                True,
        #                sofprep.envs["ROOT"] + sofprep.envs["clark"],
        #            )
        #        )

        ########################### voyager ###############################

        if self.layout.install_voyager_viral:
            success_install = sofprep.voyager_install_viruses_copy(dbname=self.organism)

            if success_install:
                self.installed_software.append(self.software_install_string("voyager"))

                self.utilities.add_software(
                    self.utilities.software_item(
                        "voyager",
                        sofprep.dbs["voyager"]["db"],
                        "viral",
                        True,
                        sofprep.envs["ROOT"] + sofprep.envs["voyager"],
                    )
                )

        ############################ metaphlan ###############################

        if self.layout.install_metaphlan:
            success_install = sofprep.install_metaphlan()
            if success_install:
                self.installed_software.append(
                    self.software_install_string("metaphlan")
                )

                self.utilities.add_software(
                    self.utilities.software_item(
                        "metaphlan",
                        sofprep.dbs["metaphlan"]["db"],
                        "default",
                        True,
                        sofprep.envs["ROOT"] + sofprep.envs["metaphlan"],
                    )
                )

        ########################## kraken2 ###############################

        if self.layout.install_kraken2:
            success_install = sofprep.kraken2_two_strategies_install(
                dbname=self.organism
            )
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
                self.installed_software.append(self.software_install_string("kraken2"))

                self.utilities.add_software(
                    self.utilities.software_item(
                        "kraken2",
                        sofprep.dbs["kraken2"]["db"],
                        "default",
                        True,
                        sofprep.envs["ROOT"] + sofprep.envs["kraken2"],
                    )
                )

        if self.layout.install_kraken2_eupathdb46:
            success_install = sofprep.kraken2_download_install(dbname="eupathdb46")
            if success_install:
                self.installed_software.append(self.software_install_string("kraken2"))

                self.utilities.add_software(
                    self.utilities.software_item(
                        "kraken2",
                        sofprep.dbs["kraken2"]["db"],
                        "eupathdb46",
                        True,
                        sofprep.envs["ROOT"] + sofprep.envs["kraken2"],
                    )
                )

        if self.layout.install_kraken2_bacteria:
            success_install = sofprep.kraken2_install(dbname="bacteria")
            krlib = f"kraken2-bacteria-library.fna.gz"
            if os.path.isfile(sofprep.dbs["kraken2"]["fasta"]):
                os.system(
                    f"mv {sofprep.dbs['kraken2']['fasta']} {prepdl.seqdir}{krlib}"
                )
            else:
                if not os.path.isfile(f"{prepdl.seqdir}{krlib}"):
                    logging.info("kraken2 database fasta not found.")

            if os.path.isfile(f"{prepdl.seqdir}{krlib}"):
                prepdl.fastas["nuc"]["kraken2-bacteria"] = [f"{prepdl.seqdir}{krlib}"]

            if success_install:
                self.installed_software.append(self.software_install_string("kraken2"))

                self.utilities.add_software(
                    self.utilities.software_item(
                        "kraken2",
                        sofprep.dbs["kraken2"]["db"],
                        "bacteria",
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

        # if self.layout.install_krakenuniq_fungi:
        #    success_install = sofprep.kuniq_install(dbname="fungi")
        #    if success_install:
        #        self.installed_software.append(
        #            self.software_install_string("krakenuniq")
        #        )
        #        self.utilities.add_software(
        #            self.utilities.software_item(
        #                "krakenuniq",
        #                sofprep.dbs["krakenuniq"]["db"],
        #                "fungi",
        #                True,
        #                sofprep.envs["ROOT"] + sofprep.envs["krakenuniq"],
        #            )
        #        )

        # install viral specific databases
        if self.layout.install_kaiju:
            sofprep.kaiju_viral_install()
            self.installed_software.append(self.software_install_string("kaiju"))

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
        self,
        prepdl: setup_dl,
        sofprep: setup_install,
    ):
        logging.info("install prepped")

        for fname, fpath in prepdl.fastas["filter"].items():

            if self.layout.install_bwa_filter:
                bwa_install = sofprep.bwa_install(dbname=fname, reference=fpath)
                self.installed_software.append(
                    self.software_install_string("bwa-filter")
                )
                self.utilities.add_software(
                    self.utilities.software_item(
                        "bwa-filter",
                        sofprep.dbs["bwa"]["fasta"],
                        fname,
                        bwa_install,
                        sofprep.envs["ROOT"] + sofprep.envs["bwa"],
                        tag="filter",
                    )
                )
        for fname, fpath in prepdl.fastas["host"].items():

            if self.layout.install_bwa_host:
                bwa_install = sofprep.bwa_install(dbname=fname, reference=fpath)
                self.installed_software.append(self.software_install_string("bwa"))
                self.utilities.add_software(
                    self.utilities.software_item(
                        "bwa",
                        sofprep.dbs["bwa"]["fasta"],
                        fname,
                        bwa_install,
                        sofprep.envs["ROOT"] + sofprep.envs["bwa"],
                        tag="host",
                    )
                )

            if self.layout.install_bowtie2_depletion:
                bowtie2_install = sofprep.bowtie2_index(dbname=fname, reference=fpath)
                if bowtie2_install:
                    self.installed_software.append(
                        self.software_install_string("bowtie2")
                    )
                    self.utilities.add_software(
                        self.utilities.software_item(
                            "bowtie2",
                            sofprep.dbs["bowtie2"]["db"],
                            fname,
                            bowtie2_install,
                            sofprep.envs["ROOT"] + sofprep.envs["bowtie2"],
                            tag="host",
                        )
                    )

            self.utilities.add_software(
                self.utilities.software_item(
                    "minimap2",
                    fpath,
                    fname,
                    True,
                    sofprep.envs["ROOT"] + sofprep.envs["bwa"],
                    tag="host",
                )
            )

        # install prot databases using local files.
        for fname, fpath in prepdl.fastas["prot"].items():
            if self.layout.install_diamond:
                install_success = sofprep.diamond_install(dbname=fname, db=fpath)

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
                        reference=fpath,
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
            for fpath in fd_list:
                # skip if file size > 15GB
                if os.path.getsize(fpath) > 15000000000:
                    continue

                self.utilities.add_software(
                    self.utilities.software_item(
                        "minimap2",
                        fpath,
                        fpath,
                        True,
                        sofprep.envs["ROOT"] + sofprep.envs["bwa"],
                    )
                )

                # if self.layout.install_bowtie2_depletion:
                #    bowtie2_install = sofprep.bowtie2_index(
                #        dbname=fname, reference=fpath
                #    )
                #    if bowtie2_install:
                #        self.installed_software.append(
                #            self.software_install_string("bowtie2")
                #        )
                #        self.utilities.add_software(
                #            self.utilities.software_item(
                #                "bowtie2",
                #                sofprep.dbs["bowtie2"]["db"],
                #                fname,
                #                True,
                #                sofprep.envs["ROOT"] + sofprep.envs["bowtie2"],
                #            )
                #        )

                if self.layout.install_fastviromeexplorer:
                    install_success = sofprep.fve_install(
                        reference=fpath,
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
                                sofprep.envs["ROOT"]
                                + sofprep.envs["fastviromeexplorer"],
                            )
                        )

                if self.layout.install_blast:
                    if fname in ["refseq", "centrifuge_bacteria"]:
                        install_success = sofprep.blast_install(
                            reference=fpath,
                            dbname=fname,
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
                                    sofprep.envs["ROOT"] + sofprep.envs["blast"],
                                )
                            )

                # if nanopore:
                #    if self.layout.install_desamba:
                #        install_success = sofprep.deSAMBA_install(
                #            reference=fpath,
                #            dbname=fname,
                #        )

    #
    #        if install_success:
    #            self.installed_software.append(
    #                self.software_install_string("desamba")
    #            )
    #
    #            self.utilities.add_software(
    #                self.utilities.software_item(
    #                    "desamba",
    #                    sofprep.dbs["desamba"]["db"],
    #                    fname,
    #                    True,
    #                    sofprep.envs["ROOT"] + sofprep.envs["desamba"],
    #                )
    #            )

    def setup_soft(self):
        if self.seqdl or self.soft:
            self.utilities.reset_tables()

            # repdl = self.setup_dir(self.INSTALL_PARAMS)
            self.prep_dl()
            logging.info("Downloading databases and software")
            self.dl_metadata_prot()

            if self.soft:

                sofprep = self.db_generate_intrinsic(self.wdir)

                self.db_generate_external(
                    self.wdir,
                    sofprep,
                )

                self.wdir.fastas["nuc"] = process_nuc_fasta_dict(
                    self.wdir.fastas["nuc"], max_file_size=4000000000
                )

                if self.layout.check_index_files:
                    self.wdir.index_nuc_fasta_files()

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
        envprep = self.env_manager
        envprep.prep_dir()

        envprep.install_deployment_software()


if __name__ == "__main__":
    from install_scripts.install_source import ENVS_PARAMS
    from install_scripts.install_source import INSTALL_PARAMS

    logging.basicConfig(
        format="%(asctime)s - %(message)s",
        datefmt="%d-%b-%y %H:%M:%S",
        level=logging.INFO,
    )

    utility_repo = Utility_Repository(db_path="./", install_type="local")
    env_manager = env_install(ENVS_PARAMS)
    dl_manager = setup_dl(INSTALL_PARAMS)
    install_manager = setup_install(
        INSTALL_PARAMS, taxdump="", test=False, organism="viral"
    )

    metagen_prep = main_setup(env_manager, dl_manager, install_manager, utility_repo)
    metagen_prep.user_input()

    metagen_prep.setup_envs_conda()
    metagen_prep.setup_envs()
    metagen_prep.setup_dir()
    metagen_prep.setup_soft()
