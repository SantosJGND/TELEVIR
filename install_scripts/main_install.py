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
from install_scripts.install_config import TelevirLayout
from install_scripts.load_sources import get_db_entry


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
            request_path = self.wdir.fastas.get("nuc", {}).get("requests", [""])[0]
            req_category, req_name = self.layout.DATABASE_NAMES.get("install_request_sequences", ("taxonomy", "requests"))
            req_entry = get_db_entry(req_category, req_name)
            req_desc = req_entry.get("description") if req_entry else None
            self.utilities.add_database(
                self.utilities.database_item(
                    f"{req_category}/{req_name}",
                    request_path,
                    request_success,
                    description=req_desc,
                )
            )

        if self.layout.install_refseq_prot:
            success_refprot = self.wdir.refseq_prot_dl()
            if success_refprot:
                self.installed_databases.append(
                    self.database_install_string("refseq_prot")
                )

            db_ver = self.wdir.db_versions.get("refseq_prot", {})
            refseq_prot_path = self.wdir.fastas.get("prot", {}).get("refseq_prot", "")
            db_category, db_name = self.layout.DATABASE_NAMES.get("install_refseq_prot", ("protein", "refseq_prot"))
            db_entry = get_db_entry(db_category, db_name)
            db_desc = db_entry.get("description") if db_entry else None
            self.utilities.add_database(
                self.utilities.database_item(
                    f"{db_category}/{db_name}",
                    refseq_prot_path,
                    success_refprot,
                    version=db_ver.get("version"),
                    source_url=db_ver.get("source_url"),
                    file_mod_date=db_ver.get("file_mod_date"),
                    description=db_desc,
                )
            )

        if self.layout.install_refseq_gen:
            success_refnuc = self.wdir.refseq_gen_dl()
            if success_refnuc:
                self.installed_databases.append(
                    self.database_install_string("refseq_gen")
                )

            db_ver = self.wdir.db_versions.get("refseq", {})
            refseq_gen_path = self.wdir.fastas.get("nuc", {}).get("refseq", [""])[0]
            db_category, db_name = self.layout.DATABASE_NAMES.get("install_refseq_gen", ("refseq", "refseq_gen"))
            db_entry = get_db_entry(db_category, db_name)
            db_desc = db_entry.get("description") if db_entry else None
            self.utilities.add_database(
                self.utilities.database_item(
                    f"{db_category}/{db_name}",
                    refseq_gen_path,
                    success_refnuc,
                    version=db_ver.get("version"),
                    source_url=db_ver.get("source_url"),
                    file_mod_date=db_ver.get("file_mod_date"),
                    description=db_desc,
                )
            )

        if self.layout.install_swissprot:
            swissprot_dl = self.wdir.swissprot_dl()
            if swissprot_dl:
                self.installed_databases.append(
                    self.database_install_string("swissprot")
                )

            swissprot_path = self.wdir.fastas.get("prot", {}).get("swissprot", "")
            db_ver = self.wdir.db_versions.get("swissprot", {})
            db_category, db_name = self.layout.DATABASE_NAMES.get("install_swissprot", ("protein", "swissprot"))
            db_entry = get_db_entry(db_category, db_name)
            db_desc = db_entry.get("description") if db_entry else None
            self.utilities.add_database(
                self.utilities.database_item(
                    f"{db_category}/{db_name}",
                    swissprot_path,
                    swissprot_dl,
                    version=db_ver.get("version"),
                    source_url=db_ver.get("source_url"),
                    file_mod_date=db_ver.get("file_mod_date"),
                    description=db_desc,
                )
            )

        if self.layout.install_refseq_16s:
            success_16s = self.wdir.refseq_16s_dl()
            if success_16s:
                self.installed_databases.append(
                    self.database_install_string("refseq_16s")
                )

            refseq_16s_path = self.wdir.fastas.get("filter", {}).get("refseq_16s", "")
            db_ver = self.wdir.db_versions.get("refseq_16s", {})
            db_category, db_name = self.layout.DATABASE_NAMES.get("install_refseq_16s", ("ribosomal_rna", "refseq_16s"))
            db_entry = get_db_entry(db_category, db_name)
            db_desc = db_entry.get("description") if db_entry else None
            self.utilities.add_database(
                self.utilities.database_item(
                    f"{db_category}/{db_name}",
                    refseq_16s_path,
                    success_16s,
                    version=db_ver.get("version"),
                    source_url=db_ver.get("source_url"),
                    file_mod_date=db_ver.get("file_mod_date"),
                    description=db_desc,
                )
            )

        if self.layout.install_ribo16s:
            success_16s = self.wdir.silva_16s_dl(fname="silva_16s")
            if success_16s:
                self.installed_databases.append(
                    self.database_install_string("silva_16s")
                )

            ribo16s_path = self.wdir.fastas.get("filter", {}).get("silva_16s", "")
            db_ver = self.wdir.db_versions.get("silva_16s", {})
            db_category, db_name = self.layout.DATABASE_NAMES.get("install_ribo16s", ("ribosomal_rna", "silva_16s"))
            db_entry = get_db_entry(db_category, db_name)
            db_desc = db_entry.get("description") if db_entry else None
            self.utilities.add_database(
                self.utilities.database_item(
                    f"{db_category}/{db_name}",
                    ribo16s_path,
                    success_16s,
                    version=db_ver.get("version"),
                    source_url=db_ver.get("source_url"),
                    file_mod_date=db_ver.get("file_mod_date"),
                    description=db_desc,
                )
            )

        for host_name in self.layout.HOSTS_TO_INSTALL:
            success_install = self.wdir.download_host(host_name)
            if success_install:
                self.installed_databases.append(self.database_install_string( host_name))

            db_ver = self.wdir.db_versions.get(host_name, {})
            host_path = self.wdir.fastas.get("host", {}).get(host_name, "")
            self.utilities.add_database(
                self.utilities.database_item(
                    name=host_name,
                    path=host_path,
                    installed=success_install,
                    software='host',
                    version=db_ver.get("version"),
                    source_url=db_ver.get("source_url"),
                    file_mod_date=db_ver.get("file_mod_date"),
                )
            )

        if self.layout.install_virosaurus:
            success_virosaurus = self.wdir.virosaurus_dl()
            if success_virosaurus:
                self.installed_databases.append(
                    self.database_install_string("virosaurus")
                )

            virosaurus_path = self.wdir.fastas.get("nuc", {}).get("virosaurus", [""])[0]
            db_ver = self.wdir.db_versions.get("virosaurus", {})
            db_category, db_name = self.layout.DATABASE_NAMES.get("install_virosaurus", ("protein", "virosaurus"))
            db_entry = get_db_entry(db_category, db_name)
            db_desc = db_entry.get("description") if db_entry else None
            self.utilities.add_database(
                self.utilities.database_item(
                    f"{db_category}/{db_name}",
                    virosaurus_path,
                    success_virosaurus,
                    version=db_ver.get("version"),
                    source_url=db_ver.get("source_url"),
                    file_mod_date=db_ver.get("file_mod_date"),
                    description=db_desc,
                )
            )

        if self.layout.install_rvdb:
            success_rvdb = self.wdir.RVDB_dl()
            if success_rvdb:
                self.installed_databases.append(
                    self.database_install_string("rvdb")
                )

            rvdb_path = self.wdir.fastas.get("prot", {}).get("rvdb", "")
            db_ver = self.wdir.db_versions.get("rvdb", {})
            db_category, db_name = self.layout.DATABASE_NAMES.get("install_rvdb", ("protein", "rvdb"))
            db_entry = get_db_entry(db_category, db_name)
            db_desc = db_entry.get("description") if db_entry else None
            self.utilities.add_database(
                self.utilities.database_item(
                    f"{db_category}/{db_name}",
                    rvdb_path,
                    success_rvdb,
                    version=db_ver.get("version"),
                    source_url=db_ver.get("source_url"),
                    file_mod_date=db_ver.get("file_mod_date"),
                    description=db_desc,
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

            sw_name, sw_tag = self.layout.SOFTWARE_NAMES.get("install_centrifuge", ("centrifuge", "viral"))
            self.utilities.add_software(
                self.utilities.software_item(
                    sw_name,
                    sofprep.dbs["centrifuge"]["db"],
                    sw_tag,
                    install_success,
                    sofprep.envs["ROOT"] + sofprep.envs["centrifuge"],
                )
            )

            # Also register as database
            db_cat, db_name = self.layout.DATABASE_NAMES.get("install_centrifuge", ("centrifuge", "viral"))
            db_entry = get_db_entry(db_cat, db_name)
            db_desc = db_entry.get("description") if db_entry else None
            self.utilities.add_database(
                self.utilities.database_item(
                    f"{db_cat}/{db_name}",
                    sofprep.dbs["centrifuge"]["db"],
                    install_success,
                    software=db_cat,
                    description=db_desc,
                )
            )

        if self.layout.install_centrifuge_bacteria:
            install_success = sofprep.centrifuge_download_install(
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

            sw_name, sw_tag = self.layout.SOFTWARE_NAMES.get("install_centrifuge_bacteria", ("centrifuge", "bacteria"))
            self.utilities.add_software(
                self.utilities.software_item(
                    sw_name,
                    sofprep.dbs["centrifuge"]["db"],
                    sw_tag,
                    install_success,
                    sofprep.envs["ROOT"] + sofprep.envs["centrifuge"],
                )
            )

            # Also register as database
            db_cat, db_name = self.layout.DATABASE_NAMES.get("install_centrifuge_bacteria", ("centrifuge", "bacteria"))
            db_entry = get_db_entry(db_cat, db_name)
            db_desc = db_entry.get("description") if db_entry else None
            self.utilities.add_database(
                self.utilities.database_item(
                    f"{db_cat}/{db_name}",
                    sofprep.dbs["centrifuge"]["db"],
                    install_success,
                    software=db_cat,
                    description=db_desc,
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

            sw_name, sw_tag = self.layout.SOFTWARE_NAMES.get("install_voyager_viral", ("voyager", "viral"))
            self.utilities.add_software(
                self.utilities.software_item(
                    sw_name,
                    sofprep.dbs["voyager"]["db"],
                    sw_tag,
                    success_install,
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

            sw_name, sw_tag = self.layout.SOFTWARE_NAMES.get("install_metaphlan", ("metaphlan", "default"))
            self.utilities.add_software(
                self.utilities.software_item(
                    sw_name,
                    sofprep.dbs["metaphlan"]["db"],
                    sw_tag,
                    success_install,
                    sofprep.envs["ROOT"] + sofprep.envs["metaphlan"],
                )
            )

        ########################## kraken2 ###############################

        if self.layout.install_kraken2:
            success_install = sofprep.kraken2_download_install(dbname=self.organism)
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

            kraken_ver = sofprep.dbs.get("kraken2", {}).get("version", "")
            sw_name, sw_tag = self.layout.SOFTWARE_NAMES.get("install_kraken2", ("kraken2", "viral"))
            self.utilities.add_software(
                self.utilities.software_item(
                    sw_name,
                    sofprep.dbs["kraken2"]["db"],
                    sw_tag,
                    success_install,
                    sofprep.envs["ROOT"] + sofprep.envs["kraken2"],
                    db_version=kraken_ver,
                )
            )

            # Also register as database
            db_cat, db_name = self.layout.DATABASE_NAMES.get("install_kraken2", ("kraken2", "viral"))
            db_entry = get_db_entry(db_cat, db_name)
            db_desc = db_entry.get("description") if db_entry else None
            self.utilities.add_database(
                self.utilities.database_item(
                    f"{db_cat}/{db_name}",
                    sofprep.dbs["kraken2"]["db"],
                    success_install,
                    software=db_cat,
                    version=kraken_ver,
                    description=db_desc,
                )
            )

        if self.layout.install_kraken2_eupathdb46:
            success_install = sofprep.kraken2_download_install(dbname="eupathdb46")
            if success_install:
                self.installed_software.append(self.software_install_string("kraken2"))

            sw_name, sw_tag = self.layout.SOFTWARE_NAMES.get("install_kraken2_eupathdb46", ("kraken2", "eupathdb46"))
            self.utilities.add_software(
                self.utilities.software_item(
                    sw_name,
                    sofprep.dbs["kraken2"]["db"],
                    sw_tag,
                    success_install,
                    sofprep.envs["ROOT"] + sofprep.envs["kraken2"],
                )
            )

            # Also register as database
            db_cat, db_name = self.layout.DATABASE_NAMES.get("install_kraken2_eupathdb46", ("kraken2", "eupathdb46"))
            db_entry = get_db_entry(db_cat, db_name)
            db_desc = db_entry.get("description") if db_entry else None
            self.utilities.add_database(
                self.utilities.database_item(
                    f"{db_cat}/{db_name}",
                    sofprep.dbs["kraken2"]["db"],
                    success_install,
                    software=db_cat,
                    description=db_desc,
                )
            )

        if self.layout.install_kraken2_bacteria:
            success_install = sofprep.kraken2_download_install(dbname="bacteria")
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

            sw_name, sw_tag = self.layout.SOFTWARE_NAMES.get("install_kraken2_bacteria", ("kraken2", "bacteria"))
            self.utilities.add_software(
                self.utilities.software_item(
                    sw_name,
                    sofprep.dbs["kraken2"]["db"],
                    sw_tag,
                    success_install,
                    sofprep.envs["ROOT"] + sofprep.envs["kraken2"],
                )
            )

            # Also register as database
            db_cat, db_name = self.layout.DATABASE_NAMES.get("install_kraken2_bacteria", ("kraken2", "bacteria"))
            db_entry = get_db_entry(db_cat, db_name)
            db_desc = db_entry.get("description") if db_entry else None
            self.utilities.add_database(
                self.utilities.database_item(
                    f"{db_cat}/{db_name}",
                    sofprep.dbs["kraken2"]["db"],
                    success_install,
                    software=db_cat,
                    description=db_desc,
                )
            )

        ########################## krakenuniq ###############################
        if self.layout.install_krakenuniq:
            success_install = sofprep.kuniq_install(dbname=self.organism)
            if success_install:
                self.installed_software.append(
                    self.software_install_string("krakenuniq")
                )
            sw_name, sw_tag = self.layout.SOFTWARE_NAMES.get("install_krakenuniq", ("krakenuniq", "default"))
            self.utilities.add_software(
                self.utilities.software_item(
                    sw_name,
                    sofprep.dbs["krakenuniq"]["db"],
                    sw_tag,
                    success_install,
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
            success_install = sofprep.kaiju_dl_install(dbname="viral")
            if success_install:
                self.installed_software.append(self.software_install_string("kaiju"))

            sw_name, sw_tag = self.layout.SOFTWARE_NAMES.get("install_kaiju", ("kaiju", "viral"))
            self.utilities.add_software(
                self.utilities.software_item(
                    sw_name,
                    sofprep.dbs["kaiju"]["db"],
                    sw_tag,
                    success_install,
                    sofprep.envs["ROOT"] + sofprep.envs["kaiju"],
                )
            )

            # Also register as database
            db_cat, db_name = self.layout.DATABASE_NAMES.get("install_kaiju", ("kaiju", "viral"))
            db_entry = get_db_entry(db_cat, db_name)
            db_desc = db_entry.get("description") if db_entry else None
            self.utilities.add_database(
                self.utilities.database_item(
                    f"{db_cat}/{db_name}",
                    sofprep.dbs["kaiju"]["db"],
                    success_install,
                    software=db_cat,
                    description=db_desc,
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

            minimap2_install = sofprep.minimap2_install(dbname=fname, reference=fpath)
            self.utilities.add_software(
                self.utilities.software_item(
                    "minimap2",
                    fpath,
                    fname,
                    minimap2_install,
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

                sw_name, sw_tag = self.layout.SOFTWARE_NAMES.get("install_diamond", ("diamond", fname))
                self.utilities.add_software(
                    self.utilities.software_item(
                        sw_name,
                        sofprep.dbs["diamond"]["db"],
                        sw_tag,
                        install_success,
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

                    sw_name, sw_tag = self.layout.SOFTWARE_NAMES.get("install_blast", ("blast", "genome"))
                    self.utilities.add_software(
                        self.utilities.software_item(
                            sw_name,
                            sofprep.dbs["blast"]["db"],
                            sw_tag,
                            success_install,
                            sofprep.envs["ROOT"] + sofprep.envs["blast"],
                        )
                    )

        # install nuc databases using local files.
        for fname, fd_list in prepdl.fastas["nuc"].items():
            for fpath in fd_list:
                if os.path.getsize(fpath) > 15000000000:
                    continue

                minimap2_install = sofprep.minimap2_install(dbname=fname, reference=fpath)
                self.utilities.add_software(
                    self.utilities.software_item(
                        "minimap2",
                        fpath,
                        fpath,
                        minimap2_install,
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

                    sw_name, sw_tag = self.layout.SOFTWARE_NAMES.get("install_fastviromeexplorer", ("fastviromeexplorer", fname))
                    self.utilities.add_software(
                        self.utilities.software_item(
                            sw_name,
                            sofprep.dbs["fastviromeexplorer"]["db"],
                            sw_tag,
                            install_success,
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

                        sw_name, sw_tag = self.layout.SOFTWARE_NAMES.get("install_blast", ("blast", "genome"))
                        self.utilities.add_software(
                            self.utilities.software_item(
                                sw_name,
                                sofprep.dbs["blast"]["db"],
                                sw_tag,
                                install_success,
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
            logging.info("Downloading databases and software")
            # repdl = self.setup_dir(self.INSTALL_PARAMS)
            self.prep_dl()
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

            #if self.wdir.update:
            #    self._reregister_on_update()
            #
            #if self.setup_install_class.test:
            #    self._reregister_on_update()
#
    def _reregister_on_update(self):
        """
        Re-register all installed databases and software after an update.
        This ensures the SQLite database has current version info.
        """
        logging.info("Re-registering all installed databases and software after update...")

        self.utilities.reset_tables()

        for db_name, db_ver in self.wdir.db_versions.items():
            fasta_path = None

            if db_name in self.wdir.fastas.get("prot", {}):
                fasta_path = self.wdir.fastas["prot"].get(db_name)
            elif db_name in self.wdir.fastas.get("nuc", {}):
                fasta_path = self.wdir.fastas["nuc"].get(db_name)
            elif db_name in self.wdir.fastas.get("host", {}):
                fasta_path = self.wdir.fastas["host"].get(db_name)
            elif db_name in self.wdir.fastas.get("filter", {}):
                fasta_path = self.wdir.fastas["filter"].get(db_name)

            print(fasta_path)
            if fasta_path and isinstance(fasta_path, list) is False:
                fasta_path = [fasta_path]

            for path in fasta_path:
                if os.path.isfile(path):
                    try:
                        self.utilities.add_database(
                            self.utilities.database_item(
                                db_name,
                                path,
                                True,
                                version=db_ver.get("version"),
                                source_url=db_ver.get("source_url"),
                                file_mod_date=db_ver.get("file_mod_date"),
                        )
                    )
                        logging.info(f"Re-registered database: {db_name}")
                    except Exception as e:
                        logging.info(f"Failed to re-register database {db_name}: {e}")

        for db_name, db_info in self.setup_install_class.dbs.items():
            if isinstance(db_info, dict) and "db" in db_info:
                try:
                    db_version = db_info.get("version", "")
                    self.utilities.add_software(
                        self.utilities.software_item(
                            db_name,
                            db_info["db"],
                            db_info.get("dbname", "default"),
                            True,
                            self.setup_install_class.envs["ROOT"] + self.setup_install_class.envs.get(db_name, ""),
                            db_version=db_version,
                        )
                    )
                    logging.info(f"Re-registered software: {db_name}")
                except Exception as e:
                    logging.info(f"Failed to re-register software {db_name}: {e}")

        self.utilities.dump_database(self.wdir.home)
        self.utilities.dump_software(self.wdir.home)
        logging.info("Re-registration complete.")

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
