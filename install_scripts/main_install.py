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
from install_scripts.load_sources import (
    get_db_entry,
    get_refseq_entry,
    get_prebuilt_index_path,
    get_prebuilt_indices_tool,
    get_prebuilt_indices,
    is_db_enabled,
    get_install_types,
    get_db_key,
    get_download_method,
    list_enabled_databases,
    list_all_db_configs,
)


# dbs key used by setup_install.install_metaphlan() when called with no dbname.
# Must match that method's default.
METAPHLAN_DEFAULT_DB = "mpa_vJan25_CHOCOPhlAnSGB_202503"


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

        if self.layout.install_clark:
           envprep.clark_install()

        if self.layout.install_fastviromeexplorer:
            envprep.fve_install()

        if self.layout.install_voyager_viral:
            envprep.voyager_install()


        envprep.centrifuge_install()

        envprep.m2m_install()

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
        Databases are now configured via sources.yaml 'install' and 'install_types' fields.
        Dynamically iterates over all enabled databases in configured categories.
        """

        ## DOWNLOAD HOSTS
        if self.layout.install_hosts:
            for host in self.layout.HOSTS_TO_INSTALL:
                _success = self.wdir.download_host(host)
                logging.info(f"host {host} download success: {_success}")

        # Process databases from sources.yaml
        download_categories = ["ribosomal_rna", "protein", "nucleotide", "refseq"]
        db_configs = list_all_db_configs(download_categories)

        for category, db_list in db_configs.items():
            if category == "refseq":
                self._download_refseq_databases_dynamic(db_list)
            else:
                self._download_category_databases(category, db_list)

        # Download taxdump if enabled in sources.yaml
        self._download_taxdump()

        # Process request sequences from sources.yaml
        self.wdir.install_requests()

    def _download_category_databases(self, category: str, db_configs: list):
        """
        Download databases for a category using config from sources.yaml.
        
        Args:
            category: Database category (e.g., 'ribosomal_rna', 'protein')
            db_configs: List of config dicts from list_all_db_configs
        """
        for config in db_configs:
            db_name = config['name']
            db_key = config['db_key']
            entry = config['entry']
            
            # Download the file using filter download method
            success = self.wdir.dl_filter_file(entry, fname=db_name, db_key=db_key)
            logging.info(f"{db_name} download success: {success}")
            
            # Get file path
            db_path = self.wdir.fastas.get(db_key, {}).get(db_name, [""])[0]
            db_ver = self.wdir.db_versions.get(db_name, {})
            db_desc = entry.get("description", "")
            
            # Register database
            self.utilities.add_database(
                self.utilities.database_item(
                    name=db_name,
                    path=db_path,
                    installed=success,
                    software=category,
                    db_type=db_key,
                    version=db_ver.get("version"),
                    source_url=db_ver.get("source_url"),
                    file_mod_date=db_ver.get("file_mod_date"),
                    description=db_desc,
                )
            )

    def _download_refseq_databases_dynamic(self, db_configs: list):
        """Download RefSeq databases using config from sources.yaml.
        
        Args:
            db_configs: List of config dicts for refseq databases
        """
        for config in db_configs:
            db_name = config['name']
            db_key = config['db_key']
            entry = config['entry']
            organism = config.get('refseq_organism')
            refseq_type = config.get('refseq_type', 'protein')
            
            # Get refseq entry
            refseq_entry = get_refseq_entry(organism, refseq_type)
            if not refseq_entry:
                logging.warning(f"No refseq entry for {organism}/{refseq_type}")
                continue
            
            url = refseq_entry.get("url", "")
            filename = refseq_entry.get("file_pattern", "")
            
            if refseq_type == "protein":
                success = self.wdir.refseq_prot_dl(url, filename, db_name)
            else:
                success = self.wdir.refseq_gen_dl(url, filename, db_name)
                
            
            if success:
                self.installed_databases.append(
                    self.database_install_string(db_name)
                )

            db_ver = self.wdir.db_versions.get(db_name, {})
            db_path = self.wdir.fastas.get(db_key, {}).get(db_name, [""])[0]
            db_desc = refseq_entry.get("description")
            
            self.utilities.add_database(
                self.utilities.database_item(
                    name=db_name,
                    path=db_path,
                    installed=success,
                    software="refseq",
                    db_type=db_key,
                    version=db_ver.get("version"),
                    source_url=db_ver.get("source_url"),
                    file_mod_date=db_ver.get("file_mod_date"),
                    description=db_desc,
                )
            )

    def _should_build_index(self, db_key: str, db_name: str, index_type: str) -> bool:
        """Check if an index should be built for a database.
        
        Args:
            db_key: The database key (filter, prot, nuc, host)
            db_name: The database name
            index_type: The index type to build (bwa, diamond, etc.)
            
        Returns:
            True if the index should be built based on sources.yaml config
        """
        category_map = {
            "filter": "ribosomal_rna",
            "prot": "protein",
            "nuc": "nucleotide",
            "host": "host_genomes",
        }
        category = category_map.get(db_key)
        if not category:
            return False
        
        install_types = get_install_types(category, db_name)
        return index_type in install_types or index_type in ["bwa", "diamond", "minimap2", "fastviromeexplorer"] and db_key in install_types

    def _download_taxdump(self):
        """Download taxdump if enabled in sources.yaml and not present."""
        taxdump_entry = get_db_entry('taxonomy', 'taxdump')
        if not taxdump_entry:
            logging.info("Taxdump entry not found in sources.yaml, skipping")
            return
        
        if not taxdump_entry.get('install', False):
            logging.info("Taxdump installation disabled in sources.yaml, skipping")
            return
        
        if self.wdir.check_taxdump_extracted():
            logging.info("Taxdump already extracted, skipping download")
            return
        
        taxdump_url = taxdump_entry.get('url')
        taxdump_path = os.path.join(self.wdir.metadir, "taxdump.tar.gz")
        
        success = self.wdir.dl_taxdump(taxdump_url, taxdump_path)
        if success:
            self.wdir.extract_taxdump(taxdump_path)
            logging.info("Taxdump downloaded and extracted successfully")
        else:
            logging.error("Failed to download taxdump")

    def dl_metadata_prot(self):
        """
        generate metadata files for prot databases.
        """
        self.wdir.prot_metadata()

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

        ######################### centrifuge (download) ###############################

        for db_name, entry in list_enabled_databases('centrifuge'):
            # key under which setup_install registered this db
            cent_key = self.organism if db_name == "viral" else db_name
            if db_name in ["viral"]:
                install_success = sofprep.centrifuge_install(dbname=self.organism)
                fasta_key = f"centrifuge_{db_name}"
            else:
                install_success = sofprep.centrifuge_download_install(
                    dbname=db_name, threads="16"
                )
                
                fasta_key = f"centrifuge_{db_name}"
            
            if install_success:
                self.installed_software.append(
                    self.software_install_string("centrifuge")
                )
            else:
                continue
            

            cent_entry = sofprep.get_dbs_entry("centrifuge", cent_key)
            centlib = f"refseq-{self.organism}.dust.fna.gz"

            if os.path.isfile(cent_entry.get("fasta", "")):
                os.system(
                    f"mv {cent_entry['fasta']} {prepdl.seqdir}{centlib}"
                )
            else:
                if not os.path.isfile(f"{prepdl.seqdir}{centlib}"):
                    logging.info(f"centrifuge {db_name} fasta database not found.")

            if os.path.isfile(f"{prepdl.seqdir}{centlib}"):
                prepdl.fastas["nuc"][fasta_key] = [
                    f"{prepdl.seqdir}{centlib}"
                ]



            sw_name, sw_tag = ("centrifuge", db_name)
            self.utilities.add_software(
                self.utilities.software_item(
                    sw_name,
                    cent_entry.get("db", ""),
                    sw_tag,
                    install_success,
                    sofprep.envs["ROOT"] + sofprep.envs["centrifuge"],
                    binary_name="centrifuge",
                )
            )

            # Also register as database
            db_cat = "centrifuge"
            db_entry = get_db_entry(db_cat, db_name)
            db_desc = db_entry.get("description") if db_entry else None
            self.utilities.add_database(
                self.utilities.database_item(
                    name=db_name,
                    path=cent_entry.get("db", ""),
                    installed=install_success,
                    software=db_cat,
                    db_type="filter",
                    description=db_desc,
                )
            )

        ######################### centrifuge (pre-built custom) ###############################
        # Install custom pre-built indices from PREBUILT_CENTRIFUGE_INDICES list
        
        for dbname in get_prebuilt_indices_tool("centrifuge"):

            install_success = sofprep.centrifuge_register_prebuilt(dbname=dbname)

            if not install_success:
                logging.info(f"Centrifuge/{dbname} not installed (prebuilt not found)")
                sofprep.dbs.setdefault("centrifuge", {})[dbname] = {
                    "status": "not_installed",
                    "db": "",
                }

            if install_success:
                self.installed_software.append(
                    self.software_install_string("centrifuge")
                )

            sw_name, sw_tag = ("centrifuge", dbname)
            self.utilities.add_software(
                self.utilities.software_item(
                    sw_name,
                    sofprep.get_dbs_entry("centrifuge", dbname).get("db", ""),
                    sw_tag,
                    install_success,
                    sofprep.envs["ROOT"] + sofprep.envs["centrifuge"],
                    binary_name="centrifuge",
                )
            )

            # Also register as database
            db_cat = "centrifuge"
            db_entry = get_db_entry(db_cat, dbname)
            db_desc = db_entry.get("description") if db_entry else None
            self.utilities.add_database(
                self.utilities.database_item(
                    name=dbname,
                    path=sofprep.get_dbs_entry("centrifuge", dbname).get("db", ""),
                    installed=install_success,
                    software=db_cat,
                    db_type="filter",
                    description=db_desc,
                )
            )

        ########################### clark ##################################

        if self.layout.install_clark:
            success_install = sofprep.clark_install(dbname=self.organism)

            if success_install:
                self.installed_software.append(self.software_install_string("clark"))

            sw_name, sw_tag = self.layout.SOFTWARE_NAMES.get("install_clark_viral", ("clark", "viral"))
            self.utilities.add_software(
                self.utilities.software_item(
                    sw_name,
                    sofprep.get_dbs_entry("clark", self.organism).get("db", ""),
                    sw_tag,
                    success_install,
                    sofprep.envs["ROOT"] + sofprep.envs["clark"],
                    binary_name="clark",
                )
            )

            # Also register as database
            db_cat = "clark"
            db_entry = get_db_entry(db_cat, dbname)
            db_desc = db_entry.get("description") if db_entry else None
            self.utilities.add_database(
                self.utilities.database_item(
                    name=dbname,
                    path=sofprep.get_dbs_entry("clark", self.organism).get("db", ""),
                    installed=install_success,
                    software=db_cat,
                    db_type="filter",
                    description=db_desc,
                )
            )


        ########################### voyager ###############################

        if self.layout.install_voyager_viral:
            success_install = sofprep.voyager_install_viruses_copy(dbname=self.organism)

            if success_install:
                self.installed_software.append(self.software_install_string("voyager"))

            sw_name, sw_tag = self.layout.SOFTWARE_NAMES.get("install_voyager_viral", ("voyager", "viral"))
            self.utilities.add_software(
                self.utilities.software_item(
                    sw_name,
                    sofprep.get_dbs_entry("voyager", self.organism).get("db", ""),
                    sw_tag,
                    success_install,
                    sofprep.envs["ROOT"] + sofprep.envs["voyager"],
                    binary_name="voyager",
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
                    sofprep.get_dbs_entry("metaphlan", METAPHLAN_DEFAULT_DB).get("db", ""),
                    sw_tag,
                    success_install,
                    sofprep.envs["ROOT"] + sofprep.envs["metaphlan"],
                    binary_name="metaphlan",
                )
            )

        ########################## kraken2 (download) ###############################

        for db_name, entry in list_enabled_databases('kraken2'):
            # Map db_name to install params
            kr_key = self.organism if db_name == "viral" else db_name
            if db_name == "viral":
                success_install = sofprep.kraken2_download_install(dbname=self.organism)
                krlib = f"kraken2-{self.organism}-library.fna.gz"
                fasta_key = "kraken2"
            elif db_name == "standard":
                success_install = sofprep.kraken2_download_install(dbname=db_name)
                krlib = f"kraken2-bacteria-library.fna.gz"
                fasta_key = "kraken2-bacteria"
            else:
                success_install = sofprep.kraken2_download_install(dbname=db_name)
                krlib = f"kraken2-{db_name}-library.fna.gz"
                fasta_key = f"kraken2-{db_name}"

            kr_entry = sofprep.get_dbs_entry("kraken2", kr_key)
            if os.path.isfile(kr_entry.get("fasta", "")):
                os.system(
                    f"mv {kr_entry['fasta']} {prepdl.seqdir}{krlib}"
                )
            else:
                if not os.path.isfile(f"{prepdl.seqdir}{krlib}"):
                    logging.info(f"kraken2 {db_name} database fasta not found.")

            if os.path.isfile(f"{prepdl.seqdir}{krlib}"):
                prepdl.fastas["nuc"][fasta_key] = [f"{prepdl.seqdir}{krlib}"]

            if success_install:
                self.installed_software.append(self.software_install_string("kraken2"))

            kraken_ver = kr_entry.get("version", "")
            sw_name, sw_tag = ("kraken2", db_name)
            self.utilities.add_software(
                self.utilities.software_item(
                    sw_name,
                    kr_entry.get("db", ""),
                    sw_tag,
                    success_install,
                    sofprep.envs["ROOT"] + sofprep.envs["kraken2"],
                    db_version=kraken_ver,
                    binary_name="kraken2",
                )
            )

            # Also register as database
            db_cat = "kraken2"
            db_entry = get_db_entry(db_cat, db_name)
            db_desc = db_entry.get("description") if db_entry else None
            self.utilities.add_database(
                self.utilities.database_item(
                    name=db_name,
                    path=kr_entry.get("db", ""),
                    installed=success_install,
                    software=db_cat,
                    db_type="filter",
                    version=kraken_ver,
                    description=db_desc,
                )
            )

        ########################## kraken2 (pre-built custom) ###############################
        # Install custom pre-built indices from PREBUILT_KRAKEN2_INDICES list
        for dbname in get_prebuilt_indices_tool("kraken2"):
            success_install = sofprep.kraken2_register_prebuilt(dbname=dbname)

            if not success_install:
                logging.info(f"Kraken2/{dbname} not installed (prebuilt not found)")
                sofprep.dbs.setdefault("kraken2", {})[dbname] = {
                    "status": "not_installed",
                    "db": "",
                }
            else:
                krlib = f"kraken2-{dbname}-library.fna.gz"
                pb_fasta = sofprep.get_dbs_entry("kraken2", dbname).get("fasta", "")
                if os.path.isfile(pb_fasta):
                    os.system(
                        f"mv {pb_fasta} {prepdl.seqdir}{krlib}"
                    )
                else:
                    if not os.path.isfile(f"{prepdl.seqdir}{krlib}"):
                        logging.info("kraken2 database fasta not found.")

                if os.path.isfile(f"{prepdl.seqdir}{krlib}"):
                    prepdl.fastas["nuc"][f"kraken2-{dbname}"] = [f"{prepdl.seqdir}{krlib}"]

                self.installed_software.append(self.software_install_string("kraken2"))

            pb_entry = sofprep.get_dbs_entry("kraken2", dbname)
            kraken_ver = pb_entry.get("version", "")
            sw_name, sw_tag = ("kraken2", dbname)
            self.utilities.add_software(
                self.utilities.software_item(
                    sw_name,
                    pb_entry.get("db", ""),
                    sw_tag,
                    success_install,
                    sofprep.envs["ROOT"] + sofprep.envs["kraken2"],
                    db_version=kraken_ver,
                    binary_name="kraken2",
                )
            )

            # Also register as database
            db_cat = "kraken2"
            db_entry = get_db_entry(db_cat, dbname)
            db_desc = db_entry.get("description") if db_entry else None
            self.utilities.add_database(
                self.utilities.database_item(
                    name=dbname,
                    path=pb_entry.get("db", ""),
                    installed=success_install,
                    software=db_cat,
                    db_type="filter",
                    version=kraken_ver,
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
            db_entry = get_db_entry("krakenuniq", self.organism)
            db_desc = db_entry.get("description") if db_entry else None
            self.utilities.add_software(
                self.utilities.software_item(
                    sw_name,
                    sofprep.get_dbs_entry("krakenuniq", self.organism).get("db", ""),
                    sw_tag,
                    success_install,
                    sofprep.envs["ROOT"] + sofprep.envs["krakenuniq"],
                    binary_name="krakenuniq",
                )
            )

            self.utilities.add_database(
                self.utilities.database_item(
                    name=sw_name,
                    path=sofprep.get_dbs_entry("krakenuniq", self.organism).get("db", ""),
                    installed=success_install,
                    software=sw_name,
                    db_type=self.organism,
                    description=db_desc,
                )
            )


        ########################## kaiju (download) ###############################

        for db_name, entry in list_enabled_databases('kaiju'):
            success_install = sofprep.kaiju_dl_install(dbname=db_name)

            if success_install:
                self.installed_software.append(self.software_install_string("kaiju"))

            kaiju_entry = sofprep.get_dbs_entry("kaiju", db_name)
            kaiju_ver = kaiju_entry.get("version", "")
            sw_name, sw_tag = ("kaiju", db_name)
            self.utilities.add_software(
                self.utilities.software_item(
                    sw_name,
                    kaiju_entry.get("db", ""),
                    sw_tag,
                    success_install,
                    sofprep.envs["ROOT"] + sofprep.envs["kaiju"],
                    db_version=kaiju_ver,
                    binary_name="kaiju",
                )
            )

            # Also register as database
            db_cat = "kaiju"
            db_entry = get_db_entry(db_cat, db_name)
            db_desc = db_entry.get("description") if db_entry else None
            self.utilities.add_database(
                self.utilities.database_item(
                    name=db_name,
                    path=kaiju_entry.get("db", ""),
                    installed=success_install,
                    software=db_cat,
                    db_type="filter",
                    version=kaiju_ver,
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
        """
        Build indices from downloaded databases based on sources.yaml install_types.
        """
        logging.info("install prepped")

        # Build indices for "filter" databases (16S rRNA, etc.)
        print(prepdl.fastas)
        for fname, fd_list in prepdl.fastas["filter"].items():
            print(f"Building index for filter database: {fname}")
            print(fd_list)
            for fpath in fd_list:
                if self._should_build_index("filter", fname, "bwa"):
                    bwa_install = sofprep.bwa_install(dbname=fname, reference=fpath)
                    self.installed_software.append(
                        self.software_install_string("bwa-filter")
                    )
                    self.utilities.add_software(
                        self.utilities.software_item(
                            "bwa-filter",
                            sofprep.get_dbs_entry("bwa", fname).get("fasta", ""),
                            fname,
                            bwa_install,
                            sofprep.envs["ROOT"] + sofprep.envs["bwa"],
                            tag="filter",
                            binary_name="bwa",
                        )
                    )
                    db_cat, _db_name = self.layout.DATABASE_NAMES.get("install_bwa_filter", ("bwa", fname))
                    db_entry = get_db_entry(db_cat, _db_name)
                    db_desc = db_entry.get("description") if db_entry else None
                    self.utilities.add_database(
                        self.utilities.database_item(
                            name=fname,
                            path=sofprep.get_dbs_entry("bwa", fname).get("fasta", ""),
                            installed=bwa_install,
                            software=db_cat,
                            db_type="filter",
                            description=db_desc,
                        )
                    )

        for fname, fd_list in prepdl.fastas["host"].items():

            for fpath in fd_list:   
                if self.layout.install_bwa_host:
                    bwa_install = sofprep.bwa_install(dbname=fname, reference=fpath)
                    self.installed_software.append(self.software_install_string("bwa"))
                    self.utilities.add_software(
                        self.utilities.software_item(
                            "bwa",
                            sofprep.get_dbs_entry("bwa", fname).get("fasta", ""),
                            fname,
                            bwa_install,
                            sofprep.envs["ROOT"] + sofprep.envs["bwa"],
                            tag="host",
                            binary_name="bwa",
                        )
                    )
                    self.utilities.add_database(
                        self.utilities.database_item(
                            name=fname,
                            path=sofprep.get_dbs_entry("bwa", fname).get("fasta", ""),
                            installed=bwa_install,
                            software="bwa",
                            db_type="host",
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
                                sofprep.get_dbs_entry("bowtie2", fname).get("db", ""),
                                fname,
                                bowtie2_install,
                                sofprep.envs["ROOT"] + sofprep.envs["bowtie2"],
                                tag="host",
                                binary_name="bowtie2",
                            )
                        )
                        self.utilities.add_database(
                            self.utilities.database_item(
                                name=fname,
                                path=sofprep.get_dbs_entry("bowtie2", fname).get("db", ""),
                                installed=bowtie2_install,
                                software="bowtie2",
                                db_type="host",
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
                        binary_name="minimap2",
                    )
                )
                self.utilities.add_database(
                    self.utilities.database_item(
                        name=fname,
                        path=fpath,
                        installed=minimap2_install,
                        software="minimap2",
                        db_type="host",
                    )
                )
        # install prot databases using local files.
        for fname, fd_list in prepdl.fastas["prot"].items():
            for fpath in fd_list:
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
                            sofprep.get_dbs_entry("diamond", fname).get("db", ""),
                            sw_tag,
                            install_success,
                            sofprep.envs["ROOT"] + sofprep.envs["diamond"],
                            binary_name="diamond",
                        )
                    )
                    self.utilities.add_database(
                        self.utilities.database_item(
                            name=sw_tag,
                            path=sofprep.get_dbs_entry("diamond", fname).get("db", ""),
                            installed=install_success,
                            software="diamond",
                            db_type="filter",
                        )
                    )

                if self.layout.install_blast:

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
                            sofprep.get_dbs_entry("blast", f"refseq_{self.organism}_prot").get("db", ""),
                            sw_tag,
                            success_install,
                            sofprep.envs["ROOT"] + sofprep.envs["blast"],
                            binary_name="blastn",
                        )
                    )

        # install nuc databases using local files.
        for fname, fd_list in prepdl.fastas["nuc"].items():
            print(f"Processing nuc database: {fname}")
            print(f"File paths: {fd_list}")
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
                        binary_name="minimap2",
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
                #                sofprep.get_dbs_entry("bowtie2", fname).get("db", ""),
                #                fname,
                #                True,
                #                sofprep.envs["ROOT"] + sofprep.envs["bowtie2"],
                #            )
                #        )

                if self.layout.install_fastviromeexplorer and len(fd_list) == 1:
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
                            sofprep.get_dbs_entry("fastviromeexplorer", fname).get("db", ""),
                            sw_tag,
                            install_success,
                            sofprep.envs["ROOT"]
                            + sofprep.envs["fastviromeexplorer"],
                            binary_name="FastViromeExplorer",
                        )
                    )

                    self.utilities.add_database(
                        self.utilities.database_item(
                            name=sw_tag,
                            path=sofprep.get_dbs_entry("fastviromeexplorer", fname).get("db", ""),
                            installed=install_success,
                            software="fastviromeexplorer",
                            db_type="viral",
                        )
                    )

                if self.layout.install_blast and len(fd_list) == 1:
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
                            sofprep.get_dbs_entry("blast", fname).get("db", ""),
                            sw_tag,
                            install_success,
                            sofprep.envs["ROOT"] + sofprep.envs["blast"],
                        )
                    )

                    self.utilities.add_database(
                        self.utilities.database_item(
                            name=sw_tag,
                            path=sofprep.get_dbs_entry("blast", fname).get("db", ""),
                            installed=install_success,
                            software="blastn",
                            db_type="viral",
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

                self.wdir.fastas["nuc"] = process_nuc_fasta_dict(
                    self.wdir.fastas["nuc"], max_file_size=4000000000
                )

                self.db_generate_external(
                    self.wdir,
                    sofprep,
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
                                name=db_name,
                                path=path,
                                installed=True,
                                version=db_ver.get("version"),
                                source_url=db_ver.get("source_url"),
                                file_mod_date=db_ver.get("file_mod_date"),
                        )
                    )
                        logging.info(f"Re-registered database: {db_name}")
                    except Exception as e:
                        logging.info(f"Failed to re-register database {db_name}: {e}")

        for tool, entries in self.setup_install_class.dbs.items():
            for db_name, db_info in entries.items():
                if isinstance(db_info, dict) and "db" in db_info:
                    try:
                        db_version = db_info.get("version", "")
                        self.utilities.add_software(
                            self.utilities.software_item(
                                tool,
                                db_info["db"],
                                db_name,
                                True,
                                self.setup_install_class.envs["ROOT"] + self.setup_install_class.envs.get(tool, ""),
                                db_version=db_version,
                            )
                        )
                        logging.info(f"Re-registered software: {tool}/{db_name}")
                    except Exception as e:
                        logging.info(f"Failed to re-register software {tool}/{db_name}: {e}")

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
    env_manager = env_install(ENVS_PARAMS, INSTALL_PARAMS)
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
