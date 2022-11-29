#!/usr/bin/python3
import logging
import os
import sys
from distutils.command.install import install
from re import I


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
        # "virsorter",
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
    installed_dbs = []

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
        self.env_install_class = env_install
        self.setup_install_class = setup_install
        self.install_config = install_config

        if not pdir:
            pdir = os.getcwd()
        if pdir[-1] != "/":
            pdir += "/"
        self.pdir = pdir

        self.setup_config()
        self.utilities = repository(db_path=self.wdir.home, install_type=install_type)

    def setup_config(self):
        if self.install_config == "full":

            try:
                from install_scripts.config import Televir_Layout_full
            except ModuleNotFoundError as e:
                print("check install config.py exists")
                sys.exit()

            self.layout = Televir_Layout_full()

        if self.install_config == "minimal":

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
        envprep = self.env_install_class()
        envprep.prep_dir(ENVS_PARAMS)
        envprep.conda_install()

    def env_prepare(self, ENVS_PARAMS):
        """install environments described in ENVS_PARAMS script.

        :param ENVS_PARAMS: dictionary of environment parameters
        :type ENVS_PARAMS: dict
        :return: None
        """
        envprep = self.env_install_class()

        envprep.prep_dir(ENVS_PARAMS)
        # envprep.rabbitqc_install()
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
                self.installed_dbs.append("requests")
                self.utilities.add_database(
                    self.utilities.database_item(
                        "requests",
                        self.wdir.fastas["nuc"]["requests"],
                        True,
                    )
                )

        if self.layout.install_refseq_prot:
            success_refprot = self.wdir.refseq_prot_dl()
            if success_refprot:
                self.installed_dbs.append("refseq_prot")

                self.utilities.add_database(
                    self.utilities.database_item(
                        "refseq_prot",
                        self.wdir.fastas["prot"]["refseq"],
                        True,
                    )
                )

        if self.layout.install_refseq_gen:
            success_refnuc = self.wdir.refseq_gen_dl()
            if success_refnuc:
                self.installed_dbs.append("refseq_gen")

                self.utilities.add_database(
                    self.utilities.database_item(
                        "refseq_gen",
                        self.wdir.fastas["nuc"]["refseq"],
                        True,
                    )
                )

        if self.layout.install_swissprot:
            swissprot_dl = self.wdir.swissprot_dl()
            if swissprot_dl:
                self.installed_dbs.append("swissprot")

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
                self.installed_dbs.append("hg38")
                self.utilities.add_database(
                    self.utilities.database_item(
                        "hg38",
                        self.wdir.fastas["host"]["hg38"],
                        True,
                    )
                )

        if self.organism == "viral":

            if self.layout.install_virosaurus:
                success_virosaurus = self.wdir.virosaurus_dl()
                if success_virosaurus:
                    self.installed_dbs.append("virosaurus")

                    self.utilities.add_database(
                        self.utilities.database_item(
                            "virosaurus",
                            self.wdir.fastas["nuc"]["virosaurus"],
                            True,
                        )
                    )

            if self.layout.install_rvdb:
                success_rvdb = self.wdir.RVDB_dl()
                if success_rvdb:
                    self.installed_dbs.append("rvdb")

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

    def db_generate(
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

        ## install databases using organism name only.
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
                ] = f"{prepdl.seqdir}{centlib}"  # add to fastas dict

                if install_success:
                    self.installed_software.append("centrifuge")

                    self.utilities.add_software(
                        self.utilities.software_item(
                            "centrifuge",
                            sofprep.dbs["centrifuge"]["db"],
                            "default",
                            True,
                            sofprep.envs["ROOT"] + sofprep.envs["centrifuge"],
                        )
                    )

        ########################## clark ##################################
        if self.layout.install_clark:
            success_install = sofprep.install_clark(dbname=self.organism)
            if success_install:
                self.installed_software.append("clark")

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
            success_install = sofprep.kraken2_install(dbname=self.organism)
            krlib = f"kraken2-{self.organism}-library.fna.gz"
            if os.path.isfile(sofprep.dbs["kraken2"]["fasta"]):
                os.system(
                    f"mv {sofprep.dbs['kraken2']['fasta']} {prepdl.seqdir}{krlib}"
                )
            else:
                if not os.path.isfile(f"{prepdl.seqdir}{krlib}"):
                    logging.info("kraken2 database not found.")

            if os.path.isfile(f"{prepdl.seqdir}{krlib}"):
                prepdl.fastas["nuc"]["kraken2"] = f"{prepdl.seqdir}{krlib}"

                if success_install:

                    self.installed_software.append("kraken2")

                    self.utilities.add_software(
                        self.utilities.software_item(
                            "kraken2",
                            sofprep.dbs["kraken2"]["db"],
                            "default",
                            True,
                            sofprep.envs["ROOT"] + sofprep.envs["kraken2"],
                        )
                    )

        ########################## krakenuniq ###############################
        if self.layout.install_krakenuniq:
            success_install = sofprep.kuniq_install(dbname=self.organism)
            if success_install:
                self.installed_software.append("krakenuniq")
                self.utilities.add_software(
                    self.utilities.software_item(
                        "krakenuniq",
                        sofprep.dbs["krakenuniq"]["db"],
                        "default",
                        True,
                        sofprep.envs["ROOT"] + sofprep.envs["krakenuniq"],
                    )
                )

        ### install viral specific databases
        if self.organism == "viral":

            if self.layout.install_kaiju:
                sofprep.kaiju_viral_install()
                self.installed_software.append("kaiju")

                self.utilities.add_software(
                    self.utilities.software_item(
                        "kaiju",
                        sofprep.dbs["kaiju"]["db"],
                        "default",
                        True,
                        sofprep.envs["ROOT"] + sofprep.envs["kaiju"],
                    )
                )

            if self.layout.install_virsorter:
                success_install = sofprep.virsorter_install()
                if success_install:
                    self.installed_software.append("virsorter")

                    self.utilities.add_software(
                        self.utilities.software_item(
                            "virsorter",
                            sofprep.dbdir + "virsorter",
                            "default",
                            True,
                            sofprep.envs["ROOT"] + sofprep.envs["virsorter"],
                        )
                    )

        ### install host dbs.

        for fname, fpath in prepdl.fastas["host"].items():
            bwa_install = sofprep.bwa_install(dbname=fname, reference=fpath)
            if bwa_install:
                self.installed_software.append("bwa")
                self.utilities.add_software(
                    self.utilities.software_item(
                        "bwa",
                        sofprep.dbs["bwa"]["db"],
                        fname,
                        True,
                        sofprep.envs["ROOT"] + sofprep.envs["bwa"],
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

        ### install prot databases using local files.
        for fname, fdb in prepdl.fastas["prot"].items():

            if self.layout.install_diamond:
                install_success = sofprep.diamond_install(dbname=fname, db=fdb)

                if install_success:
                    self.installed_software.append("diamond")

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
                        self.installed_software.append("blast")

                        self.utilities.add_software(
                            self.utilities.software_item(
                                "blastp",
                                sofprep.dbs["blast"]["db"],
                                f"refseq_{self.organism}_prot",
                                True,
                                sofprep.envs["ROOT"] + sofprep.envs["blast"],
                            )
                        )

        ### install nuc databases using local files.
        for fname, fdb in prepdl.fastas["nuc"].items():

            if self.layout.install_fastviromeexplorer:
                install_success = sofprep.fve_install(
                    reference=fdb,
                    dbname=fname,
                    virus_list=sofprep.metadir + f"{fname}-list.txt",
                    list_create=True,
                )
                if install_success:
                    self.installed_software.append("fastviromeexplorer")

                    self.utilities.add_software(
                        self.utilities.software_item(
                            "fastviromeexplorer",
                            sofprep.dbs["fastviromeexplorer"]["db"],
                            fname,
                            True,
                            sofprep.envs["ROOT"] + sofprep.envs["fastviromeexplorer"],
                        )
                    )

            if self.layout.install_blast:

                if fname == "refseq":
                    install_success = sofprep.blast_install(
                        reference=prepdl.fastas["nuc"]["refseq"],
                        dbname=f"refseq_{self.organism}_genome",
                        nuc=True,
                        taxid_map=sofprep.metadir + "acc2taxid.nuc.map",
                        args="-parse_seqids",
                        title=f"refseq {self.organism} genome",
                    )

                    if install_success:
                        self.installed_software.append("blast")

                        self.utilities.add_software(
                            self.utilities.software_item(
                                "blastn",
                                sofprep.dbs["blast"]["db"],
                                f"refseq_{self.organism}_genome",
                                True,
                                sofprep.envs["ROOT"] + sofprep.envs["blast"],
                            )
                        )

            if nanopore:

                if self.layout.install_desamba:
                    install_success = sofprep.deSAMBA_install(
                        reference=fdb,
                        dbname=fname,
                    )

                    if install_success:
                        self.installed_software.append("desamba")

                        self.utilities.add_software(
                            self.utilities.software_item(
                                "desamba",
                                sofprep.dbs["desamba"]["db"],
                                fname,
                                True,
                                sofprep.envs["ROOT"] + sofprep.envs["desamba"],
                            )
                        )

    def setup_soft(self):

        if self.seqdl or self.soft:

            # repdl = self.setup_dir(self.INSTALL_PARAMS)
            self.prep_dl()
            logging.info("Downloading databases and software")
            self.dl_metadata_prot()

            if self.soft:

                self.db_generate(
                    self.INSTALL_PARAMS,
                    self.wdir,
                    nanopore=self.nanopore,
                    taxdump=self.taxdump,
                )

                self.dl_metadata_nuc()


if __name__ == "__main__":

    from modules.db_install import setup_dl, setup_install
    from modules.env_install import env_install
    from modules.utility_manager import Utility_Repository

    logging.basicConfig(
        format="%(asctime)s - %(message)s",
        datefmt="%d-%b-%y %H:%M:%S",
        level=logging.INFO,
    )

    metagen_prep = main_setup(env_install, setup_dl, setup_install, Utility_Repository)
    metagen_prep.user_input()

    metagen_prep.setup_envs_conda()
    metagen_prep.setup_envs()
    metagen_prep.setup_dir()
    metagen_prep.setup_soft()
