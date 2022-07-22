#!/usr/bin/python3
import logging
import os
import sys


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


class main_setup:
    """
    prepare metagenomics run environments, databases.
    """

    def __init__(
        self,
        env_install,
        setup_dl,
        setup_install,
        pdir="",
        ENVS_PARAMS="",
        INSTALL_PARAMS="",
    ) -> None:
        if not ENVS_PARAMS:
            from install_config import ENVS_PARAMS
        if not INSTALL_PARAMS:
            from install_config import INSTALL_PARAMS

        self.ENVS_PARAMS = ENVS_PARAMS
        self.INSTALL_PARAMS = INSTALL_PARAMS
        self.wdir = setup_dl(INSTALL_PARAMS)
        self.env_install_class = env_install
        self.setup_install_class = setup_install

        if not pdir:
            pdir = os.getcwd()
        if pdir[-1] != "/":
            pdir += "/"
        self.pdir = pdir

    def user_input(self):
        args = get_args_install()
        self.envs = args.envs
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
        self.envs = envs
        self.seqdl = seqdl
        self.soft = soft
        self.nanopore = nanopore
        self.taxdump = taxdump
        self.setup_install_class.test = test
        self.organism = organism
        self.wdir.organism = organism

    def env_prepare(self, ENVS_PARAMS):
        """install environments described in ENVS_PARAMS script.

        :param ENVS_PARAMS: dictionary of environment parameters
        :type ENVS_PARAMS: dict
        :return: None
        """
        envprep = self.env_install_class()
        envprep.prep_dir(ENVS_PARAMS)
        envprep.conda_install()
        envprep.rabbitqc_install()
        envprep.flye_install()
        envprep.clark_install()
        envprep.fve_install()

        if self.nanopore:
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

    def setup_envs(self):

        if self.envs:
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
        self.wdir.refseq_dl()
        self.wdir.swissprot_dl()

        if self.organism == "viral":
            self.wdir.virosaurus_dl()
            self.wdir.RVDB_dl()

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
        sofprep.centrifuge_install(dbname=self.organism)
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

        ########################## clark ##################################
        # sofprep.clark_install(dbname=self.organism)

        ########################## kraken2 ###############################
        sofprep.kraken2_install(dbname=self.organism)
        krlib = f"kraken2-{self.organism}-library.fna.gz"
        if os.path.isfile(sofprep.dbs["kraken2"]["fasta"]):
            os.system(f"mv {sofprep.dbs['kraken2']['fasta']} {prepdl.seqdir}{krlib}")
        else:
            if not os.path.isfile(f"{prepdl.seqdir}{krlib}"):
                logging.info("kraken2 database not found.")

        if os.path.isfile(f"{prepdl.seqdir}{krlib}"):
            prepdl.fastas["nuc"]["kraken2"] = f"{prepdl.seqdir}{krlib}"

        ########################## krakenuniq ###############################
        sofprep.kuniq_install(dbname=self.organism)

        ### install viral specific databases
        if self.organism == "viral":
            sofprep.kaiju_viral_install()
            sofprep.virsorter_install()

        ### install prot databases using local files.
        for fname, fdb in prepdl.fastas["prot"].items():

            sofprep.diamond_install(dbname=fname, db=fdb)

            if fname == "refseq":

                sofprep.blast_install(
                    reference=fdb,
                    dbname=f"refseq_{self.organism}_prot",
                    nuc=False,
                    taxid_map=sofprep.metadir + "acc2taxid.prot.map",
                    args="-parse_seqids",
                    title=f"refseq {self.organism} prot",
                )

        ### install nuc databases using local files.
        for fname, fdb in prepdl.fastas["nuc"].items():
            sofprep.fve_install(
                reference=fdb,
                dbname=fname,
                virus_list=sofprep.metadir + f"{fname}-list.txt",
                list_create=True,
            )

            if fname == "refseq":
                sofprep.blast_install(
                    reference=prepdl.fastas["nuc"]["refseq"],
                    dbname=f"refseq_{self.organism}_genome",
                    nuc=True,
                    taxid_map=sofprep.metadir + "acc2taxid.nuc.map",
                    args="-parse_seqids",
                    title=f"refseq {self.organism} genome",
                )

            if nanopore:

                sofprep.deSAMBA_install(
                    reference=fdb,
                    dbname=fname,
                )

    def setup_soft(self):

        if self.seqdl or self.soft:

            # repdl = self.setup_dir(self.INSTALL_PARAMS)
            self.prep_dl()
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

    logging.basicConfig(
        format="%(asctime)s - %(message)s",
        datefmt="%d-%b-%y %H:%M:%S",
        level=logging.INFO,
    )

    metagen_prep = main_setup(env_install, setup_dl, setup_install)
    metagen_prep.user_input()
    metagen_prep.setup_envs()
    metagen_prep.setup_dir()
    metagen_prep.setup_soft()
