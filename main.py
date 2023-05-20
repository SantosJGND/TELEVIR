#!/usr/bin/python


def get_args():
    """
    get user defined arguments.
    """

    try:
        import argparse

        parser = argparse.ArgumentParser(description="parse arguments")
        parser.add_argument(
            "--docker", action="store_true", default=False, help="docker installation"
        )
        parser.add_argument(
            "--config",
            type=str,
            required=False,
            default="config.py",
            help="docker installation",
        )
        parser.add_argument(
            "--test",
            action="store_true",
            default=False,
            help="test software installation",
        )

        parser.add_argument(
            "--envs",
            action="store_true",
            default=False,
            help="Install software environments",
        )

        parser.add_argument(
            "--setup_conda",
            action="store_true",
            default=False,
            help="Install software conda environments",
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
            "--deploy",
            action="store_true",
            default=False,
            help="setup software deployment",
        )

        parser.add_argument(
            "--partial",
            action="store_true",
            default=False,
            help="Install software databases",
        )

        args = parser.parse_args()

    except TypeError as e:
        print("check report args")
        print(e)

    return args


def main():

    from install_scripts.install_source import ENVS_PARAMS, INSTALL_PARAMS

    CWD = os.getcwd()
    args = get_args()
    ###
    if not args.partial:
        args.setup_conda = True
        args.envs = True
        args.soft = True
        args.seqdl = True
        args.deploy = True

    if args.docker:
        SOURCE = "/opt/conda/etc/profile.d/conda.sh"
        SOURCE_DIR = "/opt/conda/"
        HOME = "/televir/mngs_benchmark/"
        ENVDIR = "/televir/mngs_benchmark/mngs_environments/"
        DEPLOYMENT_DIR = "/televir/mngs_benchmark/app/"
        DATABASE = "insaflu_db"
        TECH = "nanopore"
        TAXDUMP = "/opt/taxdump.tar.gz"
        ORGANISM = "viral"
        INSTALL_CONFIG = "minimal"
        INSTALL_TYPE = "docker"

        INSTALL_PARAMS["HOME"] = HOME
        INSTALL_PARAMS["ENVSDIR"]["SOURCE"] = SOURCE
        INSTALL_PARAMS["ENVSDIR"]["ROOT"] = ENVDIR

        ENVS_PARAMS["SOURCE"] = SOURCE
        ENVS_PARAMS["ENVSDIR"] = ENVDIR
        ENVS_PARAMS["YMLDIR"] = CWD + "/install_scripts/yaml/"
        ENV = CWD + ".venv"

    else:

        try:
            mainconf = __import__(os.path.splitext(args.config)[0])
        except ImportError as e:
            logging.info(f"failed to import config file {args.config}")
            sys.exit(1)

        DEPLOYMENT_DIR = mainconf.DEPLOYMENT_DIR
        ENVDIR = mainconf.ENVDIR
        HOME = mainconf.HOME
        SOURCE = mainconf.SOURCE
        SOURCE_DIR = mainconf.SOURCE_DIR
        TAXDUMP = mainconf.TAXDUMP
        TECH = mainconf.TECH
        ORGANISM = mainconf.ORGANISM
        ENV = mainconf.ENV
        INSTALL_CONFIG = mainconf.INSTALL_CONFIG
        INSTALL_TYPE = "local"
        DATABASE = mainconf.DATABASE

        INSTALL_PARAMS["HOME"] = HOME
        INSTALL_PARAMS["ENVSDIR"]["SOURCE"] = SOURCE
        INSTALL_PARAMS["ENVSDIR"]["ROOT"] = ENVDIR

        ENVS_PARAMS["SOURCE"] = SOURCE
        ENVS_PARAMS["ENVSDIR"] = ENVDIR
        ENVS_PARAMS["YMLDIR"] = CWD + "/install_scripts/yaml/"

    ###
    from install_scripts.main_install import main_setup
    from install_scripts.modules.db_install import setup_dl, setup_install
    from install_scripts.modules.env_install import env_install
    from install_scripts.modules.utility_manager import Utility_Repository

    metagen_prep = main_setup(
        env_install,
        setup_dl,
        setup_install,
        Utility_Repository,
        ENVS_PARAMS=ENVS_PARAMS,
        INSTALL_PARAMS=INSTALL_PARAMS,
        pdir=CWD + "/install_scripts/",
        install_config=INSTALL_CONFIG,
        install_type=INSTALL_TYPE,
    )

    metagen_prep.object_input(
        envs=args.envs,
        seqdl=args.seqdl,
        soft=args.soft,
        nanopore=(TECH == "nanopore"),
        taxdump=TAXDUMP,
        organism=ORGANISM,
        test=args.test,
    )
    print(args.setup_conda)
    if args.setup_conda:
        metagen_prep.setup_envs_conda()

    if args.envs:
        metagen_prep.setup_envs()

    os.system(
        f"cp install_scripts/bin/* {ENVDIR + INSTALL_PARAMS['ENVSDIR']['centrifuge']}/bin/"
    )

    metagen_prep.setup_dir()

    if metagen_prep.seqdl or metagen_prep.soft:
        os.system(f"cp install_scripts/metadata/* {metagen_prep.wdir.metadir}")

    metagen_prep.setup_soft()
    ############
    ############

    if args.deploy:
        from deployment_scripts.main_deployment_setup import main_deploy_prep

        # begin by installing televir deployment requirements
        metagen_prep.setup_deploy()

        deploy_prep = main_deploy_prep(pdir=CWD + "/deployment_scripts/")
        deploy_prep.object_input(
            DEPLOYMENT_DIR,
            ENVDIR,
            metagen_prep.wdir.dbdir,
            metagen_prep.wdir.seqdir,
            metagen_prep.wdir.metadir,
            SOURCE_DIR,
            tech=TECH,
            docker_home=HOME,
        )

        deploy_prep.dir_prep()
        deploy_prep.export(DATABASE)


if __name__ == "__main__":
    import logging
    import os
    import sys

    logging.basicConfig(
        format="%(asctime)s - %(message)s",
        datefmt="%d-%b-%y %H:%M:%S",
        level=logging.INFO,
    )
    main()
