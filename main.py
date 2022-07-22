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

        args = parser.parse_args()

    except TypeError as e:
        print("check report args")
        print(e)

    return args


def main():
    from install_scripts.install_config import ENVS_PARAMS, INSTALL_PARAMS

    CWD = os.getcwd()
    args = get_args()
    ###

    if args.docker:
        SOURCE = "/opt/conda/etc/profile.d/conda.sh"
        SOURCE_DIR = "/opt/conda/"
        HOME = "/usr/local/mngs_benchmark/"
        ENVDIR = "/usr/local/mngs_benchmark/mngs_environments/"
        DEPLOYMENT_DIR = "/usr/local/mngs_benchmark/deployment/"
        TECH = "nanopore"
        TAXDUMP = ""
        ORGANISM = "viral"

        INSTALL_PARAMS["HOME"] = HOME
        INSTALL_PARAMS["ENVSDIR"]["SOURCE"] = SOURCE
        INSTALL_PARAMS["ENVSDIR"]["ROOT"] = ENVDIR

        ENVS_PARAMS["SOURCE"] = SOURCE
        ENVS_PARAMS["ENVSDIR"] = ENVDIR
        ENVS_PARAMS["YMLDIR"] = CWD + "/install_scripts/yaml/"

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

    metagen_prep = main_setup(
        env_install,
        setup_dl,
        setup_install,
        ENVS_PARAMS=ENVS_PARAMS,
        INSTALL_PARAMS=INSTALL_PARAMS,
        pdir=CWD + "/install_scripts/",
    )
    metagen_prep.object_input(
        envs=True,
        seqdl=True,
        soft=True,
        nanopore=(TECH == "nanopore"),
        taxdump=TAXDUMP,
        organism=ORGANISM,
        test=args.test,
    )
    metagen_prep.setup_envs()

    os.system(
        f"cp install_scripts/bin/* {ENVDIR + INSTALL_PARAMS['ENVSDIR']['centrifuge']}/bin/"
    )

    metagen_prep.setup_dir()
    os.system(f"cp install_scripts/metadata/* {metagen_prep.wdir.metadir}")

    metagen_prep.setup_soft()
    ############
    ############
    from deployment_scripts.main_deployment_setup import main_deploy_prep

    deploy_prep = main_deploy_prep(pdir=CWD + "/deployment_scripts/")
    deploy_prep.object_input(
        DEPLOYMENT_DIR,
        ENVDIR,
        metagen_prep.wdir.dbdir,
        metagen_prep.wdir.seqdir,
        metagen_prep.wdir.metadir,
        SOURCE_DIR,
        tech=TECH,
    )

    deploy_prep.dir_prep()
    deploy_prep.export()


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
