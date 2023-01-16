#!/usr/bin/python
import os
import shutil
from distutils.dir_util import copy_tree


def get_args_deployment():
    """
    get user defined arguments.
    """
    try:
        import argparse

        parser = argparse.ArgumentParser(description="parse arguments")
        parser.add_argument(
            "--envdir",
            "-e",
            required=True,
            type=str,
            help="Environments root directory",
        ),

        parser.add_argument(
            "--dbmain",
            required=True,
            type=str,
            help="software db root directory",
        ),

        parser.add_argument(
            "--fdir",
            required=True,
            type=str,
            help="Fasta db root directory",
        ),

        parser.add_argument(
            "--mdir",
            required=True,
            type=str,
            help="metadata directory",
        ),

        parser.add_argument(
            "--source",
            type=str,
            required=True,
            help="source conda bash file.",
        )

        parser.add_argument(
            "--dir", "-d", required=True, type=str, help="deployment directory"
        ),

        parser.add_argument(
            "--technology",
            type=str,
            default="illumina",
            help="sequencing technology, illumina or nanopore [default= nanopore]",
        )

        args = parser.parse_args()

    except TypeError as e:
        print("check report args")
        print(e)

    return args


class main_deploy_prep:
    def __init__(self, pdir="") -> None:
        if not pdir:
            pdir = os.getcwd()

        if pdir[-1] != "/":
            pdir += "/"

        self.pdir = pdir
        self.bindir = self.pdir + "scripts/"
        self.module = self.pdir + "modules/metaruns_class.py"
        self.django_dir = self.pdir + "INSaFLU/"

    def user_input(self):
        args = get_args_deployment()
        self.envd = args.envdir
        self.dir = args.dir
        self.dbd = args.dbmain
        self.fmain = args.fdir
        self.tech = args.technology
        self.source = args.source
        self.metad = args.metad
        self.paramf = args.paramf

    def object_input(
        self,
        depdir,
        envd,
        dbdir,
        fdir,
        metad,
        source,
        tech="illumina",
        paramf="",
        docker_home="",
    ):
        self.envd = envd
        self.dir = depdir
        self.dbd = dbdir
        self.fmain = fdir
        self.tech = tech
        self.source = source
        self.metad = metad
        self.paramf = paramf
        self.docker_home = docker_home

    def read_available_software(self):
        """
        read available software from params file.
        """
        import sys

        sys.path.append(self.pdir)
        from utility_manager import installed_utilities

        self.software = soft

    def dir_prep(self):

        if not os.path.isdir(self.dir):
            os.makedirs(self.dir, exist_ok=True)

        print("copying files")
        print("dir: ", self.dir)
        copy_tree(self.django_dir, self.dir)

        self.app_dir = os.path.join(self.dir, "pathogen_identification") + "/"
        # os.makedirs(self.dir)

    def export(self):
        # self.paramf = self.pdir + f"params_files/params_{self.tech}.py"
        self.paramf = os.path.join(self.dir + "constants") + "/constants.py"
        env_file = self.dir + ".env_model"

        self.mainsh = self.pdir + f"main/main_{self.tech}.sh"
        new_params = self.app_dir + "televir_deploy_parameters.py"
        test_params_ont_json = os.path.join(self.dir, "product") + "/ont_params.json"
        test_params_illumina_json = (
            os.path.join(self.dir, "product") + "/illumina_params.json"
        )

        mods_dict = {
            "$ENVDIR": self.envd,
            "$DBDIR": self.dbd,
            "$FDIR": self.fmain,
            "$BINDIR": self.pdir + "scripts/",
            "$METADIR": self.metad,
            "$INSTALL_HOME": self.docker_home,
        }

        for tag, repl in mods_dict.items():
            if repl[-1] != "/":
                repl += "/"
            os.system(f"sed -i 's#{tag}#{repl}#g' {self.paramf}")
            os.system(f"sed -i 's#{tag}#{repl}#g' {env_file}")
            os.system(f"sed -i 's#{tag}#{repl}#g' {new_params}")
            # os.system(f"sed -i 's#{tag}#{repl}#g' {test_params_ont_json}")
            # os.system(f"sed -i 's#{tag}#{repl}#g' {test_params_illumina_json}")

        os.system(f"sed -i 's#$SOURCE#{self.source}#g' {new_params}")
        os.system(f"sed -i 's#$SOURCE#{self.source}#g' {self.paramf}")

        os.system(f"cp {self.pdir}metadata/taxid2desc.tsv {self.metad}")
        os.system(f"cp {self.pdir}README.md {self.fmain}")
        shutil.copy(env_file, self.dir + ".env")


if __name__ == "__main__":

    deploy_prep = main_deploy_prep()
    deploy_prep.user_input()
    deploy_prep.dir_prep()
    deploy_prep.export()
