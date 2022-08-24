#!/usr/bin/python
import os
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
        self.mainpy = self.pdir + "main/main.py"
        self.django_dir = self.pdir + "metagen_view/"

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
        self, depdir, envd, dbdir, fdir, metad, source, tech="illumina", paramf=""
    ):
        self.envd = envd
        self.dir = depdir
        self.dbd = dbdir
        self.fmain = fdir
        self.tech = tech
        self.source = source
        self.metad = metad
        self.paramf = paramf

    def dir_prep(self):

        if not os.path.isdir(self.dir):
            os.makedirs(self.dir, exist_ok=True)

        print("copying files")
        print("dir: ", self.dir)
        copy_tree(self.django_dir, self.dir)

        self.app_dir = os.path.join(self.dir, "pathogen_detection") + "/"
        # os.makedirs(self.dir)

    def export(self):
        self.paramf = self.pdir + f"params_files/params_{self.tech}.py"
        self.mainsh = self.pdir + f"main/main_{self.tech}.sh"
        new_params = self.app_dir + "params.py"
        test_params_ont_json = os.path.join(self.dir, "product") + "/ont_params.json"
        test_params_illumina_json = (
            os.path.join(self.dir, "product") + "/illumina_params.json"
        )

        os.system(f"cp {self.paramf} {new_params}")
        os.system(f"cp {self.mainpy} {self.app_dir}")

        mods_dict = {
            "$ENVDIR": self.envd,
            "$DBDIR": self.dbd,
            "$FDIR": self.fmain,
            "$BINDIR": self.pdir + "scripts/",
            "$METADIR": self.metad,
        }

        for tag, repl in mods_dict.items():
            if repl[-1] != "/":
                repl += "/"
            os.system(f"sed -i 's#{tag}#{repl}#g' {new_params}")
            os.system(f"sed -i 's#{tag}#{repl}#g' {test_params_ont_json}")
            os.system(f"sed -i 's#{tag}#{repl}#g' {test_params_illumina_json}")

        os.system(f"sed -i 's#$SOURCE#{self.source}#g' {new_params}")

        os.system(f"cp {self.pdir}metadata/taxid2desc.tsv {self.metad}")
        os.system(f"cp {self.pdir}README.md {self.fmain}")


if __name__ == "__main__":

    deploy_prep = main_deploy_prep()
    deploy_prep.user_input()
    deploy_prep.dir_prep()
    deploy_prep.export()
