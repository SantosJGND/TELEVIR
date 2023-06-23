#!/usr/bin/python

import os
import shutil
import subprocess


class env_install:
    def __init__(self, ENVDICT) -> None:
        self.envsdir = ENVDICT["ENVSDIR"]
        self.ymld = ENVDICT["YMLDIR"]
        self.envs = ENVDICT["ENVS"]
        self.git = ENVDICT["GIT"]
        self.tar = ENVDICT["TAR"]
        self.source = ENVDICT["SOURCE"]
        self.bin = ENVDICT["BIN"]

    def prep_dir(self):
        os.makedirs(self.envsdir, exist_ok=True)

    def conda_install(self, force_install=False):
        for ndir, yml in self.envs.items():
            ymlf = self.ymld + yml
            cdir = self.envsdir + ndir

            if "yml" not in ymlf:
                continue

            force = False
            if os.path.isdir(cdir):
                force = True

            command = [
                "conda",
                "env",
                "create",
                "-p",
                cdir,
                "--file",
                ymlf,
            ]

            if force:
                command.append("--force")

            if force and not force_install:
                continue

            print("installing environment %s" % cdir)

            subprocess.run(command)

    def fve_install(self, force_install=False):
        """
        FastViromeExplorer installation.
        """

        # GIT CLONE
        soft = "FastViromeExplorer/fve"
        sdir = soft.split("/")[0]
        git = self.git[soft]

        CWD = os.getcwd()
        os.chdir(self.envsdir)

        idir = git.split("/")[-1].strip(".git")
        exists = os.path.isdir(sdir)
        if not exists:
            os.mkdir(idir)
        os.chdir(idir)
        exists = os.path.isdir(idir)

        command = ["git", "clone", git]

        if not exists or force_install:
            subprocess.run(command)

            # JAVA INSTALLATION.
            tmpsh = "fve_install.sh"

            bash_lines = [
                "#!/bin/bash",
                f"source {self.source}",
                f"conda activate ./fve",
                f"cd {idir}",
                "javac -d bin src/*.java",
                "conda deactivate",
            ]

            os.system("touch " + tmpsh)
            with open(tmpsh, "w") as f:
                for l in bash_lines:
                    os.system('echo "{}" >> {}'.format(l, tmpsh))
            #                f.write("/n".join(bash_lines))

            subprocess.run(["chmod", "+x", tmpsh])
            subprocess.call(f"./{tmpsh}")

            os.system("rm " + tmpsh)
        os.chdir(CWD)

    def rabbitqc_install(self, force_install=False):
        """RabbitQC install"""

        soft = "preprocess/RabbitQC"
        sdir = os.path.join(self.envsdir, soft.split("/")[0])

        try:
            git = self.git[soft]
        except KeyError:
            print("No git repo for %s" % soft)
            return

        CWD = os.getcwd()
        os.chdir(self.envsdir)

        idir = os.path.join(sdir, git.split("/")[-1].strip(".git"))
        exists = os.path.isdir(idir)
        command = ["git", "clone", git]

        if not exists or force_install:
            os.makedirs(sdir)
            subprocess.run(command)
            os.chdir(idir)
            subprocess.run(["make"])

        os.chdir(CWD)

    def clark_install(self, force_install=False):
        """Clark install"""

        soft = "classification/Clark"
        sdir = os.path.join(self.envsdir, soft.split("/")[0])

        try:
            git = self.tar[soft]
        except KeyError:
            print("No git repo for %s" % soft)
            return

        CWD = os.getcwd()
        os.chdir(self.envsdir)

        idir = git.split("/")[-1].strip(".tar.gz").replace("CLARK", "CLARKSC")
        exists = os.path.isdir(soft)
        command = ["wget", git]

        if not exists or force_install:
            os.makedirs(soft, exist_ok=True)

            subprocess.run(command)

            subprocess.run(
                [
                    "tar",
                    "-xzvf",
                    git.split("/")[-1],
                ]
            )

            if os.path.isdir(soft):
                shutil.rmtree(soft)

            os.rename(
                idir,
                soft,
            )

            os.chdir(soft)
            subprocess.run(["sh", "install.sh"])

        os.chdir(CWD)

    def flye_install(self, force_install=False):
        """Flye install"""

        soft = "assembly/Flye"
        sdir = os.path.join(self.envsdir, soft.split("/")[0])

        try:
            git = self.git[soft]
        except KeyError:
            print("No git repo for %s" % soft)
            return

        CWD = os.getcwd()
        os.chdir(self.envsdir)

        idir = os.path.join(sdir, git.split("/")[-1].strip(".git"))
        exists = os.path.isdir(idir)
        command = ["git", "clone", git]

        if not exists or force_install:
            os.chdir(soft.split("/")[0])
            subprocess.run(command)
            os.chdir(idir)
            subprocess.run(["make"])

        os.chdir(CWD)

    def jellyfish_get(self):
        pid = "jellyfish"
        bin = os.path.join(
            self.envsdir,
            self.bin[pid],
            "bin",
            pid,
        )

        if os.path.isfile(bin):
            return bin

        return None

    def deSAMBA_install(self, force_install=False):
        """
        deSAMBA installation.
        """

        # GIT CLONE
        soft = "classm_lc/deSAMBA"
        sdir = self.envsdir + soft.split("/")[0]

        try:
            git = self.git[soft]
        except KeyError:
            print("No git repo for %s" % soft)
            return

        CWD = os.getcwd()
        os.chdir(self.envsdir)

        idir = sdir + "/" + git.split("/")[-1].strip(".git")
        exists = os.path.isdir(sdir)
        if not exists:
            os.mkdir(sdir)

        os.chdir(sdir)
        exists = os.path.isdir(idir)

        command = ["git", "clone", git]
        if not exists or force_install:
            subprocess.run(command)

            # INSTALLATION.
            ##
            tmpsh = "deSAMBA_install.sh"

            bash_lines = [
                "#!/bin/bash",
                f"cd {idir}",
                "cd ./src",
                "make -j 4",
                "cd ..",
                "mkdir bin",
                "cp ./src/deSAMBA ./bin",
            ]

            os.system("touch " + tmpsh)
            with open(tmpsh, "w") as f:
                for l in bash_lines:
                    os.system('echo "{}" >> {}'.format(l, tmpsh))

            subprocess.run(["chmod", "+x", tmpsh])
            subprocess.run(["chmod", "+x", idir + "/build"])
            subprocess.run(["chmod", "+x", idir + "/build-index"])
            print("Running deSAMBA installation script...")
            subprocess.call(f"./{tmpsh}")

            jelly_bin = self.jellyfish_get()
            if jelly_bin is not None:
                os.system(
                    "sed -i 's#./bin/jellyfish#{}#g' {}".format(
                        jelly_bin, os.path.join(idir, "build-index")
                    )
                )
            os.system("rm " + tmpsh)

        os.chdir(CWD)

    def trimmomatic_insaflu_install(self):
        """
        Trimmomatic installation.
        """

        # GIT CLONE
        soft = "trimmomatic"
        sdir = self.envsdir + soft.split("/")[0] + "/"

        CWD = os.getcwd()
        os.chdir(self.envsdir)

        ln_target = sdir + "classes/trimmomatic.jar"
        exists = os.path.exists(ln_target)
        print(ln_target, exists)
        if exists:
            return True

        os.makedirs(sdir, exist_ok=True)
        os.makedirs(sdir + "classes", exist_ok=True)
        os.makedirs(sdir + "adapters", exist_ok=True)

        os.chdir(sdir)
        zip_source = "http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip"
        zip_name = os.path.basename(zip_source).strip(".zip")

        get_cmd = f"wget -O trimmomatic-0.39.zip {zip_source}"
        unzip_cmd = "unzip trimmomatic-0.39.zip"
        rm_cmd = "rm trimmomatic-0.39.zip"

        os.system(get_cmd)
        os.system(unzip_cmd)
        os.system(rm_cmd)
        print((f"{sdir}{zip_name}/trimmomatic-0.39.jar"))

        print(os.path.exists(f"{sdir}{zip_name}/trimmomatic-0.39.jar"))

        ln_target = sdir + "classes/trimmomatic.jar"
        ln_cmd = f"ln -s {sdir}{zip_name}/trimmomatic-0.39.jar " + ln_target
        ln_targets_cmd = f"ln -s {sdir}{zip_name}/adapters/* " + sdir + "adapters/"
        print(ln_targets_cmd)

        os.system(ln_cmd)
        os.system(ln_targets_cmd)
        os.chdir(CWD)
        return True

    def fastqc_insaflu_install(self):
        """
        FastQC installation.
        """

        # GIT CLONE
        soft = "FastQC"
        sdir = self.envsdir + soft.split("/")[0] + "/"

        CWD = os.getcwd()
        os.chdir(self.envsdir)

        # exists = os.path.isdir(sdir)
        # if not exists:
        #    os.mkdir(sdir)

        zip_source = "https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip"
        zip_name = os.path.basename(zip_source).strip(".zip")
        exec_target = sdir + "/fastqc"

        exists = os.path.exists(exec_target)
        if exists:
            return True

        get_cmd = f"wget -O fastqc_v0.11.9.zip {zip_source}"
        unzip_cmd = "unzip fastqc_v0.11.9.zip"
        rm_cmd = "rm fastqc_v0.11.9.zip"

        os.system(get_cmd)
        os.system(unzip_cmd)
        os.system(rm_cmd)

        os.system(f"chmod +x {exec_target}")

        os.chdir(CWD)
        return True

    def rabbitqc_insaflu_install(self):
        # GIT CLONE
        soft = "RabbitQC"
        sdir = self.envsdir + soft.split("/")[0] + "/"

        CWD = os.getcwd()
        os.chdir(self.envsdir)

        zip_source = "https://github.com/ZekunYin/RabbitQC/archive/v0.0.1.zip"
        zip_name = os.path.basename(zip_source).strip(".zip")
        exec_target = sdir + "rabbit_qc"

        if os.path.exists(exec_target):
            return True

        try:
            wget_cmd = f"wget -O v0.0.1.zip {zip_source}"
            unzip_cmd = "unzip v0.0.1.zip"
            rm_cmd = "rm v0.0.1.zip"

            os.system(wget_cmd)
            os.system(unzip_cmd)
            os.system(rm_cmd)

            os.rename("RabbitQC-0.0.1", "RabbitQC")
            os.chdir(sdir)

            sed_make_cmd = (
                "sed 's/ -static//' Makefile > temp.txt; mv temp.txt Makefile"
            )
            os.system(sed_make_cmd)

            make_cmd = "make"
            os.system(make_cmd)

            return True

        except:
            return False

    def install_deployment_software(self):
        """Install deployment software: fastqc, trimmomatic, abricate, etc.
        INSAFLU specific.
        """

        success_fastqc = self.fastqc_insaflu_install()
        success_trimmomatic = self.trimmomatic_insaflu_install()
        success_rabbitqc = self.rabbitqc_insaflu_install()
