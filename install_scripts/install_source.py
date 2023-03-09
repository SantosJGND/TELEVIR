#!/usr/bin/python
# DIRS

# HOME = "/home/artic/Desktop/METAGEN/db_second_test/"
# SOURCE = "~/miniconda3/etc/profile.d/conda.sh"

#
# Installation directions.
# HOME : directory to install reference databases, sequences.
# ENVSDIR : Path to software binaries.
# SOURCE : conda.sh path.
# BINDIR : path to software bin not installed the conda way.
#

INSTALL_PARAMS = {
    "HOME": "/home/artic/Desktop/METAGEN/db_second_test/",
    "ENVSDIR": {
        "SOURCE": "~/miniconda3/etc/profile.d/conda.sh",
        "ROOT": "/home/artic/Desktop/METAGEN/db_second_test/mngs_environments/",
        "centrifuge": "hostDepletion/hostdep_env",
        "minimap2": "hostDepletion/hostdep_env",
        "jellyfish": "hostDepletion/hostdep_env",
        "diamond": "hostDepletion/hostdep_env",
        "kaiju": "hostDepletion/hostdep_env",
        "krakenuniq": "hostDepletion/hostdep_env",
        "blast": "hostDepletion/hostdep_env",
        "kraken2": "kraken2/kraken_env",
        "fastviromeexplorer": "FastViromeExplorer/FastViromeExplorer",
        "kallisto": "FastViromeExplorer/fve",
        "virsorter": "hostDepletion/vs2",
        "desamba": "classm_lc/deSAMBA",
        "clark": "classification/Clark",
        "bwa": "remap/remap",
        "bowtie2": "remap/remap",
    },
    "BINDIR": {"deSAMBA": "classm_lc/deSAMBA"},
    "REQUEST_REFERENCES": {
        "ACCID": "",
    },
}

#
# Directions to environment yaml and requirement files.
# SOURCE: conda source BASH script.
# ENVSDIR: directory where environments are to be installed in.
# YAMLDIR: requirements files directory.
# ENVS: environment directories (downstream of ENVSDIR) and requirement file names (in YAMLDIR).
#

ENVS_PARAMS = {
    "SOURCE": "~/miniconda3/etc/profile.d/conda.sh",
    "ENVSDIR": "/home/artic/Desktop/METAGEN/db_second_test/mngs_environments/",
    "YMLDIR": "yaml/",
    "ENVS": {
        "hostDepletion/hostdep_env": "HD.yml",
        "kraken2/kraken_env": "Krk2.yml",
        "hostDepletion/vs2": "virsorter.yml",
        "FastViromeExplorer/fve": "fve.yml",
        "assembly/assembly": "assemb.yml",
        "assembly/assemb_lc": "assemb_lc.yml",
        "preprocess/preproc": "prep.yml",
        "remap/remap": "remap.yml",
        "remap/Renv": "Renv.yml",
        "Pyenv/pyenv": "pyenv.yml",
        "classm_lc/venvlc": "venvlc.txt",
        # "plots/pyGenomeTracks": "pyGenomeTracks.yml",
    },
    "GIT": {
        "FastViromeExplorer/fve": "https://github.com/saima-tithi/FastViromeExplorer.git",
        "classm_lc/deSAMBA": "https://github.com/hitbc/deSAMBA.git",
        "preprocess/RabbitQC": "https://github.com/ZekunYin/RabbitQC.git",
        "assembly/Flye": "https://github.com/fenderglass/Flye.git",
        "trimmomatic": "https://github.com/usadellab/Trimmomatic.git",
    },
    "TAR": {
        "classification/Clark": "http://clark.cs.ucr.edu/Download/CLARKV1.2.6.1.tar.gz",
    },
    "PIP": {"classm_lc/venvlc": "venvlc.txt"},
    "BIN": {
        "jellyfish": "hostDepletion/hostdep_env",
    },
}
