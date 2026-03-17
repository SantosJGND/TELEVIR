#!/usr/bin/python
# DIRS

import os

def get_env_var(key, default=""):
    """Get environment variable with fallback to default"""
    return os.environ.get(key, default)


INSTALL_HOME = get_env_var("INSTALL_HOME", "/opt/televir")
ENVDIR = get_env_var("ENVDIR", "/opt/televir/environments")
SOURCE = get_env_var("SOURCE", "/opt/conda/etc/profile.d/conda.sh")
TAXDUMP = get_env_var("TAXDUMP", "/opt/taxdump.tar.gz")
REQUEST_SEQ_FILE = get_env_var("REQUEST_SEQ_FILE", "")
INSTALL_UPDATE = get_env_var("UPDATE", "false").lower() == "true"


INSTALL_PARAMS = {
    "HOME": INSTALL_HOME,
    "ENVSDIR": {
        "SOURCE": SOURCE,
        "ROOT": ENVDIR,
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
        "metaphlan": "Metaphlan/metaphlan",
        "voyager": "classification/Voyager",
        "bwa": "remap/remap",
        "bowtie2": "remap/remap",
        "prinseq": "preprocess/prinseq",
        "entrez_direct": "entrez_direct",
        "create_report": "Pyenv/igv_reports",
    },
    "BINDIR": {"deSAMBA": "classm_lc/deSAMBA"},
    "REQUEST_REFERENCES": {
        "ACCID": "",
        "FILE": REQUEST_SEQ_FILE,
    },
}

ENVS_PARAMS = {
    "SOURCE": SOURCE,
    "ENVSDIR": ENVDIR,
    "YMLDIR": "yaml/",
    "ENVS": {
        "hostDepletion/hostdep_env": "HD.yml",
        "kraken2/kraken_env": "Krk2.yml",
        "hostDepletion/vs2": "virsorter.yml",
        "FastViromeExplorer/fve": "fve.yml",
        "Metaphlan/metaphlan": "metaphlan.yml",
        "preprocess/preproc": "prep.yml",
        "remap/remap": "remap.yml",
        "remap/Renv": "Renv.yml",
        "Pyenv/pyenv": "pyenv.yml",
        "Pyenv/igv_reports": "igv_reports.yml",
        "classm_lc/venvlc": "venvlc.txt",
        "preprocess/prinseq": "prinseq.yml",
        "entrez_direct": "entrez_direct.yml",
    },
    "GIT": {
        "FastViromeExplorer/fve": "https://github.com/saima-tithi/FastViromeExplorer.git",
        "classm_lc/deSAMBA": "https://github.com/hitbc/deSAMBA.git",
        "preprocess/RabbitQC": "https://github.com/ZekunYin/RabbitQC.git",
        "trimmomatic": "https://github.com/usadellab/Trimmomatic.git",
    },
    "TAR": {
        "classification/Clark": "http://clark.cs.ucr.edu/Download/CLARKV1.2.6.1.tar.gz",
    },
    "PIP": {"classm_lc/venvlc": "venvlc.txt"},
    "BIN": {
        "jellyfish": "hostDepletion/hostdep_env",
    },
    "SOURCE": SOURCE,
    "BIN": ENVDIR,
}
