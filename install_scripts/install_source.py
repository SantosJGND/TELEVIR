#!/usr/bin/python
"""
Installation source configuration for TELE-Vir.

Centralizes environment parameters, Git repositories, and TAR archive URLs.
All software URLs are loaded from sources.yaml via load_sources module.

For new code, prefer importing directly from load_sources:
    from load_sources import get_db_url, get_software_url, get_git_url
"""

import os
from load_sources import get_loader, get_db_url, get_software_url, get_git_url, get_software_entry, list_software

def get_env_var(key, default=""):
    """Get environment variable with fallback to default"""
    return os.environ.get(key, default)


INSTALL_HOME = get_env_var("INSTALL_HOME", "/opt/televir")
ENVDIR = get_env_var("ENVDIR", "/opt/televir/environments")
SOURCE = get_env_var("SOURCE", "/opt/conda/etc/profile.d/conda.sh")
TAXDUMP = get_env_var("TAXDUMP", "/opt/taxdump.tar.gz")
REQUEST_SEQ_FILE = get_env_var("REQUEST_SEQ_FILE", "")
INSTALL_UPDATE = get_env_var("UPDATE", "false").lower() == "true"


def _load_software_from_sources():
    """Load software URLs from sources.yaml.
    
    Uses software keys from INSTALL_PARAMS['ENVSDIR'] to look up URLs in sources.yaml.
    Raises ValueError if required software is not found in sources.yaml.
    
    Returns:
        tuple: (git_dict, tar_dict)
    """
    software = list_software()
    archives = software.get('archives', {})
    git_repos = software.get('git_repos', {})
    
    git = {}
    tar = {}
    
    # Mapping from sources.yaml keys to env_install.py keys
    # sources.yaml name -> (env_install_key, type)
    software_mapping = {
        'fastviromeexplorer': ('FastViromeExplorer/fve', 'git'),
        'desamba': ('classm_lc/deSAMBA', 'git'),
        'rabbitqc': ('preprocess/RabbitQC', 'git'),
        'rabbitqc_git': ('preprocess/RabbitQC', 'git'),
        'trimmomatic_git': ('trimmomatic', 'git'),
        'trimmomatic': ('trimmomatic', 'archive'),
        'clark': ('classification/Clark', 'archive'),
        'voyager': ('classification/Voyager', 'archive'),
        'fastqc': ('fastqc', 'archive'),
    }
    
    for yaml_name, (env_key, src_type) in software_mapping.items():
        if src_type == 'archive':
            if yaml_name in archives:
                entry = archives[yaml_name]
                if entry and entry.get('url'):
                    tar[env_key] = entry['url']
        elif src_type == 'git':
            if yaml_name in git_repos:
                entry = git_repos[yaml_name]
                if entry and entry.get('url'):
                    git[env_key] = entry['url']
    
    return git, tar


# Load software from sources.yaml
_SOFTWARE_GIT, _SOFTWARE_TAR = _load_software_from_sources()


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
    "GIT": _SOFTWARE_GIT,
    "TAR": _SOFTWARE_TAR,
    "PIP": {"classm_lc/venvlc": "venvlc.txt"},
    "BIN": {
        "jellyfish": "hostDepletion/hostdep_env",
    },
    "SOURCE": SOURCE,
    "BIN": ENVDIR,
}


def get_git_url_for(key):
    """Get Git URL for a software package.
    
    Args:
        key: Software key (e.g., 'FastViromeExplorer/fve', 'preprocess/RabbitQC')
        
    Returns:
        Git repository URL or None
    """
    if key in ENVS_PARAMS.get('GIT', {}):
        return ENVS_PARAMS['GIT'][key]
    
    return get_git_url(key)


def get_tar_url_for(key):
    """Get TAR archive URL for a software package.
    
    Args:
        key: Software key (e.g., 'classification/Clark')
        
    Returns:
        TAR archive URL or None
    """
    if key in ENVS_PARAMS.get('TAR', {}):
        return ENVS_PARAMS['TAR'][key]
    
    return get_software_url(key)


def get_all_sources():
    """Get all source configurations."""
    return {
        'git': ENVS_PARAMS.get('GIT', {}),
        'tar': ENVS_PARAMS.get('TAR', {}),
        'sources_yaml': get_loader().sources,
    }
