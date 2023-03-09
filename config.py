#!/usr/bin/python
##
# replace the variables bellow with paths in your system. Beware of permissions.

SOURCE = (
    # path to miniconda3 executable.
    "/home/bioinf/miniconda3/etc/profile.d/conda.sh"
)
SOURCE_DIR = "/home/bioinf/miniconda3/"  # path to miniconda3 parent directory.
# /home/bioinf/Desktop/METAGEN/depo
# path to the directory where the databases will be stored.
HOME = "/home/xpto/INSaFLU/data/televir/"
# path to the directory where the environments will be stored.
ENVDIR = "/home/bioinf/Desktop/CODE/ENVS/televir_envs/"
# path to the directory where the deployment will be stored.
DEPLOYMENT_DIR = "/home/bioinf/Desktop/CODE/ARGUS_INSAFLU/"
# technology used to generate the database and deployment files. options: "illumina" or "nanopore".
TECH = "nanopore"
TAXDUMP = (
    # path to the taxdump.tar.gz file.
    "/home/bioinf/Desktop/METAGEN/taxdump.tar.gz"
)
ORGANISM = "viral"  # organism name.
ENV = "/home/bioinf/Desktop/CODE/TELEVIR/.venv"
INSTALL_CONFIG = "docker"
DATABASE = "insaflu_db"
