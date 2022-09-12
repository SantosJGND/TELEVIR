#!/usr/bin/python
##
## replace the variables bellow with paths in your system. Beware of permissions.
INSTALL_USER_CONF = "/televir/mngs_benchmark/config.py"

SOURCE = (
    "/home/bioinf/miniconda3/etc/profile.d/conda.sh"  # path to miniconda3 executable.
)
SOURCE_DIR = "/home/bioinf/miniconda3/"  # path to miniconda3 parent directory.
# /home/bioinf/Desktop/METAGEN/depo
HOME = "/home/bioinf/Desktop/METAGEN/db_third_test/"  # path to the directory where the databases will be stored.
ENVDIR = "/home/bioinf/Desktop/METAGEN/db_third_test/mngs_environments/"  # path to the directory where the environments will be stored.
DEPLOYMENT_DIR = "/home/bioinf/Desktop/CODE/ARGUS_nanopore/"  # path to the directory where the deployment will be stored.
TECH = "nanopore"  # technology used to generate the database and deployment files. options: "illumina" or "nanopore".
TAXDUMP = (
    "/home/bioinf/Desktop/METAGEN/taxdump.tar.gz"  # path to the taxdump.tar.gz file.
)
ORGANISM = "viral"  # organism name.
ENV = "/home/bioinf/Desktop/CODE/TELEVIR/.venv"
INSTALL_CONFIG = "minimal"
