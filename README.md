[![License: GPL v2](https://img.shields.io/badge/License-GPL%20v2-blue.svg)](https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html)

# TELE-Vir

TELE-Vir is a tool for the identification of viral sequences in metagenomic data. It is part of the INSaFLU project (https://insaflu.insa.pt/), a free bioinformatics platform for the laboratory surveillance of emerging and re-emerging viral pathogens. Originally developed for the identification of influenza viruses, it can be used for the identification and characterization of any viral sequence in metagenomic data.

INSaFLU-TELEVir is available as a web application at [https://insaflu.insa.pt](https://insaflu.insa.pt).

Documentation (latest) for each INSaFLU-TELEVir module is provided at [http://insaflu.readthedocs.io/](http://insaflu.readthedocs.io/).

## Synopsis

TELE-Vir is a pipeline that uses a combination of tools to identify viral sequences in metagenomic data.

It is composed of six main steps:

1. Quality control of the raw reads.
2. Viral Enrichment [optional].
3. Host Depletion [optional].
4. Assembly of the reads [optional].
5. Identification of the viral sequences.
   - Using reads.
   - using contigs (if assembled).
   - Using reference databases.
6. Remapping of the viral sequences to the reference genome.
   - This step produces the final output of the pipeline: a set of summary statistics and visualizations of the results.

This repository contains the code for the installation of the environments, software, databases and metadata used by the pipeline.

## Features

The pipeline is composed of several modules that can be used independently. The main features of the pipeline are:

- instalation general and software specific environments.
- installation of databases.
- generation of metadata files from the databases and NCBI.

Software available, by module:

Quality Control:

- FastQC
- NanoFilt
- Trimmomatic

Viral Enrichment:

- Centrifuge
- Kraken2

Host Depletion:

- BWA
- Minimap2

Assembly:

- SPAdes
- Raven
- Flye

Read Classification:

- Centrifuge
- FastViromeExplorer
- Kraken2
- Krakenuniq
- Kaiju

Contig Classification:

- Minimap2

Remapping:

- Snippy
- Bowtie2
- Minimap2

NOTE: Not all the modules are available for all the software. For example, the assembly module is only available for SPAdes and Raven.

## How to cite

If you use TELE-Vir in your research, please cite the following paper:

## Authors

João Dourado Santos, Miguel Pinheiro, Daniel Sobral, Joana Isidro, Miguel Pinto, João Paulo Gomes and Vítor Borges

## Funding

INSaFLU development is being co-funded by the European Commission on behalf of OneHealth EJP TELE-Vir project.
https://onehealthejp.eu/jrp-tele-vir/

## Acknowledgments

We thank the Infraestrutura Nacional de Computação Distribuída (INCD) for providing computational resources for testing the INSaFLU platform. INCD was funded by FCT and FEDER under the project 22153-01/SAICT/2016

# Installation

For an easy and rapid installation using docker see here [https://github.com/INSaFLU/docker](https://github.com/INSaFLU/docker).

This installation is oriented for Ubuntu Server 18.04.

### Dependencies

Install MiniConda:

```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```

Install and activate the `mngsbench_install.yml` environment provided in the root directory.

```
conda env create -f mngsbench_install.yml
conda activate mngsbench_install
```

#### Configuration

In the config.py file in the root directory, make sure that SOURCE points to an existing conda.sh file and SOURCE_DIR points to the miniconda installation directoiry. Then set the following variables:

HOME: the root directory for sequence and software database installation.
ENVSDIR: the directory where environments will be installed (typically within $HOME, but not necessarily).
TECH: the sequencing technology used. Options are: illumina, nanopore.
TAXDUMP: the path to the taxdump.tar.gz file. This should be the argument to a local instance of ncbi's taxdump.tar.gz. if given, the software will use this to rescue corruped files on download for those software that require it. Recommended when running locally.
ORGANISM: the organism to be used. Options are: viral.
ENV: the environment to be used. Options are: mngsbench.

Verify that these values correspond to those in the ENVS_PARAMS dictionary for the corresponding keys.

#### Deployment

The main_install.py script allows for four boolean tags:

- --envs: installs environments;
- --seqdl: downloads reference sequence databases;
- --soft: generates software databases.
- --nanopore: if given, will also install software specific to 3d generation sequencing technologies.

finally, main_install also accepts the argument `--taxdump`. This should be the argument to a local instance of ncbi's taxdump.tar.gz. if given, the software will use this to rescue corruped files on download for those software that require it. Recommended when running locally.

To install only the enrionments and databases required for the pipeline to run, run:

```
python main.py --envs --seqdl
```

To install all the software and databases, run:

```
python main.py
```

NOTE: when choosing `--nanopore`, installation of the software deSAMBA requires the installation of some dependencies using sudo. Verify that you have root priviledges.

Run
`python main_install.py -h `
for detail.
