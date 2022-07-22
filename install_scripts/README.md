## Metagenomic Prep

Scripts to install metagenomic databases for the software used in benchmarking
the metagenomics installation pipeline. 

### Instructions: 

#### Environmnent

Install and activate the `mngsbench_install.yml` environment provided in the root directory. 

####  Configuration

In the INSTALL_PARAMS dictionary, provide ROOT as the root directory for sequence and software database installation. Within the INSTALL_PARAMS > ENVSDIR dictionary, make sure that SOURCE points to an existing conda sh file and ROOT points to the directory where environments will be installed (typically within $HOME, but not necessarily). 

Verify that these values correspond to those in the ENVS_PARAMS dictionary for the corresponding keys. 

#### Deployment

The main_install.py script allows for four boolean tags:

- --envs: installs environments; 
- --seqdl: downloads reference sequence databases; 
- --soft: generates software databases. 
- --nanopore: if given, will also install software specific to 3d generation sequencing technologies. 

finally, main_install also accepts the argument `--taxdump`. This should be the argument to a local instance of ncbi's taxdump.tar.gz. if given, the software will use this to rescue corruped files on download for those software that require it. Recommended when running locally. 

**warning** when choosing `--nanopore`, installation of the software deSAMBA requires the installation of some dependencies using sudo. Verify that you have root priviledges.

Run
`
python main_install.py -h 
`
for detail. 