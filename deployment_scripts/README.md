## NGS Metagenomics Benchmark Deployment


### Config

All params for the benchmarking run are contained in the `params.py` file. Please verify that this script corresponds to the goals of your run. Namely, the SOFTWARE dictionary, which controls which software is to be combined at each step of the pipeline. 

**abbreviations** : QC = Quality Control; HD= Host Depletion; CLASSM= Read Classification; 

Dictionaries below this one correspoond to software arguments and may be of interest too. 

### Setting up 

To setup a deployment directory:

i. make sure that you have: a $REF_DB directory with the necessary software dabases; a $REF_FASTA directory with the necessary reference sequence databases (see params/params*.py files); a $METADIR directory with the necessary acc2taxid and taxid2desc files. If none of this makes sense, refer to the install_scripts directory in this repository.

ii. miniconda installated in your system. 

iii. An $ENVSDIR directory with all the necessary environments for your software (see params/params*.py, or head to the install_scripts directory).

Finally, provide this information to the main_deployment_setup.py script, together with the technology of your data. Use `python main_deployment_setup.py -h` for details.

### Deployment

Change directory to the deployment directory setup following the above procedure. 

Deployment of one or more metagenomic runs is done through the `main.py` function. This script accepts the following arguments: 

- `--fofn` or `-f`: a "file of files" containing one path to a fastq per line. Note that the maximum is two lines for pair-end reads specifically. Nanopore reads must be concatenated previously into a single file. 

- `--dir` or `-d`: a directory containing multiple .fofn files. if given, applies pipeline to all files in directory with extention `.fofn`.

- `--odir` or `-o`: output directory. defaults to `run_DD-MM-YYYY`.

- `nsup` and `--nlow` : number of combinations of parameters / software to run upstream and downstream of assembly step. Default to `0`, which runs all possible combinations. Set to a smaller number for testing purposes. 

### output. 

Inside the `odir` directory, one directory is created for each fofn file provided. 

Inside each fofn subdirectory: 
    
- One run* directory per *pre-assembly* combination of parameters. 

- One *output/* directory. All output reports are stored in here. Directory Hierarchy = `pre-assembly combination` > `post-assemly combination`

