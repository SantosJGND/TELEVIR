#!/usr/bin/python3
################################# DIRS

SOURCE = {
    "ENVSDIR": "$ENVDIR",
    "CONDA": "$SOURCE",
    "DBDIR_MAIN": "$DBDIR",
    "REF_FASTA": "$FDIR",
    "METAD": "$METADIR",
    "BIN": "$BINDIR",
}


DIRS = {
    "CLEAND": "reads/clean/",
    "FILTD": "reads/hd_filtered/",
    "HDOUT": "host_depletion/output/",
    "ASOUT": "assembly/output/",
    "CLASSD": "classification/",
    "REMD": "remap/",
    "LOGD": "logs/",
    "OUTD": "output/",
}

################################## MODULES

ACTIONS = {
    "CLEAN": False,
    "QCONTROL": True,
    "ENRICH": True,
    "DEPLETE": False,
    "ASSEMBLE": True,
    "CLASSIFY": True,
    "REMAP": True,
    "PHAGE_DEPL": True,
    "VIRSORT": False,
}

################################## SOFTWARE

SOFTWARE = {
    "QC": ["trimmomatic"],  # trimmomatic
    "HD": [
        "centrifuge",
    ],  # "virmet", "dvf", "minimap2", "centrifuge", "kaiju", "kuniq", "kraken2"],
    "ASSEMBLY_SOFT": ["spades"],  # "flye", "raven"],
    "ASSEMBLE_CLASS": ["diamond", "kraken2"],
    "CLASSM": ["kuniq", "kraken2"],  #
    "REMAP_SOFT": ["snippy"],  # snippy, rematch, bowtie, minimap-rem
}

################################## PARAMS

ARGS_ENRICH = {
    "trimmomatic": {
        "TRIM_ARGS": ["SLIDINGWINDOW:5:20 LEADING:3 TRAILING:3 TOPHRED33: MINLEN:35"],
        "TRIM_THREADS": [8],
        "FQ_THREADS": [8],
    },
    "kuniq": {
        "KUNIQ_ARGS": [
            "--threads 4 --gzip-compressed --fastq-input --hll-precision 12",
        ]
    },
    "kraken2": {"KRAKEN_ARGS": ["--threads 4 --gzip-compressed  --confidence .5 "]},
    "centrifuge": {"CENTRIFUGE_ARGS": ["-p 4 --time -k 3 --out-fmt sam"]},
    "diamond": {
        "DIAMOND_ARGS": [
            "-p 5 --top 5 -e 0.01 --id 65 --query-cover 50 --sensitive",
        ],
        "DIAMOND_DB": ["swissprot"],
    },
    "kaiju": {"KAIJU_ARGS": ["-z 4 -e 5 -s 65 -X -v"]},
    "minimap2": {
        "MINIMAP_ARGS": ["-t 4"],
        "MINIMAP_QCF": [20],  # filter for alignment score
        "MINIMAP_AX": ["sr"],
        "MINIMAP_DB": [
            "refseq_viral.genome.fna.gz",
        ],  ##NCBIrs_ViralCG.dustmasked.fna.gz, virosaurus90_vertebrate-20200330.fas.gz
    },
    "fve": {"FVE_ARGS": ["-cr 0.3 -co .1 -cn 5"], "FVE_DB": ["refseq", "virosaurus"]},
}


############################  ASSEMBLY PARAMS

ARGS_ASS = {
    "spades": {
        "SPADES_ARGS": [
            "--meta -t 8 --phred-offset 33 --only-assembler",
        ],
        "ASSEMBLY_LTRIM": [50],
    },
    "velvet": {
        "VELVET_K": [51],
        "VELVET_OPTS": ["-s 53 -e 65 -t 3"],
        "VELVET_ARGS": ["-exp_cov auto -cov_cutoff auto"],
        "VELVET_FILES": ["-fastq.gz"],
        "ASSEMBLY_LTRIM": [50],
    },
    "flye": {
        "FLYE_ARGS": [""],
    },
    "minimap-asm": {
        "MINIMAP_ARGS": ["-cx asm10 -t 4"],
        "MINIMAP_QCF": [20],  # filter for alignment score
        "MINIMAP_DB": [
            "refseq_viral.genome.fna.gz",
            "virosaurus90_vertebrate-20200330.fas.gz",
        ],
    },
}


############################# CLASSIFICATION PARAMS


ARGS_CLASS = {
    "Trim": {
        "TRIM_ARGS": ["SLIDINGWINDOW:5:20 LEADING:3 TRAILING:3 TOPHRED33: MINLEN:35"],
        "TRIM_THREADS": [4],
        "FQ_THREADS": [4],
    },
    "kuniq": {"KUNIQ_ARGS": ["--threads 4 --gzip-compressed --fastq-input"]},
    "kraken2": {"KRAKEN_ARGS": ["--threads 4 --gzip-compressed --confidence .5"]},
    "centrifuge": {"CENTRIFUGE_ARGS": ["-p 4 --time -k 3 --out-fmt sam"]},
    "diamond": {
        "DIAMOND_ARGS": [
            "-p 4 --top 5 -e 0.01 --id 60 --query-cover 60 --sensitive",
        ],
        "DIAMOND_DB": ["rvdb"],
    },
    "kaiju": {"KAIJU_ARGS": ["-z 3 -e 5 -s 65 -x -v"]},
    "minimap2": {
        "MINIMAP_ARGS": ["-t 4"],
        "MINIMAP_QCF": [20],  # filter for alignment score
        "MINIMAP_AX": ["map-ont"],
        "MINIMAP_DB": [
            "virosaurus90_vertebrate-20200330.fas.gz",
        ],  ##NCBIrs_ViralCG.dustmasked.fna.gz
    },
    "fve": {"FVE_ARGS": ["-cr 0.3 -co .1 -cn 5"], "FVE_DB": ["refseq", "virosaurus"]},
}


############################# REMAP PARAMS

ARGS_REMAP = {
    "snippy": {
        "SNIPPY_ARGS": ["--mapqual 60 --mincov 10"],
        "SNIPPY_RESOURCES": ["--cpus 3 --ram 8"],
        "REMAP_REF": ["refseq_viral.genome.fna.gz"],
        "MIN_COVERAGE_DEPTH": [2],
        "REF_KEEP": ["15"],
    },
    "rematch": {
        "REMATCH_RESOURCES": ["-j 4"],
        "REMATCH_ARGS": ["--minCovCall 10 --minCovPresence 5 --reportSequenceCoverage"],
        "REMATCH_ENV": ["/home/artic/Desktop/mngs_environments/remap/ReMatCh/ReMatCh/"],
        "REMAP_REF": ["refseq_viral.genome.fna.gz"],
        "MIN_COVERAGE_DEPTH": [2],
        "REF_KEEP": ["15"],
    },
    "bowtie": {
        "BOWTIE_RESOURCES": ["--threads 4"],
        "BOWTIE_ARGS": ["--sensitive-local"],
        "REMAP_REF": ["refseq_viral.genome.fna.gz"],
        "MIN_COVERAGE_DEPTH": [2],
        "REF_KEEP": ["15"],
    },
    "minimap-rem": {
        "MINIMAP_ARGS": ["-t 4"],
        "MINIMAP_AX": ["map-ont"],
        "MINIMAP_DB": ["refseq_viral.genome.fna.gz"],
        "REMAP_REF": [
            "refseq_viral.genome.fna.gz",
        ],
        "MIN_COVERAGE_DEPTH": [2],
        "REF_KEEP": ["15"],
    },
}
