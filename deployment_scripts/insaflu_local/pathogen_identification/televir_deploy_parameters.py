from settings.constants_settings import ConstantsSettings as CS


class Params_Illumina:

    DATA_TYPE = "illumina"

    ################################## CONSTANT PARAMS

    CONSTANTS = {
        "minimum_coverage_threshold": 2,
        "max_output_number": 15,
        "taxid_limit": 12,
        "sift_query": "phage",
        "assembly_contig_min_length": 300,
    }

    ################################## SOFTWARE

    SOFTWARE = {
        CS.PIPELINE_NAME_read_quality_analysis: ["trimmomatic"],  # trimmomatic
        CS.PIPELINE_NAME_viral_enrichment: [
            "centrifuge",
            # "kraken2",
            # "fastviromeexplorer",
        ],  # "virmet", "dvf", "minimap2", "centrifuge", "kaiju", "kuniq", "kraken2", "fve", "virmet"],
        CS.PIPELINE_NAME_assembly: ["spades"],  # "flye", "raven"],
        CS.PIPELINE_NAME_contig_classification: [
            "blast",
            # "minimap2_asm",
        ],  # ["blast", "minimap-asm", "diamond", "virsorter"],
        CS.PIPELINE_NAME_read_classification: [
            # "kraken2",
            # "fastviromeexplorer",
            "krakenuniq",
            # "clark",
            "kaiju",
            "centrifuge",
            # "diamond",
        ],
        CS.PIPELINE_NAME_remapping: ["snippy"],  # snippy, rematch, bowtie, minimap-rem
    }

    ################################## PARAMS
    ARGS_QC = {
        "trimmomatic": {
            "TRIM_ARGS": [
                "SLIDINGWINDOW:5:20 LEADING:3 TRAILING:3 TOPHRED33: MINLEN:35"
            ],
            "TRIM_THREADS": [8],
            "FQ_THREADS": [8],
        },
        "rabbitqc": {
            "RABBIT_ARGS": [
                "-t 3 -f 3 -g --poly_g_min_len 10 -x --poly_x_min_len 10 -w 4"
            ],
        },
    }

    ARGS_ENRICH = {
        "clark": {
            "CLARK_ARGS": ["-m 2 --light"],
            "CLARK_THREADS": [8],
            "CLARK_DB": ["clark/viruses"],
        },
        "krakenuniq": {
            "KRAKENUNIQ_ARGS": [
                "--quick --min-hits 1",
            ],
            "KRAKENUNIQ_THREADS": [8],
            "KRAKENUNIQ_DB": ["kuniq/viral/"],
        },
        "kraken2": {
            "KRAKEN_ARGS": ["--confidence 0.3 --quick"],
            "KRAKEN_DB": ["kraken2/viral"],
            "KRAKEN_THREADS": [8],
        },
        "centrifuge": {
            "CENTRIFUGE_ARGS": ["-k 5 --min-hitlen 16"],
            "CENTRIFUGE_THREADS": [4],
            "CENTRIFUGE_DB": ["centrifuge/viral/index"],
        },
        "diamond": {
            "DIAMOND_ARGS": [
                "--top 5 -e 0.01 --id 40 --query-cover 40 --fast",
                "--top 5 -e 0.01 --id 40 --query-cover 40 --sensitive",
            ],
            "DIAMOND_DB": ["diamond/swissprot"],
            "DIAMOND_THREADS": [4],
        },
        "kaiju": {
            "KAIJU_ARGS": ["-e 5 -s 50 -x -v"],
            "KAIJU_DB": ["kaiju/viral/kaiju_db_viruses.fmi"],
            "KAIJU_THREADS": [4],
        },
        "minimap2": {
            "MINIMAP_ARGS": ["sr"],
            "MINIMAP_QCF": [20],  # filter for alignment score
            "MINIMAP_AX": ["sr"],
            "MINIMAP_DB": [
                "refseq_viral.genome.fna.gz",
                "virosaurus90_vertebrate-20200330.fas.gz",
            ],  ##NCBIrs_ViralCG.dustmasked.fna.gz, virosaurus90_vertebrate-20200330.fas.gz
        },
        "fastviromeexplorer": {
            "FVE_ARGS": ["-cr 0.1 -co .2 -cn 5"],
            "FVE_DB": [
                "fve/refseq/refseq.idx",
            ],
            "FVE_THREADS": [4],
        },
    }

    ############################  ASSEMBLY PARAMS

    ARGS_ASS = {
        "spades": {
            "SPADES_ARGS": [
                "--meta --phred-offset 33 --only-assembler",
            ],
            "ASSEMBLY_LTRIM": [50],
        },
        "velvet": {
            "VELVET_K": [51],
            "VELVET_OPTS": ["-s 53 -e 65"],
            "VELVET_ARGS": ["-exp_cov auto -cov_cutoff auto"],
            "VELVET_FILES": ["-fastq.gz"],
            "ASSEMBLY_LTRIM": [50],
        },
        "flye": {
            "FLYE_ARGS": [""],
        },
        "minimap-asm": {
            "MINIMAP_ARGS": ["-cx asm10"],
            "MINIMAP_QCF": [20],  # filter for alignment score
            "MINIMAP_DB": [
                "refseq_viral.genome.fna.gz",
                "virosaurus90_vertebrate-20200330.fas.gz",
            ],
        },
    }

    ############################# CLASSIFICATION PARAMS

    ARGS_CLASS = {
        "clark": {
            "CLARK_ARGS": ["-m 0", "-m 1", "-m 3"],
            "CLARK_THREADS": [8],
            "CLARK_DB": ["clark/viruses"],
        },
        "krakenuniq": {
            "KUNIQ_ARGS": [
                "--exact",
                "--exact --hll-precision 14",
                "--hll-precision 14",
            ],
            "KUNIQ_THREADS": [8],
            "KUNIQ_DB": ["kuniq/viral"],
        },
        "kraken2": {
            "KRAKEN_ARGS": ["--confidence .7", "--confidence 0.5"],
            "KRAKEN_DB": ["kraken2/viral"],
            "KRAKEN_THREADS": [4],
        },
        "centrifuge": {
            "CENTRIFUGE_ARGS": ["-k 3 --min-hitlen 22"],
            "CENTRIFUGE_THREADS": [4],
            "CENTRIFUGE_DB": ["centrifuge/viral/index"],
        },
        "blast": {
            "BLAST_ARGS": ["-evalue 1e-5 -max_target_seqs 5"],
            "BLAST_DB": ["blast/NUC/refseq_viral_genome"],
            "BLAST_THREADS": [4],
        },
        "blastp": {
            "BLAST_ARGS": ["-evalue 1e-5 -max_target_seqs 1"],
            "BLAST_DB": ["blast/PROT/refseq_viral_prot"],
            "BLAST_THREADS": [4],
        },
        "diamond": {
            "DIAMOND_ARGS": [
                "--top 5 -e 0.01 --id 40 --query-cover 40 --fast",
                "--top 3 -e 0.01 --id 50 --query-cover 50 --sensitive",
            ],
            "DIAMOND_DB": ["diamond/rvdb"],
            "DIAMOND_THREADS": [4],
        },
        "kaiju": {
            "KAIJU_ARGS": ['-e 3 -s 70 -X -v -a "mem"', '-e 3 -s 70 -X -v -a "mem"'],
            "KAIJU_DB": ["kaiju/viral/kaiju_db_viruses.fmi"],
            "KAIJU_THREADS": [4],
        },
        "minimap2": {
            "MINIMAP_ARGS": [""],
            "MINIMAP_QCF": [20],  # filter for alignment score
            "MINIMAP_AX": ["-sr"],
            "MINIMAP_DB": [
                "refseq_viral.genome.fna.gz",
                "virosaurus90_vertebrate-20200330.fas.gz",
            ],  ##NCBIrs_ViralCG.dustmasked.fna.gz
        },
        "minimap2_asm": {
            "MINIMAP_ARGS": [""],
            "MINIMAP_QCF": [20],  # filter for alignment score
            "MINIMAP_DB": [
                "refseq_viral.genome.fna.gz",
                "virosaurus90_vertebrate-20200330.fas.gz",
            ],
        },
        "fastviromeexplorer": {
            "FVE_ARGS": ["-cr 0.5 -co .5 -cn 10", "-cr 0.7 -co .7 -cn 15"],
            "FVE_DB": [
                "fve/refseq/refseq.idx",
            ],
            "FVE_THREADS": [4],
        },
    }

    ############################# REMAP PARAMS

    ARGS_REMAP = {
        "snippy": {
            "SNIPPY_ARGS": ["--mapqual 20 --mincov 10"],
            "SNIPPY_RESOURCES": ["--cpus 3 --ram 8"],
            "REMAP_REF": ["refseq_viral.genome.fna.gz"],
        },
        "rematch": {
            "REMATCH_RESOURCES": ["-j 4"],
            "REMATCH_ARGS": [
                "--minCovCall 10 --minCovPresence 5 --reportSequenceCoverage"
            ],
            "REMATCH_ENV": [
                "/home/artic/Desktop/mngs_environments/remap/ReMatCh/ReMatCh/"
            ],
            "REMAP_REF": ["refseq_viral.genome.fna.gz"],
        },
        "bowtie": {
            "BOWTIE_RESOURCES": ["--threads 4"],
            "BOWTIE_ARGS": ["--sensitive-local"],
            "REMAP_REF": ["refseq_viral.genome.fna.gz"],
        },
        "minimap-rem": {
            "MINIMAP_ARGS": ["-t 4"],
            "MINIMAP_AX": ["asm10"],
            "MINIMAP_DB": ["refseq_viral.genome.fna.gz"],
            "REMAP_REF": [
                "refseq_viral.genome.fna.gz",
            ],
            "MINIMAP_QCF": [0],
        },
    }


class Params_Nanopore:
    DATA_TYPE = "nanopore"

    CONSTANTS = {
        "minimum_coverage_threshold": 1,
        "max_output_number": 15,
        "taxid_limit": 12,
        "sift_query": "phage",
        "assembly_contig_min_length": 500,
    }

    ################################## SOFTWARE

    SOFTWARE = {
        CS.PIPELINE_NAME_read_quality_analysis: ["nanofilt"],  # "nanofilt", trimmomatic
        CS.PIPELINE_NAME_viral_enrichment: [
            # "",
            # "kaiju",
            # "krakenuniq",
            "centrifuge",
        ],
        CS.PIPELINE_NAME_assembly: ["raven"],  # spades, velvet,
        CS.PIPELINE_NAME_contig_classification: [
            # "blast",
            "minimap2",
        ],  # , "minimap-asm"],
        CS.PIPELINE_NAME_read_classification: [
            "krakenuniq",  #
            # "fastviromeexplorer",
            # "clark",  #
            # "desamba",  #
            # "kaiju",  #
            # "centrifuge",  #
            # "diamond",
        ],
        CS.PIPELINE_NAME_remapping: [
            "minimap2"
        ],  # snippy, rematch, bowtie, minimap-rem
    }

    ################################## PARAMS

    ARGS_QC = {
        "trimmomatic": {
            "TRIM_ARGS": [
                "SLIDINGWINDOW:5:20 LEADING:3 TRAILING:3 TOPHRED33: MINLEN:35"
            ],
            "TRIM_THREADS": [8],
            "FQ_THREADS": [8],
        },
        "nanofilt": {
            "NANOFILT_ARGS": [
                "-q 8 -l 50 --headcrop 30 --tailcrop 30 --maxlength 50000"
            ]
        },
        "rabbitqc": {
            "RABBIT_ARGS": [
                "-t 3 -f 3 -g --poly_g_min_len 10 -x --poly_x_min_len 10 -w 4"
            ],
        },
    }

    ARGS_ENRICH = {
        "clark": {
            "CLARK_ARGS": ["-m 2 --light"],
            "CLARK_THREADS": [8],
            "CLARK_DB": ["clark/viruses"],
        },
        "krakenuniq": {
            "KRAKENUNIQ_ARGS": [
                "--quick --min-hits 1",
            ],
            "KRAKENUNIQ_THREADS": [8],
            "KRAKENUNIQ_DB": ["kuniq/viral/"],
        },
        "kraken2": {
            "KRAKEN_ARGS": ["--confidence 0.3 --quick"],
            "KRAKEN_DB": ["kraken2/viral"],
            "KRAKEN_THREADS": [8],
        },
        "centrifuge": {
            "CENTRIFUGE_ARGS": ["-k 5 --min-hitlen 16"],
            "CENTRIFUGE_THREADS": [4],
            "CENTRIFUGE_DB": ["centrifuge/viral/index"],
        },
        "diamond": {
            "DIAMOND_ARGS": [
                "--top 5 -e 0.01 --id 40 --query-cover 40 --fast",
            ],
            "DIAMOND_DB": ["diamond/swissprot"],
            "DIAMOND_THREADS": [4],
        },
        "kaiju": {
            "KAIJU_ARGS": ["-e 5 -s 50 -x -v"],
            "KAIJU_DB": ["kaiju/viral/kaiju_db_viruses.fmi"],
            "KAIJU_THREADS": [4],
        },
        "minimap2_ONT": {
            "MINIMAP_ARGS": [""],
            "MINIMAP_QCF": [20],  # filter for alignment score
            "MINIMAP_DB": [
                "refseq_viral.genome.fna.gz",
                "virosaurus90_vertebrate-20200330.fas.gz",
            ],
        },
        "fastviromeexplorer": {
            "FVE_ARGS": ["-cr 0.1 -co .2 -cn 5"],
            "FVE_DB": [
                "fve/refseq/refseq.idx",
            ],
            "FVE_THREADS": [4],
        },
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
            "FLYE_ARGS": ["--threads 4 --scaffold --nano-raw"],
            "ASSEMBLY_LTRIM": [50],
        },
        "raven": {
            "RAVEN_ARGS": [
                "-p2 -f0.001 -w5 -k15",
            ],
            "ASSEMBLY_LTRIM": [50],
        },
        "miniasm": {
            "MINIASM_ARGS": ["-m 80 -i 0.05 -s 1000  -e 0.01"],
            "ASSEMBLY_LTRIM": [50],
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
        "clark": {
            "CLARK_ARGS": ["-m 1 --light"],
            "CLARK_THREADS": [8],
            "CLARK_DB": ["clark/viruses"],
        },
        "krakenuniq": {
            "KRAKENUNIQ_ARGS": [
                "--exact",
                "--exact --hll-precision 14",
                "--hll-precision 14",
            ],
            "KRAKENUNIQ_THREADS": [8],
            "KRAKENUNIQ_DB": ["kuniq/viral"],
        },
        "kraken2": {
            "KRAKEN_ARGS": ["--confidence .7", "--confidence 0.5"],
            "KRAKEN_DB": ["kraken2/viral"],
            "KRAKEN_THREADS": [4],
        },
        "centrifuge": {
            "CENTRIFUGE_ARGS": ["--min-hitlen 22"],
            "CENTRIFUGE_THREADS": [4],
            "CENTRIFUGE_DB": ["centrifuge/viral/index"],
        },
        "blast": {
            "BLAST_ARGS": ["-evalue 1e-5 -max_target_seqs 5"],
            "BLAST_DB": ["blast/NUC/refseq_viral_genome"],
            "BLAST_THREADS": [4],
        },
        "blastp": {
            "BLAST_ARGS": ["-evalue 1e-5 -max_target_seqs 1"],
            "BLAST_DB": ["blast/PROT/refseq_viral_prot"],
            "BLAST_THREADS": [4],
        },
        "diamond": {
            "DIAMOND_ARGS": [
                "--top 3 -e 0.01 --id 70 --query-cover 60 --algo 0 --very-sensitive --long-reads",
            ],
            "DIAMOND_DB": ["diamond/refseq"],
            "DIAMOND_THREADS": [4],
        },
        "kaiju": {
            "KAIJU_ARGS": ['-e 3 -s 70 -X -v -a "mem"'],
            "KAIJU_DB": ["kaiju/viral/kaiju_db_viruses.fmi"],
            "KAIJU_THREADS": [4],
        },
        "minimap2_ONT": {
            "MINIMAP_ARGS": [""],
            "MINIMAP_DB": [
                "refseq_viral.genome.fna.gz",
                "virosaurus90_vertebrate-20200330.fas.gz",
            ],
        },
        "minimap2": {
            "MINIMAP2_ARGS": [""],
            "MINIMAP2_QCF": [20],  # filter for alignment score
            "MINIMAP2_DB": [
                "refseq_viral.genome.fna.gz",
                "virosaurus90_vertebrate-20200330.fas.gz",
            ],
        },
        "fastviromeexplorer": {
            "FVE_ARGS": ["-cr 0.5 -co .5 -cn 10", "-cr 0.7 -co .7 -cn 15"],
            "FVE_DB": [
                "fve/refseq/refseq.idx",
                "fve/virosaurus/virosaurus.idx",
            ],
            "FVE_THREADS": [4],
        },
        "desamba": {
            "DESAMBA_ARGS": [""],
            "DESAMBA_DB": [
                "desamba/refseq",
                "desamba/virosaurus",
            ],
            "DESAMBA_QFC": [20],
        },
    }

    ############################# REMAP PARAMS

    ARGS_REMAP = {
        "snippy": {
            "SNIPPY_ARGS": ["--mapqual 40 --mincov 10"],
            "SNIPPY_RESOURCES": ["--cpus 3 --ram 8"],
            "REMAP_REF": ["refseq_viral.genome.fna.gz"],
        },
        "rematch": {
            "REMATCH_RESOURCES": ["-j 4"],
            "REMATCH_ARGS": [
                "--minCovCall 10 --minCovPresence 5 --reportSequenceCoverage"
            ],
            "REMATCH_ENV": [
                "/home/artic/Desktop/mngs_environments/remap/ReMatCh/ReMatCh/"
            ],
            "REMAP_REF": ["refseq_viral.genome.fna.gz"],
        },
        "bowtie": {
            "BOWTIE_RESOURCES": ["--threads 4"],
            "BOWTIE_ARGS": ["--sensitive-local"],
            "REMAP_REF": ["refseq_viral.genome.fna.gz"],
        },
        "minimap2": {
            "MINIMAP_THREADS": ["-t 4"],
            "MINIMAP_ARGS": ["-t 4 --secondary=no"],
            "MINIMAP_DB": ["refseq_viral.genome.fna.gz"],
            "MINIMAP_AX": ["map-ont"],
            "REMAP_REF": [
                "refseq_viral.genome.fna.gz",
            ],
            "MINIMAP_QCF": [0],
        },
        "CONSTANT": {
            "MIN_COVERAGE_DEPTH": [2],
            "REF_KEEP": ["15"],
            "CX_MAP_AX": ["asm20"],
            "CX_MAP_ARGS": ["-k 15 -N 1 -t 4"],
        },
    }
