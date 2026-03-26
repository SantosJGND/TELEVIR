class TelevirLayout:
    # databases. chose at least one.
    # RefSeq viral (default)
    install_refseq_viral_prot = True
    install_refseq_viral_gen = True
    # RefSeq bacterial (optional)
    install_refseq_bacterial_prot = False
    install_refseq_bacterial_gen = False
    # Other protein databases
    install_swissprot = False
    install_rvdb = False
    install_virosaurus = True
    install_request_sequences = False
    install_ribo16s = True
    install_refseq_16s = True

    # hosts ### CHECK HOST LIBRARY FILE FOR AVAILABLE HOSTS ###
    HOSTS_TO_INSTALL = [
        "hg38",
        "sus_scrofa",
        "aedes_albopictus",
        #"gallus_gallus",
        #"oncorhynchus_mykiss",
        #"salmo_salar",
        #"bos_taurus",
        #"neogale_vison",
        #"marmota_marmota",
        #"culex_pipiens",
        #"anas_platyrhynchos",
        #"pipistrellus_kuhlii",
        #"phlebotomus_papatasi",
        #"felis_catus",
        #"canis_lupus_familiaris",
        #"cyprinus_carpio",
    ]

    # host mappers
    install_bowtie2_remap = True
    install_bowtie2_depletion = False
    install_bwa_host = True
    install_bwa_filter = True

    # classification software.
    install_metaphlan = False
    install_voyager_viral = False
    install_centrifuge = True
    install_centrifuge_bacteria = False
    install_kraken2 = True
    install_kraken2_bacteria = False
    install_kraken2_eupathdb46 = False
    install_krakenuniq = True
    install_kaiju = True
    install_diamond = True
    install_minimap2 = True
    install_fastviromeexplorer = True
    install_blast = True

    # Pre-built index configuration (list-based)
    # Add dbnames here to check for pre-built indices during installation
    # The corresponding entry must also be configured in sources.yaml under prebuilt_indices
    # Example:
    # PREBUILT_CENTRIFUGE_INDICES = ["my_custom", "another_index"]
    # PREBUILT_KRAKEN2_INDICES = ["my_custom", "custom_kraken"]
    PREBUILT_CENTRIFUGE_INDICES = []
    PREBUILT_KRAKEN2_INDICES = []
    
    # check files
    check_index_files = True

    # Database name mapping: install flag -> (category, yaml_name)
    # Used to register databases with names matching sources.yaml
    # Note: refseq_prot and refseq_gen are dynamically generated, not in sources.yaml
    DATABASE_NAMES = {
        # Reference databases
        # RefSeq viral
        "install_refseq_viral_prot": ("refseq_db", "viral_protein"),
        "install_refseq_viral_gen": ("refseq_db", "viral_genome"),
        # RefSeq bacterial
        "install_refseq_bacterial_prot": ("refseq_db", "bacterial_protein"),
        "install_refseq_bacterial_gen": ("refseq_db", "bacterial_genome"),
        # Other protein databases
        "install_swissprot": ("protein", "swissprot"),
        "install_rvdb": ("protein", "rvdb"),
        "install_virosaurus": ("nucleotide", "virosaurus"),
        "install_refseq_16s": ("ribosomal_rna", "refseq_16s"),
        "install_ribo16s": ("ribosomal_rna", "silva_16s"),
        "install_request_sequences": ("taxonomy", "requests"),
        
        # Classification indices as databases (also saved as software)
        "install_kraken2": ("kraken2", "viral"),
        "install_kraken2_bacteria": ("kraken2", "bacteria"),
        "install_kraken2_eupathdb46": ("kraken2", "eupathdb46"),
        "install_centrifuge": ("centrifuge", "viral"),
        "install_centrifuge_bacteria": ("centrifuge", "bacteria"),
        "install_kaiju": ("kaiju", "viral"),
    }

    # Software name mapping: install flag -> (software_name, tag)
    # software_name: the tool name (e.g., kraken2, centrifuge)
    # tag: specific database variant (e.g., viral, bacteria, default)
    SOFTWARE_NAMES = {
        "install_kraken2": ("kraken2", "viral"),
        "install_kraken2_eupathdb46": ("kraken2", "eupathdb46"),
        "install_centrifuge": ("centrifuge", "viral"),
        "install_metaphlan": ("metaphlan", "default"),
        "install_kaiju": ("kaiju", "viral"),
        "install_krakenuniq": ("krakenuniq", "default"),
        "install_diamond": ("diamond", "swissprot"),
        "install_voyager_viral": ("voyager", "viral"),
        "install_blast": ("blast", "genome"),
        "install_fastviromeexplorer": ("fastviromeexplorer", "viral"),
    }
