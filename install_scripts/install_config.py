class TelevirLayout:
    # databases. chose at least one.
    install_refseq_prot = True
    install_refseq_gen = True
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
    
    # check files
    check_index_files = True

    # Database name mapping: install flag -> (category, yaml_name)
    # Used to register databases with names matching sources.yaml
    DATABASE_NAMES = {
        "install_refseq_prot": ("refseq", "protein"),
        "install_refseq_gen": ("refseq", "genome"),
        "install_swissprot": ("protein", "swissprot"),
        "install_rvdb": ("protein", "rvdb"),
        "install_virosaurus": ("protein", "virosaurus"),
        "install_refseq_16s": ("ribosomal_rna", "refseq_16s"),
        "install_ribo16s": ("ribosomal_rna", "silva_16s"),
        "install_request_sequences": ("nucleotide", "requests"),
    }

    # Software name mapping: install flag -> (category, software_name, tag)
    # category: software category in sources.yaml (e.g., kraken2, centrifuge)
    # tag: specific database variant (e.g., viral, bacteria, default)
    SOFTWARE_NAMES = {
        "install_kraken2": ("kraken2", "viral"),
        "install_kraken2_bacteria": ("kraken2", "bacteria"),
        "install_kraken2_eupathdb46": ("kraken2", "eupathdb46"),
        "install_centrifuge": ("centrifuge", "viral"),
        "install_centrifuge_bacteria": ("centrifuge", "bacteria"),
        "install_metaphlan": ("metaphlan", "default"),
        "install_kaiju": ("kaiju", "viral"),
        "install_krakenuniq": ("krakenuniq", "default"),
        "install_diamond": ("diamond", "swissprot"),
        "install_voyager_viral": ("voyager", "viral"),
        "install_blast": ("blast", "genome"),
        "install_fastviromeexplorer": ("fastviromeexplorer", "viral"),
    }
