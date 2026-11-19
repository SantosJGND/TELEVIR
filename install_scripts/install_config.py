class TelevirLayout:
    # =============================================================================
    # DATABASE INSTALLATION IS NOW CONTROLLED IN sources.yaml
    # See sources.yaml 'install' and 'install_types' fields for each database entry
    # =============================================================================

    # hosts ### CHECK HOST LIBRARY FILE FOR AVAILABLE HOSTS ###
    HOSTS_TO_INSTALL = [
        "hg38",
        "sus_scrofa",
        #"aedes_albopictus",
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
    install_bowtie2_depletion = False
    install_bwa_host = True

    # classification software.
    # NOTE: centrifuge and kraken2 installations are now controlled
    # by the 'install' field in sources.yaml under databases.centrifuge
    # and databases.kraken2 sections.
    install_metaphlan = False
    install_voyager_viral = True
    install_krakenuniq = False
    # NOTE: kaiju installation is now controlled
    # by the 'install' field in sources.yaml under databases.kaiju section.
    install_diamond = True
    install_minimap2 = True
    install_fastviromeexplorer = True
    install_clark = False
    install_blast = True

    # Pre-built index configuration (list-based)
    # Add dbnames here to check for pre-built indices during installation
    # The corresponding entry must also be configured in sources.yaml under prebuilt_indices
    # Example:
    # PREBUILT_CENTRIFUGE_INDICES = ["my_custom", "another_index"]
    # PREBUILT_KRAKEN2_INDICES = ["my_custom", "custom_kraken"]
    PREBUILT_CENTRIFUGE_INDICES = ["fungi_curated"]
    PREBUILT_KRAKEN2_INDICES = []
    
    # check files
    check_index_files = True

    # Database name mapping: install flag -> (category, yaml_name)
    # Used to register databases with names matching sources.yaml
    # Note: refseq_prot and refseq_gen are dynamically generated, not in sources.yaml
    # NOTE: centrifuge and kraken2 entries are now handled by looping over
    # sources.yaml databases with install: true
    DATABASE_NAMES = {
        # Reference databases
        # RefSeq viral
        "install_refseq_viral_prot": ("refseq", "viral_protein"),
        "install_refseq_viral_gen": ("refseq", "viral_genome"),
        # RefSeq bacterial
        "install_refseq_bacterial_prot": ("refseq", "bacterial_protein"),
        "install_refseq_bacterial_gen": ("refseq", "bacterial_genome"),
        # Other protein databases
        "install_swissprot": ("protein", "swissprot"),
        "install_rvdb": ("protein", "rvdb"),
        "install_virosaurus": ("nucleotide", "virosaurus"),
        "install_refseq_16s": ("ribosomal_rna", "refseq_16s"),
        "install_ribo16s": ("ribosomal_rna", "silva_16s"),
        "install_ncbi_16s": ("ribosomal_rna", "ncbi_16s"),
        
        # NOTE: kaiju entries are now handled by looping over
        # sources.yaml databases with install: true
    }

    # Software name mapping: install flag -> (software_name, tag)
    # software_name: the tool name (e.g., kraken2, centrifuge)
    # tag: specific database variant (e.g., viral, bacteria, default)
    # NOTE: centrifuge and kraken2 entries are now handled by looping over
    # sources.yaml databases with install: true
    SOFTWARE_NAMES = {
        "install_metaphlan": ("metaphlan", "default"),
        # NOTE: kaiju entries are now handled by looping over
        # sources.yaml databases with install: true
        "install_krakenuniq": ("krakenuniq", "default"),
        "install_diamond": ("diamond", "swissprot"),
        "install_voyager_viral": ("voyager", "viral"),
        "install_clark": ("clark", "viral"),
        "install_blast": ("blast", "genome"),
        "install_fastviromeexplorer": ("fastviromeexplorer", "viral"),
    }


    @property
    def install_hosts(self):
        if self.install_bwa_host or self.install_bowtie2_depletion:
            return True

        return False