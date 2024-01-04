class Televir_Layout_docker:
    # databases. chose at least one.
    install_refseq_prot = True
    install_refseq_gen = True
    install_swissprot = False
    install_rvdb = False
    install_virosaurus = True
    install_request_sequences = False

    # hosts ### CHECK HOST LIBRARY FILE FOR AVAILABLE HOSTS ###
    HOSTS_TO_INSTALL = [
        "hg38",
        "sus_scrofa",
        "aedes_albopictus",
        "gallus_gallus",
        "oncorhynchus_mykiss",
        "salmo_salar",
        "bos_taurus",
        "neogale_vison",
        "marmota_marmota",
        "culex_pipiens",
        "anas_platyrhynchos",
        "pipistrellus_kuhlii",
        "phlebotomus_papatasi",
        "felis_catus",
        "canis_lupus_familiaris",
        "cyprinus_carpio",
    ]

    # host mappers
    install_bowtie2_remap = True
    install_bowtie2_depletion = False
    install_bwa = True

    # classification software.
    install_centrifuge = True
    install_centrifuge_bacteria = False
    install_kraken2 = True
    install_kraken2_bacteria = False
    install_krakenuniq = True
    install_krakenuniq_fungi = False
    install_kaiju = True
    install_diamond = True
    install_minimap2 = False
    install_fastviromeexplorer = False
    install_clark = False
    install_desamba = False
    install_blast = True

    # assemblers.
    install_spades = True
    install_raven = True
    install_flye = True

    # technology setup (exclusive technologies. overrides info above).
    install_illumina = True
    install_nanopore = True


class Televir_Layout_full:
    # databases. chose at least one.
    install_refseq_prot = True
    install_refseq_gen = True
    install_swissprot = False
    install_rvdb = False
    install_virosaurus = True
    install_request_sequences = True

    # hosts
    HOSTS_TO_INSTALL = [
        "hg38",
        "sus_scrofa",
        "aedes_albopictus",
        "gallus_gallus",
        "oncorhynchus_mykiss",
        "salmo_salar",
        "bos_taurus",
        "neogale_vison",
        "marmota_marmota",
        "culex_pipiens",
        "anas_platyrhynchos",
        "pipistrellus_kuhlii",
        "phlebotomus_papatasi",
        "felis_catus",
        "canis_lupus_familiaris",
        "cyprinus_carpio",
    ]
    # host mappers
    install_bowtie2_remap = True
    install_bowtie2_depletion = True
    install_bwa = True

    # classification software.
    install_centrifuge = True
    install_centrifuge_bacteria = True
    install_kraken2 = True
    install_kraken2_bacteria = True
    install_krakenuniq = True
    install_krakenuniq_fungi = False
    install_kaiju = True
    install_diamond = True
    install_minimap2 = True
    install_fastviromeexplorer = True
    install_clark = False
    install_desamba = False
    install_blast = True

    # assemblers.

    install_spades = True
    install_raven = True
    install_flye = False

    # technology setup (exclusive technologies. overrides info above).
    install_illumina = True
    install_nanopore = True


class Televir_Layout_minimal:
    # databases. chose at least one.
    install_refseq_prot = False
    install_refseq_gen = True
    install_swissprot = True
    install_rvdb = False
    install_virosaurus = True
    install_request_sequences = True

    # hosts
    HOSTS_TO_INSTALL = [
        "hg38",
        "sus_scrofa",
        "aedes_albopictus",
        "gallus_gallus",
        "oncorhynchus_mykiss",
        "salmo_salar",
        "bos_taurus",
        "neogale_vison",
        "marmota_marmota",
        "culex_pipiens",
        "anas_platyrhynchos",
        "pipistrellus_kuhlii",
        "phlebotomus_papatasi",
        "felis_catus",
        "canis_lupus_familiaris",
        "cyprinus_carpio",
    ]

    # host mappers
    install_bowtie2_remap = False
    install_bowtie2_depletion = False
    install_bwa = True

    # classification software.
    install_centrifuge = True
    install_centrifuge_bacteria = True
    install_kraken2 = True
    install_kraken2_bacteria = True
    install_krakenuniq = True
    install_krakenuniq_fungi = False
    install_kaiju = True
    install_diamond = False
    install_minimap2 = False
    install_fastviromeexplorer = True
    install_clark = True
    install_desamba = True
    install_blast = True

    # assemblers.
    install_spades = True
    install_raven = True
    install_flye = True

    # technology setup (exclusive technologies. overrides info above).
    install_illumina = True
    install_nanopore = True
