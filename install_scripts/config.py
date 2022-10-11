class Televir_Layout_full:

    # databases. chose at least one.
    install_refseq_prot = False
    install_refseq_gen = True
    install_swissprot = True
    install_rvdb = True
    install_virosaurus = True

    # classification software.
    install_centrifuge = True
    install_kraken2 = True
    install_krakenuniq = True
    install_kaiju = True
    install_diamond = True
    install_minimap2 = True
    install_fastviromeexplorer = True
    install_clark = False
    install_desamba = True
    install_blast = True

    # assemblers.

    install_spades = True
    install_raven = True
    install_flye = True

    # technology setup (exclusive technologies. overrides info above).
    install_illumina = True
    install_nanopore = True

    # additional.
    install_virsorter = True


class Televir_Layout_minimal:

    # databases. chose at least one.
    install_refseq_prot = False
    install_refseq_gen = True
    install_swissprot = True
    install_rvdb = False
    install_virosaurus = True

    # classification software.
    install_centrifuge = True
    install_kraken2 = False
    install_krakenuniq = True
    install_kaiju = True
    install_diamond = False
    install_minimap2 = False
    install_fastviromeexplorer = True
    install_clark = False
    install_desamba = True
    install_blast = True

    # assemblers.

    install_spades = True
    install_raven = True
    install_flye = True

    # technology setup (exclusive technologies. overrides info above).
    install_illumina = True
    install_nanopore = True

    # additional.
    install_virsorter = True
