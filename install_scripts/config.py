class Televir_Layout:

    # databases. chose at least one.
    install_refseq_prot = False
    install_refseq_gen = True
    install_swissprot = False
    install_rvdb = False
    install_virosaurus = True

    # classification software.
    install_centrifuge = True
    install_kraken2 = False
    install_krakenuniq = True
    install_kaiju = False
    install_diamond = False
    install_minimap2 = False
    install_fastviromeexplorer = True
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

    # additional.
    install_virsorter = True
