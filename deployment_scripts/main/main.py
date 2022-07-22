#!/user/bin/python


def get_args_deploy():
    """
    get user defined arguments.
    """
    try:
        import argparse

        parser = argparse.ArgumentParser(description="parse arguments")
        parser.add_argument(
            "--fofn",
            "-f",
            type=str,
            help="file of fastq files. one file path per line.",
        )
        parser.add_argument(
            "--fdir",
            "-d",
            type=str,
            help="dictory containing fofn files.",
        )

        parser.add_argument(
            "--nsup",
            "-n",
            type=int,
            default=1,
            help="number of parameter combinations before assembly, 0 will run all [default=0]",
        )

        parser.add_argument(
            "--nlow",
            "-m",
            type=int,
            default=1,
            help="number of parameter combinations downstrem of assembly, 0 will run all [default=0]",
        )

        parser.add_argument(
            "-o", "--odir", type=str, default="", help="Output directory."
        )

        parser.add_argument(
            "--clean",
            action="store_true",
            default=False,
            help="move output reports to final output directory, intermediate files and config files to run output directories",
        )

        parser.add_argument(
            "--fdel",
            action="store_true",
            default=False,
            help="clean output repositories, keep only report files and assembly file. Recommend for large benchmarking runs.",
        )

        parser.add_argument(
            "--ref",
            type=str,
            required=False,
            default="",
            help="reference genome for host depletion",
        )

        args = parser.parse_args()

    except TypeError as e:
        print("check report args")
        print(e)

    return args


if __name__ == "__main__":
    import os
    from datetime import date

    from scripts.metaruns_class import meta_orchestra

    args = get_args_deploy()
    ###
    #
    if not args.odir:
        args.odir = "run_" + str(date.today()) + "/"
    else:
        if args.odir[-1] != "/":
            args.odir += "/"

    args.odir = os.path.join(os.getcwd(), args.odir)
    if args.fofn:

        # if os.path.exists(os.path.join(args.odir, os.path.basename(args.fofn))):
        #    print("skipping {}".format(args.fofn))
        #    # pass
        #
        # else:
        event = meta_orchestra(args.fofn, sup=args.nsup, down=args.nlow, odir=args.odir)
        event.reference = args.ref

        event.data_qc()
        event.sup_deploy(args.fofn)
        event.low_deploy()
        event.record_runs()

        # event.clean(delete=args.clean)

    elif args.fdir:
        if args.fdir[-1] != "/":
            args.fdir += "/"

        flist = [
            args.fdir + x
            for x in os.listdir(args.fdir)
            if os.path.splitext(x)[1] == ".fofn"
        ]

        for fofn in flist:
            if os.path.exists(os.path.join(args.odir, os.path.basename(fofn))):
                print("skipping {}".format(fofn))
                # continue

            event = meta_orchestra(fofn, sup=args.nsup, down=args.nlow, odir=args.odir)
            event.reference = args.ref
            event.data_qc()
            #
            event.sup_deploy(fofn)
            event.low_deploy()
            event.record_runs()
            #
            if args.clean or args.fdel:
                event.clean(delete=args.fdel)
