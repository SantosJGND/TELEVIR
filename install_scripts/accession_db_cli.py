#!/usr/bin/env python3

import argparse
import os
import sys

try:
    from accession_store import AccessionStore
except ImportError:
    from install_scripts.modules.accession_store import AccessionStore

try:
    from tabulate import tabulate
    HAS_TABULATE = True
except ImportError:
    HAS_TABULATE = False


class AccessionDBCLI:
    def __init__(self, install_home=None):
        self.install_home = install_home or os.environ.get("INSTALL_HOME", "/opt/televir")
        self.metadir = os.path.join(self.install_home, "metadata")
        self.protein_db = os.path.join(self.metadir, "protein_accessions.db")
        self.nucleotide_db = os.path.join(self.metadir, "nucleotide_accessions.db")
        self.store = AccessionStore(self.metadir)

    def summary(self):
        return self.store.counts()

    def list_dbs(self, db_type="all"):
        return self.store.list_dbs(db_type)

    def search(self, accession=None, dbs=None, limit=50):
        return self.store.search(accession, dbs, limit)

    def coverage(self, missing_only=False):
        return self.store.coverage(missing_only)

    def export(self, dbs=None, output=None):
        if not output:
            raise ValueError("--output required")
        rows = []
        for seq_type in ("protein", "nuc"):
            rows.extend(self.store.export_records(seq_type, dbs))
        with open(output, "w") as handle:
            handle.write("accession\ttaxid\n")
            for acc, taxid in rows:
                handle.write(f"{acc}\t{taxid}\n")
        return len(rows)


def _print_table(rows, headers):
    if not rows:
        print("No results found.")
        return
    if HAS_TABULATE:
        print(tabulate(rows, headers=headers, tablefmt="grid"))
        return
    widths = [max(len(str(header)), *(len(str(row[i])) for row in rows)) for i, header in enumerate(headers)]
    print("  ".join(str(header).ljust(widths[i]) for i, header in enumerate(headers)))
    print("-" * (sum(widths) + 2 * (len(widths) - 1)))
    for row in rows:
        print("  ".join(str(value).ljust(widths[i]) for i, value in enumerate(row)))


def _print_summary(counts):
    print("\n=== Accession Database Summary ===\n")
    for label, seq_type in (("Protein", "protein"), ("Nucleotide", "nuc")):
        total = counts[f"{seq_type}_total"]
        with_taxid = counts[f"{seq_type}_with_taxid"]
        pct = (with_taxid / total * 100) if total else 0
        print(f"{label}: {total:,} records ({with_taxid:,} with taxid, {pct:.1f}%) | {counts[f'{seq_type}_size_mb']:.1f} MB")
    for entry in counts["metadata"]:
        print(
            f"  {entry['seq_type']}/{entry['db_key']}: state={entry['state']} "
            f"version={entry['version'] or '-'} fingerprint={entry['fingerprint'] or '-'}"
        )
    print(f"\nTotal Accessions: {counts['protein_total'] + counts['nuc_total']:,}")
    print(f"Total with Taxid: {counts['protein_with_taxid'] + counts['nuc_with_taxid']:,}")


def main():
    parser = argparse.ArgumentParser(
        description="TELE-Vir Accession Database CLI",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--home",
        default=None,
        help="Installation home directory (default: INSTALL_HOME env var or /opt/televir)",
    )
    subparsers = parser.add_subparsers(dest="command", help="Available commands")

    subparsers.add_parser("summary", help="Show overall summary statistics")

    list_parser = subparsers.add_parser("list", help="List databases with statistics")
    list_parser.add_argument(
        "--type",
        choices=["protein", "nuc", "all"],
        default="all",
        help="Database type to list (default: all)",
    )

    search_parser = subparsers.add_parser("search", help="Search accessions")
    search_parser.add_argument("accession", nargs="?", help="Accession to search for")
    search_parser.add_argument("--dbs", help="Filter by database name")
    search_parser.add_argument("--limit", type=int, default=50, help="Limit results (default: 50)")

    coverage_parser = subparsers.add_parser("coverage", help="Show taxid coverage report")
    coverage_parser.add_argument(
        "--missing", action="store_true", help="Show only databases with missing taxids"
    )

    export_parser = subparsers.add_parser("export", help="Export accession data to file")
    export_parser.add_argument("--dbs", help="Filter by database name")
    export_parser.add_argument("--output", "-o", required=True, help="Output file path")

    args = parser.parse_args()
    if not args.command:
        parser.print_help()
        sys.exit(1)

    cli = AccessionDBCLI(install_home=args.home)
    if args.command == "summary":
        _print_summary(cli.summary())
    elif args.command == "list":
        _print_table(cli.list_dbs(args.type), ["Type", "Database", "Total", "With Taxid", "Coverage"])
    elif args.command == "search":
        _print_table(cli.search(args.accession, args.dbs, args.limit), ["Type", "Accession", "Description", "TaxID", "File"])
    elif args.command == "coverage":
        rows = cli.coverage(args.missing)
        if args.missing:
            _print_table([row[:3] for row in rows], ["Type", "Database", "Missing Count"])
        else:
            _print_table(rows, ["Type", "Database", "Total", "With Taxid", "Coverage"])
    elif args.command == "export":
        try:
            count = cli.export(args.dbs, args.output)
        except ValueError as exc:
            print(f"Error: {exc}", file=sys.stderr)
            sys.exit(1)
        print(f"Exported {count:,} records to {args.output}")


if __name__ == "__main__":
    main()
