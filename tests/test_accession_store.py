import gzip
import os
import sqlite3

import pytest

from install_scripts.modules.accession_store import AccessionStore


def write_fasta(path, headers):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with gzip.open(path, "wt") as handle:
        for header in headers:
            handle.write(f">{header}\nACGT\n")


@pytest.fixture
def metadir(tmp_path):
    return str(tmp_path / "metadata")


@pytest.mark.parametrize("seq_type", ["protein", "nuc"])
def test_build_complete_and_export(metadir, seq_type):
    fasta = os.path.join(metadir, f"{seq_type}.fa.gz")
    write_fasta(fasta, ["acc1 description", "acc2 description"])
    store = AccessionStore(metadir)
    assert store.needs_rebuild(seq_type, "test_db", [fasta], "v1")

    store.begin_db(seq_type, "test_db", [fasta], version="v1", source_url="https://example.org")
    if seq_type == "protein":
        records = [
            ("test_db", "acc1", "desc1", "acc1", None, fasta),
            ("test_db", "acc2", "desc2", "acc2", None, fasta),
        ]
    else:
        records = [
            ("test_db", "acc1", "acc1", None, fasta),
            ("test_db", "acc2", "acc2", None, fasta),
        ]
    store.insert_batch(seq_type, records)
    store.update_taxids(seq_type, [("test_db", "acc1", 123), ("test_db", "acc2", 456)])
    store.complete_db(seq_type, "test_db")

    assert not store.needs_rebuild(seq_type, "test_db", [fasta], "v1")
    assert store.search("acc1", "test_db")[0][3] == 123
    assert store.list_dbs(seq_type)[0][2] == 2
    assert store.counts()[f"{seq_type}_with_taxid"] == 2

    outfile = os.path.join(metadir, f"{seq_type}.tsv")
    assert store.export_tsv(seq_type, outfile) is True
    with open(outfile) as handle:
        lines = handle.read().splitlines()
    assert lines[0].startswith("acc\ttaxid")
    assert any(line.startswith("acc1\t123") for line in lines)


def test_staging_is_not_visible_until_complete(metadir):
    fasta = os.path.join(metadir, "db.fa.gz")
    write_fasta(fasta, ["acc1"])
    store = AccessionStore(metadir)
    store.begin_db("nuc", "db", [fasta])
    store.insert_batch("nuc", [("db", "acc1", "acc1", None, fasta)])

    assert store.search("acc1") == []
    store.complete_db("nuc", "db")
    assert len(store.search("acc1")) == 1


def test_fingerprint_and_version_detect_stale_data(metadir):
    fasta = os.path.join(metadir, "db.fa.gz")
    write_fasta(fasta, ["acc1"])
    store = AccessionStore(metadir)
    store.begin_db("nuc", "db", [fasta], version="v1")
    store.insert_batch("nuc", [("db", "acc1", "acc1", None, fasta)])
    store.complete_db("nuc", "db")
    assert not store.needs_rebuild("nuc", "db", [fasta], "v1")

    assert store.needs_rebuild("nuc", "db", [fasta], "v2")
    write_fasta(fasta, ["acc1 with a longer header"])
    assert store.needs_rebuild("nuc", "db", [fasta], "v1")


def test_register_legacy_marks_rows_stale(metadir):
    fasta = os.path.join(metadir, "old.fa.gz")
    write_fasta(fasta, ["acc1"])
    store = AccessionStore(metadir)
    store.ensure_schema()
    with sqlite3.connect(store.path("nuc")) as conn:
        conn.execute(
            """
            INSERT INTO nucleotide_accessions (dbs, acc, acc_in_file, taxid, file)
            VALUES ('old_db', 'acc1', 'acc1', 7, ?)
            """,
            (fasta,),
        )

    store.register_legacy({"old_db": {"version": "v1"}})
    metadata = [entry for entry in store.metadata() if entry["db_key"] == "old_db"]
    assert metadata[0]["state"] == "imported"
    assert store.needs_rebuild("nuc", "old_db", [fasta], "v1")

    store.begin_db("nuc", "old_db", [fasta], version="v1")
    store.insert_batch("nuc", [("old_db", "acc1", "acc1", None, fasta)])
    store.update_taxids("nuc", [("old_db", "acc1", 9)])
    store.complete_db("nuc", "old_db")
    assert store.search("acc1")[0][3] == 9
    assert not store.needs_rebuild("nuc", "old_db", [fasta], "v1")


def test_scratch_helpers(metadir):
    store = AccessionStore(metadir)
    scratch = os.path.join(metadir, "matches_1.db")
    conn = store.create_scratch(scratch)
    conn.execute("INSERT INTO matches VALUES ('db', 'acc1', 42)")
    conn.commit()
    conn.close()

    assert store.scratch_has_matches(scratch)
    assert store.read_scratch_matches(scratch, "db") == [("acc1", 42)]
    store.remove_scratch(scratch)
    assert not store.scratch_has_matches(scratch)


def test_export_records_and_missing_coverage(metadir):
    fasta = os.path.join(metadir, "db.fa.gz")
    write_fasta(fasta, ["acc1", "acc2"])
    store = AccessionStore(metadir)
    store.begin_db("protein", "db_a", [fasta])
    store.insert_batch("protein", [
        ("db_a", "acc1", "", "acc1", None, fasta),
        ("db_a", "acc2", "", "acc2", None, fasta),
    ])
    store.update_taxids("protein", [("db_a", "acc1", 1)])
    store.complete_db("protein", "db_a")

    assert store.export_records("protein", "db_a") == [("acc1", 1)]
    assert store.export_records("protein", "other") == []
    assert store.coverage(missing_only=True) == [["protein", "db_a", 1, 0, 0]]
