import os

import pytest

from install_scripts.accession_db_cli import AccessionDBCLI
from install_scripts.modules.accession_store import AccessionStore


@pytest.fixture
def cli(tmp_path):
    metadir = str(tmp_path / "metadata")
    store = AccessionStore(metadir)

    store.begin_db("protein", "prot_db", [], version="v1")
    store.insert_batch("protein", [
        ("prot_db", "p1", "desc1", "p1", None, "p.fa.gz"),
        ("prot_db", "p2", "desc2", "p2", None, "p.fa.gz"),
    ])
    store.update_taxids("protein", [("prot_db", "p1", 11), ("prot_db", "p2", 12)])
    store.complete_db("protein", "prot_db")

    store.begin_db("nuc", "db_a", [], version="v1")
    store.insert_batch("nuc", [
        ("db_a", "n1", "n1", None, "n.fa.gz"),
        ("db_a", "n2", "n2", None, "n.fa.gz"),
    ])
    store.update_taxids("nuc", [("db_a", "n1", 21), ("db_a", "n2", 22)])
    store.complete_db("nuc", "db_a")

    store.begin_db("nuc", "db_b", [], version="v1")
    store.insert_batch("nuc", [("db_b", "n3", "n3", None, "n.fa.gz")])
    store.update_taxids("nuc", [("db_b", "n3", 23)])
    store.complete_db("nuc", "db_b")

    return AccessionDBCLI(install_home=str(tmp_path))


def test_search_returns_rows(cli):
    rows = cli.search("n1")
    assert rows == [("nuc", "n1", "n1", 21, "n.fa.gz")]

    rows = cli.search(accession="n1", dbs="db_a")
    assert len(rows) == 1
    assert cli.search(accession="n1", dbs="db_b") == []
    assert cli.search("missing") == []


def test_search_dbs_filter(cli):
    rows = cli.search(dbs="db_a")
    assert {row[2] for row in rows} == {"n1", "n2"}
    assert all(row[1] in {"n1", "n2"} for row in rows)


def test_summary_and_coverage(cli):
    summary = cli.summary()
    assert summary["protein_total"] == 2
    assert summary["nuc_total"] == 3
    assert summary["nuc_with_taxid"] == 3

    coverage = cli.coverage()
    assert any(row[:2] == ["protein", "prot_db"] for row in coverage)
    assert any(row[:2] == ["nuc", "db_a"] for row in coverage)
    assert cli.coverage(missing_only=True) == []


def test_export_filters_nucleotide_by_dbs(cli, tmp_path):
    output = os.path.join(str(tmp_path), "export.tsv")
    count = cli.export(dbs="db_b", output=output)
    assert count == 1
    with open(output) as handle:
        assert handle.read().splitlines() == ["accession\ttaxid", "n3\t23"]


def test_export_requires_output(cli):
    with pytest.raises(ValueError):
        cli.export(dbs="db_a")
