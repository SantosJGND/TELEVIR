#!/usr/bin/env python3

import csv
import hashlib
import json
import os
import sqlite3
from contextlib import closing

PROTEIN_DB = "protein_accessions.db"
NUCLEOTIDE_DB = "nucleotide_accessions.db"

TABLES = {
    "protein": "protein_accessions",
    "nuc": "nucleotide_accessions",
}

STAGE_TABLES = {
    "protein": "protein_accessions_stage",
    "nuc": "nucleotide_accessions_stage",
}

COLUMNS = {
    "protein": ("dbs", "acc", "description", "acc_in_file", "taxid", "file"),
    "nuc": ("dbs", "acc", "acc_in_file", "taxid", "file"),
}

TABLE_SQL = {
    "protein": """
        CREATE TABLE IF NOT EXISTS protein_accessions (
            dbs TEXT,
            acc TEXT,
            description TEXT,
            acc_in_file TEXT,
            taxid INTEGER,
            file TEXT,
            PRIMARY KEY (dbs, acc)
        )
    """,
    "nuc": """
        CREATE TABLE IF NOT EXISTS nucleotide_accessions (
            dbs TEXT,
            acc TEXT,
            acc_in_file TEXT,
            taxid INTEGER,
            file TEXT,
            PRIMARY KEY (dbs, acc)
        )
    """,
}

INDEX_SQL = {
    "protein": (
        "CREATE INDEX IF NOT EXISTS idx_protein_acc ON protein_accessions(acc)",
        "CREATE INDEX IF NOT EXISTS idx_protein_dbs ON protein_accessions(dbs)",
        "CREATE INDEX IF NOT EXISTS idx_protein_desc ON protein_accessions(description)",
    ),
    "nuc": (
        "CREATE INDEX IF NOT EXISTS idx_nuc_acc ON nucleotide_accessions(acc)",
        "CREATE INDEX IF NOT EXISTS idx_nuc_dbs ON nucleotide_accessions(dbs)",
    ),
}

METADATA_SQL = """
    CREATE TABLE IF NOT EXISTS accession_metadata (
        db_key TEXT PRIMARY KEY,
        seq_type TEXT NOT NULL,
        fingerprint TEXT,
        version TEXT,
        source_url TEXT,
        state TEXT NOT NULL,
        row_count INTEGER DEFAULT 0,
        taxid_count INTEGER DEFAULT 0,
        created_at TEXT,
        updated_at TEXT
    )
"""

SOURCE_SQL = """
    CREATE TABLE IF NOT EXISTS accession_source_file (
        db_key TEXT,
        path TEXT,
        size INTEGER,
        mtime_ns INTEGER,
        PRIMARY KEY (db_key, path)
    )
"""

BATCH_SIZE = 100000


def _now():
    from datetime import datetime

    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


class AccessionStore:
    def __init__(self, metadir, batch_size=BATCH_SIZE):
        self.metadir = metadir
        self.batch_size = batch_size
        self.db_paths = {
            "protein": os.path.join(metadir, PROTEIN_DB),
            "nuc": os.path.join(metadir, NUCLEOTIDE_DB),
        }
        self.stage_sql = {
            seq_type: sql.replace(
                f"CREATE TABLE IF NOT EXISTS {TABLES[seq_type]}",
                f"CREATE TABLE {STAGE_TABLES[seq_type]}",
            )
            for seq_type, sql in TABLE_SQL.items()
        }

    def path(self, seq_type):
        return self.db_paths[seq_type]

    def stage_table(self, seq_type):
        return STAGE_TABLES[seq_type]

    def table(self, seq_type):
        return TABLES[seq_type]

    def connect(self, seq_type):
        conn = sqlite3.connect(self.path(seq_type), timeout=30)
        conn.execute("PRAGMA busy_timeout=30000")
        try:
            conn.execute("PRAGMA journal_mode=WAL")
        except sqlite3.OperationalError:
            pass
        conn.execute("PRAGMA synchronous=NORMAL")
        return conn

    def ensure_schema(self):
        os.makedirs(self.metadir, exist_ok=True)
        for seq_type in TABLES:
            with closing(self.connect(seq_type)) as conn:
                conn.execute(TABLE_SQL[seq_type])
                for statement in INDEX_SQL[seq_type]:
                    conn.execute(statement)
                conn.execute(METADATA_SQL)
                conn.execute(SOURCE_SQL)
                conn.commit()

    def fingerprint(self, files, version=None):
        entries = []
        for path in files or []:
            try:
                stat = os.stat(path)
                entries.append([str(path), stat.st_size, stat.st_mtime_ns])
            except OSError:
                entries.append([str(path), None, None])
        payload = json.dumps({"files": sorted(entries), "version": version or ""}, sort_keys=True)
        return hashlib.sha256(payload.encode("utf-8")).hexdigest()[:20]

    def register_legacy(self, versions=None):
        versions = versions or {}
        self.ensure_schema()
        for seq_type, table in TABLES.items():
            with closing(self.connect(seq_type)) as conn:
                db_keys = [row[0] for row in conn.execute(f"SELECT DISTINCT dbs FROM {table}")]
                for db_key in db_keys:
                    exists = conn.execute(
                        "SELECT 1 FROM accession_metadata WHERE db_key = ?", (db_key,)
                    ).fetchone()
                    if exists:
                        continue
                    files = [
                        row[0]
                        for row in conn.execute(
                            f"SELECT DISTINCT file FROM {table} WHERE dbs = ? AND file IS NOT NULL",
                            (db_key,),
                        )
                    ]
                    counts = conn.execute(
                        f"SELECT COUNT(*), SUM(CASE WHEN taxid IS NOT NULL THEN 1 ELSE 0 END) FROM {table} WHERE dbs = ?",
                        (db_key,),
                    ).fetchone()
                    version = versions.get(db_key, {}).get("version")
                    source_url = versions.get(db_key, {}).get("source_url")
                    conn.execute(
                        """
                        INSERT INTO accession_metadata
                            (db_key, seq_type, fingerprint, version, source_url, state,
                             row_count, taxid_count, created_at, updated_at)
                        VALUES (?, ?, ?, ?, ?, 'imported', ?, ?, ?, ?)
                        """,
                        (
                            db_key,
                            seq_type,
                            self.fingerprint(files, version),
                            version,
                            source_url,
                            counts[0] or 0,
                            counts[1] or 0,
                            _now(),
                            _now(),
                        ),
                    )
                conn.commit()

    def metadata(self):
        self.ensure_schema()
        rows = []
        for seq_type in TABLES:
            if not os.path.exists(self.path(seq_type)):
                continue
            with closing(self.connect(seq_type)) as conn:
                for row in conn.execute(
                    """
                    SELECT db_key, seq_type, fingerprint, version, source_url, state,
                           row_count, taxid_count, updated_at
                    FROM accession_metadata ORDER BY db_key
                    """
                ):
                    rows.append(dict(zip(
                        ("db_key", "seq_type", "fingerprint", "version", "source_url",
                         "state", "row_count", "taxid_count", "updated_at"),
                        row,
                    )))
        return rows

    def needs_rebuild(self, seq_type, db_key, files, version=None):
        fingerprint = self.fingerprint(files, version)
        if not os.path.exists(self.path(seq_type)):
            return True
        with closing(self.connect(seq_type)) as conn:
            row = conn.execute(
                "SELECT state, fingerprint FROM accession_metadata WHERE db_key = ?",
                (db_key,),
            ).fetchone()
        if row is None:
            return True
        return row[0] != "complete" or row[1] != fingerprint

    def begin_db(self, seq_type, db_key, files, version=None, source_url=None):
        self.ensure_schema()
        fingerprint = self.fingerprint(files, version)
        stage = STAGE_TABLES[seq_type]
        with closing(self.connect(seq_type)) as conn:
            conn.execute("BEGIN IMMEDIATE")
            conn.execute(f"DROP TABLE IF EXISTS {stage}")
            conn.execute(self.stage_sql[seq_type])
            conn.execute("DELETE FROM accession_source_file WHERE db_key = ?", (db_key,))
            conn.executemany(
                "INSERT INTO accession_source_file (db_key, path, size, mtime_ns) VALUES (?, ?, ?, ?)",
                [
                    (db_key, str(path), *self._file_stats(path))
                    for path in files or []
                ],
            )
            conn.execute("DELETE FROM accession_metadata WHERE db_key = ?", (db_key,))
            conn.execute(
                """
                INSERT INTO accession_metadata
                    (db_key, seq_type, fingerprint, version, source_url, state,
                     row_count, taxid_count, created_at, updated_at)
                VALUES (?, ?, ?, ?, ?, 'staging', 0, 0, ?, ?)
                """,
                (db_key, seq_type, fingerprint, version, source_url, _now(), _now()),
            )
            conn.commit()
        return fingerprint

    def insert_batch(self, seq_type, records):
        records = list(records)
        if not records:
            return
        stage = STAGE_TABLES[seq_type]
        columns = ", ".join(COLUMNS[seq_type])
        placeholders = ", ".join("?" for _ in COLUMNS[seq_type])
        with closing(self.connect(seq_type)) as conn:
            for start in range(0, len(records), self.batch_size):
                conn.executemany(
                    f"INSERT OR IGNORE INTO {stage} ({columns}) VALUES ({placeholders})",
                    records[start:start + self.batch_size],
                )
                conn.commit()

    def stage_accessions(self, seq_type, db_key):
        stage = STAGE_TABLES[seq_type]
        with closing(self.connect(seq_type)) as conn:
            return {
                row[0]
                for row in conn.execute(
                    f"SELECT acc FROM {stage} WHERE dbs = ?", (db_key,)
                )
            }

    def stage_keys(self, seq_type):
        stage = STAGE_TABLES[seq_type]
        with closing(self.connect(seq_type)) as conn:
            return [row[0] for row in conn.execute(f"SELECT DISTINCT dbs FROM {stage}")]

    def update_taxids(self, seq_type, matches):
        matches = list(matches)
        if not matches:
            return
        stage = STAGE_TABLES[seq_type]
        with closing(self.connect(seq_type)) as conn:
            conn.executemany(
                f"UPDATE {stage} SET taxid = ? WHERE dbs = ? AND acc = ?",
                [(taxid, db_key, acc) for db_key, acc, taxid in matches],
            )
            conn.commit()

    def complete_db(self, seq_type, db_key):
        table = TABLES[seq_type]
        stage = STAGE_TABLES[seq_type]
        columns = ", ".join(COLUMNS[seq_type])
        with closing(self.connect(seq_type)) as conn:
            conn.execute("BEGIN IMMEDIATE")
            counts = conn.execute(
                f"SELECT COUNT(*), SUM(CASE WHEN taxid IS NOT NULL THEN 1 ELSE 0 END) FROM {stage} WHERE dbs = ?",
                (db_key,),
            ).fetchone()
            conn.execute(f"DELETE FROM {table} WHERE dbs = ?", (db_key,))
            conn.execute(
                f"INSERT INTO {table} ({columns}) SELECT {columns} FROM {stage} WHERE dbs = ?",
                (db_key,),
            )
            conn.execute(f"DROP TABLE {stage}")
            conn.execute(
                """
                UPDATE accession_metadata
                SET state = 'complete', row_count = ?, taxid_count = ?, updated_at = ?
                WHERE db_key = ?
                """,
                (counts[0] or 0, counts[1] or 0, _now(), db_key),
            )
            conn.commit()

    def abort_db(self, seq_type):
        stage = STAGE_TABLES[seq_type]
        with closing(self.connect(seq_type)) as conn:
            conn.execute(f"DROP TABLE IF EXISTS {stage}")
            conn.commit()

    def create_scratch(self, path):
        directory = os.path.dirname(os.path.abspath(path))
        os.makedirs(directory, exist_ok=True)
        conn = sqlite3.connect(path, timeout=30)
        conn.execute("DROP TABLE IF EXISTS matches")
        conn.execute("CREATE TABLE matches (dbs TEXT, acc TEXT, taxid INTEGER)")
        return conn

    def scratch_has_matches(self, path):
        if not os.path.exists(path):
            return False
        try:
            with closing(sqlite3.connect(path)) as conn:
                has_table = conn.execute(
                    "SELECT 1 FROM sqlite_master WHERE type='table' AND name='matches'"
                ).fetchone()
                if not has_table:
                    return False
                return conn.execute("SELECT COUNT(*) FROM matches").fetchone()[0] > 0
        except sqlite3.DatabaseError:
            return False

    def read_scratch_matches(self, path, db_key):
        with closing(sqlite3.connect(path)) as conn:
            return conn.execute(
                "SELECT acc, taxid FROM matches WHERE dbs = ?", (db_key,)
            ).fetchall()

    def remove_scratch(self, path):
        if os.path.exists(path):
            os.remove(path)

    def export_tsv(self, seq_type, outfile, force=False):
        if not force and os.path.exists(outfile):
            return False
        if not os.path.exists(self.path(seq_type)):
            return False
        table = TABLES[seq_type]
        if seq_type == "protein":
            header = ("acc", "taxid")
            query = f"SELECT acc, taxid FROM {table} WHERE taxid IS NOT NULL"
        else:
            header = ("acc", "taxid", "file", "acc_in_file")
            query = f"SELECT acc, taxid, file, acc_in_file FROM {table}"
        with closing(self.connect(seq_type)) as conn:
            rows = conn.execute(query).fetchall()
        os.makedirs(os.path.dirname(os.path.abspath(outfile)), exist_ok=True)
        with open(outfile, "w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(header)
            for row in rows:
                writer.writerow(["" if value is None else str(value) for value in row])
        return True

    def counts(self):
        result = {
            "protein_total": 0,
            "protein_with_taxid": 0,
            "protein_size_mb": 0.0,
            "nuc_total": 0,
            "nuc_with_taxid": 0,
            "nuc_size_mb": 0.0,
            "metadata": self.metadata(),
        }
        for seq_type in TABLES:
            if not os.path.exists(self.path(seq_type)):
                continue
            with closing(self.connect(seq_type)) as conn:
                counts = conn.execute(
                    f"SELECT COUNT(*), SUM(CASE WHEN taxid IS NOT NULL THEN 1 ELSE 0 END) FROM {TABLES[seq_type]}"
                ).fetchone()
            result[f"{seq_type}_total"] = counts[0] or 0
            result[f"{seq_type}_with_taxid"] = counts[1] or 0
            result[f"{seq_type}_size_mb"] = os.path.getsize(self.path(seq_type)) / (1024 * 1024)
        return result

    def list_dbs(self, db_type="all"):
        results = []
        for seq_type, table in TABLES.items():
            if db_type not in ("all", seq_type):
                continue
            if not os.path.exists(self.path(seq_type)):
                continue
            with closing(self.connect(seq_type)) as conn:
                rows = conn.execute(
                    f"""
                    SELECT dbs, COUNT(*) AS total,
                           SUM(CASE WHEN taxid IS NOT NULL THEN 1 ELSE 0 END) AS with_taxid
                    FROM {table} GROUP BY dbs ORDER BY dbs
                    """
                ).fetchall()
            label = "protein" if seq_type == "protein" else "nuc"
            for db_key, total, with_taxid in rows:
                pct = (with_taxid / total * 100) if total else 0
                results.append([label, db_key, total, with_taxid, f"{pct:.1f}%"])
        return results

    def coverage(self, missing_only=False):
        if not missing_only:
            return self.list_dbs("all")
        results = []
        for seq_type, table in TABLES.items():
            if not os.path.exists(self.path(seq_type)):
                continue
            with closing(self.connect(seq_type)) as conn:
                rows = conn.execute(
                    f"""
                    SELECT dbs, COUNT(*) AS missing FROM {table}
                    WHERE taxid IS NULL GROUP BY dbs ORDER BY missing DESC
                    """
                ).fetchall()
            label = "protein" if seq_type == "protein" else "nuc"
            for db_key, missing in rows:
                results.append([label, db_key, missing, 0, 0])
        return results

    def search(self, accession=None, dbs=None, limit=50):
        results = []
        clauses = []
        params = []
        if accession:
            clauses.append("acc LIKE ?")
            params.append(f"%{accession}%")
        if dbs:
            clauses.append("dbs = ?")
            params.append(dbs)
        where = f" WHERE {' AND '.join(clauses)}" if clauses else ""
        for seq_type, table in TABLES.items():
            if not os.path.exists(self.path(seq_type)):
                continue
            if seq_type == "protein":
                columns = "dbs, acc, description, taxid, file"
            else:
                columns = "dbs, acc, acc_in_file, taxid, file"
            with closing(self.connect(seq_type)) as conn:
                rows = conn.execute(
                    f"SELECT {columns} FROM {table}{where} ORDER BY acc LIMIT ?",
                    (*params, limit),
                ).fetchall()
            label = "protein" if seq_type == "protein" else "nuc"
            for row in rows:
                results.append((label, *row[1:]))
        return results

    def export_records(self, seq_type, dbs=None):
        if not os.path.exists(self.path(seq_type)):
            return []
        query = f"SELECT acc, taxid FROM {TABLES[seq_type]} WHERE taxid IS NOT NULL"
        params = ()
        if dbs:
            query += " AND dbs = ?"
            params = (dbs,)
        with closing(self.connect(seq_type)) as conn:
            return conn.execute(query, params).fetchall()

    @staticmethod
    def _file_stats(path):
        try:
            stat = os.stat(path)
            return stat.st_size, stat.st_mtime_ns
        except OSError:
            return None, None
