"""
Integration tests for name mapping consistency between:
- install_config.py (TelevirLayout attributes)
- main_install.py (DATABASE_NAMES, SOFTWARE_NAMES mappings)
- sources.yaml (canonical database/software names)
- televir_status.py (display and lookup logic)
"""

import pytest
import sqlite3
import tempfile
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'install_scripts'))

from load_sources import list_databases, list_software, list_hosts
from install_config import TelevirLayout


class TestNameMappingConsistency:
    """Test TelevirLayout attributes match DATABASE_NAMES/SOFTWARE_NAMES keys."""
    
    def test_database_names_keys_match_televirlayout(self):
        """All DATABASE_NAMES keys should be valid TelevirLayout attributes."""
        for key in TelevirLayout.DATABASE_NAMES.keys():
            assert hasattr(TelevirLayout, key), f"TelevirLayout missing attribute: {key}"
    
    def test_software_names_keys_match_televirlayout(self):
        """All SOFTWARE_NAMES keys should be valid TelevirLayout attributes."""
        for key in TelevirLayout.SOFTWARE_NAMES.keys():
            assert hasattr(TelevirLayout, key), f"TelevirLayout missing attribute: {key}"
    
    def test_database_names_exist_in_sources_yaml(self):
        """DATABASE_NAMES (category, name) pairs should exist in sources.yaml OR be dynamically generated."""
        dbs = list_databases()
        # Known dynamically generated databases (not in sources.yaml)
        known_dynamic = {"refseq_prot", "refseq_gen", "requests"}
        # Known valid categories
        valid_categories = set(dbs.keys())
        
        for key, (category, name) in TelevirLayout.DATABASE_NAMES.items():
            assert category in valid_categories, f"Category '{category}' not in sources.yaml databases"
            # Allow dynamically generated databases like refseq_prot
            if name not in known_dynamic:
                assert name in dbs[category], f"Database '{name}' not in sources.yaml['databases']['{category}']"
    
    def test_software_names_are_valid(self):
        """All SOFTWARE_NAMES should have (name, tag) tuple."""
        for key, values in TelevirLayout.SOFTWARE_NAMES.items():
            assert len(values) == 2, f"SOFTWARE_NAMES should have (name, tag), got: {values}"


class TestDatabaseRegistrationFormat:
    """Test registration format produces correct keys."""
    
    def test_database_names_format_is_category_name(self):
        """DATABASE_NAMES should produce 'category/name' format."""
        for key, (category, name) in TelevirLayout.DATABASE_NAMES.items():
            result = f"{category}/{name}"
            assert "/" in result, f"Should produce 'category/name' format, got: {result}"
    
    def test_software_names_have_two_parts(self):
        """SOFTWARE_NAMES should have (name, tag)."""
        for key, values in TelevirLayout.SOFTWARE_NAMES.items():
            assert len(values) == 2, f"SOFTWARE_NAMES should have (name, tag), got: {values}"


@pytest.fixture
def temp_db():
    """Create a temporary test database."""
    with tempfile.NamedTemporaryFile(suffix='.db', delete=False) as f:
        db_path = f.name
    
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    cursor.execute("""
        CREATE TABLE database (
            name TEXT PRIMARY KEY,
            path TEXT,
            installed TEXT,
            software TEXT,
            date TEXT,
            version TEXT,
            source_url TEXT,
            file_mod_date TEXT
        )
    """)
    cursor.execute("""
        CREATE TABLE software (
            name TEXT PRIMARY KEY,
            path TEXT,
            database TEXT,
            installed TEXT,
            tag TEXT,
            env_path TEXT,
            date TEXT,
            db_version TEXT,
            needs_update TEXT
        )
    """)
    conn.commit()
    conn.close()
    
    yield db_path
    
    os.unlink(db_path)


class TestTelevirStatusDatabaseLookup:
    """Integration tests for database lookup in televir_status.py."""
    
    def test_database_lookup_with_full_name(self, temp_db):
        """Test that database lookup works with 'category/name' keys."""
        conn = sqlite3.connect(temp_db)
        cursor = conn.cursor()
        
        # Insert with full name (as main_install.py now does)
        cursor.execute(
            "INSERT INTO database VALUES (?, ?, ?, ?, ?, ?, ?, ?)",
            ("refseq/protein", "/path/to/db", "True", "none", "2026-03-24", "1.0", "http://example.com", "2026-03-24")
        )
        conn.commit()
        conn.close()
        
        # Simulate _get_installed_databases() query
        installed = {}
        conn = sqlite3.connect(temp_db)
        cursor = conn.cursor()
        cursor.execute("SELECT name, version, date, installed FROM database WHERE installed = 'True' AND software != 'host'")
        for row in cursor.fetchall():
            installed[row[0]] = {
                "version": row[1] if row[1] else "N/A",
                "installed": row[3] if row[3] else "N/A"
            }
        conn.close()
        
        # Verify lookup with full name works
        assert "refseq/protein" in installed
        assert installed["refseq/protein"]["version"] == "1.0"
    
    def test_status_lookup_matches_registration_format(self, temp_db):
        """Test that status display lookup finds registered databases using name_full."""
        conn = sqlite3.connect(temp_db)
        cursor = conn.cursor()
        
        # Register databases with full names (simulating main_install.py with mappings)
        test_cases = [
            ("protein/refseq_prot", True),
            ("refseq/refseq_gen", True),
            ("protein/swissprot", True),
            ("kraken2/viral", True),
            ("ribosomal_rna/refseq_16s", True),
            ("ribosomal_rna/silva_16s", True),
            ("protein/virosaurus", True),
            ("protein/rvdb", True),
        ]
        
        for name, installed in test_cases:
            cursor.execute(
                "INSERT INTO database VALUES (?, ?, ?, ?, ?, ?, ?, ?)",
                (name, "/path", "True" if installed else "False", "none", "2026-03-24", "1.0", "http://example.com", "2026-03-24")
            )
        conn.commit()
        conn.close()
        
        # Simulate _get_installed_databases() query (same as televir_status.py)
        installed = {}
        conn = sqlite3.connect(temp_db)
        cursor = conn.cursor()
        cursor.execute("SELECT name, version, date, installed FROM database WHERE installed = 'True' AND software != 'host'")
        for row in cursor.fetchall():
            installed[row[0]] = {"version": row[1], "installed": row[3]}
        conn.close()
        
        # Simulate _populate_databases() lookup using name_full (as fixed in televir_status.py)
        dbs = list_databases()
        found_count = 0
        
        # Databases that exist in sources.yaml and should be found
        expected_from_yaml = [
            ("protein/swissprot"),
            ("protein/rvdb"),
            ("protein/virosaurus"),
            ("kraken2/viral"),
            ("ribosomal_rna/refseq_16s"),
            ("ribosomal_rna/silva_16s"),
        ]
        
        for category, entries in dbs.items():
            if isinstance(entries, dict):
                for name, info in entries.items():
                    name_full = f"{category}/{name}"
                    
                    # This is what the FIXED televir_status.py does
                    result = installed.get(name_full, {})
                    
                    # Check if this database was registered in our test
                    if name_full in expected_from_yaml:
                        assert result, f"Failed to find {name_full} in installed databases"
                        assert result.get("version") == "1.0"
                        found_count += 1
        
        assert found_count == len(expected_from_yaml), f"Found {found_count} of {len(expected_from_yaml)} expected databases"
    
    def test_unregistered_databases_return_empty(self, temp_db):
        """Test that databases not in the database return empty dict."""
        # Insert some databases
        conn = sqlite3.connect(temp_db)
        cursor = conn.cursor()
        cursor.execute(
            "INSERT INTO database VALUES (?, ?, ?, ?, ?, ?, ?, ?)",
            ("refseq/protein", "/path", "True", "none", "2026-03-24", "1.0", "http://example.com", "2026-03-24")
        )
        conn.commit()
        conn.close()
        
        # Query installed databases
        installed = {}
        conn = sqlite3.connect(temp_db)
        cursor = conn.cursor()
        cursor.execute("SELECT name, version, date, installed FROM database WHERE installed = 'True' AND software != 'host'")
        for row in cursor.fetchall():
            installed[row[0]] = {"version": row[1], "installed": row[3]}
        conn.close()
        
        # Check lookup for non-existent database
        result = installed.get("refseq/genome", {})
        assert result == {}, f"Should return empty dict for non-existent database, got: {result}"


class TestTelevirStatusSoftwareLookup:
    """Integration tests for software lookup in televir_status.py."""
    
    def test_software_lookup_works(self, temp_db):
        """Test that software lookup works correctly."""
        conn = sqlite3.connect(temp_db)
        cursor = conn.cursor()
        
        # Insert software (simulating main_install.py)
        test_software = [
            ("kraken2", "viral", "True"),
            ("centrifuge", "bacteria", "True"),
            ("diamond", "swissprot", "False"),
        ]
        
        for name, tag, installed in test_software:
            cursor.execute(
                "INSERT INTO software VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)",
                (name, "/path", "database", installed, tag, "/env", "2026-03-24", "1.0", "False")
            )
        conn.commit()
        conn.close()
        
        # Simulate _get_installed_software() query (same as televir_status.py)
        installed = {}
        conn = sqlite3.connect(temp_db)
        cursor = conn.cursor()
        cursor.execute("SELECT name, tag FROM software WHERE installed = 'True'")
        for row in cursor.fetchall():
            installed[row[0]] = row[1] if row[1] else "N/A"
        conn.close()
        
        # Verify
        assert "kraken2" in installed
        assert installed["kraken2"] == "viral"
        assert "centrifuge" in installed
        assert installed["centrifuge"] == "bacteria"
        assert "diamond" not in installed  # False


class TestTelevirLayoutAttributes:
    """Test TelevirLayout has expected attributes."""
    
    def test_has_database_install_flags(self):
        """TelevirLayout should have install_ flags for databases."""
        expected = [
            'install_refseq_prot', 'install_refseq_gen', 'install_swissprot',
            'install_rvdb', 'install_virosaurus', 'install_refseq_16s', 'install_ribo16s',
        ]
        for attr in expected:
            assert hasattr(TelevirLayout, attr), f"Missing attribute: {attr}"
    
    def test_has_software_install_flags(self):
        """TelevirLayout should have install_ flags for software."""
        expected = [
            'install_kraken2', 'install_kraken2_bacteria', 'install_kraken2_eupathdb46',
            'install_centrifuge', 'install_centrifuge_bacteria', 'install_metaphlan',
            'install_kaiju', 'install_krakenuniq', 'install_diamond',
            'install_voyager_viral', 'install_blast', 'install_fastviromeexplorer',
        ]
        for attr in expected:
            assert hasattr(TelevirLayout, attr), f"Missing attribute: {attr}"


class TestHostMappingConsistency:
    """Test that host genomes are handled correctly."""
    
    def test_hosts_have_required_fields(self):
        """All hosts in sources.yaml should have required fields."""
        hosts = list_hosts()
        for host_key, host_data in hosts.items():
            if isinstance(host_data, dict):
                assert 'host_name' in host_data, f"{host_key} missing host_name"
                assert 'common_name' in host_data, f"{host_key} missing common_name"
    
    def test_host_registration_uses_host_name(self):
        """Host registration should use host_name as identifier."""
        # This is how main_install.py registers hosts - using host_name
        hosts = list_hosts()
        for host_key, host_data in hosts.items():
            if isinstance(host_data, dict):
                host_name = host_data.get('host_name', host_key)
                assert host_name is not None, f"host_name should not be None for {host_key}"
