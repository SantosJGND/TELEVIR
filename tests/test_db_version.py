"""
Integration tests for db_install.py using sources.yaml.
Tests that database URLs and versions are correctly loaded from sources.yaml.
"""

import pytest
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'install_scripts'))

from load_sources import (
    get_db_url,
    get_db_version,
    get_db_entry,
    list_databases,
    extract_version_from_string,
)


class TestSourcesYamlKeysMatchDbInstall:
    """Test sources.yaml uses same keys as db_install.py expects."""

    def test_kraken2_has_required_keys(self):
        """Verify kraken2 has all keys db_install.py needs."""
        dbs = list_databases('kraken2')
        required_keys = ['viral', 'standard', 'bacteria', 'ribo16s', 'eupathdb46']
        for key in required_keys:
            assert key in dbs, f"Missing key: {key}"

    def test_kaiju_has_required_keys(self):
        """Verify kaiju has all keys db_install.py needs."""
        dbs = list_databases('kaiju')
        required_keys = ['viral', 'fungi', 'bacteria']
        for key in required_keys:
            assert key in dbs, f"Missing key: {key}"

    def test_centrifuge_has_required_keys(self):
        """Verify centrifuge has all keys db_install.py needs."""
        dbs = list_databases('centrifuge')
        required_keys = ['viral', 'bacteria']
        for key in required_keys:
            assert key in dbs, f"Missing key: {key}"


class TestDbEntryFunction:
    """Test get_db_entry returns complete database info."""

    def test_get_db_entry_returns_all_fields(self):
        """Verify get_db_entry returns url, version, file, description."""
        entry = get_db_entry('kraken2', 'viral')
        assert entry is not None
        assert 'url' in entry
        assert 'version' in entry
        assert 'file' in entry
        assert 'description' in entry

    def test_get_db_entry_returns_none_for_missing(self):
        """Verify get_db_entry returns None for missing entry."""
        entry = get_db_entry('nonexistent', 'viral')
        assert entry is None


class TestVersionExtraction:
    """Test version extraction returns valid values."""

    def test_kraken2_versions_are_not_none(self):
        """Verify Kraken2 versions are retrieved."""
        assert get_db_version('kraken2', 'viral') is not None
        assert get_db_version('kraken2', 'standard') is not None
        assert get_db_version('kraken2', 'bacteria') is not None
        assert get_db_version('kraken2', 'ribo16s') is not None
        assert get_db_version('kraken2', 'eupathdb46') is not None

    def test_kaiju_versions_are_not_none(self):
        """Verify Kaiju versions are retrieved."""
        assert get_db_version('kaiju', 'viral') is not None
        assert get_db_version('kaiju', 'fungi') is not None
        assert get_db_version('kaiju', 'bacteria') is not None

    def test_centrifuge_versions_are_not_none(self):
        """Verify Centrifuge versions are retrieved."""
        assert get_db_version('centrifuge', 'viral') is not None
        assert get_db_version('centrifuge', 'bacteria') is not None


class TestExtractVersionFromString:
    """Unit tests for version extraction from strings."""

    def test_extract_YYYYMMDD_from_filename(self):
        """Test extraction from filename like k2_viral_20250402.tar.gz"""
        result = extract_version_from_string('k2_viral_20250402.tar.gz')
        assert result == '20250402'

    def test_extract_YYYYMMDD_from_url(self):
        """Test extraction from URL."""
        result = extract_version_from_string(
            'https://example.com/kraken/k2_viral_20250402.tar.gz'
        )
        assert result == '20250402'

    def test_extract_YYYY_MM_DD(self):
        """Test extraction from date format YYYY-MM-DD."""
        result = extract_version_from_string(
            'kaiju_db_viruses_2024-08-15.tgz'
        )
        assert result == '2024-08-15'

    def test_case_insensitive_extraction(self):
        """Test case-insensitive extraction."""
        result = extract_version_from_string('K2_VIRAL_20250402.TAR.GZ')
        assert result == '20250402'

    def test_extract_vJan25_format(self):
        """Test extraction from vJan25 format."""
        result = extract_version_from_string(
            'mpa_vJan25_CHOCOPhlAnSGB_202503.tar'
        )
        assert result == '202503'

    def test_extract_returns_none_for_no_match(self):
        """Test returns None when no pattern matches."""
        result = extract_version_from_string('some_file.tar.gz')
        assert result is None

    def test_extract_from_underscore_date_format(self):
        """Test extraction from _YYYYMMDD. pattern."""
        result = extract_version_from_string('k2_standard_16gb_20230314.tar.gz')
        assert result == '20230314'
