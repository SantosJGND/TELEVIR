"""
Pytest configuration and shared fixtures for TELE-Vir source tests.
"""

import os
import sys
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'install_scripts'))

@pytest.fixture
def loader():
    """Provide a SourceLoader instance."""
    from load_sources import get_loader, reload
    
    reload()
    return get_loader()


@pytest.fixture
def sources():
    """Provide all sources dict."""
    from load_sources import sources
    return sources


@pytest.fixture
def sample_db_config():
    """Sample database configuration."""
    return {
        'url': 'https://example.com/test.tar.gz',
        'description': 'Test database',
        'file': 'test.tar.gz',
        'version': '1.0'
    }


@pytest.fixture
def sample_host_config():
    """Sample host genome configuration."""
    return {
        'host_name': 'test_organism',
        'common_name': 'test organism',
        'host': 'ftp.example.com',
        'path': '/genomes/test/',
        'file': 'test_genome.fna.gz'
    }


@pytest.fixture
def clean_loader():
    """Provide a fresh loader instance for each test."""
    from load_sources import SourceLoader
    
    SourceLoader._instance = None
    SourceLoader._sources = None
    
    yield SourceLoader()
    
    SourceLoader._instance = None
    SourceLoader._sources = None
