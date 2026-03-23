"""
Unit tests for load_sources.py - SourceLoader and accessor functions.
"""

import pytest
from load_sources import (
    get_loader, get_db_url, get_host_config, get_software_url,
    get_git_url, get_ncbi_eutils_url, list_databases, list_hosts,
    list_software, reload, get, sources, version, last_updated
)


class TestSourceLoader:
    """Tests for SourceLoader class."""

    def test_loader_initialization(self, loader):
        """Test SourceLoader initializes successfully."""
        assert loader is not None
        assert loader.sources is not None

    def test_sources_is_dict(self, loader):
        """Test that sources returns a dictionary."""
        assert isinstance(loader.sources, dict)

    def test_version_info(self, loader):
        """Test version information is available."""
        assert loader.version is not None
        assert isinstance(loader.version, str)
        assert len(loader.version) > 0

    def test_last_updated_info(self, loader):
        """Test last_updated information is available."""
        assert loader.last_updated is not None
        assert isinstance(loader.last_updated, str)

    def test_get_nested_value(self, loader):
        """Test nested value retrieval with get()."""
        value = loader.get('databases', 'kraken2', 'viral', 'url')
        assert value is not None
        assert isinstance(value, str)
        assert 'kraken' in value.lower()

    def test_get_missing_key_returns_default(self, loader):
        """Test get() returns default for missing keys."""
        result = loader.get('nonexistent', 'key', default='default_value')
        assert result == 'default_value'

    def test_singleton_pattern(self, clean_loader):
        """Test that get_loader returns the same instance."""
        loader1 = get_loader()
        loader2 = get_loader()
        assert loader1 is loader2


class TestDatabaseURLs:
    """Tests for database URL retrieval."""

    def test_get_db_url_kraken2_viral(self):
        """Test getting Kraken2 viral database URL."""
        url = get_db_url('kraken2', 'viral')
        assert url is not None
        assert 'kraken' in url.lower()
        assert 'viral' in url.lower()
        assert url.startswith('http')

    def test_get_db_url_kraken2_standard(self):
        """Test getting Kraken2 standard database URL."""
        url = get_db_url('kraken2', 'standard')
        assert url is not None
        assert 'kraken' in url.lower()
        assert 'standard' in url.lower()

    def test_get_db_url_silva_16s(self):
        """Test getting SILVA 16S database URL."""
        url = get_db_url('ribosomal_rna', 'silva_16s')
        assert url is not None
        assert 'silva' in url.lower()

    def test_get_db_url_uniref(self):
        """Test getting UniRef database URL."""
        url = get_db_url('protein', 'uniref90')
        assert url is not None
        assert 'uniprot' in url.lower() or 'uniref' in url.lower()

    def test_get_db_url_missing_category(self):
        """Test get_db_url returns None for missing category."""
        url = get_db_url('nonexistent_category', 'viral')
        assert url is None

    def test_get_db_url_missing_name(self):
        """Test get_db_url returns None for missing database name."""
        url = get_db_url('kraken2', 'nonexistent_db')
        assert url is None


class TestHostConfigs:
    """Tests for host genome configuration retrieval."""

    def test_get_host_config_homo_sapiens(self):
        """Test getting human host configuration."""
        config = get_host_config('homo_sapiens')
        assert config is not None
        assert isinstance(config, dict)
        assert config['host_name'] == 'hg38'
        assert config['common_name'] == 'human'
        assert 'ftp' in config['host']
        assert 'genomes' in config['path']
        assert config['file'].endswith('.fna.gz')

    def test_get_host_config_by_host_name(self):
        """Test getting host config by host_name (hg38)."""
        config = get_host_config('hg38')
        assert config is not None
        assert config['host_name'] == 'hg38'
        assert config['common_name'] == 'human'

    def test_get_host_config_by_common_name(self):
        """Test getting host config by common_name."""
        config = get_host_config('human')
        assert config is not None
        assert config['common_name'] == 'human'

    def test_get_host_config_sus_scrofa(self):
        """Test getting pig host configuration."""
        config = get_host_config('sus_scrofa')
        assert config is not None
        assert config['host_name'] == 'sus_scrofa'
        assert config['common_name'] == 'pig'

    def test_get_host_config_misspelled(self):
        """Test get_host_config returns None for invalid key."""
        config = get_host_config('invalid_host_name_xyz')
        assert config is None

    def test_host_config_has_required_fields(self):
        """Test that all host configs have required fields."""
        hosts = list_hosts()
        for host_key, host_data in hosts.items():
            if isinstance(host_data, dict) and 'host_name' in host_data:
                assert 'host_name' in host_data
                assert 'common_name' in host_data
                assert 'host' in host_data
                assert 'path' in host_data
                assert 'file' in host_data


class TestSoftwareURLs:
    """Tests for software URL retrieval."""

    def test_get_software_url_fastqc(self):
        """Test getting FastQC URL."""
        url = get_software_url('fastqc')
        assert url is not None
        assert 'fastqc' in url.lower()

    def test_get_software_url_trimmomatic(self):
        """Test getting Trimmomatic URL."""
        url = get_software_url('trimmomatic')
        assert url is not None
        assert 'trimmomatic' in url.lower()

    def test_get_software_url_voyager(self):
        """Test getting Voyager URL."""
        url = get_software_url('voyager')
        assert url is not None
        assert 'voyager' in url.lower()

    def test_get_software_url_missing(self):
        """Test get_software_url returns None for missing software."""
        url = get_software_url('nonexistent_software_xyz')
        assert url is None

    def test_get_git_url_fastviromeexplorer(self):
        """Test getting FastViromeExplorer Git URL."""
        url = get_git_url('fastviromeexplorer')
        assert url is not None
        assert 'github' in url.lower()
        assert 'fastviromeexplorer' in url.lower()

    def test_get_git_url_desamba(self):
        """Test getting deSAMBA Git URL."""
        url = get_git_url('desamba')
        assert url is not None
        assert 'github' in url.lower()
        assert 'desamba' in url.lower() or 'hitbc' in url.lower()


class TestNCBIEutils:
    """Tests for NCBI E-utilities URLs."""

    def test_get_ncbi_eutils_esearch(self):
        """Test getting NCBI ESearch URL."""
        url = get_ncbi_eutils_url('esearch')
        assert url is not None
        assert 'eutils.ncbi.nlm.nih.gov' in url
        assert 'esearch' in url

    def test_get_ncbi_eutils_efetch(self):
        """Test getting NCBI EFetch URL."""
        url = get_ncbi_eutils_url('efetch')
        assert url is not None
        assert 'eutils.ncbi.nlm.nih.gov' in url
        assert 'efetch' in url

    def test_get_ncbi_eutils_elink(self):
        """Test getting NCBI ELink URL."""
        url = get_ncbi_eutils_url('elink')
        assert url is not None
        assert 'eutils.ncbi.nlm.nih.gov' in url
        assert 'elink' in url


class TestListFunctions:
    """Tests for list_* functions."""

    def test_list_databases(self):
        """Test list_databases returns all databases."""
        dbs = list_databases()
        assert isinstance(dbs, dict)
        assert len(dbs) > 0
        assert 'kraken2' in dbs
        assert 'kaiju' in dbs
        assert 'protein' in dbs

    def test_list_databases_category(self):
        """Test list_databases with specific category."""
        kraken2_dbs = list_databases('kraken2')
        assert isinstance(kraken2_dbs, dict)
        assert 'viral' in kraken2_dbs
        assert 'standard' in kraken2_dbs

    def test_list_hosts(self):
        """Test list_hosts returns all hosts."""
        hosts = list_hosts()
        assert isinstance(hosts, dict)
        assert len(hosts) > 10
        assert 'homo_sapiens' in hosts
        assert 'sus_scrofa' in hosts

    def test_list_software(self):
        """Test list_software returns all software."""
        software = list_software()
        assert isinstance(software, dict)
        assert len(software) > 0
        assert 'archives' in software
        assert 'git_repos' in software


class TestModuleLevelAccessors:
    """Tests for module-level convenience accessors."""

    def test_sources_function(self):
        """Test sources function returns dict."""
        from load_sources import sources as sources_func
        s = sources_func()
        assert isinstance(s, dict)
        assert 'databases' in s

    def test_version_function(self):
        """Test version function returns string."""
        from load_sources import version as ver_func
        v = ver_func()
        assert isinstance(v, str)

    def test_last_updated_function(self):
        """Test last_updated function returns string."""
        from load_sources import last_updated as lu_func
        l = lu_func()
        assert isinstance(l, str)

    def test_get_function(self):
        """Test get() function works."""
        result = get('databases')
        assert isinstance(result, dict)


class TestReload:
    """Tests for reload functionality."""

    def test_reload_returns_loader(self):
        """Test reload returns a SourceLoader instance."""
        loader = reload()
        assert loader is not None
        assert isinstance(loader, type(get_loader()))


class TestURLFormatValidation:
    """Tests to validate URL formats in sources."""

    def test_all_db_urls_are_strings(self):
        """Test all database URLs are valid string URLs."""
        dbs = list_databases()
        for category, entries in dbs.items():
            if isinstance(entries, dict):
                for name, info in entries.items():
                    if isinstance(info, dict) and 'url' in info:
                        url = info['url']
                        assert isinstance(url, str), f"URL for {category}/{name} is not a string"
                        assert url.startswith(('http://', 'https://', 'ftp://')), \
                            f"URL for {category}/{name} doesn't start with valid protocol"

    def test_host_ftp_paths_valid(self):
        """Test all host FTP paths are non-empty strings."""
        hosts = list_hosts()
        for host_key, host_data in hosts.items():
            if isinstance(host_data, dict) and 'path' in host_data:
                path = host_data['path']
                assert isinstance(path, str), f"Path for {host_key} is not a string"
                assert len(path) > 0, f"Path for {host_key} is empty"
                assert 'genomes' in path, f"Path for {host_key} should contain 'genomes'"

    def test_host_files_have_valid_extension(self):
        """Test all host genome files have .gz extension."""
        hosts = list_hosts()
        for host_key, host_data in hosts.items():
            if isinstance(host_data, dict) and 'file' in host_data:
                file_name = host_data['file']
                assert file_name.endswith('.gz'), f"File for {host_key} should end with .gz"
