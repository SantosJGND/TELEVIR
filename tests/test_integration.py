"""
Integration tests for centralized sources configuration.
Tests interaction between load_sources, host_library, and install_source modules.
"""

import pytest
import re
import os


class TestIntegrationSources:
    """Integration tests across modules."""

    def test_load_sources_and_host_library_consistent(self):
        """Test that load_sources and host_library return consistent data."""
        from load_sources import get_host_config
        from host_library import HomoSapiens
        
        config = get_host_config('homo_sapiens')
        human = HomoSapiens()
        
        assert config['host_name'] == human.host_name
        assert config['common_name'] == human.common_name
        assert config['host'] == human.remote_host
        assert config['path'] == human.remote_path
        assert config['file'] == human.remote_filename

    def test_all_hosts_consistent_between_modules(self):
        """Test all hosts have consistent data across modules."""
        from load_sources import get_host_config
        from host_library import get_host
        
        host_keys = [
            'homo_sapiens', 'sus_scrofa', 'canis_lupus_familiaris',
            'felis_catus', 'bos_taurus', 'gallus_gallus'
        ]
        
        for key in host_keys:
            config = get_host_config(key)
            instance = get_host(config['host_name'])
            
            assert config['host_name'] == instance.host_name
            assert config['common_name'] == instance.common_name
            assert config['host'] == instance.remote_host

    def test_install_source_compatible_with_load_sources(self):
        """Test install_source module is compatible with load_sources."""
        from install_source import ENVS_PARAMS, get_git_url_for, get_tar_url_for
        from load_sources import get_git_url, get_software_url
        
        assert 'GIT' in ENVS_PARAMS
        assert 'TAR' in ENVS_PARAMS
        
        fve_url = get_git_url_for('FastViromeExplorer/fve')
        assert fve_url is not None
        assert 'github' in fve_url
        
        clark_url = get_tar_url_for('classification/Clark')
        assert clark_url is not None

    def test_backward_compatibility_envs_params(self):
        """Test ENVS_PARAMS maintains backward compatibility."""
        from install_source import ENVS_PARAMS
        
        assert 'ENVS' in ENVS_PARAMS
        assert 'GIT' in ENVS_PARAMS
        assert 'TAR' in ENVS_PARAMS
        assert 'ENVSDIR' in ENVS_PARAMS
        
        assert 'hostDepletion/hostdep_env' in ENVS_PARAMS['ENVS']
        assert 'kraken2/kraken_env' in ENVS_PARAMS['ENVS']

    def test_backward_compatibility_install_params(self):
        """Test INSTALL_PARAMS maintains backward compatibility."""
        from install_source import INSTALL_PARAMS
        
        assert 'HOME' in INSTALL_PARAMS
        assert 'ENVSDIR' in INSTALL_PARAMS
        assert 'REQUEST_REFERENCES' in INSTALL_PARAMS
        
        assert 'ENVSDIR' in INSTALL_PARAMS
        assert 'centrifuge' in INSTALL_PARAMS['ENVSDIR']
        assert 'kraken2' in INSTALL_PARAMS['ENVSDIR']


class TestDataCompleteness:
    """Tests for data completeness across all sources."""

    def test_all_kraken2_databases_have_urls(self):
        """Test all Kraken2 database entries have URLs."""
        from load_sources import list_databases
        
        kraken2_dbs = list_databases('kraken2')
        for db_name, db_info in kraken2_dbs.items():
            if isinstance(db_info, dict):
                assert 'url' in db_info, f"Kraken2 {db_name} missing URL"
                assert db_info['url'].startswith('http')

    def test_all_kaiju_databases_have_urls(self):
        """Test all Kaiju database entries have URLs."""
        from load_sources import list_databases
        
        kaiju_dbs = list_databases('kaiju')
        for db_name, db_info in kaiju_dbs.items():
            if isinstance(db_info, dict):
                assert 'url' in db_info, f"Kaiju {db_name} missing URL"

    def test_all_software_have_urls(self):
        """Test all software entries have URLs."""
        from load_sources import list_software
        
        archives = list_software().get('archives', {})
        for sw_name, sw_info in archives.items():
            if isinstance(sw_info, dict):
                assert 'url' in sw_info, f"Software {sw_name} missing URL"

    def test_all_host_genomes_have_complete_config(self):
        """Test all host genomes have complete FTP configuration."""
        from load_sources import list_hosts
        
        hosts = list_hosts()
        for host_key, host_data in hosts.items():
            if isinstance(host_data, dict) and 'host_name' in host_data:
                required_fields = ['host_name', 'common_name', 'host', 'path', 'file']
                for field in required_fields:
                    assert field in host_data, f"Host {host_key} missing {field}"
                    assert len(host_data[field]) > 0, f"Host {host_key} {field} is empty"


class TestURLValidation:
    """Tests for URL format validation across all sources."""

    def test_all_http_urls_valid_format(self):
        """Test all HTTP URLs match expected format."""
        from load_sources import list_databases, list_software
        
        http_pattern = re.compile(r'^https?://[^\s]+$')
        
        for category, entries in list_databases().items():
            if isinstance(entries, dict):
                for name, info in entries.items():
                    if isinstance(info, dict) and 'url' in info:
                        url = info['url']
                        if url.startswith('http'):
                            assert http_pattern.match(url), f"Invalid URL format: {url}"

    def test_ftp_urls_valid_format(self):
        """Test all FTP URLs match expected format."""
        from load_sources import list_hosts
        
        ftp_pattern = re.compile(r'^ftp://[^\s]+$')
        
        for host_key, host_data in list_hosts().items():
            if isinstance(host_data, dict) and 'host' in host_data:
                assert host_data['host'] == 'ftp.ncbi.nlm.nih.gov'

    def test_github_urls_are_valid(self):
        """Test all GitHub URLs are valid."""
        from load_sources import get_git_url
        
        git_urls = [
            get_git_url('fastviromeexplorer'),
            get_git_url('desamba'),
        ]
        
        for url in git_urls:
            if url:
                assert 'github.com' in url or 'github' in url.lower()


class TestConfigurationManagement:
    """Tests for configuration management functionality."""

    def test_sources_yaml_can_be_parsed(self):
        """Test sources.yaml can be parsed as YAML."""
        import yaml
        
        yaml_path = os.path.join(
            os.path.dirname(__file__), 
            '..', 
            'sources.yaml'
        )
        
        with open(yaml_path, 'r') as f:
            data = yaml.safe_load(f)
        
        assert data is not None
        assert 'databases' in data
        assert 'software' in data
        assert 'host_genomes' in data

    def test_version_and_metadata_present(self):
        """Test version and metadata fields are present."""
        import yaml
        
        yaml_path = os.path.join(
            os.path.dirname(__file__), 
            '..', 
            'sources.yaml'
        )
        
        with open(yaml_path, 'r') as f:
            data = yaml.safe_load(f)
        
        assert 'version' in data
        assert 'last_updated' in data
        assert data['version'] == '1.0'

    def test_sources_config_path_detection(self):
        """Test SourceLoader can find sources.yaml."""
        from load_sources import get_loader
        
        loader = get_loader()
        assert loader._config_path is not None
        assert os.path.exists(loader._config_path)


class TestErrorHandling:
    """Tests for error handling."""

    def test_get_db_url_missing_category_handled(self):
        """Test get_db_url handles missing category gracefully."""
        from load_sources import get_db_url
        
        result = get_db_url('nonexistent_category_xyz', 'viral')
        assert result is None

    def test_get_db_url_missing_db_handled(self):
        """Test get_db_url handles missing database gracefully."""
        from load_sources import get_db_url
        
        result = get_db_url('kraken2', 'nonexistent_db_xyz')
        assert result is None

    def test_get_host_config_invalid_key_handled(self):
        """Test get_host_config handles invalid key gracefully."""
        from load_sources import get_host_config
        
        result = get_host_config('invalid_host_name_xyz')
        assert result is None

    def test_get_software_url_missing_handled(self):
        """Test get_software_url handles missing software gracefully."""
        from load_sources import get_software_url
        
        result = get_software_url('nonexistent_software_xyz')
        assert result is None
