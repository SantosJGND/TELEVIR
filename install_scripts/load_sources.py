#!/usr/bin/python3
"""
Centralized source configuration loader for TELE-Vir.

Provides easy access to all download sources organized by category:
- databases: Kraken2, Kaiju, UniRef, SILVA, etc.
- software: Git repos, archive downloads
- host_genomes: NCBI RefSeq host genome locations
- infrastructure: Conda, tools, NCBI E-utilities

Usage:
    from load_sources import sources, get_db_url, get_host_config
    
    # Get a full section
    kraken2_dbs = sources['databases']['kraken2']
    
    # Get specific URL
    viral_db = get_db_url('kraken2', 'viral')
    
    # Get host genome config
    human = get_host_config('homo_sapiens')
"""

import os
import yaml
from pathlib import Path
from typing import Optional, Dict, Any, Union

try:
    import yaml
    Loader = yaml.CLoader if hasattr(yaml, 'CLoader') else yaml.Loader
    Dumper = yaml.CDumper if hasattr(yaml, 'CDumper') else yaml.Dumper
except ImportError:
    import yaml
    Loader = yaml.Loader
    Dumper = yaml.Dumper


class SourceLoader:
    """Centralized source configuration loader."""
    
    _instance = None
    _sources = None
    _config_path = None
    
    def __new__(cls, config_path: Optional[str] = None):
        if cls._instance is None:
            cls._instance = super().__new__(cls)
            cls._instance._initialized = False
        return cls._instance
    
    def __init__(self, config_path: Optional[str] = None):
        if self._initialized and self._sources is not None:
            return
            
        if config_path is None:
            config_path = self._find_config_path()
        
        self._config_path = config_path
        self._load_sources()
        self._initialized = True
    
    def _find_config_path(self) -> str:
        """Find sources.yaml in common locations."""
        possible_paths = [
            os.path.join(os.path.dirname(__file__), 'sources.yaml'),
            os.path.join(os.path.dirname(os.path.dirname(__file__)), 'sources.yaml'),
            os.path.join(os.getcwd(), 'sources.yaml'),
            '/opt/televir/sources.yaml',
        ]
        
        for path in possible_paths:
            if os.path.exists(path):
                return path
        
        return possible_paths[0]
    
    def _load_sources(self):
        """Load sources from YAML file."""
        try:
            with open(self._config_path, 'r') as f:
                self._sources = yaml.load(f, Loader=Loader)
        except FileNotFoundError:
            raise FileNotFoundError(
                f"Configuration file not found at: {self._config_path}\n"
                f"Please ensure sources.yaml exists in the project root or install_scripts directory."
            )
        except yaml.YAMLError as e:
            raise yaml.YAMLError(f"Error parsing sources.yaml: {e}")
    
    @property
    def sources(self) -> Dict[str, Any]:
        """Get all sources."""
        return self._sources
    
    @property
    def version(self) -> str:
        """Get config version."""
        return self._sources.get('version', 'unknown')
    
    @property
    def last_updated(self) -> str:
        """Get last update date."""
        return self._sources.get('last_updated', 'unknown')
    
    def get(self, *keys: str, default: Any = None) -> Any:
        """Get nested value from sources."""
        value = self._sources
        for key in keys:
            if isinstance(value, dict):
                value = value.get(key)
                if value is None:
                    return default
            else:
                return default
        return value
    
    def get_refseq_entry(self, organism: str, db_type: str) -> Optional[Dict]:
        """Get RefSeq entry for specific organism and type (protein/genome).
        
        Args:
            organism: 'viral' or 'bacterial'
            db_type: 'protein' or 'genome'
            
        Returns:
            Dict with url, description, file_pattern or None
        """
        return self.get('databases', 'refseq', organism, db_type)
    
    def get_db_url(self, category: str, name: str) -> Optional[str]:
        """Get database URL by category and name.
        
        Args:
            category: Database category (e.g., 'kraken2', 'kaiju', 'protein')
            name: Database name within category (e.g., 'viral', 'uniref90')
            
        Returns:
            URL string or None if not found
        """
        db_section = self._sources.get('databases', {}).get(category, {})
        if isinstance(db_section, dict):
            db_entry = db_section.get(name, {})
            if isinstance(db_entry, dict):
                return db_entry.get('url')
        return None
    
    def get_db_version(self, category: str, name: str) -> Optional[str]:
        """Get database version from sources.yaml.
        
        First checks 'version' field, then extracts from 'file' or 'url' string.
        
        Args:
            category: Database category (e.g., 'kraken2', 'kaiju', 'metaphlan')
            name: Database name within category
            
        Returns:
            Version string or None if not found
        """
        entry = self.get_db_entry(category, name)
        if not entry:
            return None
        
        version = entry.get('version')
        if version:
            return version
        
        file_str = entry.get('file', '') or entry.get('url', '')
        return self.extract_version_from_string(file_str)
    
    def get_db_entry(self, category: str, name: str) -> Optional[Dict]:
        """Get full database entry from sources.yaml.
        
        Args:
            category: Database category
            name: Database name
            
        Returns:
            Dict with url, version, file, description or None if not found
        """
        db_section = self._sources.get('databases', {}).get(category, {})
        if isinstance(db_section, dict):
            entry = db_section.get(name, {})
            if entry and isinstance(entry, dict):
                return entry
        return None
    
    def extract_version_from_string(self, s: str) -> Optional[str]:
        """Extract version from filename or URL string (case-insensitive).
        
        Patterns tried:
            _YYYYMMDD or _YYYYMMDD. -> 20250402
            YYYY-MM-DD -> 2024-08-15
            vYYYYxN -> vJan25
            
        Args:
            s: String to extract version from (filename or URL)
            
        Returns:
            Extracted version string or None
        """
        import re
        if not s:
            return None
        
        s_lower = s.lower()
        patterns = [
            (r'_(\d{8})(?:\.|_|$)', 1),  # _20250402. or _20250402_ or _20250402end
            (r'_(\d{6})(?:\.|_|$)', 1),  # _202504. or _202504_ or _202504end
            (r'\.(\d{8})\.', 1),    # .20250402.
            (r'(\d{4}-\d{2}-\d{2})', 1),  # 2024-08-15
            (r'_v([a-z]+\d{4})', 1),  # _vJan25 (month prefix)
            (r'(\d{8})(?:\.|_|$)', 1),  # 20250402. or 20250402_ or 20250402end
        ]
        
        for pattern, group in patterns:
            match = re.search(pattern, s_lower)
            if match:
                return match.group(1)
        return None
    
    def get_host_config(self, host_key: str) -> Optional[Dict[str, str]]:
        """Get host genome configuration.
        
        Args:
            host_key: Host identifier (e.g., 'homo_sapiens', 'hg38')
            
        Returns:
            Dict with host_name, common_name, host, path, file or None
        """
        hosts = self._sources.get('host_genomes', {})
        
        if host_key in hosts:
            config = hosts[host_key]
            if isinstance(config, dict):
                return config
            return None
        
        for host_data in hosts.values():
            if isinstance(host_data, dict):
                if host_data.get('host_name') == host_key or host_data.get('common_name') == host_key:
                    return host_data
        
        return None
    
    def get_software_url(self, name: str, repo_type: str = 'archives') -> Optional[str]:
        """Get software download URL.
        
        Args:
            name: Software name (e.g., 'clark', 'fastqc')
            repo_type: 'archives' or 'git_repos'
            
        Returns:
            URL string or None if not found
        """
        software = self._sources.get('software', {}).get(repo_type, {})
        if isinstance(software, dict):
            entry = software.get(name, {})
            if isinstance(entry, dict):
                return entry.get('url')
        return None
    
    def get_git_url(self, name: str) -> Optional[str]:
        """Get Git repository URL."""
        return self.get_software_url(name, 'git_repos')
    
    def get_infrastructure_url(self, category: str, name: str) -> Optional[str]:
        """Get infrastructure tool URL."""
        infra = self._sources.get('infrastructure', {}).get(category, {})
        if isinstance(infra, dict):
            entry = infra.get(name, {})
            if isinstance(entry, dict):
                return entry.get('url')
        return None
    
    def get_ncbi_eutils_url(self, endpoint: str) -> Optional[str]:
        """Get NCBI E-utilities endpoint URL."""
        return self.get_infrastructure_url('ncbi_eutils', endpoint)
    
    def list_databases(self, category: Optional[str] = None) -> Union[Dict, Dict[str, Dict]]:
        """List available databases.
        
        Args:
            category: Optional specific category to list
            
        Returns:
            Dict of databases with their metadata
        """
        dbs = self._sources.get('databases', {})
        if category:
            return dbs.get(category, {})
        return dbs
    
    def list_hosts(self) -> Dict[str, Dict]:
        """List all available host genomes."""
        return self._sources.get('host_genomes', {})
    
    def list_software(self) -> Dict[str, Dict]:
        """List all available software."""
        return self._sources.get('software', {})


_loader = None

def get_loader(config_path: Optional[str] = None) -> SourceLoader:
    """Get singleton SourceLoader instance."""
    global _loader
    if _loader is None:
        _loader = SourceLoader(config_path)
    return _loader


def reload(config_path: Optional[str] = None):
    """Reload configuration from file."""
    global _loader
    _loader = SourceLoader(config_path)
    return _loader


def sources() -> Dict[str, Any]:
    """Get all sources (convenience accessor)."""
    return get_loader().sources


def version() -> str:
    """Get config version."""
    return get_loader().version


def last_updated() -> str:
    """Get last update date."""
    return get_loader().last_updated


def get(*keys: str, default: Any = None) -> Any:
    """Get nested value from sources."""
    return get_loader().get(*keys, default=default)

def get_db_url(category: str, name: str) -> Optional[str]:
    """Get database URL."""
    return get_loader().get_db_url(category, name)

def get_db_version(category: str, name: str) -> Optional[str]:
    """Get database version from sources.yaml.
    
    First checks 'version' field, then extracts from 'file' or 'url' string.
    """
    return get_loader().get_db_version(category, name)

def get_db_entry(category: str, name: str) -> Optional[Dict]:
    """Get full database entry from sources.yaml."""
    return get_loader().get_db_entry(category, name)

def get_refseq_entry(organism: str, db_type: str) -> Optional[Dict]:
    """Get RefSeq entry for specific organism and type (protein/genome).
    
    Args:
        organism: 'viral' or 'bacterial'
        db_type: 'protein' or 'genome'
        
    Returns:
        Dict with url, description, file_pattern or None
    """
    return get_loader().get_refseq_entry(organism, db_type)

def extract_version_from_string(s: str) -> Optional[str]:
    """Extract version from filename or URL string (case-insensitive)."""

    temp = get_loader().extract_version_from_string(s)
    print(temp)
    return temp


def get_host_config(host_key: str) -> Optional[Dict[str, str]]:
    """Get host genome configuration."""
    return get_loader().get_host_config(host_key)

def get_software_url(name: str, repo_type: str = 'archives') -> Optional[str]:
    """Get software download URL."""
    return get_loader().get_software_url(name, repo_type)

def get_software_entry(name: str) -> Optional[Dict]:
    """Get software entry from sources.yaml (checks archives and git_repos).
    
    Args:
        name: Software name to look up
        
    Returns:
        Dict with url, description, file, branch or None if not found
    """
    software = get_loader().sources.get('software', {})
    
    archives = software.get('archives', {})
    if name in archives:
        return archives[name]
    
    git_repos = software.get('git_repos', {})
    if name in git_repos:
        return git_repos[name]
    
    return None

def get_git_url(name: str) -> Optional[str]:
    """Get Git repository URL."""
    return get_loader().get_git_url(name)

def get_ncbi_eutils_url(endpoint: str) -> Optional[str]:
    """Get NCBI E-utilities endpoint URL."""
    return get_loader().get_ncbi_eutils_url(endpoint)

def list_databases(category: Optional[str] = None):
    """List available databases."""
    return get_loader().list_databases(category)

def list_hosts():
    """List all host genomes."""
    return get_loader().list_hosts()

def list_software():
    """List all software."""
    return get_loader().list_software()


class LazySource:
    """Lazy loader for individual source entries."""
    
    def __init__(self, loader: SourceLoader, *keys: str):
        self._loader = loader
        self._keys = keys
        self._value = None
    
    def __str__(self) -> str:
        if self._value is None:
            self._value = self._loader.get(*self._keys)
        return str(self._value) if self._value else ""
    
    def __repr__(self) -> str:
        return f"LazySource({', '.join(repr(k) for k in self._keys)})"
    
    def get(self) -> Any:
        """Get the actual value."""
        if self._value is None:
            self._value = self._loader.get(*self._keys)
        return self._value


if __name__ == '__main__':
    print("TELE-Vir Source Configuration")
    print("=" * 50)
    print(f"Version: {version()}")
    print(f"Last Updated: {last_updated()}")
    print()
    
    print("Available Categories:")
    print("-" * 30)
    for section in sources().keys():
        if section not in ('version', 'last_updated'):
            print(f"  - {section}")
    print()
    
    print("Example Queries:")
    print("-" * 30)
    print(f"  Kraken2 viral URL: {get_db_url('kraken2', 'viral')}")
    print(f"  Human host config: {get_host_config('homo_sapiens')}")
    print(f"  FastQC URL: {get_software_url('fastqc')}")
    print(f"  SILVA 16S URL: {get_db_url('ribosomal_rna', 'silva_16s')}")
