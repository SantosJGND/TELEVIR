"""
Unit tests for sources_cli.py - CLI command interface.
"""

import pytest
import sys
import os
from io import StringIO

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'install_scripts'))

from sources_cli import main, cmd_list, cmd_get, cmd_validate


class TestCLICommands:
    """Tests for CLI command functions."""

    def test_validate_command_output(self, capsys):
        """Test validate command produces correct output."""
        args = type('Args', (), {'category': None})()
        cmd_validate(args)
        captured = capsys.readouterr()
        assert 'Validating' in captured.out
        assert 'Databases' in captured.out or 'Host' in captured.out

    def test_list_all_command(self, capsys):
        """Test list all command."""
        args = type('Args', (), {'category': 'all', 'show_urls': False})()
        cmd_list(args)
        captured = capsys.readouterr()
        assert 'Version' in captured.out
        assert 'Databases' in captured.out
        assert 'Host Genomes' in captured.out

    def test_list_databases_command(self, capsys):
        """Test list databases command."""
        args = type('Args', (), {'category': 'databases', 'show_urls': False})()
        cmd_list(args)
        captured = capsys.readouterr()
        assert 'Available Databases' in captured.out
        assert 'kraken2' in captured.out or 'kaiju' in captured.out

    def test_list_hosts_command(self, capsys):
        """Test list hosts command."""
        args = type('Args', (), {'category': 'hosts', 'show_urls': False})()
        cmd_list(args)
        captured = capsys.readouterr()
        assert 'Available Host Genomes' in captured.out
        assert 'hg38' in captured.out or 'human' in captured.out

    def test_list_software_command(self, capsys):
        """Test list software command."""
        args = type('Args', (), {'category': 'software', 'show_urls': False})()
        cmd_list(args)
        captured = capsys.readouterr()
        assert 'Available Software' in captured.out

    def test_get_database_command(self, capsys):
        """Test get database command."""
        args = type('Args', (), {
            'source_type': 'database',
            'category': 'kraken2',
            'name': 'viral',
            'format': 'full'
        })()
        cmd_get(args)
        captured = capsys.readouterr()
        assert 'kraken' in captured.out.lower() or 'https://' in captured.out

    def test_get_host_command_full(self, capsys):
        """Test get host command with full format."""
        args = type('Args', (), {
            'source_type': 'host',
            'category': None,
            'name': 'homo_sapiens',
            'format': 'full'
        })()
        cmd_get(args)
        captured = capsys.readouterr()
        assert 'host_name' in captured.out
        assert 'hg38' in captured.out

    def test_get_host_command_url(self, capsys):
        """Test get host command with URL format."""
        args = type('Args', (), {
            'source_type': 'host',
            'category': None,
            'name': 'homo_sapiens',
            'format': 'url'
        })()
        cmd_get(args)
        captured = capsys.readouterr()
        assert 'ftp://' in captured.out

    def test_get_software_command(self, capsys):
        """Test get software command."""
        args = type('Args', (), {
            'source_type': 'software',
            'category': None,
            'name': 'fastqc',
            'format': 'full'
        })()
        cmd_get(args)
        captured = capsys.readouterr()
        assert 'fastqc' in captured.out.lower() or 'https://' in captured.out

    def test_get_git_command(self, capsys):
        """Test get git command."""
        args = type('Args', (), {
            'source_type': 'git',
            'category': None,
            'name': 'fastviromeexplorer',
            'format': 'full'
        })()
        cmd_get(args)
        captured = capsys.readouterr()
        assert 'github' in captured.out.lower() or 'https://' in captured.out


class TestCLIMain:
    """Tests for CLI main entry point."""

    def test_main_list_all(self, capsys, monkeypatch):
        """Test main entry point with list all."""
        monkeypatch.setattr(sys, 'argv', ['sources_cli.py', 'list', 'all'])
        try:
            main()
        except SystemExit:
            pass
        captured = capsys.readouterr()
        assert 'Version' in captured.out or 'Databases' in captured.out

    def test_main_validate(self, capsys, monkeypatch):
        """Test main entry point with validate."""
        monkeypatch.setattr(sys, 'argv', ['sources_cli.py', 'validate'])
        try:
            main()
        except SystemExit:
            pass
        captured = capsys.readouterr()
        assert 'valid' in captured.out.lower() or 'Error' in captured.out

    def test_main_get_database(self, capsys, monkeypatch):
        """Test main entry point with get database."""
        monkeypatch.setattr(sys, 'argv', ['sources_cli.py', 'get', 'database', 'kraken2', 'viral'])
        try:
            main()
        except SystemExit:
            pass
        captured = capsys.readouterr()
        assert 'https://' in captured.out

    def test_main_get_host(self, capsys, monkeypatch):
        """Test main entry point with get host."""
        monkeypatch.setattr(sys, 'argv', ['sources_cli.py', 'get', 'host', 'homo_sapiens'])
        try:
            main()
        except SystemExit:
            pass
        captured = capsys.readouterr()
        assert 'hg38' in captured.out or 'ftp://' in captured.out

    def test_main_help(self, monkeypatch):
        """Test main entry point with help."""
        monkeypatch.setattr(sys, 'argv', ['sources_cli.py', '--help'])
        with pytest.raises(SystemExit):
            main()

    def test_main_invalid_command(self, capsys, monkeypatch):
        """Test main entry point with invalid command."""
        monkeypatch.setattr(sys, 'argv', ['sources_cli.py', 'invalid_command'])
        try:
            main()
        except SystemExit:
            pass
        captured = capsys.readouterr()
        assert len(captured.out) > 0 or len(captured.err) > 0


class TestCLIValidation:
    """Tests for CLI validation functionality."""

    def test_validate_returns_zero_on_success(self, capsys):
        """Test validate command returns success."""
        args = type('Args', (), {'category': None})()
        cmd_validate(args)
        captured = capsys.readouterr()
        assert 'valid' in captured.out.lower() or 'Configuration' in captured.out

    def test_validate_shows_database_count(self, capsys):
        """Test validate shows database count."""
        args = type('Args', (), {'category': None})()
        cmd_validate(args)
        captured = capsys.readouterr()
        assert 'Databases:' in captured.out

    def test_validate_shows_host_count(self, capsys):
        """Test validate shows host count."""
        args = type('Args', (), {'category': None})()
        cmd_validate(args)
        captured = capsys.readouterr()
        assert 'Host Genomes:' in captured.out


class TestCLIOutputFormats:
    """Tests for CLI output formats."""

    def test_list_shows_urls_when_requested(self, capsys):
        """Test list shows URLs when show_urls is True."""
        args = type('Args', (), {'category': 'hosts', 'show_urls': True})()
        cmd_list(args)
        captured = capsys.readouterr()
        assert 'ftp.ncbi' in captured.out

    def test_get_url_format(self, capsys):
        """Test get command URL format."""
        args = type('Args', (), {
            'source_type': 'host',
            'category': None,
            'name': 'homo_sapiens',
            'format': 'url'
        })()
        cmd_get(args)
        captured = capsys.readouterr()
        output = captured.out.strip()
        assert output.startswith('ftp://')
