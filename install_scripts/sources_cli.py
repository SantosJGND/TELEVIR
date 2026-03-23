#!/usr/bin/env python3
"""
Sources management CLI for TELE-Vir.

Usage:
    python sources_cli.py list databases
    python sources_cli.py list hosts
    python sources_cli.py get kraken2 viral
    python sources_cli.py get host homo_sapiens
    python sources_cli.py validate
    python sources_cli.py diff <other_yaml>
"""

import argparse
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from load_sources import (
    get_loader, get_db_url, get_host_config, get_software_url,
    get_git_url, list_databases, list_hosts, list_software
)


def cmd_list(args):
    """List sources by category."""
    if args.category == 'databases':
        dbs = list_databases()
        print("Available Databases:")
        print("-" * 60)
        for cat, entries in dbs.items():
            if isinstance(entries, dict):
                print(f"\n  [{cat}]")
                for name, info in entries.items():
                    if isinstance(info, dict):
                        url = info.get('url', 'N/A')
                        desc = info.get('description', '')
                        print(f"    {name}: {desc}")
                        if args.show_urls:
                            print(f"      URL: {url}")
    elif args.category == 'hosts':
        hosts = list_hosts()
        print("Available Host Genomes:")
        print("-" * 60)
        for name, info in hosts.items():
            if isinstance(info, dict) and 'host_name' in info:
                print(f"  {info['host_name']} ({info['common_name']})")
                if args.show_urls:
                    print(f"    Host: {info['host']}")
                    print(f"    Path: {info['path']}")
    elif args.category == 'software':
        software = list_software()
        print("Available Software:")
        print("-" * 60)
        for repo_type, entries in software.items():
            if isinstance(entries, dict):
                print(f"\n  [{repo_type}]")
                for name, info in entries.items():
                    if isinstance(info, dict):
                        desc = info.get('description', '')
                        print(f"    {name}: {desc}")
                        if args.show_urls:
                            print(f"      URL: {info.get('url', 'N/A')}")
    elif args.category == 'all':
        cmd_list_types(args)


def cmd_list_types(args):
    """List all source types."""
    loader = get_loader()
    print("TELE-Vir Source Configuration")
    print("=" * 60)
    print(f"Version: {loader.version}")
    print(f"Last Updated: {loader.last_updated}")
    print()
    
    for cat in ['databases', 'software', 'host_genomes']:
        if cat == 'databases':
            dbs = list_databases()
            count = sum(len(v) for v in dbs.values() if isinstance(v, dict))
            print(f"Databases: {count} entries")
        elif cat == 'software':
            sw = list_software()
            count = sum(len(v) for v in sw.values() if isinstance(v, dict))
            print(f"Software: {count} entries")
        elif cat == 'host_genomes':
            hosts = list_hosts()
            count = sum(1 for v in hosts.values() if isinstance(v, dict) and 'host_name' in v)
            print(f"Host Genomes: {count} entries")
    print()


def cmd_get(args):
    """Get specific source URL or config."""
    if args.source_type == 'database':
        url = get_db_url(args.category, args.name)
        if url:
            print(url)
        else:
            print(f"Database not found: {args.category}/{args.name}", file=sys.stderr)
            sys.exit(1)
    elif args.source_type == 'host':
        config = get_host_config(args.name)
        if config:
            if args.format == 'url':
                print(f"ftp://{config['host']}{config['path']}{config['file']}")
            else:
                for k, v in config.items():
                    print(f"{k}: {v}")
        else:
            print(f"Host not found: {args.name}", file=sys.stderr)
            sys.exit(1)
    elif args.source_type == 'software':
        url = get_software_url(args.name)
        if url:
            print(url)
        else:
            print(f"Software not found: {args.name}", file=sys.stderr)
            sys.exit(1)
    elif args.source_type == 'git':
        url = get_git_url(args.name)
        if url:
            print(url)
        else:
            print(f"Git repo not found: {args.name}", file=sys.stderr)
            sys.exit(1)


def cmd_validate(args):
    """Validate sources configuration."""
    print("Validating sources configuration...")
    print("-" * 60)
    
    errors = []
    loader = get_loader()
    
    if not loader.sources:
        errors.append("Failed to load sources.yaml")
    
    dbs = list_databases()
    if not dbs:
        errors.append("No databases found")
    else:
        print(f"  Databases: {sum(len(v) for v in dbs.values() if isinstance(v, dict))} entries")
    
    hosts = list_hosts()
    host_count = sum(1 for v in hosts.values() if isinstance(v, dict) and 'host_name' in v)
    if host_count == 0:
        errors.append("No host genomes found")
    else:
        print(f"  Host Genomes: {host_count} entries")
    
    software = list_software()
    if software:
        print(f"  Software: {sum(len(v) for v in software.values() if isinstance(v, dict))} entries")
    
    print()
    if errors:
        print("ERRORS:")
        for e in errors:
            print(f"  - {e}")
        sys.exit(1)
    else:
        print("Configuration is valid!")


def cmd_diff(args):
    """Compare two source configurations."""
    import yaml
    
    if not os.path.exists(args.other_yaml):
        print(f"File not found: {args.other_yaml}", file=sys.stderr)
        sys.exit(1)
    
    with open(args.other_yaml) as f:
        other = yaml.safe_load(f)
    
    current = get_loader().sources
    
    print("Comparing source configurations...")
    print(f"  Current: sources.yaml")
    print(f"  Other:   {args.other_yaml}")
    print()
    
    for section in ['databases', 'software', 'host_genomes']:
        if section in current:
            c_items = current[section]
            o_items = other.get(section, {})
            
            if isinstance(c_items, dict) and isinstance(o_items, dict):
                c_keys = set(k for k, v in c_items.items() if isinstance(v, dict))
                o_keys = set(k for k, v in o_items.items() if isinstance(v, dict))
                
                added = c_keys - o_keys
                removed = o_keys - c_keys
                
                if added:
                    print(f"[{section}] Added categories: {added}")
                if removed:
                    print(f"[{section}] Removed categories: {removed}")


def main():
    parser = argparse.ArgumentParser(
        description="TELE-Vir source configuration management",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    subparsers = parser.add_subparsers(dest='command', help='Available commands')
    
    list_parser = subparsers.add_parser('list', help='List sources by category')
    list_parser.add_argument('category', choices=['databases', 'hosts', 'software', 'all'])
    list_parser.add_argument('-u', '--show-urls', action='store_true', help='Show URLs')
    
    get_parser = subparsers.add_parser('get', help='Get specific source')
    get_parser.add_argument('source_type', choices=['database', 'host', 'software', 'git'])
    get_parser.add_argument('category', nargs='?', help='Database category (for database type)')
    get_parser.add_argument('name', help='Source name')
    get_parser.add_argument('-f', '--format', choices=['url', 'full'], default='full', help='Output format')
    
    subparsers.add_parser('validate', help='Validate configuration')
    
    diff_parser = subparsers.add_parser('diff', help='Compare configurations')
    diff_parser.add_argument('other_yaml', help='Path to other YAML file')
    
    args = parser.parse_args()
    
    if args.command == 'list':
        cmd_list(args)
    elif args.command == 'get':
        cmd_get(args)
    elif args.command == 'validate':
        cmd_validate(args)
    elif args.command == 'diff':
        cmd_diff(args)
    else:
        parser.print_help()


if __name__ == '__main__':
    main()
