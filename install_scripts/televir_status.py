#!/usr/bin/env python3
"""
TELE-Vir Status GUI Application

Displays all available databases, host genomes, and software from sources.yaml
with their installation status from utility_local.db.
"""

import tkinter as tk
from tkinter import ttk
import os
import sys
import sqlite3
from pathlib import Path

sys.path.insert(0, os.path.dirname(__file__))

from load_sources import list_databases, list_software, list_hosts
from config import DATABASE_FILENAME
from install_source import ENVS_PARAMS


class TelevirStatusApp:
    def __init__(self, root):
        self.root = root
        self.root.title("TELE-Vir Status")
        self.root.geometry("1000x700")
        
        self.install_home = os.environ.get('INSTALL_HOME', '/opt/televir')
        self.db_path = os.path.join(self.install_home, DATABASE_FILENAME)
        self.env_root = os.path.join(self.install_home, "envs")
        
        # Mapping from archive/git names to their install directories
        self.ARCHIVE_INSTALL_PATHS = {
            "clark": "classification/Clark",
            "voyager": "classification/Voyager/voyager-cli",
            "trimmomatic": "trimmomatic",
            "fastqc": "fastqc",
            "rabbitqc": "RabbitQC",
        }
        
        self.GIT_INSTALL_PATHS = {
            "fastviromeexplorer": "FastViromeExplorer",
            "desamba": "classm_lc/deSAMBA",
            "rabbitqc_git": "RabbitQC",
            "trimmomatic_git": "trimmomatic",
        }
        
        self._create_widgets()
        
        self.db_tree.tag_configure("installed", background="#3EBE0B")
        self.hosts_tree.tag_configure("installed", background="#3EBE0B")
        self.soft_tree.tag_configure("installed", background="#3EBE0B")

    def _create_widgets(self):
        title = tk.Label(self.root, text="TELE-Vir Status", font=("Arial", 18, "bold"))
        title.pack(pady=10)
        
        db_frame = tk.LabelFrame(self.root, text="Databases", font=("Arial", 12, "bold"))
        db_frame.pack(fill="both", expand=True, padx=10, pady=5)
        
        self.db_tree = ttk.Treeview(
            db_frame, 
            columns=("name", "category", "description", "file", "available", "installed", "version", "date"), 
            show="headings"
        )
        self.db_tree.heading("name", text="Name")
        self.db_tree.heading("category", text="Category")
        self.db_tree.heading("description", text="Description")
        self.db_tree.heading("file", text="File")
        self.db_tree.heading("available", text="Available")
        self.db_tree.heading("installed", text="Installed")
        self.db_tree.heading("version", text="Version")
        self.db_tree.heading("date", text="Date")
        self.db_tree.column("name", width=150)
        self.db_tree.column("category", width=100)
        self.db_tree.column("description", width=200)
        self.db_tree.column("file", width=150)
        self.db_tree.column("available", width=70)
        self.db_tree.column("installed", width=70)
        self.db_tree.column("version", width=80)
        self.db_tree.column("date", width=90)
        self.db_tree.pack(fill="both", expand=True, padx=5, pady=5)
        
        hosts_frame = tk.LabelFrame(self.root, text="Host Genomes", font=("Arial", 12, "bold"))
        hosts_frame.pack(fill="both", expand=True, padx=10, pady=5)
        
        self.hosts_tree = ttk.Treeview(
            hosts_frame, 
            columns=("name", "common_name", "available", "installed", "filename"), 
            show="headings"
        )
        self.hosts_tree.heading("name", text="Name")
        self.hosts_tree.heading("common_name", text="Common Name")
        self.hosts_tree.heading("available", text="Available")
        self.hosts_tree.heading("installed", text="Installed")
        self.hosts_tree.heading("filename", text="Filename")
        self.hosts_tree.column("name", width=120)
        self.hosts_tree.column("common_name", width=100)
        self.hosts_tree.column("available", width=80)
        self.hosts_tree.column("installed", width=80)
        self.hosts_tree.column("filename", width=150)
        self.hosts_tree.pack(fill="both", expand=True, padx=5, pady=5)
        
        soft_frame = tk.LabelFrame(self.root, text="Software", font=("Arial", 12, "bold"))
        soft_frame.pack(fill="both", expand=True, padx=10, pady=5)
        
        self.soft_tree = ttk.Treeview(
            soft_frame, 
            columns=("name", "type", "available", "installed", "tag"), 
            show="headings"
        )
        self.soft_tree.heading("name", text="Name")
        self.soft_tree.heading("type", text="Type")
        self.soft_tree.heading("available", text="Available")
        self.soft_tree.heading("installed", text="Installed")
        self.soft_tree.heading("tag", text="Tag")
        self.soft_tree.column("name", width=200)
        self.soft_tree.column("type", width=100)
        self.soft_tree.column("available", width=80)
        self.soft_tree.column("installed", width=80)
        self.soft_tree.column("tag", width=100)
        self.soft_tree.pack(fill="both", expand=True, padx=5, pady=5)
        
        check_btn = tk.Button(
            self.root, 
            text="Check Status", 
            command=self.check_status, 
            bg="#4CAF50", 
            fg="white", 
            font=("Arial", 12, "bold"),
            padx=20,
            pady=5
        )
        check_btn.pack(side="left", padx=5, pady=10)
        
        export_btn = tk.Button(
            self.root, 
            text="Export to TSV", 
            command=self.export_tsv, 
            bg="#2196F3", 
            fg="white", 
            font=("Arial", 12, "bold"),
            padx=20,
            pady=5
        )
        export_btn.pack(side="left", padx=5, pady=10)
    
    def check_status(self):
        self._populate_databases()
        self._populate_hosts()
        self._populate_software()
    
    def _populate_databases(self):
        for item in self.db_tree.get_children():
            self.db_tree.delete(item)
        
        dbs = list_databases()

        installed_dbs = self._get_installed_databases()

        for category, entries in dbs.items():
            if isinstance(entries, dict):
                for name, info in entries.items():
                    name_full = f"{category}/{name}"
                    available = "✓" if info else "✗"
                    description = installed_dbs.get(name_full, {}).get("description", "N/A")
                    file_info = installed_dbs.get(name_full, {}).get("path", "N/A")
                    version = installed_dbs.get(name_full, {}).get("version", "N/A")
                    installed = installed_dbs.get(name_full, {}).get("installed", "N/A")
                    date = installed_dbs.get(name_full, {}).get("date", "N/A")
                    tag = "installed" if installed == "True" else "not_installed"
                    self.db_tree.insert("", "end", values=(name_full, category, description, file_info, available, installed, version, date), tags=(tag,))
    
    def _populate_hosts(self):
        for item in self.hosts_tree.get_children():
            self.hosts_tree.delete(item)
        
        hosts = list_hosts()
        installed_hosts = self._get_installed_hosts()
        
        for key, info in hosts.items():
            if isinstance(info, dict):
                host_name = info.get('host_name', key)
                common = info.get('common_name', 'N/A')
                available = "✓" if info else "✗"
                filename = "N/A"
                installed = "✗"
                for display_name, fn in installed_hosts.items():
                    if host_name.lower() in display_name.lower():
                        filename = fn
                        installed = "✓"
                        break
                tag = "installed" if installed == "✓" else "not_installed"
                self.hosts_tree.insert("", "end", values=(host_name, common, available, installed, filename), tags=(tag,))
    
    def _populate_software(self):
        for item in self.soft_tree.get_children():
            self.soft_tree.delete(item)
        
        soft = list_software()
        installed_soft = self._get_installed_software()
        
        for name, info in soft.get('archives', {}).items():
            if info and isinstance(info, dict):
                available = "✓" if self._check_archive_installed(name, info) else "✗"
                tag = installed_soft.get(name, "N/A")
                installed = "✓" if tag != "N/A" else "✗"
                self.soft_tree.insert("", "end", values=(name, "archive", available, installed, tag), tags=("installed" if installed == "✓" else "not_installed",))
        
        for name, info in soft.get('git_repos', {}).items():
            if info and isinstance(info, dict):
                available = "✓" if self._check_git_installed(name, info) else "✗"
                tag = installed_soft.get(name, "N/A")
                installed = "✓" if tag != "N/A" else "✗"
                self.soft_tree.insert("", "end", values=(name, "git", available, installed, tag), tags=("installed" if installed == "✓" else "not_installed",))
        
        for name, info in soft.get('conda_tools', {}).items():
            if info and isinstance(info, dict):
                yaml_file = info.get('yaml', '')
                binary_name = info.get('binary', name)
                available = "✓" if self._check_conda_binary(yaml_file, binary_name) else "✗"
                tag = installed_soft.get(name, "N/A")
                installed = "✓" if tag != "N/A" else "✗"
                self.soft_tree.insert("", "end", values=(name, "conda", available, installed, tag), tags=("installed" if installed == "✓" else "not_installed",))
    
    def _get_installed_databases(self):
        """Query installed databases (non-host) with version."""
        installed = {}
        if not os.path.exists(self.db_path):
            return installed
        
        try:
            conn = sqlite3.connect(self.db_path)
            cursor = conn.cursor()
            cursor.execute(
                "SELECT name, db_category, db_name, version, date, installed, description, path FROM database WHERE installed = 'True' AND db_type != 'host'"
            )
            for row in cursor.fetchall():
                name, db_category, db_name = row[0], row[1], row[2]
                display_name = f"{db_category}/{db_name}" if db_category and db_name else name
                installed[display_name] = {
                    "version": row[3] if row[3] else "N/A",
                    "date": row[4] if row[4] else "N/A",
                    "installed": row[5] if row[5] else "N/A",
                    "description": row[6] if row[6] else "N/A",
                    "path": row[7] if row[7] else "N/A"
                }
            conn.close()
        except sqlite3.OperationalError:
            cursor.execute(
                "SELECT name, version, date, installed, description, path FROM database WHERE installed = 'True' AND software != 'host'"
            )
            for row in cursor.fetchall():
                installed[row[0]] = {
                    "version": row[1] if row[1] else "N/A",
                    "date": row[2] if row[2] else "N/A",
                    "installed": row[3] if row[3] else "N/A",
                    "description": row[4] if row[4] else "N/A",
                    "path": row[5] if row[5] else "N/A"
                }
            conn.close()
        except Exception as e:
            print(f"Error reading database: {e}")
        
        return installed

    def _get_all_databases(self):
        installed = {}

        if not os.path.exists(self.db_path):
            return {}

        try:
            conn = sqlite3.connect(self.db_path)
            cursor = conn.cursor()
            cursor.execute("SELECT * FROM database")
            for row in cursor.fetchall():
                installed[row[0]] = row
            conn.close()
        except Exception as e:
            print(f"Error reading database: {e}")

        return installed

    def _get_installed_hosts(self):
        """Query installed host genomes with filename."""
        installed = {}
        if not os.path.exists(self.db_path):
            return installed
        
        try:
            conn = sqlite3.connect(self.db_path)
            cursor = conn.cursor()
            cursor.execute(
                "SELECT name, db_category, db_name, path FROM database WHERE installed = 'True' AND db_type = 'host'"
            )
            for row in cursor.fetchall():
                name, db_category, db_name = row[0], row[1], row[2]
                display_name = f"{db_category}/{db_name}" if db_category and db_name else name
                filename = os.path.basename(row[3]) if row[3] else "N/A"
                installed[display_name] = filename
            conn.close()
        except sqlite3.OperationalError:
            cursor.execute(
                "SELECT name, path FROM database WHERE installed = 'True' AND software = 'host'"
            )
            for row in cursor.fetchall():
                filename = os.path.basename(row[1]) if row[1] else "N/A"
                installed[row[0]] = filename
            conn.close()
        except Exception as e:
            print(f"Error reading database: {e}")
        
        return installed
    
    def _get_installed_software(self):
        """Query installed software with tag."""
        installed = {}
        if not os.path.exists(self.db_path):
            return installed
        
        try:
            conn = sqlite3.connect(self.db_path)
            cursor = conn.cursor()
            cursor.execute("SELECT name, tag FROM software WHERE installed = 'True'")
            for row in cursor.fetchall():
                installed[row[0]] = row[1] if row[1] else "N/A"
            conn.close()
        except Exception as e:
            print(f"Error reading database: {e}")
        
        return installed

    def _get_env_path_from_yaml(self, yaml_file: str) -> str:
        """Get environment path from yaml file using ENVS_PARAMS mapping."""
        # ENVS maps env_path -> yaml_file, so we need to reverse lookup
        envs_map = ENVS_PARAMS.get("ENVS", {})
        for env_path, yml in envs_map.items():
            if yml == yaml_file:
                return env_path
        return ""

    def _check_conda_binary(self, yaml_file: str, binary_name: str) -> bool:
        """Check if binary exists in conda environment."""
        env_path = self._get_env_path_from_yaml(yaml_file)
        if not env_path:
            return False
        binary_path = os.path.join(self.env_root, env_path, "bin", binary_name)
        return os.path.isfile(binary_path)

    def _check_archive_installed(self, name: str, soft_info: dict) -> bool:
        """Check if archive-based software is installed."""
        install_path = self.ARCHIVE_INSTALL_PATHS.get(name, "")
        if not install_path:
            return False
        full_path = os.path.join(self.env_root, install_path)
        return os.path.exists(full_path)

    def _check_git_installed(self, name: str, soft_info: dict) -> bool:
        """Check if git repo is cloned."""
        install_path = self.GIT_INSTALL_PATHS.get(name, "")
        if not install_path:
            return False
        full_path = os.path.join(self.env_root, install_path)
        return os.path.exists(full_path)

    def export_tsv(self):
        from tkinter import filedialog
        
        filepath = filedialog.asksaveasfilename(
            defaultextension=".tsv",
            filetypes=[("TSV files", "*.tsv"), ("All files", "*.*")],
            initialfile="televir_status.tsv"
        )
        if not filepath:
            return
        
        with open(filepath, 'w') as f:
            f.write("# TELE-Vir Status Export\n\n")
            
            f.write("# Databases\n")
            f.write("Name\tCategory\tDescription\tFile\tAvailable\tInstalled\tVersion\tDate\n")
            for row in self.db_tree.get():
                f.write("\t".join(str(v) for v in row) + "\n")
            
            f.write("\n# Host Genomes\n")
            f.write("Name\tCommon Name\tAvailable\tInstalled\tFilename\n")
            for row in self.hosts_tree.get():
                f.write("\t".join(str(v) for v in row) + "\n")
            
            f.write("\n# Software\n")
            f.write("Name\tType\tAvailable\tInstalled\tTag\n")
            for row in self.soft_tree.get():
                f.write("\t".join(str(v) for v in row) + "\n")


def main():
    root = tk.Tk()
    app = TelevirStatusApp(root)
    root.mainloop()


if __name__ == "__main__":
    main()
