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


class TelevirStatusApp:
    def __init__(self, root):
        self.root = root
        self.root.title("TELE-Vir Status")
        self.root.geometry("1000x700")
        
        data_home = os.environ.get('TELEVIR_DATA_HOME', '/home/xpto/INSaFLU/data/televir')
        self.db_path = os.path.join(data_home, 'utility_docker.safe.db')
        
        self._create_widgets()
    
    def _create_widgets(self):
        title = tk.Label(self.root, text="TELE-Vir Status", font=("Arial", 18, "bold"))
        title.pack(pady=10)
        
        db_frame = tk.LabelFrame(self.root, text="Databases", font=("Arial", 12, "bold"))
        db_frame.pack(fill="both", expand=True, padx=10, pady=5)
        
        self.db_tree = ttk.Treeview(
            db_frame, 
            columns=("name", "category", "available", "installed", "version"), 
            show="headings"
        )
        self.db_tree.heading("name", text="Name")
        self.db_tree.heading("category", text="Category")
        self.db_tree.heading("available", text="Available")
        self.db_tree.heading("installed", text="Installed")
        self.db_tree.heading("version", text="Version")
        self.db_tree.column("name", width=200)
        self.db_tree.column("category", width=120)
        self.db_tree.column("available", width=80)
        self.db_tree.column("installed", width=80)
        self.db_tree.column("version", width=100)
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
        check_btn.pack(pady=10)
    
    def check_status(self):
        self._populate_databases()
        self._populate_hosts()
        self._populate_software()
    
    def _populate_databases(self):
        for item in self.db_tree.get_children():
            self.db_tree.delete(item)
        
        dbs = list_databases()
        print(dbs)

        all_dbs = self._get_all_databases()
        print(all_dbs)
        installed_dbs = self._get_installed_databases()

        print(installed_dbs)
        print("####")
        for category, entries in dbs.items():
            if isinstance(entries, dict):
                for name, info in entries.items():
                    name_full = f"{category}/{name}"
                    available = "✓" if info else "✗"
                    version = installed_dbs.get(name, {}).get("version", "N/A")
                    installed = installed_dbs.get(name, {}).get("installed", "N/A")
                    self.db_tree.insert("", "end", values=(name_full, category, available, installed, version))
    
    def _populate_hosts(self):
        for item in self.hosts_tree.get_children():
            self.hosts_tree.delete(item)
        
        hosts = list_hosts()
        installed_hosts = self._get_installed_hosts()
        
        for key, info in hosts.items():
            if isinstance(info, dict):
                name = info.get('host_name', key)
                common = info.get('common_name', 'N/A')
                available = "✓" if info else "✗"
                filename = installed_hosts.get(name, "N/A")
                installed = "✓" if filename != "N/A" else "✗"
                self.hosts_tree.insert("", "end", values=(name, common, available, installed, filename))
    
    def _populate_software(self):
        for item in self.soft_tree.get_children():
            self.soft_tree.delete(item)
        
        soft = list_software()
        installed_soft = self._get_installed_software()
        
        for name, info in soft.get('archives', {}).items():
            if info and isinstance(info, dict):
                available = "✓"
                tag = installed_soft.get(name, "N/A")
                installed = "✓" if tag != "N/A" else "✗"
                self.soft_tree.insert("", "end", values=(name, "archive", available, installed, tag))
        
        for name, info in soft.get('git_repos', {}).items():
            if info and isinstance(info, dict):
                available = "✓"
                tag = installed_soft.get(name, "N/A")
                installed = "✓" if tag != "N/A" else "✗"
                self.soft_tree.insert("", "end", values=(name, "git", available, installed, tag))
    
    def _get_installed_databases(self):
        """Query installed databases (non-host) with version."""
        installed = {}
        if not os.path.exists(self.db_path):
            return installed
        
        try:
            conn = sqlite3.connect(self.db_path)
            cursor = conn.cursor()
            cursor.execute(
                "SELECT name, version, date, installed FROM database WHERE installed = 'True' AND software != 'host'"
            )
            for row in cursor.fetchall():
                installed[row[0]] = {
                    "version": row[1] if row[1] else "N/A",
                    "date": row[2] if row[2] else "N/A",
                    "installed": row[3] if row[3] else "N/A"
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


def main():
    root = tk.Tk()
    app = TelevirStatusApp(root)
    root.mainloop()


if __name__ == "__main__":
    main()
