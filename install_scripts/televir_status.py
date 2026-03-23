#!/usr/bin/env python3
"""
TELE-Vir Status GUI Application

Displays all available databases and software from sources.yaml
with their installation status from utility_local.db.
"""

import tkinter as tk
from tkinter import ttk
import os
import sys
import sqlite3
from pathlib import Path

sys.path.insert(0, os.path.dirname(__file__))

from load_sources import list_databases, list_software


class TelevirStatusApp:
    def __init__(self, root):
        self.root = root
        self.root.title("TELE-Vir Status")
        self.root.geometry("900x600")
        
        self.install_home = os.environ.get('INSTALL_HOME', '/opt/televir')
        self.db_path = os.path.join(self.install_home, 'utility_local.db')
        
        self._create_widgets()
    
    def _create_widgets(self):
        title = tk.Label(self.root, text="TELE-Vir Status", font=("Arial", 18, "bold"))
        title.pack(pady=10)
        
        db_frame = tk.LabelFrame(self.root, text="Databases", font=("Arial", 12, "bold"))
        db_frame.pack(fill="both", expand=True, padx=10, pady=5)
        
        self.db_tree = ttk.Treeview(
            db_frame, 
            columns=("name", "category", "available", "installed"), 
            show="headings"
        )
        self.db_tree.heading("name", text="Name")
        self.db_tree.heading("category", text="Category")
        self.db_tree.heading("available", text="Available")
        self.db_tree.heading("installed", text="Installed")
        self.db_tree.column("name", width=250)
        self.db_tree.column("category", width=150)
        self.db_tree.column("available", width=80)
        self.db_tree.column("installed", width=80)
        self.db_tree.pack(fill="both", expand=True, padx=5, pady=5)
        
        soft_frame = tk.LabelFrame(self.root, text="Software", font=("Arial", 12, "bold"))
        soft_frame.pack(fill="both", expand=True, padx=10, pady=5)
        
        self.soft_tree = ttk.Treeview(
            soft_frame, 
            columns=("name", "type", "available", "installed"), 
            show="headings"
        )
        self.soft_tree.heading("name", text="Name")
        self.soft_tree.heading("type", text="Type")
        self.soft_tree.heading("available", text="Available")
        self.soft_tree.heading("installed", text="Installed")
        self.soft_tree.column("name", width=250)
        self.soft_tree.column("type", width=150)
        self.soft_tree.column("available", width=80)
        self.soft_tree.column("installed", width=80)
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
                    installed = "✓" if name_full in installed_dbs else "✗"
                    self.db_tree.insert("", "end", values=(name_full, category, available, installed))
    
    def _populate_software(self):
        for item in self.soft_tree.get_children():
            self.soft_tree.delete(item)
        
        soft = list_software()
        installed_soft = self._get_installed_software()
        
        for name, info in soft.get('archives', {}).items():
            if info and isinstance(info, dict):
                available = "✓"
                installed = "✓" if name in installed_soft else "✗"
                self.soft_tree.insert("", "end", values=(name, "archive", available, installed))
        
        for name, info in soft.get('git_repos', {}).items():
            if info and isinstance(info, dict):
                available = "✓"
                installed = "✓" if name in installed_soft else "✗"
                self.soft_tree.insert("", "end", values=(name, "git", available, installed))
    
    def _get_installed_databases(self):
        installed = set()
        if not os.path.exists(self.db_path):
            return installed
        
        try:
            conn = sqlite3.connect(self.db_path)
            cursor = conn.cursor()
            cursor.execute("SELECT name FROM database WHERE installed = 1")
            for row in cursor.fetchall():
                installed.add(row[0])
            conn.close()
        except Exception as e:
            print(f"Error reading database: {e}")
        
        return installed
    
    def _get_installed_software(self):
        installed = set()
        if not os.path.exists(self.db_path):
            return installed
        
        try:
            conn = sqlite3.connect(self.db_path)
            cursor = conn.cursor()
            cursor.execute("SELECT name FROM software WHERE installed = 1")
            for row in cursor.fetchall():
                installed.add(row[0])
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
