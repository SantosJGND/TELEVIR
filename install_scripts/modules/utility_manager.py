# from sqlalchemy import create_engine, Column, Integer, String, Boolean, ForeignKey

import datetime
import logging
import os
from abc import abstractmethod
from typing import List, Optional
import sys 
sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))

from install_scripts.config import DATABASE_FILENAME

from sqlalchemy import (
    Boolean,
    Column,
    Integer,
    MetaData,
    String,
    Table,
    create_engine,
    text,
)


class SoftwareItem:
    def __init__(
        self, name, path, database, installed, env_path, 
        tag: str = "undefined", db_version: Optional[str] = None, 
        needs_update: bool = False, binary_name: Optional[str] = None
    ) -> None:
        self.name = name
        self.path = path
        self.database = database
        self.installed = installed
        self.env_path = env_path
        self.date = datetime.datetime.now().strftime("%Y-%m-%d")
        self.tag = tag
        self.db_version = db_version
        self.needs_update = needs_update
        self.binary_name = binary_name

    def __repr__(self) -> str:
        return f"({self.name}, {self.path}, {self.database}, {self.installed}, {self.env_path})"


class DatabaseItem:
    def __init__(self, name, path, installed, software: str = "none",
                 version: Optional[str] = None, source_url: Optional[str] = None, 
                 file_mod_date: Optional[str] = None, description: Optional[str] = None,
                 db_category: Optional[str] = None, db_name: Optional[str] = None,
                 db_type: Optional[str] = None) -> None:
        self.name = name
        self.path = path
        self.installed = installed
        self.software = software
        self.date = datetime.datetime.now().strftime("%Y-%m-%d")
        self.version = version
        self.source_url = source_url
        self.file_mod_date = file_mod_date
        self.description = description
        self._parse_name_fields(name, db_category, db_name, db_type)

    def _parse_name_fields(self, name: str, db_category: Optional[str], 
                           db_name: Optional[str], db_type: Optional[str]):
        if '/' in name and db_category is None:
            parts = name.split('/', 1)
            self.db_category = parts[0]
            self.db_name = parts[1]
        else:
            self.db_name = db_name if db_name else name
            self.db_category = db_category if db_category else self.software
        self.db_type = db_type if db_type else self.software

    def __repr__(self) -> str:
        return f"({self.name}, {self.path}, {self.installed})"


class Utility_Repository:
    database_item = DatabaseItem
    software_item = SoftwareItem
    dbtype_local: str = "sqlite"
    SOFTWARE_TABLE_NAME: str = "software"
    DATABASE_TABLE_NAME: str = "database"

    tables: list = [SOFTWARE_TABLE_NAME, DATABASE_TABLE_NAME]

    def __init__(self, db_path="", install_type="local", file_prefix="utility") -> None:
        print(f"Initializing Utility Repository: {db_path}, {install_type}, {file_prefix}")
        self.db_path = db_path
        self.metadata = MetaData()
        self.engine_name = DATABASE_FILENAME
        self.engine_filepath = os.path.join(*self.db_path.split("/"), self.engine_name)

        self.setup_engine(install_type)

        self.create_tables()

    def clear_existing_repo(self):
        """
        Delete the database
        """

        self.metadata.drop_all(self.engine)
        self.metadata.create_all(self.engine)

    def setup_engine(self, install_type):
        """
        setup the engine
        """

        # if os.path.exists(self.engine_filepath):
        #    os.remove(self.engine_filepath)

        self.engine = create_engine(f"{self.dbtype_local}:////" + self.engine_filepath)
        # if install_type == "local":
        #    self.setup_engine_local()
        # elif install_type == "docker":
        #    self.setup_engine_docker()

    def setup_engine_postgres(self):
        from decouple import config

        self.engine = create_engine(
            f"postgresql+psycopg2://{config('DB_USER')}:{config('DB_PASSWORD')}@{config('DB_HOST')}:{config('DB_PORT')}/{config('DB_NAME')}"
        )

    def create_software_table(self):
        self.software = Table(
            self.SOFTWARE_TABLE_NAME,
            self.metadata,
            Column("name", String, primary_key=True),
            Column("path", String),
            Column("database", String),
            Column("installed", Boolean),
            Column("tag", String, default="undefined"),
            Column("env_path", String),
            Column("date", String),
            Column("db_version", String),
            Column("needs_update", Boolean),
        )

        self.engine_execute(
            "CREATE TABLE IF NOT EXISTS software (name TEXT PRIMARY KEY, path TEXT, database TEXT, installed BOOLEAN, tag TEXT, env_path TEXT, date TEXT, db_version TEXT, needs_update BOOLEAN)"
        )

    def create_database_table(self):
        self.database = Table(
            self.DATABASE_TABLE_NAME,
            self.metadata,
            Column("id", Integer, primary_key=True, autoincrement=True),
            Column("name", String),
            Column("db_category", String),
            Column("db_name", String),
            Column("db_type", String),
            Column("path", String),
            Column("installed", Boolean),
            Column("software", String),
            Column("date", String),
            Column("version", String),
            Column("source_url", String),
            Column("file_mod_date", String),
            Column("description", String),
        )

        self.engine_execute(
            "CREATE TABLE IF NOT EXISTS database (id INTEGER PRIMARY KEY AUTOINCREMENT, name TEXT, db_category TEXT, db_name TEXT, db_type TEXT, path TEXT, installed BOOLEAN, software TEXT, date TEXT, version TEXT, source_url TEXT, file_mod_date TEXT, description TEXT)"
        )

    def migrate_database_table(self):
        self.engine_execute(
            "ALTER TABLE database ADD COLUMN db_category TEXT"
        )
        self.engine_execute(
            "ALTER TABLE database ADD COLUMN db_name TEXT"
        )
        self.engine_execute(
            "ALTER TABLE database ADD COLUMN db_type TEXT"
        )
        rows = self.engine_execute_return_table("SELECT name, software FROM database")
        for row in rows:
            name, software = row
            if '/' in name:
                parts = name.split('/', 1)
                db_category, db_name = parts[0], parts[1]
            else:
                db_category = software if software else name
                db_name = name
            db_type = software if software else name
            self.engine_execute(
                f"UPDATE database SET db_category='{db_category}', db_name='{db_name}', db_type='{db_type}' WHERE name='{name}'"
            )
        print("[DEBUG] Database migration completed: added db_category, db_name, db_type columns")

    def delete_tables(self):
        for table in self.tables:
            self.engine_execute(f"DROP TABLE {table}")

    def delete_table(self, table_name):
        self.engine_execute(f"DROP TABLE {table_name}")

    def clear_tables(self):
        self.clear_table(self.SOFTWARE_TABLE_NAME)
        self.clear_table(self.DATABASE_TABLE_NAME)

    def clear_table(self, table_name):
        self.engine_execute(f"DELETE FROM {table_name}")

    def print_table_schema(self, table_name: str):
        result = self.engine_execute_return_table(f"PRAGMA table_info({table_name})")
        print(result)

    def reset_tables(self):
        """
        Create the tables
        """
        self.clear_tables()

        self.metadata.create_all(self.engine)
    
    def engine_execute(self, string: str):
        sql = text(string)

        with self.engine.connect() as conn:
            result = conn.execute(sql)
            #conn.commit()

        return result

    def engine_execute_return_table(self, string: str):
        sql = text(string)

        rows = None

        with self.engine.connect() as conn:
            result = conn.execute(sql)
            #conn.commit()
            rows = result.fetchall()

        return rows

    def create_tables(self):
        """
        Create the tables
        """
        self.create_software_table()
        self.create_database_table()

        self.metadata.create_all(self.engine)


    def check_table_exists(self, table_name):
        """
        Check if a table exists in the database
        """
        find = self.engine.execute(
            f"SELECT name FROM sqlite_master WHERE type='table' AND name='{table_name}'"
        ).fetchall()

        if len(find) > 0:
            return True
        else:
            return False

    def check_tables_exists(self):
        """
        Check if the tables exist in the database
        """
        for table_name in self.tables:
            if not self.check_table_exists(table_name):
                return False
        return True


    def dump_software(self, directory: str):
        """
        Dump the software table to a tsv file
        """
        self.dump_table_tsv("software", directory)

    def dump_database(self, directory: str):
        """
        Dump the database table to a tsv file
        """

        self.dump_table_tsv("database", directory)

    def dump_table_tsv(self, table_name: str, directory: str):
        """
        Dump a table to a tsv file
        """

        if table_name not in self.tables:
            print(f"Table {table_name} not found. Available tables: {self.tables}")
            return

        if not os.path.exists(directory):
            os.makedirs(directory, exist_ok=True)

        table_rows = self.engine_execute_return_table(f"SELECT * FROM {table_name}")

        with open(os.path.join(directory, f"{table_name}.tsv"), "w") as f:
            for row in table_rows:
                f.write("\t".join([str(x) for x in row]) + "\n")


    def get_by_name(self, table_name, id):
        """
        Get a record by id from a table
        """

        return self.engine.execute(f"SELECT * FROM {table_name} WHERE name='{id}'")


    def get_list_tables(self):
        """
        Get a list of tables
        """

        find = self.engine.execute("SELECT name FROM sqlite_master WHERE type='table'")
        find = [i[0] for i in find]
        return find

    def get_list_unique_field(self, table_name, id):
        """
        Get a list of unique values in a field
        """

        find = self.engine.execute(f"SELECT DISTINCT {id} FROM {table_name}")

        find = [i[0] for i in find]
        return find

    def select_explicit_statement(
        self, table_name, field, id, filters: List[tuple] = []
    ):
        """
        select from table.
        """
        sql_statement = f"SELECT * FROM {table_name} WHERE {field}='{id}'"
        for filter in filters:
            column_name, value = filter

            if self.check_column_exists(table_name, column_name) is False:
                continue
            if value is None:
                continue

            sql_statement += f" AND {column_name}='{value}'"

        return sql_statement

    def check_column_exists(self, table_name, column_name):
        """
        Check if a column exists in a table
        """

        from sqlalchemy import inspect

        inspector = inspect(self.engine)
        columns = inspector.get_columns(table_name)
        find = any([i["name"] == column_name for i in columns])

        if find:
            return True
        else:
            return False

    def check_exists(self, table_name: str, id: str):
        """
        Check if a record exists in a table
        """

        check_list = [id]
        if "_" in id:
            check_list.append(id.split("_")[0])
        check_list = [f"'{i}'" for i in check_list]
        check_list = ",".join(check_list)

        find = self.engine_execute_return_table(
            f"SELECT * FROM {table_name} WHERE name='{id}'"
        )
        find = len(find) > 0
        if find:
            return True
        else:
            return False

    @abstractmethod
    def add_software(self, item: SoftwareItem):
        """
        Add a record to a table (uses INSERT OR REPLACE for upsert behavior)
        """
        db_version = item.db_version if item.db_version else "NULL"
        needs_update = 1 if item.needs_update else 0

        actual_installed = item.installed
        if item.binary_name:
            binary_path = os.path.join(item.env_path, "bin", item.binary_name)
            if not os.path.isfile(binary_path):
                logging.warning(f"Binary not found at {binary_path}, marking as not installed")
                actual_installed = False

        try:
            _ = self.engine_execute(
                f"INSERT OR REPLACE INTO software (name, path, database, installed, tag, env_path, date, db_version, needs_update) VALUES ('{item.name}', '{item.path}', '{item.database}', '{actual_installed}', '{item.tag}', '{item.env_path}', '{item.date}', {db_version}, {needs_update})"
            )

            verify = self.engine_execute_return_table(f"SELECT name FROM software WHERE name='{item.name}'")
            if verify:
                print(verify)
                print(f"[DEBUG] Verified software exists in DB: {item.name}")
            else:
                print(f"[DEBUG] WARNING: Software not found after insert: {item.name}")

        except Exception as e:
            print(e)
            print(
                "error adding software: delete currently existing utility_docker.db and re-run the script"
            )

    @abstractmethod
    def add_database(self, item: DatabaseItem):
        """
        Add a record to a table (uses INSERT OR REPLACE for upsert behavior)
        New schema includes db_category, db_name, db_type columns
        """
        version = f"'{item.version}'" if item.version else "NULL"
        source_url = f"'{item.source_url}'" if item.source_url else "NULL"
        file_mod_date = f"'{item.file_mod_date}'" if item.file_mod_date else "NULL"
        description = f"'{item.description}'" if item.description else "NULL"
        db_category = f"'{item.db_category}'" if item.db_category else f"'{item.name}'"
        db_name = f"'{item.db_name}'" if item.db_name else f"'{item.name}'"
        db_type = f"'{item.db_type}'" if item.db_type else f"'{item.software}'"
        print(f"Adding database: {item.name}, db_category: {item.db_category}, db_type: {item.db_type}")
        try:
            _ = self.engine_execute(
                f"INSERT OR REPLACE INTO database (name, db_category, db_name, db_type, path, installed, software, date, version, source_url, file_mod_date, description) VALUES ('{item.name}', {db_category}, {db_name}, {db_type}, '{item.path}', '{item.installed}', '{item.software}', '{item.date}', {version}, {source_url}, {file_mod_date}, {description})"
            )
            
            verify = self.engine_execute_return_table(f"SELECT name FROM database WHERE name='{item.name}'")
                
        except Exception as e:
            print(e)
            print(
                "error adding database: delete currently existing utility_docker.db and re-run the script"
            )

    def update_software(self, name: str, db_version: Optional[str] = None, needs_update: Optional[bool] = None):
        """
        Update software record by name
        """
        set_clauses = []
        if db_version is not None:
            set_clauses.append(f"db_version = '{db_version}'")
        if needs_update is not None:
            set_clauses.append(f"needs_update = {1 if needs_update else 0}")

        if set_clauses:
            try:
                self.engine_execute(
                    f"UPDATE software SET {', '.join(set_clauses)} WHERE name = '{name}'"
                )
            except Exception as e:
                print(e)

    def update_database(self, name: str, version: Optional[str] = None, 
                        source_url: Optional[str] = None, file_mod_date: Optional[str] = None):
        """
        Update database record by name
        """
        set_clauses = []
        if version is not None:
            set_clauses.append(f"version = '{version}'")
        if source_url is not None:
            set_clauses.append(f"source_url = '{source_url}'")
        if file_mod_date is not None:
            set_clauses.append(f"file_mod_date = '{file_mod_date}'")

        if set_clauses:
            try:
                self.engine_execute(
                    f"UPDATE database SET {', '.join(set_clauses)} WHERE name = '{name}'"
                )
            except Exception as e:
                print(e)
