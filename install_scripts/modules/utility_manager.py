# from sqlalchemy import create_engine, Column, Integer, String, Boolean, ForeignKey

import datetime
import os
from abc import abstractmethod

from sqlalchemy import (
    Boolean,
    Column,
    MetaData,
    String,
    Table,
    create_engine,
    text,
)


class SoftwareItem:
    def __init__(
        self, name, path, database, installed, env_path, tag: str = "undefined"
    ) -> None:
        self.name = name
        self.path = path
        self.database = database
        self.installed = installed
        self.env_path = env_path
        self.date = datetime.datetime.now().strftime("%Y-%m-%d")
        self.tag = tag

    def __repr__(self) -> str:
        return f"({self.name}, {self.path}, {self.database}, {self.installed}, {self.env_path})"


class DatabaseItem:
    def __init__(self, name, path, installed, software: str = "none") -> None:
        self.name = name
        self.path = path
        self.installed = installed
        self.software = software
        self.date = datetime.datetime.now().strftime("%Y-%m-%d")

    def __repr__(self) -> str:
        return f"({self.name}, {self.path}, {self.installed})"


class Utility_Repository:
    database_item = DatabaseItem
    software_item = SoftwareItem
    dbtype_local: str = "sqlite"

    tables: list = ["software", "database"]

    def __init__(self, db_path="", install_type="local", file_prefix="utility") -> None:
        self.db_path = db_path
        self.metadata = MetaData()
        self.engine_name = f"{file_prefix}_{install_type}.db"
        self.engine_filepath = os.path.join(*self.db_path.split("/"), self.engine_name)

        self.setup_engine(install_type)

        self.create_tables()

    def clear_existing_repo(self):
        """
        Delete the database
        """

        for table in self.tables:
            self.engine_execute(f"DROP TABLE {table}")
        # self.metadata.drop_all(self.engine)
        self.create_tables()

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
            "software",
            self.metadata,
            Column("name", String),
            Column("path", String),
            Column("database", String),
            Column("installed", Boolean),
            Column("tag", String, default="undefined"),
            Column("env_path", String),
            Column("date", String),
        )

        self.engine_execute(
            "CREATE TABLE IF NOT EXISTS software (name TEXT, path TEXT, database TEXT, installed BOOLEAN, tag TEXT, env_path TEXT, date TEXT)"
        )

    def create_database_table(self):
        self.database = Table(
            "database",
            self.metadata,
            Column("name", String),
            Column("path", String),
            Column("installed", Boolean),
            Column("software", String),
            Column("date", String),
        )

        self.engine_execute(
            "CREATE TABLE IF NOT EXISTS database (name TEXT, path TEXT, installed BOOLEAN, software TEXT, date TEXT)"
        )

    def delete_tables(self):
        self.delete_table("software")
        self.delete_table("database")

    def delete_table(self, table_name):
        self.engine_execute(f"DROP TABLE {table_name}")

    def clear_tables(self):
        self.clear_table("software")
        self.clear_table("database")

    def clear_table(self, table_name):
        self.engine_execute(f"DELETE FROM {table_name}")

    def print_table_schema(self, table_name):
        print(self.engine_execute(f"PRAGMA table_info({table_name})").fetchall())

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

            conn.commit()

        return result

    def create_tables(self):
        """
        Create the tables
        """
        self.create_software_table()
        self.create_database_table()

        self.metadata.create_all(self.engine)

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

        table_rows = self.engine_execute(f"SELECT * FROM {table_name}")

        with open(os.path.join(directory, f"{table_name}.tsv"), "w") as f:
            for row in table_rows:
                f.write("\t".join([str(x) for x in row]) + "\n")

    def get(self, table_name, id):
        """
        Get a record by id from a table
        """

        return self.engine_execute(f"SELECT * FROM {table_name} WHERE name='{id}'")

    def check_exists(self, table_name, id):
        """
        Check if a record exists in a table
        """

        find = self.engine_execute(
            f"SELECT * FROM {table_name} WHERE name='{id}'"
        ).fetchall()
        find = len(find) > 0
        if find:
            return True
        else:
            return False

    @abstractmethod
    def add_software(self, item: SoftwareItem):
        """
        Add a record to a table
        """
        # print("adding software")

        try:
            _ = self.engine_execute(
                f"INSERT INTO software (name, path, database, installed, tag, env_path, date) VALUES ('{item.name}', '{item.path}', '{item.database}', '{item.installed}', '{item.tag}', '{item.env_path}', '{item.date}')"
            )

        except Exception as e:
            print(e)
            print(
                "error adding software: delete currently existing utility_docker.db and re-run the script"
            )

    @abstractmethod
    def add_database(self, item: database_item):
        """
        Add a record to a table
        """
        try:
            _ = self.engine_execute(
                f"INSERT INTO database (name, path, installed, date) VALUES ('{item.name}', '{item.path}', '{item.installed}', '{item.date}')"
            )
        except Exception as e:
            print(e)
            print(
                "error adding database: delete currently existing utility_docker.db and re-run the script"
            )
