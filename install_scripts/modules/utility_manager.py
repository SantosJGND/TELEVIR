# from sqlalchemy import create_engine, Column, Integer, String, Boolean, ForeignKey

from abc import ABC, abstractmethod

import sqlalchemy
from sqlalchemy import (
    Boolean,
    Column,
    ForeignKey,
    Integer,
    MetaData,
    String,
    Table,
    create_engine,
)
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import mapper, sessionmaker


class software_item:
    def __init__(self, name, path, database, installed, env_path) -> None:
        self.name = name
        self.path = path
        self.database = database
        self.installed = installed
        self.env_path = env_path

    def __repr__(self) -> str:
        return f"({self.name}, {self.path}, {self.database}, {self.installed}, {self.env_path})"


class database_item:
    def __init__(self, name, path, installed) -> None:
        self.name = name
        self.path = path
        self.installed = installed

    def __repr__(self) -> str:
        return f"({self.name}, {self.path}, {self.installed})"


class Utility_Repository:

    database_item = database_item
    software_item = software_item

    def __init__(self, db_path="", dbtype="sqlite") -> None:
        self.engine = create_engine(f"{dbtype}:///{db_path}/utility.db")
        self.metadata = MetaData()
        self.create_tables()

    def create_software_table(self):

        self.software = Table(
            "software",
            self.metadata,
            Column("name", String),
            Column("path", String),
            Column("dabatabse", String),
            Column("installed", Boolean),
            Column("env_path", String),
        )

    def create_database_table(self):
        self.database = Table(
            "database",
            self.metadata,
            Column("name", String),
            Column("path", String),
            Column("installed", Boolean),
        )

    def create_tables(self):
        self.create_software_table()
        self.create_database_table()
        self.metadata.create_all(self.engine)

    def get(self, table_name, id):
        """
        Get a record by id from a table
        """

        return self.engine.execute(f"SELECT * FROM {table_name} WHERE name='{id}'")

    def check_exists(self, table_name, id):
        """
        Check if a record exists in a table
        """

        find = self.engine.execute(f"SELECT * FROM {table_name} WHERE name='{id}'")
        if find:
            return True
        else:
            return False

    def add_software(self, item):
        """
        Add a record to a table
        """
        if not self.check_exists("software", item.name):
            self.engine.execute(
                f"INSERT INTO software (name, path, database, installed, env_path) VALUES {item}"
            )

    @abstractmethod
    def add_database(self, item):
        """
        Add a record to a table
        """

        if not self.check_exists("database", item.name):
            self.engine.execute(
                f"INSERT INTO database (name, path, installed) VALUES {item}"
            )
