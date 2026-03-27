
from typing import List
from install_scripts.modules.utility_manager import Utility_Repository
import pandas as pd
import logging
import os


class Software:
    TELEVIR_GLOBAL_TYPES = [
        "alignment",
        "variant_calling",
        "quantification",
        "annotation",
    ]

    def __init__(self, name, version):
        self.name = name
        self.version = version

class ConstantsSettings:
    PIPELINE_STEPS_DB_DEPENDENT = [
        "alignment",
        "variant_calling",
        "quantification",
        "annotation",
    ]

class CS:
    PIPELINE_NAME_remapping = "remapping"
    PIPELINE_NAME_read_quality_analysis = "read_quality_analysis"
    PIPELINE_NAME_extra_qc = "extra_qc"
    PIPELINE_NAME_assembly = "assembly"
    PIPELINE_NAME_contig_classification = "contig_classification"
    PIPELINE_NAME_viral_enrichment = "viral_enrichment"
    PIPELINE_NAME_read_classification = "read_classification"
    PIPELINE_NAME_remap_filtering = "remap_filtering"

class Televir_Directory_Constants:
    conda_directory = "/opt/televir/conda"
    environments_directory = "/opt/televir/environments"

class Televir_Metadata:

    BINARIES = {
        "SOURCE": Televir_Directory_Constants.conda_directory,
        "ROOT": Televir_Directory_Constants.environments_directory,
        "software": {
            CS.PIPELINE_NAME_contig_classification: "hostDepletion/hostdep_env",
            CS.PIPELINE_NAME_read_classification: "hostDepletion/hostdep_env",
            CS.PIPELINE_NAME_viral_enrichment: "hostDepletion/hostdep_env",
            "centrifuge": "hostDepletion/hostdep_env",
            "diamond": "hostDepletion/hostdep_env",
            "kaiju": "hostDepletion/hostdep_env",
            "krakenuniq": "hostDepletion/hostdep_env",
            "blast": "hostDepletion/hostdep_env",
            "kraken2": "kraken2/kraken_env",
            "fastviromeexplorer": "FastViromeExplorer/FastViromeExplorer",
            "kallisto": "FastViromeExplorer/fve",
            "java": "FastViromeExplorer/fve",
            "virsorter": "hostDepletion/vs2",
            "desamba": "classm_lc/deSAMBA",
            "minimap-rem": "hostDepletion/hostdep_env",
            "flye": "assembly/Flye",
            "clark": "classification/Clark",
            "blastn": "hostDepletion/hostdep_env",
            "blastp": "hostDepletion/hostdep_env",
            "snippy": "remap/remap",
            "bamutil": "remap/remap",
            "msamtools": "remap/remap",
            "bcftools": "remap/remap",
            "samtools": "remap/remap",
            "bwa": "remap/remap",
            "bwa-filter": "remap/remap",
            "bwa-mem": "remap/remap",
            "bowtie2": "remap/remap",
            "bowtie2_remap": "remap/remap",
            "minimap2": "hostDepletion/hostdep_env",
            "minimap2_asm": "hostDepletion/hostdep_env",
            "minimap2_illu": "hostDepletion/hostdep_env",
            "minimap2_remap": "hostDepletion/hostdep_env",
            "minimap2_ont": "hostDepletion/hostdep_env",
            "metaphlan": "Metaphlan/metaphlan",
            "voyager": "classification",
            "snippy_pi": config("DIR_SOFTWARE") + "/snippy",
            "prinseq++": "preprocess/prinseq",
            "prinseq": "preprocess/prinseq",
            "collapsibleTree": "remap/Renv",
            "create_report": "Pyenv/igv_reports",
            "entrez_direct": "entrez_direct",
            "bgzip": "hostDepletion/hostdep_env",
        },
        CS.PIPELINE_NAME_remapping: {"default": "remap/remap"},
        CS.PIPELINE_NAME_remap_filtering: {"default": "remap/remap"},
        CS.PIPELINE_NAME_read_quality_analysis: {"default": "preprocess/preproc"},
        CS.PIPELINE_NAME_extra_qc: {"default": "preprocess/preproc"},
        CS.PIPELINE_NAME_assembly: {"default": "assembly/assembly"},
    }



class Televir_Directories:
    docker_app_directory = "/opt/televir/docker"
    local_app_directory = "/opt/televir/local"
    docker_install_directory = "/opt/televir/docker/install"


class UtilityDB:
    def __init__(self, db_path: str, install_type: str):
        self.repository = Utility_Repository(
            db_path=db_path, install_type=install_type
        )

    def query_databases(self, category: str = None, db_type: str = None) -> pd.DataFrame:
        sql_parts = ["SELECT * FROM database"]
        params = []
        
        conditions = []
        if category:
            conditions.append(f"db_category = '{category}'")
        if db_type:
            conditions.append(f"db_type = '{db_type}'")
        
        if conditions:
            sql_parts.append("WHERE " + " AND ".join(conditions))
        
        sql = " ".join(sql_parts)
        return self._execute_query(sql)

    def get_software_dbs(self, category: str, db_type: str = None) -> pd.DataFrame:
        return self.query_databases(category=category.lower(), db_type=db_type)

    def get_host_dbs(self, category: str = None) -> pd.DataFrame:
        return self.query_databases(category=category.lower() if category else None, db_type="host")

    def get_filter_dbs(self, category: str) -> pd.DataFrame:
        return self.query_databases(category=category.lower(), db_type="filter")

    def get_all_databases(self) -> pd.DataFrame:
        return self.query_databases()

    def get_unique_categories(self) -> list:
        rows = self.repository.engine_execute_return_table(
            "SELECT DISTINCT db_category FROM database WHERE db_category IS NOT NULL"
        )
        return [r[0] for r in rows if r[0]]

    def get_unique_software(self) -> list:
        rows = self.repository.engine_execute_return_table(
            "SELECT DISTINCT name FROM software WHERE name IS NOT NULL ORDER BY name"
        )
        return [r[0] for r in rows if r[0]]

    def get_software_from_database(self) -> list:
        rows = self.repository.engine_execute_return_table(
            "SELECT DISTINCT software FROM database WHERE software IS NOT NULL ORDER BY software"
        )
        return [r[0] for r in rows if r[0]]

    def _execute_query(self, sql: str) -> pd.DataFrame:
        try:
            with self.repository.engine.connect() as conn:
                result = conn.execute(sql)
                rows = result.fetchall()
                if not rows:
                    return pd.DataFrame()
                columns = result.keys()
                df = pd.DataFrame(rows, columns=columns)
                return df
        except Exception as e:
            print(f"Query error: {e}")
            return pd.DataFrame()


class Utility_Pipeline_Manager:
    """
    Takes a combined table and generates a pipeline tree.
    Combined table is a table with the following columns:
    - sample_name
    - sample_id
    - pipeline_step
    - software_name
    - parameter
    - value

    Uitility_Pipeline_Manager Comunicates with utility repository to complement the information
    with software specific installed databases. Creates a pipeline tree from the combined information.
    """

    software_name_list: list
    existing_pipeline_order: list
    combined_table: pd.DataFrame
    software_dbs_dict: dict
    technology: str
    new_variable: int
    pipeline_order: list
    pipeline_makeup: int

    def __init__(self):
        self.utility_repository = Utility_Repository(
            db_path=Televir_Directories.docker_app_directory, install_type="docker"
        )
        self.utility_db = UtilityDB(
            db_path=Televir_Directories.docker_app_directory, install_type="docker"
        )

        self.steps_db_dependant = ConstantsSettings.PIPELINE_STEPS_DB_DEPENDENT
        self.binaries = Televir_Metadata.BINARIES

        self.logger = logging.getLogger(__name__)
        if self.logger.hasHandlers():
            self.logger.handlers.clear()
        self.logger.setLevel(logging.ERROR)
        self.logger.addHandler(logging.StreamHandler())
        self.host_dbs = {}
        self.filter_dbs = {}
        self.software_dbs_dict = {}

    def check_software_is_installed(self, software_name: str) -> bool:
        """
        Check if a software is installed
        """
        software_lower = software_name.lower()
        if software_lower in self.binaries["software"].keys():
            bin_path = os.path.join(
                Televir_Directories.docker_install_directory,
                self.binaries["software"][software_lower],
                "bin",
                software_lower,
            )
            return os.path.isfile(bin_path)
        else:
            for pipeline in [
                CS.PIPELINE_NAME_remapping,
                CS.PIPELINE_NAME_read_quality_analysis,
                CS.PIPELINE_NAME_extra_qc,
                CS.PIPELINE_NAME_assembly,
            ]:
                if os.path.exists(
                    os.path.join(
                        Televir_Directories.docker_install_directory,
                        self.binaries[pipeline]["default"],
                        "bin",
                        software_lower,
                    )
                ):
                    return True

        return False

    def normalize_name(self, name: str) -> str:
        return name.lower().split('/')[0]

    def set_software_list(self, software_list):
        self.software_name_list = software_list


    ##############################
    # Software DBs
    def get_software_dbs_if_exist(
        self, software_name: str, filters: List[tuple] = []
    ) -> pd.DataFrame:
        print(f"[DEBUG] get_software_dbs_if_exist: {software_name}, filters: {filters}")
        
        db_type_filter = None
        for col, val in filters:
            if col == "software" and val == "host":
                db_type_filter = "host"
            elif col == "tag" and val == "filter":
                db_type_filter = "filter"
        
        df = self.utility_db.get_software_dbs(
            category=self.normalize_name(software_name),
            db_type=db_type_filter
        )
        
        print(f"[DEBUG] Found {len(df)} databases")
        print(df.columns.tolist() if not df.empty else "Empty")
        
        return df

    def check_tables_exist(self):
        """
        Check if the software table exist
        """
        return self.utility_repository.check_tables_exists()

    def check_software_db_available(self, software_name: str) -> bool:
        """
        Check if a software is installed in the database.
        """
        return self.utility_repository.check_exists("software", software_name.lower())

    def get_software_db_dict(self):
        software_list = self.utility_db.get_unique_software()
        print(f"[DEBUG] Unique software: {software_list}")

        self.software_dbs_dict = {
            software: self.utility_db.query_databases(category=software.lower())['path'].unique().tolist()
            for software in software_list
        }

        print(f"[DEBUG] Software DBs dict: {self.software_dbs_dict}")

    def get_host_dbs(self):
        all_host_dbs = self.utility_db.query_databases(db_type="host")
        
        hosts_dbs_dict = {}
        if not all_host_dbs.empty and 'db_category' in all_host_dbs.columns:
            for category in all_host_dbs['db_category'].unique():
                category_df = all_host_dbs[all_host_dbs['db_category'] == category]
                if len(category_df) > 0:
                    hosts_dbs_dict[category] = category_df

        try:
            import numpy as np
            for software in hosts_dbs_dict.keys():
                if 'db_name' in hosts_dbs_dict[software].columns:
                    hosts_dbs_dict[software] = hosts_dbs_dict[software].rename(
                        columns={'db_name': 'database'}
                    )
                if 'host_name' not in hosts_dbs_dict[software].columns:
                    hosts_dbs_dict[software]["host_name"] = np.nan
                    hosts_dbs_dict[software]["host_filename"] = hosts_dbs_dict[software].get("database", "")
                    hosts_dbs_dict[software]["file_str"] = hosts_dbs_dict[software].get("database", "")
        except ImportError:
            pass

        self.host_dbs = hosts_dbs_dict

    def get_filter_dbs(self):
        all_filter_dbs = self.utility_db.query_databases(db_type="filter")

        filter_dbs_dict = {}
        if not all_filter_dbs.empty and 'db_category' in all_filter_dbs.columns:
            for category in all_filter_dbs['db_category'].unique():
                category_df = all_filter_dbs[all_filter_dbs['db_category'] == category]
                if len(category_df) > 0:
                    filter_dbs_dict[category] = category_df

        filter_dbs_dict = {k: v for k, v in filter_dbs_dict.items() if len(v) > 0}
        for sof, filter_df in filter_dbs_dict.items():
            if 'path' in filter_df.columns:
                filter_df = filter_df.copy()
                filter_df["file_str"] = filter_df.apply(
                    lambda x: os.path.basename(x.path) if pd.notna(x.path) else "", axis=1
                )
                filter_dbs_dict[sof] = filter_df
        self.filter_dbs = filter_dbs_dict

    ##################################
    ##################################

    def get_from_software_db_dict(self, software_name: str, empty=None):
        if empty is None:
            empty = []
        
        possibilities = self._get_name_possibilities(software_name)

        for possibility in possibilities:
            if possibility in self.software_dbs_dict.keys():
                return self.software_dbs_dict[possibility]

        return empty

    def get_from_host_db(self, software_name: str, empty=None):
        if empty is None:
            empty = []
        possibilities = self._get_name_possibilities(software_name)

        for possibility in possibilities:
            if possibility in self.host_dbs.keys():
                host_df = self.host_dbs[possibility]
                if 'host_name' in host_df.columns:
                    try:
                        human_reference = HomoSapiens()
                        if human_reference.host_name in host_df.host_name.unique():
                            host_df = host_df[host_df.host_name == human_reference.host_name]
                            host_df = pd.concat(
                                [host_df, self.host_dbs[possibility].drop(host_df.index)]
                            )
                    except (NameError, AttributeError):
                        pass

                if 'path' in host_df.columns and 'file_str' in host_df.columns:
                    return list(
                        host_df[["path", "file_str"]].itertuples(index=False, name=None)
                    )

        return empty

    def get_from_filter_dbs(self, software_name: str, empty=None):
        if empty is None:
            empty = ["None"]
        possibilities = self._get_name_possibilities(software_name)

        for possibility in possibilities:
            if possibility in self.filter_dbs.keys():
                filter_df = self.filter_dbs[possibility]
                if 'path' in filter_df.columns and 'file_str' in filter_df.columns:
                    return list(
                        filter_df[["path", "file_str"]].itertuples(index=False, name=None)
                    )
        
        return ["None"]

    def _get_name_possibilities(self, software_name: str) -> list:
        possibilities = [software_name, software_name.lower()]
        if "_" in software_name:
            element = software_name.split("_")[0]
            possibilities.append(element)
            possibilities.append(element.lower())
        return possibilities

