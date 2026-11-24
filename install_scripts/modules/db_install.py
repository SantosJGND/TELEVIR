#!/usr/bin/python3

import datetime
import gzip
import logging
import os
import shutil
import subprocess
from contextlib import closing
from ftplib import FTP
from threading import Thread
import sys
import pandas as pd
from install_scripts.host_library import Host
from install_scripts.load_sources import get_db_url, get_db_version, get_db_file
from install_scripts.modules.accession_store import AccessionStore

from fastq_filter import file_to_fastq_records

from xopen import xopen

from decouple import config

BGZIP_BIN = "bgzip"
try:
    BGZIP_BIN = config("BGZIP_BIN")
except:
    pass


def compress_using_xopen(fq_in: str, fq_out: str):
    """
    compress using fastq_filter generator"""

    records = file_to_fastq_records(fq_in)
    with xopen(fq_out, mode="wb", threads=0, compresslevel=2) as output_h:
        for record in records:
            header = ">" + record.name + "\n"
            header = header.encode("ascii")
            output_h.write(header)
            sequence = record.sequence + "\n"
            output_h.write(sequence.encode("ascii"))


def entrez_fetch_sequence(accid, outfile):
    """return fasta from ncbi nuccore db using accid"""

    os.system(f"esearch -db nuccore -query {accid} | efetch -format fasta >> {outfile}")


def verify_file_accessible(filepath: str) -> bool:
    """
    Check that a file exists, is readable, and has non-zero size.
    Returns True if file is valid, False otherwise.
    """
    if not os.path.isfile(filepath):
        logging.error(f"File not found: {filepath}")
        return False
    if os.path.getsize(filepath) == 0:
        logging.error(f"File is empty: {filepath}")
        return False
    if not os.access(filepath, os.R_OK):
        logging.error(f"File is not readable: {filepath}")
        return False
    return True


def verify_gzip_integrity(filepath: str) -> bool:
    """
    Verify gzip file integrity by attempting to read the gzip header and
    checking the file is not truncated.
    Returns True if gzip file is complete, False if corrupted or incomplete.
    """
    if not verify_file_accessible(filepath):
        return False

    try:
        subprocess.run(
            ["gzip", "-t", filepath],
            capture_output=True,
            check=True,
        )
        return True
    except subprocess.CalledProcessError:
        logging.error(
            f"Corrupted or incomplete gzip file (unexpected end-of-file): {filepath}"
        )
        return False
    except Exception as e:
        logging.error(f"Failed to verify gzip integrity for {filepath}: {e}")
        return False


def verify_xz_integrity(filepath: str) -> bool:
    """
    Verify xz file integrity by attempting to read the xz header.
    Returns True if xz file is complete, False if corrupted or incomplete.
    """
    if not verify_file_accessible(filepath):
        return False

    try:
        subprocess.run(
            ["xz", "-t", filepath],
            capture_output=True,
            check=True,
        )
        return True
    except subprocess.CalledProcessError:
        logging.error(
            f"Corrupted or incomplete xz file (unexpected end-of-file): {filepath}"
        )
        return False
    except Exception as e:
        logging.error(f"Failed to verify xz integrity for {filepath}: {e}")
        return False


def verify_fasta_integrity(filepath: str) -> bool:
    """
    Verify FASTA file integrity by checking it can be read and has valid headers.
    Returns True if FASTA is valid, False otherwise.
    """
    if not verify_file_accessible(filepath):
        return False

    try:
        result = subprocess.run(
            ["grep", "-c", "^>", filepath],
            capture_output=True,
            text=True,
            check=True,
        )
        seq_count = int(result.stdout.strip())
        if seq_count == 0:
            logging.error(f"FASTA file has no sequences: {filepath}")
            return False
        logging.info(f"FASTA file verified: {seq_count} sequences found")
        return True
    except subprocess.CalledProcessError:
        logging.error(f"Failed to verify FASTA integrity: {filepath}")
        return False


class setup_dl:
    # Shared configuration
    BATCH_SIZE = 100000
    
    # Database filenames
    PROTEIN_DB = "protein_accessions.db"
    NUCLEOTIDE_DB = "nucleotide_accessions.db"
    
    # NCBI download directory names
    PROT_ACC2TAX_DIR = "prot.accession2taxid"
    NUCL_ACC2TAX_DIR = "nucl.accession2taxid"
    
    def __init__(
        self,
        INSTALL_PARAMS,
        organism="viral",
        home="",
        bindir="",
        test=False,
        update=False,
    ):
        if not INSTALL_PARAMS["HOME"]:
            home = os.getcwd()

        else:
            home = INSTALL_PARAMS["HOME"]

        if home[-1] != "/":
            home += "/"

        if not len(bindir):
            bindir = home + "scripts/"

        self.dbdir = home + "ref_db/"
        self.seqdir = home + "ref_fasta/"
        self.metadir = home + "metadata/"
        self.accessions = AccessionStore(self.metadir)

        self.envs = INSTALL_PARAMS["ENVSDIR"]
        self.source = self.envs["SOURCE"]
        self.requests = INSTALL_PARAMS["REQUEST_REFERENCES"]
        self.home = home
        self.bindr = bindir
        self.fastas = {"prot": {}, "nuc": {}, "host": {}, "filter": {}}
        self.meta = {}
        self.db_versions = {}
        self.test = test
        self.update = update
        self.batch_size = self.BATCH_SIZE

        self.organism = organism

    @property
    def protein_db_path(self):
        return os.path.join(self.metadir, self.PROTEIN_DB)
    
    @property
    def nucleotide_db_path(self):
        return os.path.join(self.metadir, self.NUCLEOTIDE_DB)
    
    @property
    def prot_acc2tax_dir(self):
        return os.path.join(self.metadir, self.PROT_ACC2TAX_DIR)
    
    @property
    def nucl_acc2tax_dir(self):
        return os.path.join(self.metadir, self.NUCL_ACC2TAX_DIR)

    def get_file_mod_date(self, filepath: str):
        """Get file modification date as YYYY-MM-DD"""
        if os.path.isfile(filepath):
            timestamp = os.path.getmtime(filepath)
            return datetime.datetime.fromtimestamp(timestamp).strftime("%Y-%m-%d")
        return ""

    def mkdirs(self):
        if not os.path.isdir(self.home):
            os.mkdir(self.home)
        for dr in self.dbdir, self.seqdir, self.metadir:
            if not os.path.isdir(dr):
                os.mkdir(dr)

    def verify_file_integrity(self, filepath: str, description: str = "") -> bool:
        """
        Check if file exists and verify integrity.
        If corrupted, delete and return False so caller re-downloads.
        Returns True if file is valid, False if missing or corrupted.
        """
        if not os.path.isfile(filepath):
            return False

        if filepath.endswith(".gz"):
            if verify_gzip_integrity(filepath):
                return True
            else:
                logging.warning(f"Removing corrupted {description}: {filepath}")
                os.remove(filepath)
                return False
        elif filepath.endswith(".xz"):
            if verify_xz_integrity(filepath):
                return True
            else:
                logging.warning(f"Removing corrupted {description}: {filepath}")
                os.remove(filepath)
                return False
        elif filepath.endswith(".tar.gz") or filepath.endswith(".tgz"):
            if verify_gzip_integrity(filepath):
                return True
            else:
                logging.warning(f"Removing corrupted {description}: {filepath}")
                os.remove(filepath)
                return False
        return True

    @staticmethod
    def check_fasta_bgziped(fasta_path: str):
        """
        use samtools faidx to check if fasta is bgzipped.
        check return of running samtools faidx on fasta file. if return is 0, file is bgzipped.
        if return is 1, file is not bgzipped. unzipped file and bgzip it.
        """
        if not os.path.isfile(fasta_path):
            return False
        if not os.path.isfile(fasta_path + ".fai"):
            return False
        try:
            subprocess.run(
                [
                    "samtools",
                    "faidx",
                    fasta_path,
                ]
            )
            return True
        except subprocess.CalledProcessError:
            return False

    def bgzip_file(self, filename):
        """
        bgzip file.
        :param filename:
        :return:
        """
        flname = os.path.basename(filename)
        basename = os.path.splitext(filename)[0]
        blname = os.path.basename(basename)

        if os.path.isfile(basename) and os.path.isfile(filename):
            os.remove(basename)

        file_is_bgzipped = False
        if os.path.isfile(filename):
            file_is_bgzipped = self.check_fasta_bgziped(filename)

        if not file_is_bgzipped:
            logging.info(f"gunzipping {flname}")
            subprocess.run(["gunzip", filename])

            if os.path.isfile(basename):
                if os.path.isfile(filename):
                    os.remove(filename)

            logging.info(f"bgzipping {flname}")
            subprocess.run([BGZIP_BIN, basename])

        return basename + ".gz"

    def index_nuc_fasta_files(self):
        """
        index fasta files.
        :return:
        """
        logging.info("Checking nuc fasta files for index ")
        for k, v in self.fastas["nuc"].items():

            for fl in v:
                flname = os.path.basename(fl)
                if not self.check_fasta_bgziped(fl):
                    logging.info(f"{flname} not bgzipped.")
                    self.bgzip_file(fl)
                    logging.info(f"indexing {flname} using samtools faidx")
                    subprocess.run(["samtools", "faidx", fl])

                    logging.info(f"{flname} indexed.")
                else:
                    logging.info(f"{flname} already bgzipped and indexed.")

    def dl_filter_file(self, dl_entry: dict, fname: str = None, db_key: str = "filter"):
        """
        download generic filter file
        
        Args:
            dl_entry: dict with url, filename (or file), description
            fname: optional name for fastas dict key (defaults to filename without extension)
            db_key: which fastas dict to use ("filter", "prot", "nuc")
        """
        filename = dl_entry.get("filename") or dl_entry.get("file", "")
        url = dl_entry.get("url", "")
        version = dl_entry.get("version") or dl_entry.get("version_default", "")
        
        if not filename or not url:
            logging.error("Invalid download entry.")
            return False

        if "{version}" in url and version:
            url = url.replace("{version}", version)
            filename = filename.replace("{version}", version)

        source_url = url
        filepath = self.seqdir + filename
        
        if fname is None:
            fname = os.path.splitext(filename)[0]

        if os.path.isfile(filepath):
            if self.verify_file_integrity(filepath, filename):
                self.fastas[db_key][fname] = [filepath]
                self.db_versions[fname] = {
                    "version": version or self.get_file_mod_date(filepath),
                    "source_url": source_url,
                    "file_mod_date": self.get_file_mod_date(filepath)
                }
                logging.info(f"{filename} found and verified.")
                return True
            logging.warning(f"{filename} exists but is corrupted. Re-downloading...")
        try:
            subprocess.run(
                ["wget", "--quiet", url, "-O", filepath],
                check=False,
            )
        except subprocess.CalledProcessError:
            logging.info(f"{filename} not found.")
            return False

        if self.verify_file_integrity(filepath, filename):
            self.fastas[db_key][fname] = [filepath]
            self.db_versions[fname] = {
                "version": version or self.get_file_mod_date(filepath),
                "source_url": source_url,
                "file_mod_date": self.get_file_mod_date(filepath)
            }
            logging.info(f"{filename} found and verified.")
            return True
        else:
            logging.error(f"Downloaded {filename} is corrupted.")
            return False

    def dl_taxdump(self, url: str = None, dest_path: str = None) -> bool:
        """Download taxdump.tar.gz from NCBI to metadata directory."""
        if url is None:
            url = "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz"
        
        if dest_path is None:
            dest_path = os.path.join(self.metadir, "taxdump.tar.gz")
        
        if os.path.exists(dest_path):
            logging.info("Taxdump already exists")
            return True
        
        os.makedirs(os.path.dirname(dest_path), exist_ok=True)
        
        logging.info(f"Downloading taxdump from {url}")
        try:
            subprocess.run(
                ["wget", "-q", "-O", dest_path, url],
                check=True,
            )
        except subprocess.CalledProcessError:
            logging.error("Failed to download taxdump with wget, trying curl")
            try:
                subprocess.run(
                    ["curl", "-fsSL", "-o", dest_path, url],
                    check=True,
                )
            except subprocess.CalledProcessError:
                logging.error("Failed to download taxdump")
                return False
        
        logging.info("Taxdump downloaded successfully")
        return True

    def extract_taxdump(self, taxdump_path: str = None) -> bool:
        """Extract taxdump to metadata/taxonomy/ directory."""
        if taxdump_path is None:
            taxdump_path = os.path.join(self.metadir, "taxdump.tar.gz")
        
        if not os.path.exists(taxdump_path):
            logging.error(f"Taxdump not found: {taxdump_path}")
            return False
        
        taxonomy_dir = os.path.join(self.metadir, "taxonomy")
        os.makedirs(taxonomy_dir, exist_ok=True)
        
        logging.info(f"Extracting taxdump to {taxonomy_dir}")
        try:
            subprocess.run(
                ["tar", "-xvzf", taxdump_path, "-C", taxonomy_dir],
                check=True,
            )
            logging.info(f"Extracted taxdump to {taxonomy_dir}")
            return True
        except subprocess.CalledProcessError as e:
            logging.error(f"Failed to extract taxdump: {e}")
            return False

    def check_taxdump_extracted(self) -> bool:
        """Check if taxdump has been extracted (names.dmp exists)."""
        names_dmp = os.path.join(self.metadir, "taxonomy", "names.dmp")
        return os.path.exists(names_dmp)

    def ftp_host_file(self, host, source, filename, fname):
        """
        download file from ftp host.
        :param host:
        :param source:
        :param filename:
        :param fname:
        :return:
        """
        filepath = self.seqdir + filename
        if os.path.isfile(filepath):
            if self.verify_file_integrity(filepath, filename):
                self.fastas["host"][fname] = [filepath]
                logging.info(f"{filename} found and verified.")
                return True
            logging.warning(f"{filename} exists but is corrupted. Re-downloading...")
        else:
            if self.test:
                logging.info(f"{filename} not found.")
                return False

        try:
            ftp = FTP(host)
        except Exception as e:
            logging.info(f"{fname} ftp attempt failed. Check internet connection.")
            return False

        ftp.login()
        ftp.cwd(source)
        files = ftp.nlst()
        ftp.quit()

        if filename not in files:
            logging.info(f"{filename} not found.")
            return False

        if self.test:
            logging.info(f"{filename} not found.")
            return False
        else:
            logging.info(f"{filename} not found. downloading...")
            sep = "" if source.startswith("/") else "/"

            try:
                subprocess.run(
                    [
                        "wget",
                        f"ftp://{host}{sep}{source}{filename}",
                        "-P",
                        self.seqdir,
                    ],
                    check=False,
                )
            except subprocess.CalledProcessError:
                logging.info(f"{filename} not found.")
                return False

            if self.verify_file_integrity(filepath, filename):
                self.fastas["host"][fname] = [filepath]
                return True
            else:
                logging.error(f"Downloaded {filename} is corrupted. Download failed.")
                return False

    def find_host(self, host_name):
        """
        find host in host library
        :param host_name:
        :return:
        """
        host = None
        for h in Host.__subclasses__():
            if h().host_name == host_name:
                host = h()
                break

        return host

    def get_host_common_name(self, host_name: str):
        """
        get common name of host
        :param host_name:
        :return:
        """
        host = self.find_host(host_name)
        if not host:
            logging.info(f"{host_name} not found in host library.")
            return False

        return host.common_name

    def download_host(self, host_name: str) -> bool:
        host = self.find_host(host_name)
        if not host:
            logging.info(f"{host_name} not found in host library.")
            return False

        try:
            ftp_download = self.ftp_host_file(
                host.remote_host,
                host.remote_path,
                host.remote_filename,
                host.host_name,
            )

            if ftp_download:
                fpath = self.fastas["host"].get(host_name, [None])[0]
                if fpath:
                    gcf_version = host.remote_filename.split("_")[0]
                    source_url = f"ftp://{host.remote_host}/{host.remote_path}{host.remote_filename}"
                    self.db_versions[host_name] = {
                        "version": gcf_version,
                        "source_url": source_url,
                        "file_mod_date": self.get_file_mod_date(fpath)
                    }

            return ftp_download
        except Exception as e:
            import traceback
            traceback.print_exc()
            logging.warning(f"failed to download {host_name} from ftp.")
            return False

    def get_latest_assembly(
        host, base_path, latest_assembly_dir="latest_assembly_versions"
    ):
        """
        Identify the latest assembly file in the latest_assembly_versions directory.

        :param host: FTP host address.
        :param base_path: Base path to the organism directory on the FTP server.
        :param latest_assembly_dir: Directory containing the latest assembly versions.
        :return: Full path to the latest assembly file ending with '_genomic.fna.gz'.
        """
        try:
            ftp = FTP(host)
            ftp.login()

            # Navigate to the latest_assembly_versions directory
            latest_assembly_path = os.path.join(base_path, latest_assembly_dir)
            ftp.cwd(latest_assembly_path)

            # Get the single subdirectory
            subdirectories = ftp.nlst()
            if len(subdirectories) != 1:
                raise ValueError(
                    "Expected a single subdirectory in latest_assembly_versions."
                )

            # Navigate to the subdirectory
            ftp.cwd(subdirectories[0])

            # Find the file ending with '_genomic.fna.gz'
            files = ftp.nlst()
            for file in files:
                if file.endswith("_genomic.fna.gz"):
                    ftp.quit()
                    return os.path.join(latest_assembly_path, subdirectories[0], file)

            ftp.quit()
            raise FileNotFoundError("No file ending with '_genomic.fna.gz' found.")

        except Exception as e:
            print(f"Error: {e}")
            return None


    def install_requests(self):
        files = self.requests.get("FILES", [])
        if not files:
            return False

        success = False
        for seq_entry in files:
            if not isinstance(seq_entry, dict):
                continue

            path = seq_entry.get("path", "")
            description = seq_entry.get("description", "")

            if not path or not os.path.isfile(path):
                logging.warning(f"Request sequence file not found: {path}")
                continue

            basename = os.path.basename(path)
            for ext in [".gz", ".fa", ".fasta", ".fna"]:
                if basename.endswith(ext):
                    basename = basename[: -len(ext)]
                    break
            key = f"requests_{basename}"

            dest_file = os.path.join(self.seqdir, os.path.basename(path))
            if path != dest_file and not os.path.exists(dest_file):
                shutil.copy(path, dest_file)
                logging.info(f"Copied request sequences file to {dest_file}")
            else:
                dest_file = path if os.path.exists(path) else dest_file

            if not dest_file.endswith(".gz") and os.path.exists(dest_file):
                os.system(f"{BGZIP_BIN} {dest_file}")
                dest_file = dest_file + ".gz"

            if os.path.isfile(dest_file):
                self.fastas["nuc"][key] = [dest_file]
                self.db_versions[key] = {
                    "version": self.get_file_mod_date(dest_file),
                    "source_url": path,
                    "file_mod_date": self.get_file_mod_date(dest_file),
                    "description": description,
                }
                logging.info(f"Request sequences '{key}' prepped from {path}.")
                success = True

        if self.requests.get("ACCID"):
            acc_tsv = self.requests["ACCID"]
            tempfile = os.path.join(self.metadir, "accession_requests.tsv")
            references_file = os.path.join(self.seqdir, "request_references.fa")

            if os.path.isfile(acc_tsv):
                shutil.copy(acc_tsv, tempfile)
            elif "https" in acc_tsv:
                try:
                    subprocess.run(["curl", "-o", tempfile, acc_tsv])
                except subprocess.CalledProcessError:
                    logging.error("accid request file not found")
                    return success

            accid_list = []
            with open(tempfile, "r") as tf:
                accid_list = tf.readlines()
            accid_list = [x.strip() for x in accid_list]

            for accid in accid_list:
                entrez_fetch_sequence(accid, references_file)

            if os.path.isfile(references_file) and os.path.getsize(references_file):
                os.system(f"{BGZIP_BIN} {references_file}")
                references_file = references_file + ".gz"
                self.fastas["nuc"]["requests"] = [references_file]
                self.db_versions["requests"] = {
                    "version": self.get_file_mod_date(references_file),
                    "source_url": "entrez_accid",
                    "file_mod_date": self.get_file_mod_date(references_file),
                }
                logging.info("Request sequences from ACCID prepped.")
                success = True

        return success

    def refseq_prot_dl(self, url: str, filename: str, db_key: str = "refseq_prot"):
        """
        parse and download latest refseq protein db from ncbi ftp.
        :param url: Source URL (e.g., https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/)
        :param filename: File pattern to match (e.g., viral.protein.faa.gz)
        :param db_key: Key to use in fastas/db_versions dicts
        :return:
        """
        host = "ftp.ncbi.nlm.nih.gov"
        source = url.replace(f"https://{host}/", "")
        source_url = url

        fprot = filename
        fprot_suf = os.path.splitext(fprot)[0]
        fprot_path = self.seqdir + fprot

        if os.path.isfile(fprot_path):
            if self.verify_file_integrity(fprot_path, fprot):
                self.fastas["prot"][db_key] = [fprot_path]
                self.db_versions[db_key] = {
                    "version": self.get_file_mod_date(fprot_path),
                    "source_url": source_url,
                    "file_mod_date": self.get_file_mod_date(fprot_path)
                }
                logging.info(f"{fprot_suf} found and verified.")
                return True
            logging.warning(f"{fprot_suf} exists but is corrupted. Re-downloading...")

        else:
            if self.test:
                logging.info(f"{fprot_suf} not found.")
                return False

        try:
            ftp = FTP(host)
        except:
            logging.info("refseq ftp attempt failed. Check internet connection.")
            return False

        ftp.login()
        ftp.cwd(source)
        files = ftp.nlst()
        ftp.quit()

        ext_dict = [x.split(".") for x in files]
        ext_dict = [[".".join(x), ".".join(x[-3:])] for x in ext_dict]
        #
        extset = set([x[1] for x in ext_dict])
        ext_dict = {z: [x[0] for x in ext_dict if x[1] == z] for z in extset}

        protf = [g for x, g in ext_dict.items() if "protein.faa" in x][0]


        logging.info(f"{fprot_suf} not found. downloading...")
        self.get_concat(protf, fprot_suf, host, source)
        if self.verify_file_integrity(fprot_path, fprot):
            self.fastas["prot"][db_key] = [fprot_path]
            self.db_versions[db_key] = {
                "version": self.get_file_mod_date(fprot_path),
                "source_url": source_url,
                "file_mod_date": self.get_file_mod_date(fprot_path)
            }
            return True
        else:
            logging.error(f"Downloaded {fprot_suf} is corrupted.")
            return False

    def refseq_gen_dl(self, url: str, filename: str, db_key: str = "refseq"):
        """
        parse and download latest refseq genome db from ncbi ftp.
        :param url: Source URL (e.g., https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/)
        :param filename: File pattern to match (e.g., viral.genome.fna.gz)
        :param db_key: Key to use in fastas/db_versions dicts
        :return:
        """
        host = "ftp.ncbi.nlm.nih.gov"
        source = url.replace(f"https://{host}/", "")
        source_url = url

        fnuc = filename
        fnuc_suf = os.path.splitext(fnuc)[0]
        fnuc_path = self.seqdir + fnuc

        if self.update:
            if os.path.isfile(fnuc_path):
                logging.info(f"{fnuc_suf} found, removing for Update.")
                os.remove(fnuc_path)

        if os.path.isfile(fnuc_path):
            if self.verify_file_integrity(fnuc_path, fnuc):
                self.fastas["nuc"][db_key] = [fnuc_path]
                self.db_versions[db_key] = {
                    "version": self.get_file_mod_date(fnuc_path),
                    "source_url": source_url,
                    "file_mod_date": self.get_file_mod_date(fnuc_path)
                }
                logging.info(f"{fnuc} found and verified.")
                return True
            logging.warning(f"{fnuc_suf} exists but is corrupted. Re-downloading...")
        else:
            if self.test:
                logging.info(f"{fnuc_suf} not found.")
                return False
        try:
            ftp = FTP(host)
        except:
            logging.info("refseq ftp failed. Check internet connection.")
            return False

        ftp.login()
        ftp.cwd(source)
        files = ftp.nlst()
        ftp.quit()

        ext_dict = [x.split(".") for x in files]
        ext_dict = [[".".join(x), ".".join(x[-3:])] for x in ext_dict]
        #
        extset = set([x[1] for x in ext_dict])
        ext_dict = {z: [x[0] for x in ext_dict if x[1] == z] for z in extset}

        nucf = [g for x, g in ext_dict.items() if "genomic.fna" in x][0]

        fnuc = filename
        fnuc_suf = os.path.splitext(fnuc)[0]
        fnuc_path = self.seqdir + fnuc

        logging.info(f"{fnuc_suf} not found. downloading...")
        self.get_concat(nucf, fnuc_suf, host, source)
        if self.verify_file_integrity(fnuc_path, fnuc):
            self.fastas["nuc"][db_key] = [fnuc_path]
            self.db_versions[db_key] = {
                "version": self.get_file_mod_date(fnuc_path),
                "source_url": source_url,
                "file_mod_date": self.get_file_mod_date(fnuc_path)
            }
            return True
        else:
            logging.error(f"Downloaded {fnuc_suf} is corrupted.")
            return False

    def get_concat(self, flist, outf, host, source):
        """
        download files in list and concatenate into single file. gzip that file
        :param flist: list of files
        :param outf: concatenate output filepath.
        :param host: ftp host
        :param source: ftp diectory
        :return:
        """
        for fl in flist:
            if not os.path.isfile(self.seqdir + fl):
                if self.test:
                    logging.info(f"{fl} not found.")
                else:
                    logging.info(f"{fl} not found. downloading...")
                    correctly_downloaded = 0
                    link = "https://{}/{}{}".format(host, source, fl)

                    while not correctly_downloaded:
                        subprocess.run(["wget", link, "-P", self.seqdir])
                        try:
                            with gzip.open(os.path.join(self.seqdir, fl)) as fd:
                                fd.read()
                            correctly_downloaded = 1
                        except EOFError:
                            correctly_downloaded = 0

        fls = [self.seqdir + fl for fl in sorted(flist)]
        fls = [fl for fl in fls if os.path.isfile(fl)]

        if len(fls) == 0:
            logging.info("No files found.")
            return

        with open(self.seqdir + outf, "wb") as ft:
            for fl in fls:
                try:
                    with gzip.open(fl, "rb") as inf:
                        ft.write(inf.read())
                except EOFError:
                    if os.path.isfile(fl):
                        os.remove(fl)

        os.system("rm {}".format(" ".join(fls)))
        subprocess.run([BGZIP_BIN, self.seqdir + outf])

    def nuc_metadata(self, outfile="acc2taxid.tsv"):
        self.accessions.register_legacy(versions=self.db_versions)
        pending = []
        for dbs, fl_list in self.fastas["nuc"].items():
            files = fl_list if isinstance(fl_list, list) else [fl_list]
            files = [fl for fl in files if os.path.isfile(fl)]
            if not files:
                continue
            if self.test:
                logging.info(f"Test mode: would process {dbs} from nuc fasta files")
                continue
            version = self.db_versions.get(dbs, {}).get("version")
            if self.accessions.needs_rebuild("nuc", dbs, files, version):
                pending.append((dbs, files, version))
        rebuilt = bool(pending)
        for dbs, files, version in pending:
            logging.info(f"Building nucleotide accession index for {dbs}")
            self.accessions.begin_db(
                "nuc",
                dbs,
                files,
                version=version,
                source_url=self.db_versions.get(dbs, {}).get("source_url"),
            )
            for fl in files:
                filename = os.path.basename(fl)
                if dbs == "kraken2":
                    records = self._parse_kraken2_fasta(fl, dbs, filename)
                elif dbs == "virosaurus":
                    records = self._parse_virosaurus_fasta(fl, dbs, filename)
                else:
                    records = self._parse_standard_fasta(fl, dbs, filename)
                self.accessions.insert_batch("nuc", records)
        if pending:
            self.populate_nucleotide_taxids([dbs for dbs, _, _ in pending])
            for dbs, _, _ in pending:
                self.accessions.complete_db("nuc", dbs)
        self.accessions.export_tsv("nuc", self.metadir + outfile, force=self.update or rebuilt)

    def prot_metadata(self):
        self.accessions.register_legacy(versions=self.db_versions)
        pending = []
        for dbs, fl in self.fastas["prot"].items():
            files = fl if isinstance(fl, list) else [fl]
            files = [fpath for fpath in files if os.path.isfile(fpath)]
            if not files:
                continue
            if self.test:
                logging.info(f"Test mode: would process {dbs} from protein fasta files")
                continue
            version = self.db_versions.get(dbs, {}).get("version")
            if self.accessions.needs_rebuild("protein", dbs, files, version):
                pending.append((dbs, files, version))
        rebuilt = bool(pending)
        for dbs, files, version in pending:
            logging.info(f"Building protein accession index for {dbs}")
            self.accessions.begin_db(
                "protein",
                dbs,
                files,
                version=version,
                source_url=self.db_versions.get(dbs, {}).get("source_url"),
            )
            for fpath in files:
                batch = []
                with gzip.open(fpath, "rt") as fn:
                    for line in fn:
                        if not line.startswith(">"):
                            continue
                        line = line[1:].strip()
                        parts = line.split()
                        if not parts:
                            continue
                        acc = parts[0]
                        description = ""
                        if dbs == "refseq_prot":
                            if "[" in line and "]" in line:
                                description = line[line.find("[") + 1:line.find("]")]
                        elif dbs == "rvdb":
                            if "|" in acc:
                                acc_parts = acc.split("|")
                                if len(acc_parts) >= 3:
                                    acc = acc_parts[2]
                        batch.append((dbs, acc, description, acc, None, fpath))
                        if len(batch) >= self.batch_size:
                            self.accessions.insert_batch("protein", batch)
                            batch = []
                if batch:
                    self.accessions.insert_batch("protein", batch)
        if pending:
            pending_keys = [dbs for dbs, _, _ in pending]
            self.populate_protein_taxids(pending_keys)
            for dbs in pending_keys:
                self.accessions.complete_db("protein", dbs)
        outfile = os.path.join(self.metadir, "protein_acc2taxid.tsv")
        self.accessions.export_tsv("protein", outfile, force=self.update or rebuilt)
        if os.path.isfile(outfile):
            for dbs in self.fastas["prot"].keys():
                self.meta[dbs] = outfile

    def populate_protein_taxids(self, dbs_list):
        stage = self.accessions.stage_table("protein")
        if "refseq_prot" in dbs_list:
            with closing(self.accessions.connect("protein")) as conn:
                self._populate_refseq_prot_taxids(conn, stage)
        other_dbs = [dbs for dbs in dbs_list if dbs != "refseq_prot"]
        if not other_dbs:
            return
        acc2tax_dir = self.get_prot()
        threads = [
            Thread(target=self._parse_protein_taxids_thread, args=(dci, acc2tax_dir, other_dbs))
            for dci in range(1, 11)
        ]
        for th in threads:
            th.start()
        for th in threads:
            th.join()
        logging.info("Protein taxids populated")

    def _populate_refseq_prot_taxids(self, conn, stage):
        taxid2desc_path = self.metadir + "taxid2desc.tsv"
        if not os.path.isfile(taxid2desc_path):
            logging.warning(f"taxid2desc.tsv not found at {taxid2desc_path}, generating from taxdump")
            self._generate_taxid2desc_from_taxdump(conn)
            if not os.path.isfile(taxid2desc_path):
                logging.error("Failed to generate taxid2desc.tsv")
                return
        conn.execute("CREATE TABLE IF NOT EXISTS taxid_desc (taxid INTEGER PRIMARY KEY, description TEXT)")
        conn.execute("DELETE FROM taxid_desc")
        batch = []
        with open(taxid2desc_path, "r") as f:
            f.readline()
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) >= 2:
                    try:
                        taxid = int(parts[0])
                        description = parts[1]
                        batch.append((taxid, description))
                        if len(batch) >= self.batch_size:
                            conn.executemany("INSERT OR REPLACE INTO taxid_desc (taxid, description) VALUES (?, ?)", batch)
                            conn.commit()
                            batch = []
                    except ValueError:
                        continue
        if batch:
            conn.executemany("INSERT OR REPLACE INTO taxid_desc (taxid, description) VALUES (?, ?)", batch)
            conn.commit()
        conn.execute(f"""
            UPDATE {stage}
            SET taxid = (
                SELECT t.taxid
                FROM taxid_desc t
                WHERE t.description = {stage}.description
            )
            WHERE dbs = 'refseq_prot' AND description IS NOT NULL AND description != ''
        """)
        conn.commit()
        matched = conn.execute(
            f"SELECT COUNT(*) FROM {stage} WHERE dbs = 'refseq_prot' AND taxid IS NOT NULL"
        ).fetchone()[0]
        total = conn.execute(
            f"SELECT COUNT(*) FROM {stage} WHERE dbs = 'refseq_prot'"
        ).fetchone()[0]
        logging.info(f"refseq_prot: {matched}/{total} accessions matched to taxid")
        conn.execute("DROP TABLE IF EXISTS taxid_desc")

    def _generate_taxid2desc_from_taxdump(self, conn):
        """
        Generate taxid2desc.tsv from taxdump names.dmp if available.
        """
        taxdump_dir = None
        
        for root, dirs, files in os.walk(self.metadir):
            if "names.dmp" in files:
                taxdump_dir = root
                break
        
        if not taxdump_dir:
            logging.warning("taxdump not found in metadata directory")
            return
        
        names_dmp = os.path.join(taxdump_dir, "names.dmp")
        if not os.path.isfile(names_dmp):
            logging.warning(f"names.dmp not found at {names_dmp}")
            return
        
        outfile = self.metadir + "taxid2desc.tsv"
        
        data = []
        with open(names_dmp, "r") as f:
            for line in f:
                parts = line.strip().split("|")
                if len(parts) >= 4:
                    try:
                        taxid = int(parts[0].strip())
                        name = parts[1].strip()
                        name_class = parts[3].strip()
                        
                        if name_class == "scientific name":
                            data.append({"taxid": taxid, "description": name})
                    except (ValueError, IndexError):
                        continue
        
        if data:
            df = pd.DataFrame(data)
            df = df.drop_duplicates(subset="taxid")
            df = df.sort_values("taxid")
            df.to_csv(outfile, sep="\t", index=False)
            logging.info(f"Generated taxid2desc.tsv at {outfile}")
        else:
            logging.warning("Failed to generate taxid2desc.tsv")

    def _parse_protein_taxids_thread(self, dci, acc2tax_dir, dbs_list, chunksize=int(5e5)):
        match_db = os.path.join(self.metadir, f"matches_{dci}.db")
        if not self.accessions.scratch_has_matches(match_db):
            match_conn = self.accessions.create_scratch(match_db)
            acc_sets = {
                dbs: self.accessions.stage_accessions("protein", dbs)
                for dbs in dbs_list
            }
            doc = os.path.join(acc2tax_dir, f"prot.accession2taxid.FULL.{dci}.gz")
            if os.path.isfile(doc):
                try:
                    for chunk in pd.read_csv(
                        doc, compression="gzip", sep="\t",
                        chunksize=chunksize, names=["acc", "taxid"],
                    ):
                        for dbs in dbs_list:
                            matches = chunk[chunk["acc"].isin(acc_sets[dbs])]
                            if not matches.empty:
                                match_conn.executemany(
                                    "INSERT INTO matches VALUES (?, ?, ?)",
                                    [
                                        (dbs, row["acc"], row["taxid"])
                                        for _, row in matches.iterrows()
                                    ],
                                )
                except Exception as e:
                    logging.error(f"Error processing {doc}: {e}")
            match_conn.commit()
            match_conn.close()
        for dbs in dbs_list:
            matches = self.accessions.read_scratch_matches(match_db, dbs)
            self.accessions.update_taxids(
                "protein", [(dbs, acc, taxid) for acc, taxid in matches]
            )
        self.accessions.remove_scratch(match_db)

    def get_prot(self):
        """
        download ncbi protein acc2taxid files.
        :return:
        """
        acc2tax_dir = self.prot_acc2tax_dir

        if not os.path.isdir(acc2tax_dir):
            os.makedirs(acc2tax_dir, exist_ok=True)

        for si in range(1, 11):
            file = f"https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.FULL.{si}.gz"
            filename = os.path.basename(file)
            destination_file = os.path.join(acc2tax_dir, filename)
            fexist = os.path.exists(destination_file)
            print(f"file {filename} exists: {fexist}")
            tries = 0
            while not fexist:
                try:
                    subprocess.run(
                        [
                            "wget",
                            "-P",
                            acc2tax_dir,
                            file,
                        ],
                        check=False,
                    )
                except subprocess.CalledProcessError as e:
                    print(f"failed download protein taxonomy {filename}")
                    tries += 1
                    if tries == 10:
                        logging.info(
                            f"tried downloading {filename} 10 times. check connection. exiting."
                        )
                        raise SystemExit()
                else:
                    fexist = os.path.exists(destination_file)

        return acc2tax_dir

    def get_nuc(self):
        """
        download ncbi nucleotide acc2taxid files.
        Downloads nucl_gb.accession2taxid.gz and nucl_wgs.accession2taxid.gz
        :return:
        """
        import time
        acc2tax_dir = self.nucl_acc2tax_dir

        if not os.path.isdir(acc2tax_dir):
            os.makedirs(acc2tax_dir, exist_ok=True)

        files_to_download = [
            "nucl_gb.accession2taxid.gz",
            "nucl_wgs.accession2taxid.gz"
        ]

        for filename in files_to_download:
            file_url = (
                f"https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/"
                f"accession2taxid/{filename}"
            )

            destination_file = os.path.join(acc2tax_dir, filename)

            if os.path.exists(destination_file) and os.path.getsize(destination_file) > 0:
                print(f"{filename} already exists. Skipping.")
                continue

            print(f"Downloading {filename}")

            max_tries = 10
            success = False

            for attempt in range(1, max_tries + 1):
                try:
                    result = subprocess.run(
                        [
                            "wget",
                            "--continue",
                            "--tries=2",
                            "--timeout=60",
                            "--dns-timeout=30",
                            "-P",
                            acc2tax_dir,
                            file_url,
                        ],
                        check=False,
                        capture_output=True,
                        text=True,
                    )

                    if (
                        result.returncode == 0
                        and os.path.exists(destination_file)
                        and os.path.getsize(destination_file) > 0
                    ):
                        success = True
                        print(f"Successfully downloaded {filename}")
                        break

                    logging.warning(
                        f"Attempt {attempt}/{max_tries} failed for {filename}: "
                        f"{result.stderr[-1000:]}"
                    )

                except Exception as e:
                    logging.warning(f"Download error: {e}")

                if attempt < max_tries:
                    wait_time = min(30 * attempt, 300)
                    print(f"Retrying in {wait_time} seconds...")
                    time.sleep(wait_time)

            if not success:
                logging.error(
                    f"Failed to download {filename} after {max_tries} attempts."
                )
                raise SystemExit(1)

        return acc2tax_dir


    def _parse_kraken2_fasta(self, fl, dbs, filename):
        """
        Parse kraken2 FASTA format: >accession|taxid|...
        Returns all records.
        """
        records = []
        with gzip.open(fl, "rt") as f:
            for line in f:
                if not line.startswith(">"):
                    continue
                line = line[1:].strip()
                parts = line.split("|")
                if len(parts) >= 3:
                    suffix = parts[0]
                    taxid = parts[1]
                    acc = parts[2]
                    acc_in_file = f"{suffix}|{taxid}|{acc}"
                    records.append(
                        (dbs, acc, acc_in_file, int(taxid) if taxid.isdigit() else None, filename)
                    )
        return records

    def _parse_virosaurus_fasta(self, fl, dbs, filename):
        """
        Parse virosaurus FASTA format: >accession.version:GENE:...
        Returns all records.
        """
        records = []
        with gzip.open(fl, "rt") as f:
            for line in f:
                if not line.startswith(">"):
                    continue
                line = line[1:].strip()
                parts = line.split()
                if not parts:
                    continue
                acc_full = parts[0]
                acc = acc_full.split(".")[0]
                acc_in_file = f"{acc}:{acc};"
                records.append((dbs, acc, acc_in_file, None, filename))
        return records

    def _parse_standard_fasta(self, fl, dbs, filename):
        """
        Parse standard FASTA format: >accession description
        Returns all records.
        """
        records = []
        with gzip.open(fl, "rt") as f:
            for line in f:
                if not line.startswith(">"):
                    continue
                line = line[1:].strip()
                parts = line.split()
                if not parts:
                    continue
                acc = parts[0]
                acc_in_file = acc
                records.append((dbs, acc, acc_in_file, None, filename))
        return records

    def populate_nucleotide_taxids(self, dbs_list):
        acc2tax_dir = self.get_nuc()
        files = [
            "nucl_gb.accession2taxid.gz",
            "nucl_wgs.accession2taxid.gz",
        ]
        worker_errors = []
        threads = [
            Thread(
                target=self._parse_nucleotide_taxids_thread,
                args=(filename, acc2tax_dir, dbs_list),
                kwargs={"errors": worker_errors},
            )
            for filename in files
        ]
        for th in threads:
            th.start()
        for th in threads:
            th.join()
        for filename, err in worker_errors:
            logging.error(f"nucleotide taxid worker failed for {filename}: {err}")

    def _parse_nucleotide_taxids_thread(self, filename, acc2tax_dir, dbs_list, chunksize=int(2e6), errors=None):
        match_db = os.path.join(
            self.metadir,
            f"nuc_matches_{filename.replace('.gz', '').replace('.', '_')}.db",
        )
        try:
            if not self.accessions.scratch_has_matches(match_db):
                match_conn = self.accessions.create_scratch(match_db)
                acc_sets = {
                    dbs: self.accessions.stage_accessions("nuc", dbs)
                    for dbs in dbs_list
                }
                doc = os.path.join(acc2tax_dir, filename)
                if os.path.isfile(doc):
                    for chunk in pd.read_csv(
                        doc,
                        compression="gzip",
                        sep="\t",
                        chunksize=chunksize,
                        names=["accession", "accession_version", "taxid", "gi"],
                    ):
                        for dbs in dbs_list:
                            matches = chunk[chunk["accession"].isin(acc_sets[dbs])]
                            if not matches.empty:
                                match_conn.executemany(
                                    "INSERT INTO matches VALUES (?, ?, ?)",
                                    [
                                        (dbs, row["accession"], row["taxid"])
                                        for _, row in matches.iterrows()
                                    ],
                                )
                match_conn.commit()
                match_conn.close()
            for dbs in dbs_list:
                matches = self.accessions.read_scratch_matches(match_db, dbs)
                self.accessions.update_taxids(
                    "nuc", [(dbs, acc, taxid) for acc, taxid in matches]
                )
            self.accessions.remove_scratch(match_db)
        except Exception as e:
            logging.error(f"nucleotide taxid worker {filename} raised: {e}")
            if errors is not None:
                errors.append((filename, str(e)))
            self.accessions.remove_scratch(match_db)


def untax_get(taxdump, odir, dbname, sdir="taxonomy/"):
    """
    download and unzip taxdump.
    :param odir: directory to store taxdump
    :param dbname: name of database directory to store taxdump
    """
    sdir = f"/{sdir}"
    try:
        subprocess.run(
            [
                f"tar",
                "-xvzf",
                f"{odir + dbname}{sdir}taxdump.tar.gz",
                "-C",
                f"{odir + dbname}{sdir}",
            ],
            check=True,
        )

    except subprocess.CalledProcessError:
        logging.info("failed to extract taxdump.")
        if taxdump:
            logging.info(f"getting local {taxdump}")
            os.system(f"cp {taxdump} {odir + dbname}{sdir}taxdump.tar.gz")
        else:
            logging.info("taxdump not provided. downloading using wget.")
            subprocess.run(
                [
                    "wget",
                    "-P",
                    f"{odir + dbname}{sdir}",
                    "ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz",
                ]
            )

        subprocess.run(
            [
                f"tar",
                "-xvzf",
                f"{odir + dbname}/taxonomy/taxdump.tar.gz",
                "-C",
                f"{odir + dbname}/taxonomy/",
            ],
            check=True,
        )


class setup_install(setup_dl):
    def __init__(
        self,
        INSTALL_PARAMS,
        home="",
        bindir="",
        taxdump="",
        test=False,
        update=False,
    ):
        super().__init__(INSTALL_PARAMS, home, bindir, test=test, update=update)
        self.taxdump = taxdump

        if not self.taxdump:
            logging.info(
                "taxdump not provided. will use default software \
                    download. May encounter issues. Suggest provinding."
            )

    def install_prep(self):
        """
        initializes dbs dictionary to store installation directories by name, know what you have installed.
        """

        self.dbs = {}

    def get_dbs_entry(self, tool: str, dbname: str) -> dict:
        """Return the dbs entry for (tool, dbname), or {} if not registered."""
        return self.dbs.get(tool, {}).get(dbname, {})

    def centrifuge_download_install(
            self,
            dbname="viral",
            threads="3",
            id="centrifuge",
            dbdir="centrifuge",
            dlp="wget",
    ):
        
        odir = self.dbdir + dbdir + "/"
        bin = self.envs["ROOT"] + self.envs[id] + "/bin/"
        sdir = odir + dbname + "/" + dbname
        index_file_prefix = f"{odir}{dbname}/{dbname}_index"
        

        href = get_db_url('centrifuge', dbname)
        filename = get_db_file('centrifuge', dbname)
        if not href:
            logging.error(f"No source URL found for centrifuge/{dbname}")
            return False

        if self.update:
            logging.info(f"Updating centrifuge db {dbname}.")
            if os.path.exists(f"{odir}{dbname}"):
                logging.info(f"Removing old centrifuge db {dbname} index.")
                shutil.rmtree(f"{odir}{dbname}")
        
        if os.path.isfile(index_file_prefix + ".1.cf"):
            logging.info(f"Centrifuge db {dbname} index is installed.")
            centrifuge_fasta = f"{sdir}/complete.fna.gz"
            if os.path.isfile(os.path.splitext(centrifuge_fasta)[0]):
                os.system(f"bgzip {sdir}/complete.fna")

            self.dbs.setdefault(id, {})[dbname] = {
                "dir": odir,
                "dbname": dbname,
                "fasta": "",
                "db": index_file_prefix,
                "status": "success",
            }
            return True

        else:
            if self.test:
                logging.info(f"Centrifuge db {dbname} is not installed.")
                return False
            else:
                logging.info(f"Centrifuge db {dbname} is not installed. Installing...")


        if os.path.isfile(index_file_prefix + ".1.cf"):
            logging.info(f"Centrifuge db {dbname} index is installed.")
            self.dbs.setdefault(id, {})[dbname] = {
                "dir": odir,
                "dbname": dbname,
                "fasta": "",
                "db": index_file_prefix,
                "status": "success",
            }
            return True
        
        try:

            os.makedirs(sdir, exist_ok=True)
            os.system(f"wget {href} -O {sdir}/{filename}")
            os.system(f" tar -xvzf {sdir}/{filename} -C {sdir}")
            os.system(f"rm {sdir}/{filename}")

            files_in_directory = os.listdir(sdir)
            index_files = sorted([f for f in files_in_directory if f.endswith(".cf")])
            for ix, f in enumerate(index_files):
                os.rename(os.path.join(sdir, f), f"{index_file_prefix}.{ix+1}.cf")


            if os.path.isfile(index_file_prefix + ".1.cf"):
                logging.info(f"Centrifuge db {dbname} index is installed.")
                self.dbs.setdefault(id, {})[dbname] = {
                    "dir": odir,
                    "dbname": dbname,
                    "fasta": "",
                    "db": index_file_prefix,
                    "status": "success",
                }
                return True
            else:
                logging.info(f"Centrifuge db {dbname} index is not installed.")

                return False
        except subprocess.CalledProcessError:
            logging.error(f"Error occurred while installing centrifuge db {dbname}.")
            return False

    def centrifuge_register_prebuilt(self, dbname, id = "centrifuge"):
        """Register a pre-built Centrifuge index.

        Args:
            dbname: Name of the database (e.g., 'viral', 'bacteria')

        Returns:
            True if prebuilt index found and registered, False otherwise
        """
        prebuilt_path = os.path.join("/opt/data/prebuilt/", "centrifuge", dbname)
        odir = self.dbdir + "centrifuge" + "/"
        index_file_prefix = f"{odir}{dbname}/{dbname}_index"

        #sys.path.insert(0, os.path.dirname(__file__))
        #from load_sources import get_prebuilt_index_path
        logging.info(f"############# PREBUILT NEW {index_file_prefix}")
        print(os.path.isfile(index_file_prefix + ".1.cf"))

        #prebuilt_path = get_prebuilt_index_path("centrifuge", dbname)
        logging.info(f"############# PREBUILT {prebuilt_path}")
        if prebuilt_path is None or os.path.exists(prebuilt_path) == False:
            logging.info(f"No prebuilt path configured for centrifuge/{dbname}: {prebuilt_path}")
            return False

        index_files = [f for f in os.listdir(prebuilt_path) if f.endswith(".cf")]
        logging.info(f"Found prebuilt centrifuge index files: {index_files}")
        if len(index_files) < 2:
            logging.info(f"Prebuilt centrifuge/{dbname} index files not found at {prebuilt_path}")
            return False

        dest_files = [f"{index_file_prefix}.{i+1}.cf" for i in range(len(index_files))]
        dest_files_exist = [os.path.exists(f) for f in dest_files]
        if all(dest_files_exist):
            logging.info(f"Found prebuilt centrifuge/{dbname} at {prebuilt_path}")
            if self.update:
                logging.info(f"Updating centrifuge/{dbname} index from {prebuilt_path}")
                shutil.rmtree(f"{odir}{dbname}", ignore_errors=True)
            else:

                self.dbs.setdefault(id, {})[dbname] = {
                    "dir": odir,
                    "dbname": dbname,
                    "db": index_file_prefix,
                    "status": "success",
                }
                return True


        os.makedirs(odir + dbname, exist_ok=True)
        index_files = sorted(index_files)
        for i in range(len(index_files)):
            src = os.path.join(prebuilt_path, index_files[i])
            dst = f"{index_file_prefix}.{i+1}.cf"
            if not os.path.isfile(dst):
                shutil.move(src, dst)

        self.dbs.setdefault("centrifuge", {})[dbname] = {
            "dir": prebuilt_path + "/",
            "dbname": dbname,
            "fasta": "",
            "db": index_file_prefix,
            "status": "success",
        }
        return True

    def kraken2_register_prebuilt(self, dbname="viral", id="kraken2"):
        """Register a pre-built Kraken2 index.

        Args:
            dbname: Name of the database (e.g., 'viral', 'bacteria', 'eupathdb46')

        Returns:
            True if prebuilt index found and registered, False otherwise
        """
        prebuilt_path = os.path.join("/opt/data/prebuilt/", "kraken2", dbname)
        odir = self.dbdir + "kraken2/"

        import sys
        sys.path.insert(0, os.path.dirname(__file__))
        from load_sources import get_prebuilt_index_path

        #prebuilt_path = get_prebuilt_index_path("kraken2", dbname)

        if not os.path.exists(prebuilt_path):
            logging.info(f"Prebuilt kraken2/{dbname} not found at {prebuilt_path}")
            return False

        if not os.path.isfile(prebuilt_path + "/taxo.k2d"):
            logging.info(f"Prebuilt kraken2/{dbname} index file (taxo.k2d) not found at {prebuilt_path}")
            return False

        if os.path.exists(prebuilt_path + "/taxo.k2d"):
            logging.info(f"Prebuilt kraken2/{dbname} index file (taxo.k2d) found at {prebuilt_path}")
            if self.update:
                logging.info(f"Updating kraken2/{dbname} index from {prebuilt_path}")
                shutil.rmtree(os.path.join(odir, dbname), ignore_errors=True)
            else:
                self.dbs.setdefault(id, {})[dbname] = {
                    "dir": prebuilt_path + "/",
                    "dbname": dbname,
                    "db": prebuilt_path,
                    "version": "prebuilt",
                    "source_url": "user-provided",
                    "status": "success",
                }
                return True

        logging.info(f"Found prebuilt kraken2/{dbname} at {prebuilt_path}")
        shutil.copytree(prebuilt_path, os.path.join(odir, dbname), dirs_exist_ok=True)

        self.dbs.setdefault("kraken2", {})[dbname] = {
            "dir": prebuilt_path + "/",
            "dbname": dbname,
            "db": prebuilt_path,
            "version": "prebuilt",
            "source_url": "user-provided",
            "status": "success",
        }
        return True

    def centrifuge_install(
        self,
        dbname="viral",
        threads="3",
        id="centrifuge",
        dbdir="centrifuge",
        dlp="wget",
    ):
        """
        install centrifuge.
        :param dbname: name of centrifuge db.
        :param threads: number of threads to use.
        :return:
        """
        if not dbdir:
            dbdir = id

        odir = self.dbdir + dbdir + "/"
        bin = self.envs["ROOT"] + self.envs[id] + "/"
        sdir = odir + dbname + "/" + dbname
        index_file_prefix = f"{odir}{dbname}/{dbname}_index"
        old_index_file_prefix = f"{odir}{dbname}/index"

        if self.update:
            logging.info(f"Updating centrifuge db {dbname} index.")
            if os.path.exists(f"{odir}{dbname}"):
                logging.info(f"Removing old centrifuge db {dbname} index.")
                shutil.rmtree(f"{odir}{dbname}")

        if os.path.isfile(index_file_prefix + ".1.cf"):
            logging.info(f"Centrifuge db {dbname} index is installed.")
            centrifuge_fasta = f"{sdir}/complete.fna.gz"
            if os.path.isfile(os.path.splitext(centrifuge_fasta)[0]):
                os.system(f"bgzip {sdir}/complete.fna")
                # compress_using_xopen(f"{sdir}/complete.fna", f"{sdir}/complete.fna.gz")

            self.dbs.setdefault(id, {})[dbname] = {
                "dir": odir,
                "dbname": dbname,
                "fasta": f"{sdir}/complete.fna.gz",
                "db": index_file_prefix,
            }
            return True

        elif os.path.isfile(old_index_file_prefix + ".1.cf"):

            logging.info(f"Centrifuge db {dbname} index is installed.")
            centrifuge_fasta = f"{sdir}/complete.fna.gz"
            if os.path.isfile(os.path.splitext(centrifuge_fasta)[0]):
                os.system(f"bgzip {sdir}/complete.fna")
            #                compress_using_xopen(f"{sdir}/complete.fna", f"{sdir}/complete.fna.gz")

            # create symlink to new index for files that use old index
            files_in_directory = os.listdir(odir + dbname)

            for file in files_in_directory:
                if file.startswith("index"):
                    new_filename = file.replace("index", dbname + "_index")
                    new_filepath = os.path.join(odir + dbname, new_filename)
                    old_file_path = os.path.join(odir + dbname, file)
                    os.symlink(old_file_path, new_filepath)

            self.dbs.setdefault(id, {})[dbname] = {
                "dir": odir,
                "dbname": dbname,
                "fasta": f"{sdir}/complete.fna.gz",
                "db": index_file_prefix,
            }
            return True

        else:
            if self.test:
                logging.info(f"Centrifuge db {dbname} is not installed.")
                return False
            else:
                logging.info(f"Centrifuge db {dbname} is not installed. Installing...")

        subprocess.run(["mkdir", "-p", odir])
        #
        tax_command = [
            bin + "centrifuge-download",
            "-o",
            odir + dbname,
            "-P",
            threads,
            "-g",
            dlp,
            "taxonomy",
        ]
        #
        seqmap_command = [
            bin + "centrifuge-download",
            "-o",
            odir + dbname,
            "-P",
            threads,
            "-m",
            "-d",
            dbname,
            "-g",
            dlp,
            "refseq",
        ]

        build_command = [
            bin + "centrifuge-build",
            "-o",
            odir + dbname,
            "-p",
            threads,
            "--conversion-table",
            f"{odir}{dbname}.seq2taxid.map",
            "--taxonomy-tree",
            f"{odir}{dbname}/nodes.dmp",
            "--name-table",
            f"{odir}{dbname}/names.dmp",
            f"{sdir}/complete.fna",
            index_file_prefix,
        ]

        ###
        try:
            subprocess.run(tax_command)
            os.system(" ".join(seqmap_command) + f" > {odir}{dbname}.seq2taxid.map")

            # iterate over .fna files and concatenate them into one file
            fna_files = [f for f in os.listdir(sdir) if f.endswith(".fna")]
            for f in fna_files:
                os.system(f"cat {sdir}/{f} >> {sdir}/complete.fna")
            #

            try:
                subprocess.run(build_command)
            except subprocess.CalledProcessError:
                logging.info(f"failed to build centrifuge db {dbname}")

                if not os.path.exists(
                    f"{odir}{dbname}/nodes.dmp"
                ) or not os.path.exists(f"{odir}{dbname}/names.dmp"):
                    if not self.taxdump:
                        cmd_dl_taxonomy = [
                            "wget",
                            "--no-check-certificate",
                            "-P",
                            f"{odir}{dbname}",
                            "ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz",
                        ]
                        subprocess.run(cmd_dl_taxonomy)

                    else:
                        shutil.copy(self.taxdump, f"{odir}{dbname}/taxdump.tar.gz")

                    subprocess.run(
                        [
                            f"tar",
                            "-xvzf",
                            f"{odir}{dbname}/taxdump.tar.gz",
                            "-C",
                            f"{odir}{dbname}",
                        ],
                        check=False,
                    )

                subprocess.run(build_command)

            os.system(f"bgzip {sdir}/complete.fna")
            # compress_using_xopen(f"{sdir}/complete.fna", f"{sdir}/complete.fna.gz")

            self.dbs.setdefault(id, {})[dbname] = {
                "dir": odir,
                "dbname": dbname,
                "fasta": f"{sdir}/complete.fna.gz",
                "db": index_file_prefix,
            }
            return True

        except subprocess.CalledProcessError:
            logging.warning(f"failed to download centrifuge db {dbname}")
            return False

    def voyager_install_viruses_copy(
        self,
        dbname="viral",
        id="voyager",
        dbdir="voyager",
    ):
        db_key = dbname  # caller's name; dbname is translated/mutated below

        dbname_translate = {
            "viral": "viruses",
            "bacteria": "bacteria",
            "archaea": "archaea",
            "fungi": "fungi",
        }
        dbname = dbname_translate[dbname]
        odir = self.dbdir + dbdir + "/"

        if self.update:
            logging.info(f"Updating voyager db {dbname}.")
            if os.path.exists(odir + dbname):
                logging.info(f"Removing old voyager db {dbname}.")
                shutil.rmtree(odir + dbname)

        print(f"odir: {odir}, dbname: {dbname}")
        print(os.path.isfile(os.path.join(odir, dbname, f"{dbname}.idx")))
        if os.path.isfile(os.path.join(odir, dbname, f"{dbname}.idx")):
            logging.info(f"Voyager db {dbname} is installed.")
            self.dbs.setdefault(id, {})[db_key] = {
                "dir": odir,
                "dbname": dbname,
                "db": os.path.join(odir, dbname, f"{dbname}.idx"),
            }
            return True
        else:
            if self.test:
                logging.info(f"Voyager db {dbname} is not installed.")
                return False
            else:
                logging.info(f"Voyager db {dbname} is not installed. Installing...")

        subprocess.run(["mkdir", "-p", odir])
        subprocess.run(["mkdir", "-p", odir])
        subprocess.run(["cp", "-r", "install_scripts/software/viruses.tar.gz", odir])
        subprocess.run(["tar", "-xvf", odir + f"/viruses.tar.gz", "-C", odir])
        subprocess.run(["rm", odir + f"/viruses.tar.gz"])

        self.dbs.setdefault(id, {})[db_key] = {
            "dir": odir,
            "dbname": dbname,
            "db": os.path.join(odir, dbname, f"{dbname}.idx"),
        }

        if os.path.isfile(os.path.join(odir, dbname, f"{dbname}.idx")):
            logging.info(f"Voyager db {dbname} is installed.")
            return True
        else:
            logging.info(f"Voyager db {dbname} is not installed.")
            return False

    def voyager_install_download(
        self,
        dbname="viral",
        id="voyager",
        dbdir="voyager",
    ):
        db_key = dbname  # caller's name; dbname is translated/mutated below

        dbname_translate = {
            "viral": "viruses",
            "bacteria": "bacteria",
            "archaea": "archaea",
            "fungi": "fungi",
        }

        dbname = dbname_translate[dbname]
        odir = self.dbdir + dbdir + "/"

        if dbname == "viral":
            source = "ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/voyager/viral/"
        elif dbname == "bacteria":
            source = "ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/voyager/bacterial/"

        if self.update:
            logging.info(f"Updating voyager db {dbname}.")
            if os.path.exists(odir + dbname):
                logging.info(f"Removing old voyager db {dbname}.")
                shutil.rmtree(odir + dbname)

        if os.path.isfile(os.path.join(odir, dbname, f"{dbname}.idx")):
            logging.info(f"Voyager db {dbname} is installed.")
            self.dbs.setdefault(id, {})[db_key] = {
                "dir": odir,
                "dbname": dbname,
                "db": os.path.join(odir, dbname, f"{dbname}.idx"),
            }
            return True
        else:
            if self.test:
                logging.info(f"Voyager db {dbname} is not installed.")
                return False
            else:
                logging.info(f"Voyager db {dbname} is not installed. Installing...")
        subprocess.run(["mkdir", "-p", odir])
        subprocess.run(["mkdir", "-p", odir + dbname])
        subprocess.run(["wget", "-P", odir + dbname, source])
        subprocess.run(
            ["tar", "-xvzf", odir + dbname + f"/{dbname}.tar.gz", "-C", odir + dbname]
        )
        subprocess.run(["rm", odir + dbname + f"/{dbname}.tar.gz"])

        self.dbs.setdefault(id, {})[db_key] = {
            "dir": odir,
            "dbname": dbname,
            "db": os.path.join(odir, dbname, f"{dbname}.idx"),
        }
        if os.path.isfile(os.path.join(odir, dbname, f"{dbname}.idx")):
            logging.info(f"Voyager db {dbname} is installed.")
            return True
        else:
            logging.info(f"Voyager db {dbname} is not installed.")
            return False

    def install_metaphlan(
        self,
        dbname="mpa_vJan25_CHOCOPhlAnSGB_202503",
        id="metaphlan",
        dbdir="metaphlan",
        dlp="wget",
    ):
        """
        install metaphlan database.
        """
        odir = self.dbdir + dbdir + "/"
        # dbname = dbname + ".mpa"

        if self.update:
            logging.info(f"Updating metaphlan db {dbname}.")
            if os.path.exists(odir + dbname):
                logging.info(f"Removing old metaphlan db {dbname}.")
                shutil.rmtree(odir + dbname)

        if os.path.isfile(
            odir + dbname + "/{}.pkl".format(os.path.splitext(dbname)[0])
        ):
            logging.info(f"Metaphlan db {dbname} is installed.")
            self.dbs.setdefault(id, {})[dbname] = {
                "dir": odir,
                "dbname": dbname,
                "db": f"{odir}{dbname}/{dbname}.pkl",
                "status": "success",
            }
            return True

        else:
            if self.test:
                logging.info(f"Metaphlan db {dbname} is not installed.")
                self.dbs.setdefault(id, {})[dbname] = {
                    "dir": odir,
                    "dbname": dbname,
                    "db": f"{odir}{dbname}/{dbname}.pkl",
                    "status": "test",
                }
                return False
            else:
                logging.info(f"Metaphlan db {dbname} is not installed. Installing...")

        subprocess.run(["mkdir", "-p", odir])
        subprocess.run(["mkdir", "-p", odir + dbname])

        try:
            bin = self.envs["ROOT"] + self.envs[id] + "/bin/"
            cmd = [
                f"{bin}metaphlan",
                "--install",
                "--index",
                f"{dbname}",
                "--bowtie2db",
                f"{odir}",
            ]
            cmd = " ".join(cmd)

            # subprocess.run(cmd, check=True)
            metaphlan_install_script = f"{odir + dbname}/install.sh"
            with open(metaphlan_install_script, "w") as f:
                f.write("#!/bin/bash\n")
                f.write('eval "$(conda shell.bash hook)"\n')
                f.write(f"conda activate {self.envs['ROOT'] + self.envs[id]}\n")
                f.write(f"export PATH={bin}:$PATH\n")
                f.write(cmd + "\n")
                f.write("conda deactivate\n")
            os.chmod(metaphlan_install_script, 0o755)
            subprocess.run([metaphlan_install_script], check=True)
            subprocess.run(["rm", metaphlan_install_script], check=True)

            self.dbs.setdefault(id, {})[dbname] = {
                "dir": odir,
                "dbname": dbname,
                "db": f"{odir}{dbname}/{dbname}",
                "status": "success",
            }
            return True

        except subprocess.CalledProcessError:
            logging.warning(f"failed to download metaphlan db {dbname}")
            self.dbs.setdefault(id, {})[dbname] = {
                "dir": odir,
                "dbname": dbname,
                "db": f"{odir}{dbname}/{dbname}",
                "status": "failed",
            }
            return False

    def install_metaphlan_dl(
        self,
        dbname="mpa_vFeb24_CDIFF_CHOCOPhlAnSGB_20240910",
        id="metaphlan",
        dbdir="metaphlan",
    ):
        db_key = dbname  # caller's name; dbname is translated/mutated below

        odir = self.dbdir + dbdir + "/"
        dbname = dbname + ".mpa"

        if self.update:
            logging.info(f"Updating metaphlan db {dbname}.")
            if os.path.exists(odir + dbname):
                logging.info(f"Removing old metaphlan db {dbname}.")
                shutil.rmtree(odir + dbname)

        if os.path.isfile(
            odir + dbname + "/{}.pkl".format(os.path.splitext(dbname)[0])
        ):
            logging.info(f"Metaphlan db {dbname} is installed.")
            self.dbs.setdefault(id, {})[db_key] = {
                "dir": odir,
                "dbname": dbname,
                "db": f"{odir}{dbname}/{dbname}.pkl",
            }
            return True
        else:
            if self.test:
                logging.info(f"Metaphlan db {dbname} is not installed.")
                return False
            else:
                logging.info(f"Metaphlan db {dbname} is not installed. Installing...")

        subprocess.run(["mkdir", "-p", odir])
        subprocess.run(["mkdir", "-p", odir + dbname])

        source = get_db_url('metaphlan', dbname)
        if not source:
            source = get_db_url('metaphlan', 'default')
        
        if not source:
            logging.error(f"No source URL found for metaphlan/{dbname}")
            return False

        try:

            source_file = source.split("/")[-1]
            subprocess.run(["wget", "-P", odir + dbname, source], check=True)
            subprocess.run(
                ["tar", "-xvf", odir + dbname + "/" + source_file, "-C", odir + dbname],
                check=True,
            )
            subprocess.run(["rm", odir + dbname + "/" + source_file], check=True)
            self.dbs.setdefault(id, {})[db_key] = {
                "dir": odir,
                "dbname": dbname,
                "db": f"{odir}{dbname}/{dbname}.pkl",
            }
            return True

        except subprocess.CalledProcessError:
            logging.warning(f"failed to download metaphlan db {dbname}")
            return False

    def clark_install(
        self,
        dbname="viral",
        id="clark",
        dbdir="clark",
    ):
        db_key = dbname  # caller's name; dbname is translated/mutated below
        dbname_translate = {
            "viral": "viruses",
            "bacteria": "bacteria",
            "archaea": "archaea",
            "fungi": "fungi",
        }

        dbname = dbname_translate[dbname]

        odir = self.dbdir + dbdir + "/"
        bin = self.envs["ROOT"] + self.envs[id] + "/"

        if os.path.isfile(odir + dbname + f"/.{dbname}"):
            logging.info(f"{id} db {dbname} is installed.")

            self.dbs.setdefault(id, {})[db_key] = {
                "dir": odir,
                "dbname": dbname,
                "db": f"{odir}{dbname}",
            }
            return True
        else:
            if self.test:
                logging.info(f"CLARK db {dbname} is not installed.")
                return False
            else:
                logging.info(f"CLARK db {dbname} is not installed. Installing...")

        subprocess.run(["mkdir", "-p", odir])
        ##

        lib_command = [
            bin + "set_targets.sh",
            odir + dbname,
            dbname,
            "--species",
        ]

        spaced_command = [bin + "buildSpacedDB.sh"]

        try:
            try:
                subprocess.run(lib_command)

            except subprocess.CalledProcessError as e:
                logging.error(f"CLARK db {dbname} failed to download. {e}")
                return

            try:
                subprocess.run(spaced_command)

            except subprocess.CalledProcessError as e:
                logging.error(f"CLARK db spaced DB failed to build. {e}")
                return

            self.dbs.setdefault(id, {})[db_key] = {
                "dir": odir,
                "dbname": dbname,
                "fasta": odir + dbname + "/library/" + dbname + "/library.fna.gz",
                "db": f"{odir}{dbname}",
            }
            return True
        except subprocess.CalledProcessError:
            logging.warning(f"failed to download CLARK db {dbname}")
            return False

    def kraken2_download_install(
        self,
        dbname="viral",
        id="kraken2",
        dbdir="kraken2",
    ):
        odir = self.dbdir + dbdir + "/"

        source = get_db_url('kraken2', dbname)
        kraken_version = get_db_version('kraken2', dbname)

        if not source:
            logging.error(f"No source URL found for kraken2/{dbname}")
            return False

        source_file = source.split("/")[-1]

        if self.update:
            logging.info(f"Updating Kraken2 db {dbname}.")
            if os.path.exists(odir + dbname):

                logging.info(f"Removing old Kraken2 db {dbname}.")
                shutil.rmtree(odir + dbname)

        if os.path.isfile(odir + dbname + "/taxo.k2d"):
            logging.info(f"Kraken2 db {dbname} k2d file exists. Kraken2 is installed.")

            self.dbs.setdefault(id, {})[dbname] = {
                "dir": odir,
                "dbname": dbname,
                "fasta": odir + dbname + "/library/" + dbname + "/library.fna.gz",
                "db": odir + dbname,
                "version": kraken_version,
                "source_url": source,
                "status": "success",
            }
            return True
        else:
            if self.test:
                logging.info(f"Kraken2 db {dbname} is not installed.")
                self.dbs.setdefault(id, {})[dbname] = {
                    "dir": odir,
                    "dbname": dbname,
                    "fasta": odir + dbname + "/library/" + dbname + "/library.fna.gz",
                    "db": odir + dbname,
                    "version": kraken_version,
                    "source_url": source,
                    "status": "test",
                }
                return False
            else:
                logging.info(
                    f"Kraken2 db {dbname} is not installed. Downloading from {source}."
                )

        sdir = odir + dbname + "/"
        subprocess.run(["mkdir", "-p", sdir])
        subprocess.run(["wget", "-P", sdir, source])
        subprocess.run(["tar", "-xvzf", sdir + source_file, "-C", sdir])
        subprocess.run(["rm", sdir + source_file])

        if os.path.isfile(odir + dbname + "/taxo.k2d"):
            logging.info(f"Kraken2 db {dbname} k2d file exists. Kraken2 is installed.")
            self.dbs.setdefault(id, {})[dbname] = {
                "dir": odir,
                "dbname": dbname,
                "fasta": odir + dbname + "/library/" + dbname + "/library.fna.gz",
                "db": odir + dbname,
                "version": kraken_version,
                "source_url": source,
                "status": "success",
            }
            return True
        else:
            logging.warning(f"failed to download Kraken2 db {dbname}")
            self.dbs.setdefault(id, {})[dbname] = {
                "dir": odir,
                "dbname": dbname,
                "fasta": odir + dbname + "/library/" + dbname + "/library.fna.gz",
                "db": odir + dbname,
                "version": kraken_version,
                "source_url": source,
                "status": "failed",
            }
            return False

    def kraken2_install(
        self,
        dbname="viral",
        threads="5",
        id="kraken2",
        dbdir="kraken2",
        build_args="--max-db-size 18000000000 --kmer-len 31",
        ftp=False,
    ):
        odir = self.dbdir + dbdir + "/"
        bin = self.envs["ROOT"] + self.envs[id] + "/bin/"

        if self.update:
            logging.info(f"Updating Kraken2 db {dbname}.")
            if os.path.exists(odir + dbname):
                logging.info(f"Removing old Kraken2 db {dbname}.")
                shutil.rmtree(odir + dbname)

        if os.path.isfile(odir + dbname + "/taxo.k2d"):
            logging.info(f"Kraken2 db {dbname} k2d file exists. Kraken2 is installed.")
            krk2_fasta = odir + dbname + "/library/" + dbname + "/library.fna.gz"

            if os.path.isfile(os.path.splitext(krk2_fasta)[0]):
                os.system(f"{BGZIP_BIN} " + os.path.splitext(krk2_fasta)[0])

            self.dbs.setdefault(id, {})[dbname] = {
                "dir": odir,
                "dbname": dbname,
                "fasta": odir + dbname + "/library/" + dbname + "/library.fna.gz",
                "db": odir + dbname,
            }

            return True
        else:
            if self.test:
                logging.info(f"Kraken2 db {dbname} is not installed.")
                return False
            else:
                logging.info(f"Kraken2 db {dbname} is not installed. Installing...")

        subprocess.run(["mkdir", "-p", odir])
        ##

        lib_command = [
            bin + "kraken2-build",
            "--download-library",
            dbname,
            "--db",
            odir + dbname,
            "--threads",
            threads,
        ]
        tax_command = [
            bin + "kraken2-build",
            "--download-taxonomy",
            "--db",
            odir + dbname,
            "--threads",
            threads,
        ]
        build_command = [
            bin + "kraken2-build",
            "--build",
            build_args,
            "--db",
            odir + dbname,
            "--threads",
            threads,
        ]

        if ftp:
            lib_command.append("--use-ftp")
            tax_command.append("--use-ftp")

        try:
            try:
                subprocess.run(lib_command)

            except subprocess.CalledProcessError as e:
                if not ftp:
                    logging.info(
                        f"Kraken2 db {dbname} library download failed. Attempting to download from ftp."
                    )
                    lib_command.append("--ftp")
                    subprocess.run(lib_command)

                logging.error("kraken2 library download failed command.")
            ##
            try:
                subprocess.run(tax_command)

            except subprocess.CalledProcessError as e:
                if not ftp:
                    logging.info(
                        f"Kraken2 db {dbname} taxonomy download failed. Attempting to download from ftp."
                    )
                    tax_command.append("--ftp")
                    subprocess.run(tax_command)
                logging.error("kraken2 taxonomy download failed command.")

            subprocess.call(" ".join(build_command), shell=True)
            untax_get(self.taxdump, odir, dbname)
            os.system(
                f"{BGZIP_BIN} " + odir + dbname + "/library/" + dbname + "/library.fna"
            )

            if os.path.isfile(odir + dbname + "/taxo.k2d"):
                logging.info(
                    f"Kraken2 db {dbname} k2d file exists. Kraken2 is installed."
                )
                self.dbs.setdefault(id, {})[dbname] = {
                    "dir": odir,
                    "dbname": dbname,
                    "fasta": odir + dbname + "/library/" + dbname + "/library.fna.gz",
                    "db": odir + dbname,
                }
                return True
            else:
                logging.warning(f"failed to download Kraken2 db {dbname}")
                self.dbs.setdefault(id, {})[dbname] = {
                    "dir": odir,
                    "dbname": dbname,
                    "fasta": odir + dbname + "/library/" + dbname + "/library.fna.gz",
                    "db": odir + dbname,
                }
                return False

        except subprocess.CalledProcessError:
            logging.warning(f"failed to download Kraken2 db {dbname}")
            self.dbs.setdefault(id, {})[dbname] = {
                "dir": odir,
                "dbname": dbname,
                "fasta": odir + dbname + "/library/" + dbname + "/library.fna.gz",
                "db": odir + dbname,
            }
            return False

    def kraken2_two_strategies_install(
        self,
        dbname="viral",
        threads="5",
        id="kraken2",
        dbdir="kraken2",
        build_args="--max-db-size 18000000000 --kmer-len 31",
        ftp=False,
    ):
        odir = self.dbdir + dbdir + "/"

        traditional_install = self.kraken2_install(
            dbname=dbname,
            threads=threads,
            id=id,
            dbdir=dbdir,
            build_args=build_args,
            ftp=ftp,
        )

        if traditional_install:
            return True

        # if traditional install fails, try to download and install
        # remove the dbdir and try again
        os.system("rm -rf " + odir)

        kraken2_download_install = self.kraken2_download_install(
            dbname=dbname,
            id=id,
            dbdir=dbdir,
        )

        if kraken2_download_install:
            return True

        return False

    def diamond_install(
        self, id="diamond", dbdir="diamond", dbname="swissprot", db="swissprot.gz"
    ):
        odir = self.dbdir + dbdir + "/"
        bin = self.envs["ROOT"] + self.envs[id] + "/bin/"

        if self.update:
            if os.path.exists(odir + dbname):
                logging.info(f"Removing old Diamond db {dbname}.")
                shutil.rmtree(odir + dbname)

        if os.path.isfile(odir + dbname + ".dmnd"):
            logging.info(f"Diamond db {dbname}.dmnd present. Diamond prepped.")
            self.dbs.setdefault(id, {})[dbname] = {
                "dir": odir,
                "dbname": dbname,
                "db": odir + dbname,
                "status": "success",
            }

            return True
        else:
            if self.test:
                logging.info(f"Diamond db {dbname} not installed.")
                self.dbs.setdefault(id, {})[dbname] = {
                    "dir": odir,
                    "dbname": dbname,
                    "db": odir + dbname,
                    "status": "test",
                }
                return False
            else:
                logging.info(f"Diamond db {dbname} . Installing...")

        try:
            subprocess.run(["mkdir", "-p", odir])
            command = [bin + "diamond", "makedb", "--in", db, "--db", odir + dbname]

            subprocess.call(" ".join(command), shell=True)

            self.dbs.setdefault(id, {})[dbname] = {
                "dir": odir,
                "dbname": dbname,
                "db": odir + dbname,
                "status": "success",
            }
            return True

        except subprocess.CalledProcessError:
            logging.warning(f"failed to download Diamond db {dbname}")
            self.dbs.setdefault(id, {})[dbname] = {
                "dir": odir,
                "dbname": dbname,
                "db": odir + dbname,
                "status": "failed",
            }
            return False

    def process_kuniq_files(self, dbname, odir):

        def check_map_orig_processed(map_orig_file):

            if os.path.exists(map_orig_file) is False:
                return False

            seqid = pd.read_csv(map_orig_file, sep="\t")
            if seqid.shape[1] == 4:
                if "GTDB" in seqid.columns and "description" in seqid.columns:
                    return True

            return False

        def process_map_orig(map_orig_file, out_map_orig):

            if os.path.exists(map_orig_file) is False:
                return

            seqid = pd.read_csv(map_orig_file, sep="\t", header=None)

            seqid.columns = ["refseq", "taxid", "merge"]

            def split_merge(x: str):
                if x is None:
                    return ["", ""]
                else:
                    x = x.split(" ")
                    if len(x) == 1:
                        return [x[0], ""]
                    else:
                        return [x[0], " ".join(x[1:])]

            seqid[["GTDB", "description"]] = seqid["merge"].apply(
                lambda x: pd.Series(split_merge(x))
            )

            if "merge" in seqid.columns:
                new_columns = [x for x in seqid.columns if x != "merge"]
                seqid = seqid[new_columns]

            seqid.to_csv(
                out_map_orig,
                sep="\t",
                header=True,
                index=False,
            )

        ############
        ############

        def check_map_file(map_file):
            if os.path.exists(map_file) is False:
                return False

            seqmap = pd.read_csv(map_file, sep="\t")
            if seqmap.shape[1] == 2:
                if "acc" in seqmap.columns and "protid" in seqmap.columns:
                    return True

            return False

        def process_map_file(map_file, out_map_file):
            if os.path.exists(map_file) is False:
                return

            seqmap = pd.read_csv(map_file, sep="\t")
            seqmap.columns = ["acc", "protid"]
            seqmap.to_csv(
                out_map_file,
                sep="\t",
                header=True,
                index=False,
            )

        ####################
        ####################

        map_orig_file = f"{odir + dbname}/seqid2taxid.map.orig"
        out_map_orig = f"{odir + dbname}/seqid2taxid.map.orig"

        if check_map_orig_processed(map_orig_file) is False:
            process_map_orig(map_orig_file, out_map_orig)

        ######
        map_file = f"{odir + dbname}/seqid2taxid.map"
        out_map_file = f"{self.metadir}/protein_acc2protid.tsv"

        if check_map_file(map_file) is False:
            process_map_file(map_file, out_map_file)

    def kuniq_install(
        self,
        id="krakenuniq",
        dbdir="kuniq",
        dbname="viral",
        threads="6",
        dl_args="--force --min-seq-len 300 --dust",
        build_args="--work-on-disk --jellyfish-hash-size 10M --kmer-len 31 --taxids-for-genomes --taxids-for-sequences",
    ):
        odir = self.dbdir + dbdir + "/"
        bin = self.envs["ROOT"] + self.envs[id] + "/bin/"

        if os.path.isfile(odir + dbname + "/taxDB"):
            logging.info(f"Krakenuniq {dbname} taxDB present. prepped.")
            self.process_kuniq_files(dbname, odir)
            self.dbs.setdefault(id, {})[dbname] = {
                "dir": odir,
                "dbname": dbname,
                "db": odir + dbname,
            }

            if not os.path.isfile(f"{self.metadir}/protein_acc2protid.tsv"):

                self.process_kuniq_files(dbname, odir)
                # seqmap = pd.read_csv(f"{odir + dbname}/seqid2taxid.map", sep="\t")
                # seqmap.columns = ["acc", "protid"]
                # seqmap.to_csv(
                #    f"{self.metadir}/protein_acc2protid.tsv",
                #    sep="\t",
                #    header=True,
                #    index=False,
                # )

            return True
        else:
            if self.test:
                logging.info(f"Krakenuniq {dbname} taxDB not present.")
                return False
            else:
                logging.info(
                    f"Krakenuniq {dbname} viral is not installed. Installing..."
                )

        tax_command = [
            bin + "krakenuniq-download",
            "--db",
            odir + dbname,
            "--threads",
            threads,
            "taxonomy",
        ]
        lib_command = [
            bin + "krakenuniq-download",
            "--db",
            odir + dbname,
            "--threads",
            threads,
            dl_args,
            f"refseq/{dbname}",
        ]

        build_command = [
            bin + "krakenuniq-build",
            "--db",
            odir + dbname,
            "--threads",
            threads,
            "--jellyfish-bin",
            bin + "jellyfish",
            build_args,
        ]

        try:
            subprocess.run(["mkdir", "-p", odir])

            try:
                subprocess.run(tax_command)
                untax_get(self.taxdump, odir, dbname)
                #
                subprocess.run(" ".join(lib_command), shell=True)

                subprocess.call(" ".join(build_command), shell=True)

            except subprocess.CalledProcessError:
                logging.error("failed to install krakenuniq db")

            self.process_kuniq_files(dbname, odir)

            self.dbs.setdefault(id, {})[dbname] = {"dir": odir, "dbname": dbname, "db": odir + dbname}

            return True
        except subprocess.CalledProcessError:
            logging.warning(f"failed to download Krakenuniq db {dbname}")
            return False

    def kaiju_dl_install(
        self,
        id="kaiju",
        dbdir="kaiju",
        dbname="viral",
    ):
        db_online = get_db_url('kaiju', dbname)
        
        if not db_online:
            logging.error(f"No source URL found for kaiju/{dbname}")
            return False

        odir = self.dbdir + dbdir + "/"
        bin = self.envs["ROOT"] + self.envs[id] + "/bin/"
        subdb = odir + dbname + "/"

        fmi = self._find_kaiju_fmi(subdb)
        if fmi:
            logging.info(f"Kaiju {dbname}  is installed.")
            self.dbs.setdefault(id, {})[dbname] = {
                "dir": odir,
                "dbname": dbname,
                "db": fmi,
                "status": "success",
            }
            return True
        else:
            if self.test:
                logging.info(f"Kaiju {dbname} db is not installed.")
                self.dbs.setdefault(id, {})[dbname] = {
                    "dir": subdb,
                    "dbname": dbname,
                    "db": "",
                    "status": "test",
                }
                return False
            else:
                logging.info(f"Kaiju {dbname} db is not installed. Installing...")

        try:
            subprocess.run(["mkdir", "-p", odir])

            subprocess.run(
                ["wget", "-P", subdb, db_online, "--no-check-certificate"], check=True
            )
            CWD = os.getcwd()
            os.chdir(subdb)
            subprocess.run(["tar", "-zxvf", os.path.basename(db_online)], check=True)
            subprocess.run(["rm", os.path.basename(db_online)], check=True)
            os.chdir(CWD)

            fmi = self._find_kaiju_fmi(subdb)
            if not fmi:
                logging.error(f"No .fmi index found in Kaiju db {dbname} at {subdb}")
                self.dbs.setdefault(id, {})[dbname] = {
                    "dir": subdb,
                    "dbname": dbname,
                    "db": "",
                    "status": "failed",
                }
                return False

            self.dbs.setdefault(id, {})[dbname] = {
                "dir": subdb,
                "dbname": dbname,
                "db": fmi,
                "status": "success",
            }
            return True
        except subprocess.CalledProcessError:
            logging.warning(f"failed to download Kaiju db {dbname}")
            self.dbs.setdefault(id, {})[dbname] = {
                "dir": subdb,
                "dbname": dbname,
                "db": "",
                "status": "failed",
            }
            return False

    @staticmethod
    def _find_kaiju_fmi(subdb: str) -> str:
        """Return the path of the first *.fmi file in subdb, or '' if none."""
        if not os.path.isdir(subdb):
            return ""
        fmis = sorted(f for f in os.listdir(subdb) if f.endswith(".fmi"))
        return os.path.join(subdb, fmis[0]) if fmis else ""

    def kaiju_viral_install(self, id="kaiju", dbdir="kaiju", dbname="viral"):
        odir = self.dbdir + dbdir + "/"
        bin = self.envs["ROOT"] + self.envs[id] + "/bin/"
        subdb = odir + dbname + "/"
        if dbname == "viral":
            db_online = (
                "https://kaiju.binf.ku.dk/database/kaiju_db_viruses_2021-02-24.tgz"
            )
        elif dbname == "fungi":
            db_online = "https://kaiju-idx.s3.eu-central-1.amazonaws.com/2023/kaiju_db_fungi_2023-05-26.tgz"
        file = os.path.basename(db_online)
        fmi = self._find_kaiju_fmi(subdb)
        if fmi:
            logging.info(f"Kaiju {dbname}  is installed.")
            self.dbs.setdefault(id, {})[dbname] = {
                "dir": odir,
                "dbname": dbname,
                "db": fmi,
            }
            return True
        else:
            if self.test:
                logging.info(f"Kaiju {dbname} db is not installed.")

                return False
            else:
                logging.info(f"Kaiju {dbname} db is not installed. Installing...")

        try:
            subprocess.run(["mkdir", "-p", odir])

            subprocess.run(["wget", "-P", subdb, db_online, "--no-check-certificate"])
            CWD = os.getcwd()
            os.chdir(subdb)
            subprocess.run(["tar", "-zxvf", file])
            subprocess.run(["rm", file])
            os.chdir(CWD)

            fmi = self._find_kaiju_fmi(subdb)
            if not fmi:
                logging.error(f"No .fmi index found in Kaiju db {dbname} at {subdb}")
                self.dbs.setdefault(id, {})[dbname] = {
                    "dir": subdb,
                    "dbname": dbname,
                    "db": "",
                    "status": "failed",
                }
                return False

            self.dbs.setdefault(id, {})[dbname] = {
                "dir": subdb,
                "dbname": dbname,
                "db": fmi,
            }
            return True
        except subprocess.CalledProcessError:
            logging.warning(f"failed to download Kaiju db {dbname}")
            return False

    def blast_install(
        self,
        id="blast",
        dbdir="blast",
        reference="",
        dbname="viral",
        nuc=True,
        taxid_map="",
        args="-parse_seqids",
        title="viral db",
    ):
        odir = self.dbdir + dbdir + "/"
        dbtype = "nucl"
        sdir = "NUC/"

        if not nuc:
            dbtype = "prot"
            sdir = "PROT/"
        sdir = odir + sdir

        bin = self.envs["ROOT"] + self.envs[id] + "/bin/"
        db = sdir + dbname

        if self.update:
            if os.path.exists(f"{sdir}"):
                shutil.rmtree(f"{sdir}")

        if os.path.isfile(db + f".{dbtype[0]}db"):
            logging.info(f"blast index for {dbname} is installed.")
            self.dbs.setdefault(id, {})[dbname] = {"dir": odir, "dbname": dbname, "db": db, "status": "success"}
            return True
        else:
            if self.test:
                logging.info(f"blast index for {dbname} is not installed.")
                self.dbs.setdefault(id, {})[dbname] = {"dir": odir, "dbname": dbname, "db": db, "status": "test"}
                return False
            else:
                logging.info(
                    f"blast index for {dbname} is not installed. Installing..."
                )

        try:
            subprocess.run(["mkdir", "-p", odir])

            gzipped = False
            if reference[-3:] == ".gz":
                gzipped = True
                reference_unzip = os.path.splitext(reference)[0]
                if os.path.exists(reference_unzip):
                    subprocess.run(["rm", reference_unzip])

                subprocess.run(["gunzip", reference])
                reference = reference_unzip

            commands = [
                bin + "makeblastdb",
                "-in",
                reference,
                "-out",
                db,
                "-dbtype",
                dbtype,
                "-title",
                title,
                args,
            ]

            if taxid_map:
                commands += ["-taxid_map", taxid_map]

            try:
                subprocess.run(commands)
            finally:
                if gzipped:
                    subprocess.run([BGZIP_BIN, reference])
                    reference = reference + ".gz"

            self.dbs.setdefault(id, {})[dbname] = {"dir": sdir, "dbname": dbname, "db": db, "status": "success"}

            return True

        except subprocess.CalledProcessError:
            logging.warning(f"failed to download blast db {dbname}")
            self.dbs.setdefault(id, {})[dbname] = {"dir": sdir, "dbname": dbname, "db": db, "status": "failed"}
            return False

    def bowtie2_index(
        self,
        id="bowtie2",
        dbdir="bowtie2",
        reference="",
        dbname="viral",
    ):
        odir = self.dbdir + dbdir + "/"
        bin = self.envs["ROOT"] + self.envs[id] + "/bin/"
        db = odir + dbname
        if os.path.isfile(db + ".1.bt2"):
            logging.info(f"bowtie2 index for {dbname} is installed.")
            self.dbs.setdefault(id, {})[dbname] = {"dir": odir, "dbname": dbname, "db": db, "status": "success"}
            return True
        else:
            if self.test:
                logging.info(f"bowtie2 index for {dbname} is not installed.")
                self.dbs.setdefault(id, {})[dbname] = {"dir": odir, "dbname": dbname, "db": db, "status": "test"}

                return False
            else:
                logging.info(
                    f"bowtie2 index for {dbname} is not installed. Installing..."
                )

        try:
            subprocess.run(["mkdir", "-p", odir])

            gzipped = False
            if reference[-3:] == ".gz":
                gzipped = True
                subprocess.run(["gunzip", reference])
                reference = os.path.splitext(reference)[0]

            commands = [
                bin + "bowtie2-build",
                reference,
                db,
            ]

            try:
                subprocess.run(commands)
            finally:
                if gzipped:
                    subprocess.run([BGZIP_BIN, reference])
                    reference = reference + ".gz"

            self.dbs.setdefault(id, {})[dbname] = {"dir": odir, "dbname": dbname, "db": db, "status": "success"}

            return True

        except subprocess.CalledProcessError:
            logging.warning(f"failed to download bowtie2 index {dbname}")
            self.dbs.setdefault(id, {})[dbname] = {"dir": odir, "dbname": dbname, "db": db, "status": "failed"}
            return False

    def bwa_install(
        self,
        dbname="bwa",
        url="",
        reference="",
        id="bwa",
        dbdir="bwa",
        dlp="wget",
        update=False,
    ):
        """ """

        odir = self.dbdir + dbdir + "/"
        bin = self.envs["ROOT"] + self.envs[id] + "/bin/"
        sdir = odir + dbname + "/" + dbname

        if update:
            if os.path.exists(f"{odir}{dbname}"):
                shutil.rmtree(f"{odir}{dbname}")

        if not url and not reference:
            logging.info(
                "Please provide either sequence or url for bwa install. Skipping."
            )
            return False

        if os.path.isfile(f"{odir}{dbname}/{dbname}.bwt"):
            logging.info(f"BWA db {dbname} is installed.")
            self.dbs.setdefault(id, {})[dbname] = {
                "dir": odir,
                "dbname": dbname,
                "fasta": f"{sdir}.fa",
                "db": f"{odir}{dbname}/{dbname}",
                "status": "success",
            }
            return True

        else:
            if self.test:
                logging.info(f"BWA db {dbname} is not installed.")
                self.dbs.setdefault(id, {})[dbname] = {
                    "dir": odir,
                    "dbname": dbname,
                    "fasta": f"{sdir}.fa",
                    "db": f"{odir}{dbname}/{dbname}",
                    "status": "test",
                }
                return False
            else:
                logging.info(f"BWA db {dbname} is not installed. Installing...")

        if not verify_file_accessible(reference):
            logging.error(f"BWA install failed: reference file is inaccessible: {reference}")
            self.dbs.setdefault(id, {})[dbname] = {
                "dir": odir,
                "dbname": dbname,
                "fasta": f"{sdir}.fa",
                "db": f"{odir}{dbname}/{dbname}",
                "status": "failed",
            }
            return False

        subprocess.run(["mkdir", "-p", sdir])

        gzipped = False

        if reference[-3:] == ".gz":
            gzipped = True
            if not verify_gzip_integrity(reference):
                logging.error(
                    f"BWA install failed: gzipped reference file is corrupted or incomplete: {reference}"
                )
                self.dbs.setdefault(id, {})[dbname] = {
                    "dir": odir,
                    "dbname": dbname,
                    "fasta": f"{sdir}.fa",
                    "db": f"{odir}{dbname}/{dbname}",
                    "status": "failed",
                }
                return False
            subprocess.run(["gunzip", reference])
            reference = os.path.splitext(reference)[0]
            if not verify_fasta_integrity(reference):
                logging.error(
                    f"BWA install failed: decompressed FASTA file is invalid: {reference}"
                )
                if gzipped:
                    subprocess.run([BGZIP_BIN, reference], capture_output=True)
                self.dbs.setdefault(id, {})[dbname] = {
                    "dir": odir,
                    "dbname": dbname,
                    "fasta": f"{sdir}.fa",
                    "db": f"{odir}{dbname}/{dbname}",
                    "status": "failed",
                }
                return False
            subprocess.run(["samtools", "faidx", reference])

        command = [bin + "bwa", "index", "-p", f"{odir}{dbname}/{dbname}", reference]
        # command = " ".join(command)

        try:
            subprocess.run(command)
            shutil.copy(reference, f"{odir}{dbname}/{dbname}.fa")
            self.dbs.setdefault(id, {})[dbname] = {
                "dir": odir,
                "dbname": dbname,
                "fasta": f"{sdir}.fa",
                "db": f"{odir}{dbname}/{dbname}",
                "status": "success",
            }
            return True
        except subprocess.CalledProcessError as e:
            logging.error(f"BWA index failed for {dbname} (exit code {e.returncode})")
            logging.error(f"Reference file: {reference}")
            import traceback
            traceback.print_exc()
            self.dbs.setdefault(id, {})[dbname] = {
                "dir": odir,
                "dbname": dbname,
                "fasta": f"{sdir}.fa",
                "db": f"{odir}{dbname}/{dbname}",
                "status": "failed",
            }
            return False
        except FileNotFoundError as e:
            logging.error(f"BWA install failed: required tool not found: {e}")
            self.dbs.setdefault(id, {})[dbname] = {
                "dir": odir,
                "dbname": dbname,
                "fasta": f"{sdir}.fa",
                "db": f"{odir}{dbname}/{dbname}",
                "status": "failed",
            }
            return False
        finally:
            if gzipped:
                subprocess.run([BGZIP_BIN, reference], capture_output=True)
                reference = reference + ".gz"

    def minimap2_install(self, id="minimap2", dbname="", reference=""):
        """
        Record minimap2 reference - no index needed.
        """
        odir = self.dbdir + id + "/"
        
        if os.path.isfile(reference):
            self.dbs.setdefault(id, {})[dbname] = {
                "dir": odir,
                "dbname": dbname,
                "db": reference,
                "fasta": reference,
                "status": "success",
            }
            return True
        elif self.test:
            self.dbs.setdefault(id, {})[dbname] = {
                "dir": odir,
                "dbname": dbname,
                "db": reference,
                "fasta": reference,
                "status": "test",
            }
            return False
        else:
            self.dbs.setdefault(id, {})[dbname] = {
                "dir": odir,
                "dbname": dbname,
                "db": reference,
                "fasta": reference,
                "status": "failed",
            }
            return False

    def virsorter_install(
        self, id="virsorter", dbdir="virsorter", dbname="viral", threads="4"
    ):
        """
        install virsorter
        :param id:
        :param dbdir:
        :param dbname:
        :param threads:
        :return:
        """
        odir = self.dbdir + dbdir + "/"
        bin = self.envs["ROOT"] + self.envs[id] + "/bin/"
        if os.path.isfile(odir + "Done_all_setup"):
            logging.info("Virsorter is installed.")
            return True
        else:
            if self.test:
                logging.info("Virsorter is not installed.")
                return False
            else:
                logging.info("Virsorter is not installed. Installing...")

        commands = [
            bin + "virsorter",
            "setup",
            "-d",
            odir,
            "-j",
            threads,
        ]

        tmpsh = "virsorter_install.sh"

        bash_lines = [
            "#!/bin/bash",
            f"source {self.source}",
            f"conda activate {self.envs['ROOT'] + self.envs[id]}",
            " ".join(commands),
            "conda deactivate",
        ]
        try:
            subprocess.run(["mkdir", "-p", odir])

            os.system("touch " + tmpsh)
            with open(tmpsh, "w") as f:
                for l in bash_lines:
                    os.system('echo "{}" >> {}'.format(l, tmpsh))
            #                f.write("/n".join(bash_lines))

            subprocess.run(["chmod", "+x", tmpsh])
            subprocess.call(f"./{tmpsh}")

            os.system("rm " + tmpsh)

            self.dbs.setdefault(id, {})[dbname] = {"dir": odir, "dbname": dbname}
            return True

        except subprocess.CalledProcessError:
            logging.info(f"failed to install virsorter")
            return False

    def fve_install(
        self,
        id="fastviromeexplorer",
        dbdir="fve",
        dbname="viral",
        virus_list="viruslist.txt",
        list_create=False,
        reference="",
    ):
        """
        install Fast Virome Explorer
        :param id:
        :param dbdir:
        :param dbname:
        :param threads:
        :return:
        """
        odir = self.dbdir + dbdir + "/"
        bin = self.envs["ROOT"] + self.envs[id] + "/bin/"
        subdir = odir + dbname + "/"
        fidx = subdir + dbname + ".idx"

        if self.update:
            if os.path.exists(subdir):
                shutil.rmtree(subdir)

        if os.path.isfile(fidx):

            logging.info(f"FastViromeExplorer index for {reference} is installed.")
            self.dbs.setdefault(id, {})[dbname] = {"dir": odir, "dbname": dbname, "db": fidx}
            return True
        else:
            if self.test:
                logging.info(f"FastViromeExplorer {reference} index is not installed.")
                self.dbs.setdefault(id, {})[dbname] = {"dir": odir, "dbname": dbname, "db": fidx}
                return False
            else:
                logging.info(
                    f"FastViromeExplorer {reference} index is not installed. Installing..."
                )

        try:
            subprocess.run(["mkdir", "-p", subdir])

            gzipped = False
            if reference[-3:] == ".gz":
                gzipped = True
                subprocess.run(["gunzip", reference])
                reference = os.path.splitext(reference)[0]

            genlistbin = (
                self.envs["ROOT"]
                + self.envs[id]
                + "/utility-scripts/"
                + "generateGenomeList.sh"
            )

            comm_vlist = [
                genlistbin,
                reference,
                virus_list,
            ]

            com_kallisto = [
                self.envs["ROOT"] + self.envs["kallisto"] + "/bin/kallisto",
                "index",
                "-i",
                fidx,
                reference,
                "--make-unique",
            ]

            try:
                if list_create or not os.path.exists(virus_list):
                    os.system(f"chmod +x {genlistbin}")
                    subprocess.call(" ".join(comm_vlist), shell=True)

                subprocess.run(com_kallisto)
                os.system(f"mv {virus_list} {subdir}")

            finally:
                if gzipped:
                    subprocess.run([BGZIP_BIN, reference])
                    reference = reference + ".gz"

            self.dbs.setdefault(id, {})[dbname] = {"dir": subdir, "dbname": dbname, "db": fidx}
            return True

        except subprocess.CalledProcessError:
            logging.info(f"failed to install FastViromeExplorer")
            return False

    def deSAMBA_install(
        self, id="desamba", dbdir="desamba", dbname="viral", reference=""
    ):
        """
        install virsorter
        :param id:
        :param dbdir:
        :param dbname:
        :param threads:
        :return:
        """
        odir = self.dbdir + dbdir + "/"
        sdir = odir + dbname
        bin = self.envs["ROOT"] + self.envs[id]

        if os.path.isfile(sdir + "/deSAMBA.bwt"):
            logging.info(f"deSAMBA db {dbname} is installed.")
            self.dbs.setdefault(id, {})[dbname] = {"dir": sdir, "dbname": dbname, "db": sdir}
            return True
        else:
            if self.test:
                logging.info(f"deSAMBA db {dbname} is not installed.")
                self.dbs.setdefault(id, {})[dbname] = {"dir": odir, "dbname": dbname, "db": sdir}
                return False
            else:
                logging.info(f"deSAMBA db {dbname} is not installed. Installing...")

        try:
            subprocess.run(["mkdir", "-p", odir])

            gzipped = False
            if reference[-3:] == ".gz":
                gzipped = True
                subprocess.run(["gunzip", reference])
                reference = os.path.splitext(reference)[0]

            build_command = [bin + "/build-index", reference, sdir]
            try:
                CWD = os.getcwd()
                os.chdir(bin)
                subprocess.call(" ".join(build_command), shell=True)
                os.chdir(CWD)

            finally:
                if gzipped:
                    subprocess.run([BGZIP_BIN, reference])
                    reference = reference + ".gz"

            self.dbs.setdefault(id, {})[dbname] = {"dir": odir, "dbname": dbname, "db": sdir}
            return True

        except subprocess.CalledProcessError:
            logging.info(f"failed to install deSAMBA")
            return False
