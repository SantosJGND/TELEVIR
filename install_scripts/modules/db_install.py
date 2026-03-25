#!/usr/bin/python3

import datetime
import gzip
import logging
import os
import shutil
import subprocess
from ftplib import FTP
from pathlib import Path
from random import randint
from threading import Thread

import pandas as pd
from install_scripts.host_library import Host
from install_scripts.load_sources import get_db_url, get_db_version

from fastq_filter import file_to_fastq_records

from numpy import int0
from xopen import xopen

from decouple import config

BGZIP_BIN = "bgzip"
try:
    BGZIP_BIN = config("BGZIP_BIN")
except:
    pass


def grep_sequence_identifiers(str_input, output, ignore=""):
    """
    grep sequence identifiers from fasta file.
    """

    get_pattern = f"zgrep -P '^>' {str_input}"
    filter_pattern = f"grep -v {ignore}"
    process_pattern = "sed 's/^>//; s/[ ].*$//g'"

    if ignore:
        command = "{} | {} | {}".format(get_pattern, filter_pattern, process_pattern)
    else:
        command = "{} | {}".format(get_pattern, process_pattern)

    command = command + " > {}".format(output)

    os.system(command)


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


def sed_out_after_dot(file):
    """remove everything after the dot in the file name"""
    os.system("sed -i 's/[:].*$//g' {}".format(file))


def entrez_ncbi_taxid_command(lines, tempfile, outdir, outfile):
    Path(tempfile).touch()

    with open(tempfile, "w") as ftemp:
        ftemp.write(lines)

    os.system(
        f"cat {tempfile} | epost -db nuccore | esummary -db nuccore | xtract -pattern DocumentSummary -element AccessionVersion,TaxId >> {outdir}{outfile}"
    )


def entrez_fetch_sequence(accid, outfile):
    """return fasta from ncbi nuccore db using accid"""

    os.system(f"esearch -db nuccore -query {accid} | efetch -format fasta >> {outfile}")


def entrez_ncbi_taxid(file, outdir, outfile, nmax=500):
    """get taxids from ncbi accessions in single column file. return both, using Entrez Utilities"""

    Path(f"{outdir}{outfile}").touch()

    tempdir = os.path.dirname(outfile)
    tempfile = os.path.join(tempdir, f"temp{randint(10000, 99999)}.txt")
    current_batch = 0
    lines = ""

    with open(file, "r") as f:
        line = f.readline()
        while line:
            lines += line
            current_batch += 1

            if current_batch == nmax:
                entrez_ncbi_taxid_command(lines, tempfile, outdir, outfile)
                lines = ""
                current_batch = 0

            line = f.readline()

    if lines:
        entrez_ncbi_taxid_command(lines, tempfile, outdir, outfile)

    if os.path.exists(tempfile):
        os.remove(tempfile)


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

        self.organism = organism

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

    def refseq_16s_dl(self, fname="refseq_16s"):
        """
        download 16s from https://ftp.ncbi.nlm.nih.gov/refseq/TargetedLoci/Bacteria/bacteria.16SrRNA.fna.gz
        save to fastas["filter"]["16s"]
        """
        host = "ftp.ncbi.nlm.nih.gov"
        source = "refseq/TargetedLoci/Bacteria/bacteria.16SrRNA.fna.gz"
        filename = "bacteria.16SrRNA.fna.gz"
        filepath = self.seqdir + filename
        source_url = f"https://{host}/{source}"

        if os.path.isfile(filepath):
            if self.verify_file_integrity(filepath, filename):
                self.fastas["filter"][fname] = filepath
                self.db_versions[fname] = {
                    "version": self.get_file_mod_date(filepath),
                    "source_url": source_url,
                    "file_mod_date": self.get_file_mod_date(filepath)
                }
                logging.info(f"{filename} found and verified.")
                return True
            logging.warning(f"{filename} exists but is corrupted. Re-downloading...")

        try:
            subprocess.run(
                ["wget", f"https://{host}/{source}", "-P", self.seqdir],
                check=False,
            )
        except subprocess.CalledProcessError:
            logging.info(f"{filename} not found.")
            return False

        if self.verify_file_integrity(filepath, filename):
            self.fastas["filter"][fname] = filepath
            self.db_versions[fname] = {
                "version": self.get_file_mod_date(filepath),
                "source_url": source_url,
                "file_mod_date": self.get_file_mod_date(filepath)
            }
            return True
        else:
            logging.error(f"Downloaded {filename} is corrupted.")
            return False

    def silva_16s_dl(self, fname="silva_16s"):
        """
        download silva 16s from https://www.arb-silva.de/fileadmin/silva_databases/release_138_1/Exports/SILVA_138.1_LSURef_NR99_tax_silva.fasta.gz
        save to fastas["filter"]["silva"]
        """
        host = "www.arb-silva.de"
        source = "fileadmin/silva_databases/release_138_1/Exports/"
        filename = "SILVA_138.1_LSURef_NR99_tax_silva.fasta.gz"
        filepath = self.seqdir + filename
        source_url = f"https://{host}/{source}{filename}"

        if os.path.isfile(filepath):
            if self.verify_file_integrity(filepath, filename):
                self.fastas["filter"][fname] = filepath
                self.db_versions[fname] = {
                    "version": self.get_file_mod_date(filepath),
                    "source_url": source_url,
                    "file_mod_date": self.get_file_mod_date(filepath)
                }
                logging.info(f"{filename} found and verified.")
                return True
            logging.warning(f"{filename} exists but is corrupted. Re-downloading...")

        try:
            subprocess.run(
                ["wget", f"https://{host}/{source}{filename}", "-P", self.seqdir],
                check=False,
            )
        except subprocess.CalledProcessError:
            logging.info(f"{filename} not found.")
            return False

        if self.verify_file_integrity(filepath, filename):
            self.fastas["filter"][fname] = filepath
            self.db_versions[fname] = {
                "version": self.get_file_mod_date(filepath),
                "source_url": source_url,
                "file_mod_date": self.get_file_mod_date(filepath)
            }
            return True
        else:
            logging.error(f"Downloaded {filename} is corrupted.")
            return False

    def ncbi_16s_dl(self, fname="ncbi_ribo16s"):
        """
        donwload 16s from https://ftp.ncbi.nlm.nih.gov/blast/db/16S_ribosomal_RNA.tar.gz
        save to fastas["filter"]["16s"]
        """
        host = "ftp.ncbi.nlm.nih.gov"
        source = "blast/db/16S_ribosomal_RNA.tar.gz"
        filename = "16S_ribosomal_RNA.tar.gz"
        filepath = self.seqdir + filename

        if os.path.isfile(filepath):
            if self.verify_file_integrity(filepath, filename):
                self.fastas["filter"][fname] = filepath
                logging.info(f"{filename} found and verified.")
                return True
            logging.warning(f"{filename} exists but is corrupted. Re-downloading...")

        try:
            subprocess.run(
                ["wget", f"ftp://{host}/{source}", "-P", self.seqdir],
                check=False,
            )
        except subprocess.CalledProcessError:
            logging.info(f"{filename} not found.")
            return False

        if self.verify_file_integrity(filepath, filename):
            self.fastas["filter"][fname] = filepath
            return True
        else:
            logging.error(f"Downloaded {filename} is corrupted.")
            return False

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
                self.fastas["host"][fname] = filepath
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
                self.fastas["host"][fname] = filepath
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
                fpath = self.fastas["host"].get(host_name)
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
        references_file = os.path.join(self.seqdir, "request_references.fa")

        if os.path.isfile(references_file):
            self.fastas["nuc"]["requests"] = [references_file]
            self.db_versions["requests"] = {
                "version": self.get_file_mod_date(references_file),
                "source_url": "user_provided",
                "file_mod_date": self.get_file_mod_date(references_file)
            }
            return True

        if os.path.isfile(references_file + ".gz"):
            self.fastas["nuc"]["requests"] = [references_file + ".gz"]
            self.db_versions["requests"] = {
                "version": self.get_file_mod_date(references_file + ".gz"),
                "source_url": "user_provided",
                "file_mod_date": self.get_file_mod_date(references_file + ".gz")
            }
            return True

        seq_file = self.requests.get("FILE", "")
        if seq_file and os.path.isfile(seq_file):
            logging.info(f"Using provided request sequences file: {seq_file}")
            shutil.copy(seq_file, references_file)
            if not references_file.endswith(".gz"):
                os.system(f"{BGZIP_BIN} {references_file}")
                references_file = references_file + ".gz"
            
            self.fastas["nuc"]["requests"] = [references_file]
            self.db_versions["requests"] = {
                "version": self.get_file_mod_date(references_file),
                "source_url": seq_file,
                "file_mod_date": self.get_file_mod_date(references_file)
            }
            logging.info("request sequences from file prepped.")
            return True

        if self.requests.get("ACCID"):
            acc_tsv = self.requests["ACCID"]
            tempfile = os.path.join(self.metadir, "accession_requests.tsv")
            if os.path.isfile(acc_tsv):
                shutil.copy(acc_tsv, tempfile)
            elif "https" in acc_tsv:
                try:
                    subprocess.run(["curl", "-o", tempfile, acc_tsv])

                except subprocess.CalledProcessError:
                    logging.error("accid request file not found")
                    return False

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
                "file_mod_date": self.get_file_mod_date(references_file)
            }
            logging.info("request sequences prepped.")
            return True
        else:
            return False

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
                self.fastas["prot"][db_key] = fprot_path
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
            self.fastas["prot"][db_key] = fprot_path
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

    def uniprot_dl(self, vs="90"):
        """
        download uniprot db.
        :param vs: uniprot version.
        :return:
        """
        fl = "https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref90/uniref{}.fasta.gz"
        fl = fl.format(vs)
        uniprot_file = os.path.basename(fl)
        filepath = self.seqdir + uniprot_file

        if os.path.isfile(filepath):
            if self.verify_file_integrity(filepath, uniprot_file):
                logging.info("uniref{}.fasta found and verified.".format(vs))
            else:
                logging.warning("uniref{}.fasta exists but is corrupted. Re-downloading...".format(vs))
        elif self.test:
            logging.info("uniref{}.fasta not found.".format(vs))
            return False
        else:
            logging.info("uniref{}.fasta not found. downloading...".format(vs))
            subprocess.run(["wget", fl, "-P", self.seqdir])

        if self.verify_file_integrity(filepath, uniprot_file):
            self.fastas["prot"]["uniprot"] = filepath
            return True
        else:
            logging.error("Downloaded uniref{}.fasta is corrupted.".format(vs))
            return False

    def swissprot_dl(self):
        """
        download swissprot db from ncbi.
        :return:
        """
        fl = "https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/swissprot.gz"
        basename_fl = os.path.basename(fl)
        filepath = self.seqdir + basename_fl

        if self.update:
            if os.path.isfile(filepath):
                os.remove(filepath)

        if os.path.isfile(filepath):
            if self.verify_file_integrity(filepath, basename_fl):
                logging.info("swissprot.gz found and verified.")
            else:
                logging.warning("swissprot.gz exists but is corrupted. Re-downloading...")
        elif self.test:
            logging.info("swissprot.gz not found.")
            return False
        else:
            logging.info("swissprot.gz not found. downloading...")
            subprocess.run(["wget", "--quiet", "-P", self.seqdir, fl])

        if self.verify_file_integrity(filepath, basename_fl):
            self.fastas["prot"]["swissprot"] = filepath
            self.db_versions["swissprot"] = {
                "version": self.get_file_mod_date(filepath),
                "source_url": fl,
                "file_mod_date": self.get_file_mod_date(filepath)
            }
            return True
        else:
            logging.error("Downloaded swissprot.gz is corrupted.")
            return False

    def RVDB_dl(self, vs="22.0"):
        """
        download rvdb.
        :param vs: rvdb version.
        :return:
        """
        fl = "https://rvdb-prot.pasteur.fr/files/U-RVDBv{}-prot.fasta.xz"
        fl = fl.format(vs)
        basename_fl = os.path.basename(fl).replace(".xz", ".gz")
        filepath = self.seqdir + basename_fl

        if os.path.isfile(filepath):
            if self.verify_file_integrity(filepath, basename_fl):
                logging.info("U-RVDBv{}.fasta.xz found and verified.".format(vs))
            else:
                logging.warning("U-RVDBv{}.fasta.xz exists but is corrupted. Re-downloading...".format(vs))
        elif self.test:
            logging.info("U-RVDBv{}.fasta.xz not found.".format(vs))
            return False
        else:
            logging.info("U-RVDBv{}.fasta.xz not found. downloading...".format(vs))
            xz_path = self.seqdir + os.path.basename(fl)
            subprocess.run(["wget", "-P", self.seqdir, fl])
            if os.path.isfile(xz_path):
                if self.verify_file_integrity(xz_path, os.path.basename(fl)):
                    subprocess.run(["unxz", xz_path])
                    unxz_path = os.path.splitext(xz_path)[0]
                    subprocess.run([BGZIP_BIN, unxz_path])
                    filepath = unxz_path + ".gz"
                else:
                    logging.error("Downloaded U-RVDBv{}.fasta.xz is corrupted.".format(vs))
                    return False
            else:
                logging.error("Download of U-RVDBv{}.fasta.xz failed.".format(vs))
                return False

        if self.verify_file_integrity(filepath, basename_fl):
            self.fastas["prot"]["rvdb"] = filepath
            self.db_versions["rvdb"] = {
                "version": vs,
                "source_url": fl,
                "file_mod_date": self.get_file_mod_date(filepath)
            }
            return True
        else:
            logging.error("Downloaded U-RVDBv{}.fasta is corrupted.".format(vs))
            return False

    def virosaurus_dl(self):
        """
        download virossaurus from viralzone.
        :return:
        """

        fl = "https://ftp.expasy.org/databases/viralzone/2020%5F4/virosaurus90%5Fvertebrate-20200330.fas.gz"
        basename_fl = os.path.basename(fl)
        filepath = self.seqdir + basename_fl

        if os.path.isfile(filepath) and os.path.getsize(filepath) > 100:
            if self.verify_file_integrity(filepath, basename_fl):
                logging.info("virosaurus90_vertebrate_20200330.fas found and verified.")
            else:
                logging.warning("virosaurus90_vertebrate_20200330.fas exists but is corrupted. Re-downloading...")
        elif self.test:
            logging.info("virosaurus90_vertebrate_20200330.fas not found.")
            return False
        else:
            try:
                logging.info(
                    "virosaurus90_vertebrate_20200330.fas not found. downloading..."
                )
                subprocess.run(
                    ["wget", "--no-check-certificate", fl, "-P", self.seqdir]
                )
            except subprocess.CalledProcessError:
                logging.info("virosaus download failed.")
                return False

            if not os.path.isfile(filepath) or os.path.getsize(filepath) <= 100:
                logging.info("wget failed. trying curl...")
                subprocess.run(
                    ["curl", fl, "-o", filepath]
                )

            if not os.path.isfile(filepath) or os.path.getsize(filepath) <= 100:
                logging.info("virosaurus download failed.")
                return False

        if self.verify_file_integrity(filepath, basename_fl):
            self.fastas["nuc"]["virosaurus"] = [filepath]
            self.db_versions["virosaurus"] = {
                "version": self.get_file_mod_date(filepath),
                "source_url": fl,
                "file_mod_date": self.get_file_mod_date(filepath)
            }
            return True
        else:
            logging.error("Downloaded virosaurus90_vertebrate_20200330.fas is corrupted.")
            return False

    def nuc_metadata(self, outfile="acc2taxid.tsv"):
        """
        merge accession and taxonomy info from nuc fasta files.
        """

        if self.update:
            if os.path.isfile(self.metadir + outfile):
                os.remove(self.metadir + outfile)

        if os.path.isfile(self.metadir + outfile):
            acc2tax = pd.read_csv(self.metadir + outfile, sep="\t")
            check = []
            for dbs, fl_list in self.fastas["nuc"].items():
                for fl in fl_list:
                    flb = os.path.basename(fl)
                    if flb not in acc2tax.file.values:
                        check.append(flb)

            if len(check) == 0:
                logging.info("acc2taxid.tsv found for all nuc files.")
                return
            else:
                if self.test:
                    logging.info("acc2taxid.tsv not found for {}".format(check))
                    return
                else:
                    logging.info(
                        f"acc2taxid.tsv not found for nuc files: {check}. creating.."
                    )
                    os.system(f"rm {self.metadir + outfile}")
        else:
            if self.test:
                logging.info("acc2taxid.tsv not found.")
            else:
                logging.info("acc2taxid.tsv not found. creating...")
        ###
        ###
        tax2acc = []

        for dbs, fl_list in self.fastas["nuc"].items():
            for fl in fl_list:
                temp_file = self.metadir + dbs + "_temp.tsv"

                ignore_patterns = ""
                if dbs == "virosaurus":
                    ignore_patterns = "GENE"

                grep_sequence_identifiers(fl, temp_file, ignore=ignore_patterns)

                if dbs == "kraken2":
                    dbacc = pd.read_csv(
                        temp_file, sep="|", names=["suffix", "taxid", "acc"]
                    )
                    dbacc["taxid"] = dbacc["taxid"].astype(str)
                    dbacc["file"] = os.path.basename(fl)
                    dbacc["acc_in_file"] = dbacc[["suffix", "taxid", "acc"]].agg(
                        "|".join, axis=1
                    )

                    dbacc = dbacc[["acc", "taxid", "file", "acc_in_file"]]

                    tax2acc.append(dbacc)
                    continue

                sed_out_after_dot(temp_file)

                # entrez_ncbi_taxid(temp_file, self.metadir, "nuc_tax.tsv")

                dbacc = pd.read_csv(temp_file, sep="\t", header=None)

                dbacc = dbacc.rename(columns={0: "acc"})
                dbacc["file"] = os.path.basename(fl)

                if dbs == "virosaurus":

                    def viro_acc(x):
                        acc = x.split(".")[0]
                        return f"{acc}:{acc};"

                    dbacc["acc_in_file"] = dbacc.acc.apply(viro_acc)

                else:
                    dbacc["acc_in_file"] = dbacc.acc

                tax2acc.append(dbacc)

                os.system(f"rm {temp_file}")
                # os.system("rm {}".format("{}nuc_tax.tsv".format(self.metadir)))

        tax2acc = pd.concat(tax2acc)

        tax2acc.to_csv(self.metadir + outfile, sep="\t", index=False)

    def prot_metadata(self):
        """
        get or produce accession to taxid files for each fasta in fasta.prot.
        :return: self
        """

        self.prot2taxid_rescue()
        self.parse_refseq_prot()

    def parse_refseq_prot(self):
        if "refseq_prot" not in self.fastas["prot"]:
            logging.info("refseq_prot not found.")
            return

        refseq_prot = self.fastas["prot"]["refseq_prot"]
        outfile = f"{self.metadir}refseq_prot_acc2taxid.tsv"

        if os.path.exists(outfile):
            logging.info(f"{outfile} found.")
            return

        def retrieve_within_square_brackets(string):
            return string[string.find("[") + 1 : string.find("]")]

        def retrieve_acc_string(string):
            return string.split()[0][1:]

        lines = []
        with gzip.open(refseq_prot, "rt") as f:
            for line in f:
                if line.startswith(">"):
                    acc = retrieve_acc_string(line)
                    description = retrieve_within_square_brackets(line)
                    lines.append([acc, description])

        df = pd.DataFrame(lines, columns=["acc", "description"])
        tax2description = pd.read_csv(f"{self.metadir}/taxid2desc.tsv", sep="\t")

        df = df.merge(tax2description, on="description", how="left")
        df = df[["acc", "taxid"]]
        df = df.dropna()
        df.taxid = df.taxid.astype(int)

        df.to_csv(outfile, sep="\t", index=False)

        self.meta["refseq_prot"] = outfile

    def prot2taxid_rescue(self):
        """
        parse accession to taxid files for each fasta in fasta.prot.
        """

        dict_ids = self.temp_nucmeta()

        if dict_ids:
            acc2tax_dir = self.get_prot()
            id_files = {i: [] for i in dict_ids}

            threads = [
                Thread(target=self.prot2taxid_parse, args=(dci, dict_ids, acc2tax_dir))
                for dci in range(1, 11)
            ]
            for th in threads:
                th.start()
            for th in threads:
                th.join()

            for dbi in list(id_files):
                id_files[dbi] = [
                    pd.read_csv(self.metadir + f"{dbi}_a2p_{dci}.tsv", sep="\t")
                    for dci in range(1, 11)
                ]
                fdb = pd.concat(id_files[dbi], axis=0)
                #
                fdb.to_csv(self.metadir + f"{dbi}_acc2taxid.tsv", sep="\t", index=False)
                report = pd.merge(
                    fdb, dict_ids[dbi], on="acc", how="outer", indicator=True
                )
                report = report._merge.value_counts()
                report.to_csv(
                    self.metadir + f"{dbi}_acc2taxid.merge.report",
                    sep="\t",
                    index=False,
                )
                #
                for dci in range(1, 11):
                    os.system("rm {}".format(self.metadir + f"{dbi}_a2p_{dci}.tsv"))

                self.meta[dbi] = "refseq_prot"

        logging.info(
            f"accession to taxid mapping done. You can now delete the directory {self.metadir}prot.accession2taxid/"
        )

    def prot2taxid_parse(
        self, dci: int, meta_dict: dict, acc2tax_dir: str, chunksize: int0 = 8e6
    ):
        """
        parse prot2taxid files. given dictionary of accession names, merge these with ncbi two column files.

        :param dci: index number of ncbi file to parse. 10 files in total, named 1-10
        :param meta_dict: dictionary of accession names per fasta (accession ids stored in single column pandas dfs, colname=acc)
        :param acc2tax_dir: direcctory where ncbi files are stored.
        :param id_files: dictionary of empty lists for things to be appended to. same keys as meta_dict.
        :param chunksize: chunck siwe to use in pd.read_csv. very large files.
        :return:
        """
        doc = acc2tax_dir + f"prot.accession2taxid.FULL.{dci}.gz"
        mchunks = {i: [] for i in meta_dict}
        processed = 0
        reader = pd.read_csv(
            doc, compression="gzip", sep="\t", chunksize=int(chunksize)
        )  # , iterator=True)
        for ix, docf in enumerate(reader):
            docf.columns = ["acc", "taxid"]
            for dbi, ids in meta_dict.items():
                rnv = pd.merge(left=ids, right=docf, left_on="acc", right_on="acc")
                mchunks[dbi].append(rnv)

            processed += docf.shape[0]
            #
            print(f"dci: {dci}, {processed} lines processed")

        for dbi in mchunks.keys():
            chk = pd.concat(mchunks[dbi])
            chk.to_csv(self.metadir + f"{dbi}_a2p_{dci}.tsv", sep="\t", index=False)

    def generate_main_protacc_to_taxid(self):
        """
        generates concatenated file of all protein accession to species taxid tsvs.
        """
        final_db_path = os.path.join(self.metadir, "protein_acc2taxid.tsv")
        to_concat = []
        if os.path.isfile(final_db_path):
            logging.info("protein_acc2taxid.tsv file found.")
            return

        for dbs, fl in self.fastas["prot"].items():
            outfile = self.metadir + f"{dbs}_acc2taxid.tsv"
            if os.path.isfile(outfile):
                p2t = pd.read_csv(outfile, sep="\t")
                to_concat.append(p2t)

        if to_concat:
            general_db = pd.concat(to_concat, axis=0)
            general_db.columns = ["prot_acc", "taxid"]
            general_db.drop_duplicates(subset="prot_acc")

        else:
            general_db = pd.DataFrame(columns=["prot_acc", "taxid"])

        general_db.to_csv(final_db_path, header=True, index=False, sep="\t")

    def temp_nucmeta(self):
        """
        read acc ids from fastas in self.fasta.prot.
        :return:
        """
        dict_ids = {}

        for dbs, fl in self.fastas["prot"].items():
            if dbs in ["refseq_prot"]:
                continue

            outfile = self.metadir + f"{dbs}_acc2taxid.tsv"
            if os.path.isfile(outfile):
                self.meta[dbs] = outfile
                logging.info(f"acc2taxid map file {outfile} exists, continuing.")
                continue
            else:
                if self.test:
                    logging.info(f"acc2taxid map file {outfile} not found.")
                    continue
                else:
                    logging.info(f"acc2taxid map file {outfile} not found. creating")

            kept = []
            with gzip.open(fl, "rb") as fn:
                ln = str(fn.readline(), "utf-8")
                while ln:
                    if ln[0] == ">":
                        tp = ln.split()[0][1:]

                        if dbs == "rvdb":
                            tp = tp.split("|")[2]
                        kept.append(tp)
                    else:
                        ln = str(fn.readline(), "utf-8")
                        continue
                    ln = str(fn.readline(), "utf-8")

            dict_ids[dbs] = pd.DataFrame(kept, columns=["acc"])

        return dict_ids

    def get_prot(self):
        """
        download ncbi protein acc2taxid files.
        :return:
        """
        acc2tax_dir = self.metadir + "prot.accession2taxid/"

        if not os.path.isdir(acc2tax_dir):
            os.makedirs(acc2tax_dir, exist_ok=True)

        for si in range(1, 11):
            file = f"https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.FULL.{si}.gz"
            filename = os.path.basename(file)
            destination_file = os.path.join(
                self.metadir, "prot.accession2taxid", filename
            )
            fexist = os.path.exists(destination_file)
            print(f"file {filename} exists: {fexist}")
            tries = 0
            while not fexist:
                try:
                    subprocess.run(
                        [
                            "wget",
                            "-P",
                            self.metadir + "prot.accession2taxid/",
                            f"https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.FULL.{si}.gz",
                        ]
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
                    fexist = True

        return acc2tax_dir


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

            self.dbs[id] = {
                "dir": odir,
                "dbname": dbname,
                "fasta": f"{sdir}/complete.fna.gz",
                "db": index_file_prefix,
                "status": "success",
            }
            return True

        else:
            if self.test:
                logging.info(f"Centrifuge db {dbname} is not installed.")
                self.dbs[id] = {
                    "dir": odir,
                    "dbname": dbname,
                    "fasta": "",
                    "db": index_file_prefix,
                    "status": "test",
                }
                return False
            else:
                logging.info(f"Centrifuge db {dbname} is not installed. Installing...")

        try:

            os.makedirs(sdir, exist_ok=True)
            os.system(f"wget -P {sdir} {href}")
            os.system(f"tar -xvzf {sdir}/p_compressed_2018_4_15.tar.gz -C {sdir}")
            os.system(f"rm {sdir}/p_compressed_2018_4_15.tar.gz")

            files_in_directory = os.listdir(sdir)
            index_files = [f for f in files_in_directory if f.endswith(".cf")]
            for f in index_files:
                os.rename(os.path.join(sdir, f), f"{index_file_prefix}.{f.split('.')[0][-1]}.cf")


            if os.path.isfile(index_file_prefix + ".1.cf"):
                logging.info(f"Centrifuge db {dbname} index is installed.")
                self.dbs[id] = {
                    "dir": odir,
                    "dbname": dbname,
                    "fasta": "",
                    "db": index_file_prefix,
                    "status": "success",
                }
                return True
            else:
                logging.info(f"Centrifuge db {dbname} index is not installed.")
                self.dbs[id] = {
                    "dir": odir,
                    "dbname": dbname,
                    "fasta": "",
                    "db": index_file_prefix,
                    "status": "failed",
                }
                return False
        except subprocess.CalledProcessError:
            logging.error(f"Error occurred while installing centrifuge db {dbname}.")
            self.dbs[id] = {
                "dir": odir,
                "dbname": dbname,
                "fasta": "",
                "db": index_file_prefix,
                "status": "failed",
            }
            return False

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
        bin = self.envs["ROOT"] + self.envs[id] + "/bin/"
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

            self.dbs[id] = {
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

            self.dbs[id] = {
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

            self.dbs[id] = {
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
            self.dbs[id] = {
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

        self.dbs[id] = {
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
            self.dbs[id] = {
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

        self.dbs[id] = {
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
            self.dbs[id] = {
                "dir": odir,
                "dbname": dbname,
                "db": f"{odir}{dbname}/{dbname}.pkl",
                "status": "success",
            }
            return True

        else:
            if self.test:
                logging.info(f"Metaphlan db {dbname} is not installed.")
                self.dbs[id] = {
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

            self.dbs[id] = {
                "dir": odir,
                "dbname": dbname,
                "db": f"{odir}{dbname}/{dbname}",
                "status": "success",
            }
            return True

        except subprocess.CalledProcessError:
            logging.warning(f"failed to download metaphlan db {dbname}")
            self.dbs[id] = {
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
            self.dbs[id] = {
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
            self.dbs[id] = {
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

            self.dbs[id] = {
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

            self.dbs[id] = {
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

            self.dbs[id] = {
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
                self.dbs[id] = {
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
            self.dbs[id] = {
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
            self.dbs[id] = {
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

            self.dbs[id] = {
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
                self.dbs[id] = {
                    "dir": odir,
                    "dbname": dbname,
                    "fasta": odir + dbname + "/library/" + dbname + "/library.fna.gz",
                    "db": odir + dbname,
                }
                return True
            else:
                logging.warning(f"failed to download Kraken2 db {dbname}")
                self.dbs[id] = {
                    "dir": odir,
                    "dbname": dbname,
                    "fasta": odir + dbname + "/library/" + dbname + "/library.fna.gz",
                    "db": odir + dbname,
                }
                return False

        except subprocess.CalledProcessError:
            logging.warning(f"failed to download Kraken2 db {dbname}")
            self.dbs[id] = {
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
            self.dbs[id] = {
                "dir": odir,
                "dbname": dbname,
                "db": odir + dbname,
                "status": "success",
            }

            return True
        else:
            if self.test:
                logging.info(f"Diamond db {dbname} not installed.")
                self.dbs[id] = {
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

            self.dbs[id] = {
                "dir": odir,
                "dbname": dbname,
                "db": odir + dbname,
                "status": "success",
            }
            return True

        except subprocess.CalledProcessError:
            logging.warning(f"failed to download Diamond db {dbname}")
            self.dbs[id] = {
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
            self.dbs[id] = {
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

            self.dbs[id] = {"dir": odir, "dbname": dbname, "db": odir + dbname}

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

        if os.path.isfile(subdb + "kaiju_db_viruses.fmi"):
            logging.info(f"Kaiju {dbname}  is installed.")
            self.dbs[id] = {
                "dir": odir,
                "dbname": dbname,
                "db": subdb + "kaiju_db_viruses.fmi",
                "status": "success",
            }
            return True
        else:
            if self.test:
                logging.info(f"Kaiju {dbname} db is not installed.")
                self.dbs[id] = {
                    "dir": subdb,
                    "dbname": dbname,
                    "db": subdb + "kaiju_db_viruses.fmi",
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

            self.dbs[id] = {
                "dir": subdb,
                "dbname": dbname,
                "db": subdb + "kaiju_db_viruses.fmi",
                "status": "success",
            }
            return True
        except subprocess.CalledProcessError:
            logging.warning(f"failed to download Kaiju db {dbname}")
            self.dbs[id] = {
                "dir": subdb,
                "dbname": dbname,
                "db": subdb + "kaiju_db_viruses.fmi",
                "status": "failed",
            }
            return False

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
        if os.path.isfile(subdb + "kaiju_db_viruses.fmi"):
            logging.info(f"Kaiju {dbname}  is installed.")
            self.dbs[id] = {
                "dir": odir,
                "dbname": dbname,
                "db": subdb + "kaiju_db_viruses.fmi",
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

            self.dbs[id] = {
                "dir": subdb,
                "dbname": dbname,
                "db": subdb + "kaiju_db_viruses.fmi",
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
            self.dbs[id] = {"dir": odir, "dbname": dbname, "db": db, "status": "success"}
            return True
        else:
            if self.test:
                logging.info(f"blast index for {dbname} is not installed.")
                self.dbs[id] = {"dir": odir, "dbname": dbname, "db": db, "status": "test"}
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

            self.dbs[id] = {"dir": sdir, "dbname": dbname, "db": db, "status": "success"}

            return True

        except subprocess.CalledProcessError:
            logging.warning(f"failed to download blast db {dbname}")
            self.dbs[id] = {"dir": sdir, "dbname": dbname, "db": db, "status": "failed"}
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
            self.dbs[id] = {"dir": odir, "dbname": dbname, "db": db, "status": "success"}
            return True
        else:
            if self.test:
                logging.info(f"bowtie2 index for {dbname} is not installed.")
                self.dbs[id] = {"dir": odir, "dbname": dbname, "db": db, "status": "test"}

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

            self.dbs[id] = {"dir": odir, "dbname": dbname, "db": db, "status": "success"}

            return True

        except subprocess.CalledProcessError:
            logging.warning(f"failed to download bowtie2 index {dbname}")
            self.dbs[id] = {"dir": odir, "dbname": dbname, "db": db, "status": "failed"}
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
            self.dbs[id] = {
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
                self.dbs[id] = {
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
            self.dbs[id] = {
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
                self.dbs[id] = {
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
                self.dbs[id] = {
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
            self.dbs[id] = {
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
            self.dbs[id] = {
                "dir": odir,
                "dbname": dbname,
                "fasta": f"{sdir}.fa",
                "db": f"{odir}{dbname}/{dbname}",
                "status": "failed",
            }
            return False
        except FileNotFoundError as e:
            logging.error(f"BWA install failed: required tool not found: {e}")
            self.dbs[id] = {
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
            self.dbs[id] = {
                "dir": odir,
                "dbname": dbname,
                "db": reference,
                "fasta": reference,
                "status": "success",
            }
            return True
        elif self.test:
            self.dbs[id] = {
                "dir": odir,
                "dbname": dbname,
                "db": reference,
                "fasta": reference,
                "status": "test",
            }
            return False
        else:
            self.dbs[id] = {
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

            self.dbs[id] = {"dir": odir, "dbname": dbname}
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
            self.dbs[id] = {"dir": odir, "dbname": dbname, "db": fidx}
            return True
        else:
            if self.test:
                logging.info(f"FastViromeExplorer {reference} index is not installed.")
                self.dbs[id] = {"dir": odir, "dbname": dbname, "db": fidx}
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

            self.dbs[id] = {"dir": subdir, "dbname": dbname, "db": fidx}
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
            self.dbs[id] = {"dir": sdir, "dbname": dbname, "db": sdir}
            return True
        else:
            if self.test:
                logging.info(f"deSAMBA db {dbname} is not installed.")
                self.dbs[id] = {"dir": odir, "dbname": dbname, "db": sdir}
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

            self.dbs[id] = {"dir": odir, "dbname": dbname, "db": sdir}
            return True

        except subprocess.CalledProcessError:
            logging.info(f"failed to install deSAMBA")
            return False
