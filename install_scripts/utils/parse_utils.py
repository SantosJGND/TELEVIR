#!/usr/bin/python3

import os
from typing import Iterable, List

import dnaio
import xopen
from fastq_filter import file_to_fastq_records


def fastq_records_to_files(
    records: Iterable[dnaio.Sequence],
    filepath_template: str = "test",
    compression_level: int = 2,
    max_file_size: int = 1000000000,
    margin: int = 1000000,
) -> List[str]:
    """
    split fastq records into files of max_file_size"""

    files = []
    file_index = 0
    filepath = filepath_template + f"_{file_index}" + ".fastq.gz"

    output_h = xopen.xopen(
        filepath, mode="wb", threads=0, compresslevel=compression_level
    )

    print(f"writing to {filepath}")

    for record in records:
        if os.path.getsize(filepath) > (max_file_size - margin):
            output_h.close()
            files.append(filepath)
            file_index += 1
            filepath = filepath_template + f"_{file_index}" + ".fastq.gz"
            output_h = xopen.xopen(
                filepath, mode="wb", threads=0, compresslevel=compression_level
            )
            print(f"writing to {filepath}")

        header = ">" + record.name + "\n"
        header = header.encode("ascii")
        output_h.write(header)
        sequence = record.sequence + "\n"
        output_h.write(sequence.encode("ascii"))

    output_h.close()
    files.append(filepath)
    return files


def bioinf_splitext(filepath: str):
    basename, ext = os.path.splitext(filepath)
    if ext == ".gz":
        basename, ext = os.path.splitext(basename)

    return basename, ext


def predict_files_split(
    filepath: str, max_file_size: int = 1000000000, file_template: str = "test"
):
    """
    predict files that will be created by split_file"""
    file_size = os.path.getsize(filepath)
    n_files = file_size / max_file_size
    # round up
    n_files = int(n_files) + 1
    files_predict = [file_template + f"_{i}" + ".fastq.gz" for i in range(n_files)]
    return files_predict


def check_proceed(filepath: str, file_template: str, max_file_size: int = 1000000000):
    """
    check if needed to continue"""
    files_predict = predict_files_split(
        filepath, max_file_size=max_file_size, file_template=file_template
    )

    files_exist = [
        (os.path.exists(file) and os.path.getsize(file) > 100) for file in files_predict
    ]
    if all(files_exist):
        return False
    else:
        return True


def split_file(filepath: str, max_file_size=1000000000) -> List[str]:
    """
    split file into smaller files"""

    if os.path.getsize(filepath) < max_file_size:
        return [filepath]

    else:
        records = file_to_fastq_records(filepath)
        basename, _ = bioinf_splitext(filepath)
        if check_proceed(filepath, basename, max_file_size=max_file_size):
            print(f"splitting {filepath}")
            return fastq_records_to_files(
                records, filepath_template=basename, max_file_size=max_file_size
            )

        else:
            return predict_files_split(
                filepath, max_file_size=max_file_size, file_template=basename
            )


def process_file_list(
    file_list: List[str], max_file_size: int = 1000000000
) -> List[str]:
    """
    split files in file_list into smaller files"""

    output_files = []
    for filepath in file_list:
        output_files += split_file(filepath, max_file_size=max_file_size)

    return output_files


def process_nuc_fasta_dict(
    nuc_fasta_dict: dict, max_file_size: int = 1000000000
) -> dict:
    """
    split files in file_list into smaller files"""

    output_files = {}
    for software, file_list in nuc_fasta_dict.items():
        output_files[software] = process_file_list(
            file_list, max_file_size=max_file_size
        )

    return output_files
