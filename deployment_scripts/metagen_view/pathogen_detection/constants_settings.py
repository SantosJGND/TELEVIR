"""
Ceated on 06/05/2022
@author: joao santos
"""

import json
import os

from metagen_view.settings import STATIC_URL, STATICFILES_DIRS


class ConstantsSettings:

    project_directory = "/mnt/sdc/TELEVIR/projects/"
    static_directory = STATIC_URL
    exclude_descriptions = ["phage", "Phage"]


class TestConstants:
    Test_Static_Directory = os.path.join(STATICFILES_DIRS[0], "tests")
    Test_Temp_Directory = "/mnt/sdc/TELEVIR/temp/"

    ont_params_python = os.path.join(Test_Static_Directory, "ont_params.py")
    ont_params_json = os.path.join(Test_Static_Directory, "ont_params.json")

    ont_fastq_gz_file = "ont_fastq.gz"
    ont_fastq_gz_file_path = os.path.join(Test_Static_Directory, ont_fastq_gz_file)

    illumina_params_python = os.path.join(Test_Static_Directory, "illumina_params.py")
    illumina_params_json = os.path.join(Test_Static_Directory, "illumina_params.json")

    illumina_fastq_gz_file_r1 = "illumina_r1.fq.gz"
    illumina_fastq_gz_file_r1_path = os.path.join(
        Test_Static_Directory, illumina_fastq_gz_file_r1
    )
    illumina_fastq_gz_file_r2 = "illumina_r2.fq.gz"
    illumina_fastq_gz_file_r2_path = os.path.join(
        Test_Static_Directory, illumina_fastq_gz_file_r2
    )
