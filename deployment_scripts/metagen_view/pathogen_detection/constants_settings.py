"""
Ceated on 06/05/2022
@author: joao santos
"""

import json
import os

from metagen_view.settings import STATICFILES_DIRS


class ConstantsSettings:

    project_directory = "/mnt/sdc/TELEVIR/projects/"
    static_directory = STATICFILES_DIRS[0]


class TestConstants:
    Test_Static_Directory = os.path.join(STATICFILES_DIRS[0], "tests")
    Test_Temp_Directory = "/mnt/sdc/TELEVIR/temp/"

    ont_params_python = os.path.join(Test_Static_Directory, "ont_params.py")
    ont_params_json = os.path.join(Test_Static_Directory, "ont_params.json")

    ont_fastq_gz_file = "ont_fastq.gz"
    ont_fastq_gz_file_path = os.path.join(Test_Static_Directory, ont_fastq_gz_file)
