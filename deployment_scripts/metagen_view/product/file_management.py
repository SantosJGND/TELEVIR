import os
import shutil
import subprocess

from metagen_view import settings
from metagen_view.settings import MEDIA_ROOT

from product.constants_settings import ConstantsSettings
from product.models import Fastq_Input


class Ephemeral_Project_Manager:

    constants = ConstantsSettings()

    root = os.path.join(MEDIA_ROOT)

    def handle_uploaded_file(self, f):
        print(f.name)
        os.makedirs(self.root, exist_ok=True)
        with open(os.path.join(self.root, f.name), "wb+") as destination:
            for chunk in f.chunks():
                destination.write(chunk)

        return os.path.join(self.root, f.name)

    def submit_job(self, fastq_input: Fastq_Input):

        job_file_name = f"{fastq_input.project_name}_{fastq_input.pk}.sh"
        job_file_path = os.path.join(self.constants.job_directory, job_file_name)
        base_dir = settings.BASE_DIR
        python_bin = self.constants.django_env + "/bin/python"
        log_file = os.path.join(
            self.constants.job_directory,
            f"{fastq_input.project_name}_{fastq_input.pk}.log",
        )

        job_file = open(os.path.join(self.constants.job_directory, job_file_name), "w+")
        job_file.write(f"#!/bin/bash\n")

        job_file.write(f"\n")
        job_file.write(f"# submitted : {fastq_input.date_created}\n")
        job_file.write(f"# project : {fastq_input.project_name}\n")
        job_file.write(f"# fastq_input : {fastq_input.pk}\n")
        job_file.write(f"# log_file : {log_file}\n")
        job_file.write(f"# base_dir : {base_dir}\n")
        job_file.write(f"# fastq_input_path : {fastq_input.file_r1}\n")
        job_file.write(
            f"# fastq_input_path : {fastq_input.file_r2}\n"
        ) if fastq_input.file_r2 else None
        job_file.write(f"# technology : {fastq_input.technology}\n")
        job_file.write(f"\n")

        job_file.write(
            f"{python_bin} {base_dir}/manage.py product_deploy --pk {fastq_input.pk}"
        )

        job_file.write(f"\n")

        job_file.close()
        print(job_file_path)

        subprocess.run(
            [
                "nohup",
                "sh",
                job_file_path,
                "&>",
                log_file,
                "&",
            ]
        )
        # (f"nohup sh {job_file_path} &> {log_file} &")
