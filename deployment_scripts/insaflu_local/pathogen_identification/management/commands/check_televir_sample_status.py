#####
### generate tree
#####

import sys

from django.contrib.auth.models import User
from django.core.management.base import BaseCommand
from pathogen_identification.models import ParameterSet, PIProject_Sample, Projects
from settings.models import Sample


class Command(BaseCommand):
    help = "Populates the DBS"

    def add_arguments(self, parser):
        parser.add_argument(
            "--sample",
            "-s",
            type=str,
            help="file of fastq files. one file path per line.",
        )

        parser.add_argument(
            "-p", "--project", type=str, required=True, help="Output directory."
        )

        parser.add_argument(
            "-u", "--user", type=str, required=True, help="User to assign to project."
        )

    def handle(self, *args, **options):
        ###
        #

        project_name = options["project"]
        user_name = options["user"]
        sample_name = options["sample"]

        # get user
        try:
            user = User.objects.get(username=user_name)
        except User.DoesNotExist:
            print("User does not exist.")
            sys.exit(1)
        # get project

        # check sample exists:

        try:
            sample = Sample.objects.get(name=sample_name)
        except Sample.DoesNotExist:
            print("Sample does not exist.")
            print("STATUS 0")
            sys.exit(1)

        try:
            project = Projects.objects.get(name=project_name, owner=user)
        except Projects.DoesNotExist:
            print("Project does not exist.")
            print("STATUS 0")
            sys.exit(1)

        ########################################################
        ## Possible sample status:
        # 0. not started
        # 1. running
        # 2. finished
        # 3. failed

        # check televir project sample exists:
        try:
            televir_project_sample = PIProject_Sample.objects.get(
                project=project, sample=sample
            )
        except PIProject_Sample.DoesNotExist:
            print("STATUS 0")
            sys.exit(1)

        # check if sample is running using parameter sets.

        parameter_sets = ParameterSet.objects.filter(
            project=project, sample=televir_project_sample
        )

        if parameter_sets.count() == 0:
            print("STATUS 0")
            sys.exit(1)

        for parameter_set in parameter_sets:
            if parameter_set.status in [
                ParameterSet.STATUS_RUNNING,
                ParameterSet.STATUS_QUEUED,
            ]:
                print("STATUS 1")
                sys.exit(1)

        # check if sample is finished using parameter sets.

        for parameter_set in parameter_sets:
            if parameter_set.status == ParameterSet.STATUS_FINISHED:
                print("STATUS 2")
                sys.exit(1)

        # check if sample is failed using parameter sets.

        for parameter_set in parameter_sets:
            if parameter_set.status == ParameterSet.STATUS_ERROR:
                print("STATUS 3")
                sys.exit(1)
