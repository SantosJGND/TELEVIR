import os
import sys
from datetime import date
from typing import List

from django.contrib.auth.models import User
from django.core.management.base import BaseCommand
from managing_files.models import ProcessControler
from pathogen_identification.models import PIProject_Sample, Projects
from pathogen_identification.utilities.insaflu_cli import Insaflu_Cli
from settings.constants_settings import ConstantsSettings


class Command(BaseCommand):
    help = "deploy run"

    def add_arguments(self, parser):
        parser.add_argument(
            "--user_name",
            type=str,
            help="user deploying the run (string)",
        )

        parser.add_argument(
            "--metadata",
            "-m",
            type=str,
            default=None,
            help="file of metadata. one file path per line.",
        )

        parser.add_argument(
            "--fofn",
            "-f",
            type=str,
            default=None,
            help="file of fastq files. one file path per line.",
        )

        parser.add_argument(
            "-t",
            "--technology",
            type=str,
            default="ONT",
            help="technology",
        )

    @staticmethod
    def add_samples(user: User, project: Projects, **options):
        insaflu_cli = Insaflu_Cli()
        sample_list = []

        technology = project.technology

        if options["fofn"] is not None:
            sample_fofn = options["fofn"]
            sample = insaflu_cli.create_sample_from_fofn(sample_fofn, user, technology)
            sample_list.append(sample)

        if options["metadata"] is not None:
            metadata_fofn = options["metadata"]
            samples = insaflu_cli.create_sample_from_metadata(metadata_fofn, user)
            sample_list.extend(samples)

        for sample in sample_list:
            project_sample = insaflu_cli.piproject_sample_from_sample(
                sample, project, user
            )

    def process_arguments(self, *args, **options):

        try:
            user = User.objects.get(username=options["user_name"])
        except User.DoesNotExist:
            print("User does not exist")
            sys.exit(1)

        try:
            project = Projects.objects.get(name=options["project_name"], owner=user)
        except Projects.DoesNotExist:
            print(f"Project {options['project_name']} does not exist.")

            if options["technology"] is None:
                print("Please provide a technology for the project")
                sys.exit(1)
            else:
                technology = options["technology"]
                print(
                    f"Creating project {options['project_name']} with technology {technology}"
                )

            insaflu_cli = Insaflu_Cli()

            project = insaflu_cli.create_televir_project_if_not_exists(
                options["project_name"], user, technology
            )

        return user, project

    def handle(self, *args, **options):
        ###
        #### SETUP

        insaflu_cli = Insaflu_Cli()

        technology = options["technology"]

        if technology.lower() in ["nanopore", "ont"]:
            technology = ConstantsSettings.TECHNOLOGY_minion
        elif technology.lower() in ["illumina", "pacbio", "Illumina/IonTorrent"]:
            technology = ConstantsSettings.TECHNOLOGY_illumina
        else:
            print("Technology not supported")
            sys.exit(1)

        try:
            user = User.objects.get(username=options["user_name"])
        except User.DoesNotExist:
            print("User does not exist")
            sys.exit(1)

        if options["fofn"] is not None:
            sample_fofn = options["fofn"]
            sample = insaflu_cli.create_sample_from_fofn(sample_fofn, user, technology)

        if options["metadata"] is not None:
            metadata_fofn = options["metadata"]

            samples = insaflu_cli.create_sample_from_metadata(metadata_fofn, user)
