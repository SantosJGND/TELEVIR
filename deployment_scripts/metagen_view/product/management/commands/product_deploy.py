import os
from datetime import date

from django.core.management.base import BaseCommand
from product.utils import Run_Main_from_Fastq_Input


class Command(BaseCommand):
    help = "deploy run"

    def add_arguments(self, parser):
        parser.add_argument(
            "--pk",
            type=int,
            help="fastq input pk to run",
        )

    def handle(self, *args, **options):
        ###
        #

        run = Run_Main_from_Fastq_Input(options["pk"])
        run.Submit()
