import sys

from django.core.management.base import BaseCommand, CommandError
from result_display.local_parse import Parser, Project_create


class Command(BaseCommand):
    help = "Populates the DBS"

    def handle(self, *args, **options):
        Project_create()

        self.stdout.write(self.style.SUCCESS("DBS populated"))
        sys.exit(0)
