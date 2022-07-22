import os
import sys

from django.core.management.base import BaseCommand, CommandError
from metagen_view.settings import STATICFILES_DIRS


class Command(BaseCommand):
    help = "Clears the static files"

    def handle(self, *args, **options):
        static_dirs = [
            "assemblies/",
            "classification_reports/",
            "depleted_reads/",
            "igv_files/",
            "refa_dotplots/",
            "ref_coverage/",
        ]

        for static_dir in static_dirs:
            static_dir = os.path.join(STATICFILES_DIRS[0], static_dir)
            if not os.path.exists(static_dir):
                os.makedirs(static_dir)

            for f in os.listdir(static_dir):
                os.remove(os.path.join(static_dir, f))
