import os
import sys

from django.core.management.base import BaseCommand, CommandError
from metagen_view.settings import STATICFILES_DIRS


class Command(BaseCommand):
    help = "Clears the static files"

    def handle(self, *args, **options):
        static_dirs = [
            "classification_reports/",
            "igv/",
            "plots/",
        ]

        for static_dir in static_dirs:
            static_dir = os.path.join(STATICFILES_DIRS[0], "product", static_dir)
            if not os.path.exists(static_dir):
                os.makedirs(static_dir)

            for f in os.listdir(static_dir):
                os.remove(os.path.join(static_dir, f))
