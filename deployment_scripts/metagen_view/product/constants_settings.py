"""
Ceated on 06/05/2022
@author: joao santos
"""

import configparser
import os

from metagen_view.settings import BASE_DIR, STATICFILES_DIRS


class ConstantsSettings:

    product_directory = "/mnt/sdc/TELEVIR/products/"
    static_directory = STATICFILES_DIRS[0]
    job_directory = "/mnt/sdc/TELEVIR/jobs/"

    django_env: str

    def __init__(self):

        if not os.path.exists(self.product_directory):
            os.makedirs(self.product_directory, exist_ok=True)

        if not os.path.exists(self.static_directory):
            os.makedirs(self.static_directory, exist_ok=True)

        if not os.path.exists(self.job_directory):
            os.makedirs(self.job_directory, exist_ok=True)

        self.env_file = os.path.join(BASE_DIR, ".env")

        # with open(self.env_file, "r") as f:
        #    self.django_env = f.read().strip()
