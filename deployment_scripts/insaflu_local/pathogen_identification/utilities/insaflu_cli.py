#####
### generate tree
#####
import logging
import os
import shutil
from datetime import date

from constants.meta_key_and_values import MetaKeyAndValue
from django.conf import settings
from django.contrib.auth.models import User
from extend_user.models import Profile
from pathogen_identification.models import PIProject_Sample, Projects, SoftwareTree
from settings.constants_settings import ConstantsSettings
from settings.models import Sample
from utils.process_SGE import ProcessSGE
from utils.utils import Utils


class Insaflu_Cli:
    def __init__(self):
        self.logger = logging.getLogger(__name__)

    @staticmethod
    def create_televir_project_if_not_exists(project_name, user, technology):
        try:
            project = Projects.objects.get(name=project_name)
        except Projects.DoesNotExist:
            project = Projects()
            project.name = project_name
            project.owner = user
            project.technology = technology
            project.save()
        return project

    @staticmethod
    def read_fofn(fofn):

        with open(fofn) as f:
            fastq_paths = f.read().splitlines()

        r1 = fastq_paths[0]
        if len(fastq_paths) > 1:
            r2 = fastq_paths[1]
        else:
            r2 = None

        return r1, r2

    def temp_reads_fofn(self, fofn):
        utils = Utils()
        r1, r2 = self.read_fofn(fofn)
        temp_dir = utils.get_temp_dir()

        new_r1 = os.path.join(temp_dir, os.path.basename(r1))
        new_r2 = None

        shutil.copy(r1, new_r1)

        if r2:
            new_r2 = os.path.join(temp_dir, os.path.basename(r2))
            shutil.copy(r2, new_r2)

        return new_r1, new_r2

    def temp_reads(self, r1, r2):
        utils = Utils()
        temp_dir = utils.get_temp_dir()

        new_r1 = os.path.join(temp_dir, os.path.basename(r1))
        new_r2 = None

        shutil.copy(r1, new_r1)

        if r2:
            new_r2 = os.path.join(temp_dir, os.path.basename(r2))
            shutil.copy(r2, new_r2)

        return new_r1, new_r2

    def create_sample_from_fofn(self, fofn, user: User, technology: str):

        r1, r2 = self.read_fofn(fofn)
        name = os.path.basename(fofn).split(".")[0]

        try:
            sample = Sample.objects.get(name=name, owner=user)

        except Sample.DoesNotExist:
            tmp_r1, tmp_r2 = self.temp_reads_fofn(fofn)
            sample = self.sample_save(name, user, tmp_r1, tmp_r2, technology)
            self.move_sample(sample, user)
            self.sample_preprocess(sample, user)

        return sample

    def sample_save(self, name, user, r1, r2, technology):
        utils = Utils()
        if not os.path.exists(r1):
            raise FileNotFoundError(f"File {r1} does not exist")

        sample = Sample.objects.create(
            name=name,
            owner=user,
            path_name_1=r1,
            path_name_2=r2,
            type_of_fastq=int(technology == ConstantsSettings.TECHNOLOGY_minion),
            date_of_onset=date.today(),
            date_of_collection=date.today(),
        )

        sample.is_deleted = False
        sample.is_obsolete = False
        sample.file_name_1 = utils.clean_name(os.path.basename(sample.path_name_1.name))
        sample.is_valid_1 = True
        if sample.exist_file_2():
            sample.file_name_2 = utils.clean_name(
                os.path.basename(sample.path_name_2.name)
            )
            sample.is_valid_2 = True
        else:
            sample.is_valid_2 = False
        sample.has_files = True

        sample.manual_upload()

        sample.save()

        return sample

    def move_sample(self, sample: Sample, user: User):
        utils = Utils()

        sz_file_to = os.path.join(
            getattr(settings, "MEDIA_ROOT", None),
            utils.get_path_to_fastq_file(user.id, sample.id),
            sample.file_name_1,
        )
        utils.move_file(
            os.path.join(
                getattr(settings, "MEDIA_ROOT", None), sample.path_name_1.name
            ),
            sz_file_to,
        )
        sample.path_name_1.name = os.path.join(
            utils.get_path_to_fastq_file(user.id, sample.id),
            sample.file_name_1,
        )

        if sample.exist_file_2():
            sz_file_to = os.path.join(
                getattr(settings, "MEDIA_ROOT", None),
                utils.get_path_to_fastq_file(user.id, sample.id),
                sample.file_name_2,
            )
            utils.move_file(
                os.path.join(
                    getattr(settings, "MEDIA_ROOT", None), sample.path_name_2.name
                ),
                sz_file_to,
            )
            sample.path_name_2.name = os.path.join(
                utils.get_path_to_fastq_file(user.id, sample.id),
                sample.file_name_2,
            )
        sample.save()

    def sample_preprocess(self, sample: Sample, user: User):
        from managing_files.manage_database import ManageDatabase

        process_SGE = ProcessSGE()
        job_name = None

        try:
            (job_name_wait, job_name) = user.profile.get_name_sge_seq(
                Profile.SGE_PROCESS_clean_sample, Profile.SGE_SAMPLE
            )
            if sample.is_type_fastq_gz_sequencing():  ### default is Illumina
                taskID = process_SGE.set_run_trimmomatic_species(sample, user, job_name)
            else:  ### Minion, codify with other
                taskID = process_SGE.set_run_clean_minion(sample, user, job_name)

            sample.is_ready_for_projects = True
            sample.save()

        except Exception as e:
            self.logger.error("Fail to run: ProcessSGE - " + str(e))
            return

        ## refresh sample list for this user
        if not job_name is None:
            process_SGE.set_create_sample_list_by_user(user, [job_name])
        ###

        manageDatabase = ManageDatabase()

        manageDatabase.set_sample_metakey(
            sample,
            user,
            MetaKeyAndValue.META_KEY_Queue_TaskID,
            MetaKeyAndValue.META_VALUE_Queue,
            taskID,
        )

    def piproject_sample_from_sample(
        self, sample: Sample, project: Projects, user: User
    ):

        try:
            project_sample = PIProject_Sample.objects.get(
                project=project, sample=sample
            )
        except PIProject_Sample.DoesNotExist:
            project_sample_input = sample.file_name_1
            if sample.is_valid_2:
                project_sample_input += ";" + sample.file_name_2

            project_sample = PIProject_Sample.objects.create(
                project=project,
                sample=sample,
                name=sample.name,
                input=project_sample_input,
                technology=sample.type_of_fastq,
                report="report",
            )
            project_sample.save()

        return project_sample

    def create_benchmark_softwaretree_if_not_exists(self):

        try:
            benchmark_softwaretree = SoftwareTree.objects.filter(name="default")
        except SoftwareTree.DoesNotExist:
            benchmark_softwaretree = SoftwareTree.objects.create(
                name="default",
                owner=User.objects.get(username="admin"),
            )
            benchmark_softwaretree.save()

        return benchmark_softwaretree
