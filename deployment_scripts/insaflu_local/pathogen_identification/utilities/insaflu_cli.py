#####
# generate tree
#####
import logging
import ntpath
import os
import shutil
from datetime import date
from typing import List

from constants.constants import Constants, TypeFile
from constants.meta_key_and_values import MetaKeyAndValue
from django.conf import settings
from django.contrib.auth.models import User
from extend_user.models import Profile
from managing_files.models import (MetaKey, ProcessControler, Sample,
                                   UploadFiles)
from pathogen_identification.models import (PIProject_Sample, Projects,
                                            SoftwareTree)
from settings.constants_settings import ConstantsSettings
from settings.models import Sample
from utils.parse_in_files import ParseInFiles
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

    def temp_reads(self, fofn):
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

    def create_sample_from_fofn(self, fofn, user: User, technology: str):

        r1, r2 = self.read_fofn(fofn)
        name = os.path.basename(fofn)
        name = os.path.splitext(fofn)[0]

        try:
            sample = Sample.objects.get(name=name, owner=user)

        except Sample.DoesNotExist:
            tmp_r1, tmp_r2 = self.temp_reads(fofn)
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
            type_of_fastq=int(
                technology == ConstantsSettings.TECHNOLOGY_minion),
            date_of_onset=date.today(),
            date_of_collection=date.today(),
        )

        sample.is_deleted = False
        sample.is_obsolete = False
        sample.file_name_1 = utils.clean_name(
            os.path.basename(sample.path_name_1.name))
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
                    getattr(settings, "MEDIA_ROOT",
                            None), sample.path_name_2.name
                ),
                sz_file_to,
            )
            sample.path_name_2.name = os.path.join(
                utils.get_path_to_fastq_file(user.id, sample.id),
                sample.file_name_2,
            )
        sample.save()

    def metadata_check_errors(self, metadata_full_path: str, user: User):

        if not os.path.exists(metadata_full_path):
            self.logger.info(
                "Metadata file {} could not be found".format(
                    metadata_full_path)
            )
            return False

        # Process the metadata file to check if everything is ok
        parse_in_files = self.metadata_parse(metadata_full_path, user)

        if parse_in_files.get_errors().has_errors():

            self.logger.info(
                "Errors found processing the metadata file {}".format(
                    metadata_full_path
                )
            )
            self.logger.info(str(parse_in_files.get_errors()))
            # self.logger_debug.error("Errors found processing the metadata table")
            # self.logger_debug.erro(str(parse_in_files.get_errors()))
            return False

        return True

    def metadata_parse(self, metadata_full_path: str, user: User):
        parse_in_files = ParseInFiles()
        b_test_char_encoding = False
        parse_in_files.parse_sample_files(
            metadata_full_path,
            user,
            b_test_char_encoding,
            ParseInFiles.STATE_READ_all,
        )

        return parse_in_files

    def metadata_fastqs(self, sample_list: List[Sample]):
        fastq_files_to_upload = []
        missing_fastqs = False
        for sample in sample_list:
            fastq1 = sample[0].candidate_file_name_1

            fastq_full_path = os.path.join(
                getattr(settings, "MEDIA_ROOT", None),
                Constants.DIR_PROCESSED_FILES_UPLOADS,
                fastq1,
            )
            if not os.path.exists(fastq_full_path):
                self.logger.info(
                    "Fastq file {} could not be found".format(fastq_full_path)
                )
                missing_fastqs = True
            fastq_files_to_upload.append(fastq_full_path)

            fastq2 = ""
            # Carefull this may be sensitive to spaces in the table
            if sample[0].candidate_file_name_2.strip() != "":
                fastq2 = sample[0].candidate_file_name_2
                fastq_full_path = os.path.join(
                    getattr(settings, "MEDIA_ROOT", None),
                    Constants.DIR_PROCESSED_FILES_UPLOADS,
                    fastq2,
                )
                fastq_files_to_upload.append(fastq_full_path)
                if not os.path.exists(fastq_full_path):
                    self.logger.info(
                        "Fastq file {} could not be found".format(
                            fastq_full_path)
                    )
                    missing_fastqs = True

            self.logger.info(
                "sample {} file(s) to be processed: {} {} ".format(
                    sample[0].name, fastq1, fastq2
                )
            )

        if missing_fastqs:
            self.logger.info("Fastq files are missing, cannot continue")
            return []

        self.logger.info(
            " {} samples are going to be processed".format(len(sample_list))
        )

        return fastq_files_to_upload

    def metadata_upload_prep(self, metadata_full_path: str, user: User):
        # Add the metadata file as an upload
        # may not be needed, but for consistency with website we do it
        metadata_file = os.path.basename(metadata_full_path)
        utils = Utils()

        sample_file_upload_files = UploadFiles()
        sample_file_upload_files.file_name = metadata_file
        sample_file_upload_files.is_valid = True
        sample_file_upload_files.is_processed = False
        sample_file_upload_files.is_deleted = False
        sample_file_upload_files.number_errors = 0
        sample_file_upload_files.number_files_processed = 0

        try:
            type_file = MetaKey.objects.get(
                name=TypeFile.TYPE_FILE_sample_file)
        except MetaKey.DoesNotExist:
            type_file = MetaKey()
            type_file.name = TypeFile.TYPE_FILE_sample_file
            type_file.save()

        sample_file_upload_files.type_file = type_file
        sample_file_upload_files.owner = user
        sample_file_upload_files.description = ""

        # move the files to the right place
        sz_file_to = os.path.join(
            getattr(settings, "MEDIA_ROOT", None),
            utils.get_path_upload_file(
                user.id, TypeFile.TYPE_FILE_sample_file),
            metadata_file,
        )

        # get unique file name, as the user can upload files with same name...
        sz_file_to, path_added = utils.get_unique_file(sz_file_to)

        # Add this back in the end... to "consume" the file
        utils.copy_file(metadata_full_path, sz_file_to)

        if path_added is None:
            sample_file_upload_files.path_name.name = os.path.join(
                utils.get_path_upload_file(
                    user.id, TypeFile.TYPE_FILE_sample_file),
                ntpath.basename(sz_file_to),
            )
        else:
            sample_file_upload_files.path_name.name = os.path.join(
                utils.get_path_upload_file(
                    user.id, TypeFile.TYPE_FILE_sample_file),
                path_added,
                ntpath.basename(sz_file_to),
            )

        self.logger.info(
            "{} file was processed".format(
                sample_file_upload_files.path_name.name)
        )

        return sample_file_upload_files

    def sample_file_process(self, fastq_to_upload, user: User) -> UploadFiles:
        utils = Utils()

        self.logger.info("Fastq file to upload: {}".format(fastq_to_upload))

        fastq_upload_files = UploadFiles()
        fastq_upload_files.file_name = utils.clean_name(
            os.path.basename(fastq_to_upload)
        )

        # move the files to the right place
        sz_file_to = os.path.join(
            getattr(settings, "MEDIA_ROOT", None),
            utils.get_path_upload_file(user.id, TypeFile.TYPE_FILE_fastq_gz),
            fastq_upload_files.file_name,
        )
        sz_file_to, path_added = utils.get_unique_file(
            sz_file_to
        )  # get unique file name, user can upload files with same name...
        #     utils.move_file(temp_file, sz_file_to)
        utils.copy_file(fastq_to_upload, sz_file_to)

        # test if file exists (may fail due to full disk or other error)
        if not os.path.exists(sz_file_to) and os.path.getsize(sz_file_to) > 10:
            self.logger.info(
                " Error copying file {} file to {}".format(
                    fastq_to_upload, sz_file_to)
            )
            # If we do a return here then we need an atomic transaction or a way to prevent inconsistencies...
            return False

        if path_added is None:
            fastq_upload_files.path_name.name = os.path.join(
                utils.get_path_upload_file(
                    user.id, TypeFile.TYPE_FILE_fastq_gz),
                fastq_upload_files.file_name,
            )
        else:
            fastq_upload_files.path_name.name = os.path.join(
                utils.get_path_upload_file(
                    user.id, TypeFile.TYPE_FILE_fastq_gz),
                path_added,
                fastq_upload_files.file_name,
            )

        try:
            type_file = MetaKey.objects.get(name=TypeFile.TYPE_FILE_fastq_gz)
        except MetaKey.DoesNotExist:
            type_file = MetaKey()
            type_file.name = TypeFile.TYPE_FILE_fastq_gz
            type_file.save()

        fastq_upload_files.is_valid = True
        fastq_upload_files.is_processed = False  # True when all samples are set
        fastq_upload_files.owner = user
        fastq_upload_files.type_file = type_file
        fastq_upload_files.number_files_to_process = 1
        fastq_upload_files.number_files_processed = 0
        fastq_upload_files.description = ""

        return fastq_upload_files

    def process_sample_list(self, sample_list: List[str], user: User):
        upload_list = []

        for sample in sample_list:
            sample_object = self.sample_file_process(sample, user)
            upload_list.append(sample_object)

        return upload_list

    @staticmethod
    def save_objects(objects_to_save: List[object]):
        for object in objects_to_save:
            object.save()

    @staticmethod
    def retrieve_saved_sample(file_parser: ParseInFiles, user: User):
        samples = set()
        for vect_sample in file_parser.vect_samples:

            try:
                sample = Sample.objects.get(
                    name__iexact=vect_sample[0].name, owner=user, is_deleted=False
                )
                # if exist don't add
                samples.add(sample)
                continue  # already exist
            except Sample.DoesNotExist as e:
                pass

        return samples

    def create_sample_from_metadata(self, metadata_full_path: str, user: User):

        print("metadata_full_path: ", metadata_full_path)
        print("user: ", user)
        print("errors: ", self.metadata_check_errors(metadata_full_path, user))

        if not self.metadata_check_errors(metadata_full_path, user):
            return False

        parse_in_files = self.metadata_parse(metadata_full_path, user)

        fastq_files_to_upload = self.metadata_fastqs(
            parse_in_files.get_vect_samples())
        print("fastq_files_to_upload: ", fastq_files_to_upload)

        if not fastq_files_to_upload:

            return []

        sample_file_upload_files = self.metadata_upload_prep(
            metadata_full_path, user)

        sample_file_upload_files.number_files_to_process = len(
            parse_in_files.get_vect_samples()
        )

        samples_to_save = self.process_sample_list(fastq_files_to_upload, user)
        upload_files = samples_to_save + [sample_file_upload_files]

        self.save_objects(upload_files)

        parse_in_files.create_samples(sample_file_upload_files, user)
        parse_in_files.link_files(user, False)

        self.logger.info("End")

        samples = self.retrieve_saved_sample(parse_in_files, user)

        return samples

    def sample_preprocess(self, sample: Sample, user: User):
        from managing_files.manage_database import ManageDatabase

        process_SGE = ProcessSGE()
        job_name = None

        print("sample type sequencing: ", sample.is_type_fastq_gz_sequencing())

        try:
            (job_name_wait, job_name) = user.profile.get_name_sge_seq(
                Profile.SGE_PROCESS_clean_sample, Profile.SGE_SAMPLE
            )
            if sample.is_type_fastq_gz_sequencing():  # default is Illumina
                taskID = process_SGE.set_run_trimmomatic_species(
                    sample, user, job_name)
            else:  # Minion, codify with other
                taskID = process_SGE.set_run_clean_minion(
                    sample, user, job_name)

            sample.is_ready_for_projects = True
            sample.save()

        except Exception as e:
            self.logger.error("Fail to run: ProcessSGE - " + str(e))
            return

        # refresh sample list for this user
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
            benchmark_softwaretree = SoftwareTree.objects.filter(
                name="default")
        except SoftwareTree.DoesNotExist:
            benchmark_softwaretree = SoftwareTree.objects.create(
                name="default",
                owner=User.objects.get(username="admin"),
            )
            benchmark_softwaretree.save()

        return benchmark_softwaretree
