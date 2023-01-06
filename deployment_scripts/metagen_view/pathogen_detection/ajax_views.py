##################################
## O código AJAX dentro do django
import os

from django.http import JsonResponse
from django.utils.safestring import mark_safe
from pathogen_detection.models import ReferenceMap_Main, RunMain, Sample


def show_igv(request, sample_name, run_name, reference):
    """
    get data for IGV
    """
    print("AJAX request to show_igv")

    is_ajax = request.META.get("HTTP_X_REQUESTED_WITH") == "XMLHttpRequest"
    print(request.META.get("HTTP_X_REQUESTED_WITH"))
    print("is_ajax:", is_ajax)

    data = {"is_ok": False}
    key_with_project_sample_id = "project_sample_id"
    # if key_with_project_sample_id in request.GET:
    try:
        sample = Sample.objects.get(name=sample_name)
        run = RunMain.objects.get(name=run_name, sample=sample)
        ref_map = ReferenceMap_Main.objects.get(
            reference=reference, sample=sample, run=run
        )

        path_name_bam = ref_map.bam_file_path
        path_name_bai = ref_map.bai_file_path
        path_name_reference = ref_map.fasta_file_path
        path_name_reference_index = ref_map.fai_file_path
        reference_name = ref_map.reference
        print("path_name_bam:", path_name_bam)
        print("path_name_bai:", path_name_bai)
        print("path_name_reference:", path_name_reference)
        print("path_name_reference_index:", path_name_reference_index)
        print("reference_name:", reference_name)

        # path_name_reference = (
        #    project_sample.project.reference.get_reference_fasta(
        #        TypePath.MEDIA_URL
        #    )
        # )
        # path_name_reference_index = (
        #    project_sample.project.reference.get_reference_fasta_index(
        #        TypePath.MEDIA_URL
        #    )
        # )
        #
        data["is_ok"] = True
        data["path_bam"] = mark_safe(request.build_absolute_uri(path_name_bam))
        data["path_reference"] = mark_safe(
            request.build_absolute_uri(path_name_reference)
        )
        data["path_reference_index"] = mark_safe(
            request.build_absolute_uri(path_name_reference_index)
        )
        data["reference_name"] = reference_name
        data["sample_name"] = sample_name

        #### other files
        data["bam_file_id"] = mark_safe(
            '<strong>Bam file:</strong> <a href="{}" download="{}"> {}</a>'.format(
                path_name_bam,
                os.path.basename(path_name_bam),
                os.path.basename(path_name_bam),
            )
        )
        data["bai_file_id"] = mark_safe(
            '<strong>Bai file:</strong> <a href="{}" download="{}"> {}</a>'.format(
                path_name_bai,
                os.path.basename(path_name_bai),
                os.path.basename(path_name_bai),
            )
        )
        data["reference_id"] = mark_safe(
            '<strong>Reference:</strong> <a href="{}" download="{}"> {}</a>'.format(
                path_name_reference,
                os.path.basename(path_name_reference),
                os.path.basename(path_name_reference),
            )
        )
        data["reference_index_id"] = mark_safe(
            '<strong>Ref. index:</strong> <a href="{}" download="{}"> {}</a>'.format(
                path_name_reference_index,
                os.path.basename(path_name_reference_index),
                os.path.basename(path_name_reference_index),
            )
        )
    except ReferenceMap_Main.DoesNotExist as e:
        pass
    return JsonResponse(data)
