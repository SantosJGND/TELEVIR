import mimetypes
import os
from typing import Final

from django import forms
from django.forms.models import model_to_dict
from django.http import Http404, HttpResponseNotFound
from django.http.response import HttpResponse
from django.shortcuts import render
from django.utils.safestring import mark_safe
from django.views.generic import ListView

# Create your views here.
from metagen_view.settings import STATICFILES_DIRS

from result_display import igv_app
from result_display.models import (
    QC_REPORT,
    ContigClassification,
    FinalReport,
    Projects,
    ReadClassification,
    ReferenceContigs,
    ReferenceMap_Main,
    RunAssembly,
    RunDetail,
    RunMain,
    RunRemapMain,
    Sample,
    SampleQC,
)
from result_display.tables import ContigTable, RunMainTable, SampleQCTable, SampleTable

# from result_display.forms import IGVform


class IGVform(forms.Form):
    sample_name = forms.CharField(max_length=100)
    run_name = forms.CharField(max_length=100)
    reference = forms.CharField(max_length=100)
    unique_id = forms.CharField(max_length=100)
    project_name = forms.CharField(max_length=100)


class download_form(forms.Form):
    file_path = forms.CharField(max_length=300)

    class Meta:

        widgets = {
            "myfield": forms.TextInput(
                attrs={"style": "border-color:darkgoldenrod; border-radius: 10px;"}
            ),
        }


################################################################


class Project_page(ListView):
    """Project page"""

    template_name = "result_display/projects_view.html"

    def get_queryset(self):
        return None

    def get_context_data(self, **kwargs):
        context = super(Project_page, self).get_context_data(**kwargs)

        projects = Projects.objects.all()
        context["projects"] = projects

        return context


def MainPage(request, project_name):
    """
    home page
    """

    template_name = "result_display/main_page.html"

    try:
        samples = Sample.objects.filter(project__name=project_name)
    except Sample.DoesNotExist:
        samples = parser.get_samples()

    samples = SampleTable(samples)
    context = {}
    context["samples"] = samples
    context["project_name"] = project_name

    return render(
        request,
        template_name,
        context,
    )


def Sample_main(requesdst, project_name, sample_name):
    """
    sample main page
    """
    template_name = "result_display/sample_main.html"

    try:
        runs = RunMain.objects.filter(
            sample__name=sample_name, project__name=project_name
        )
    except RunMain.DoesNotExist:
        runs = None

    try:
        sampleqc = SampleQC.objects.filter(sample__name=sample_name)

    except SampleQC.DoesNotExist:
        sampleqc = None

    sampleqc_table = SampleQCTable(sampleqc)
    runs = RunMainTable(runs)

    return render(
        requesdst,
        template_name,
        {
            "runs": runs,
            "qc_table": sampleqc_table,
            "sampleqc": [list(sampleqc.values())][0][0],
            "name": sample_name,
            "project_main": True,
            "project_name": project_name,
        },
    )


def sample_QCall(requestdst, project):
    """
    sample main page
    """
    template_name = "fastqc_html/fastqc_all.html"

    try:
        sampleqc = SampleQC.objects.all()

    except SampleQC.DoesNotExist:
        sampleqc = None

    return render(
        requestdst,
        template_name,
        {"all_reports": sampleqc, "project": project},
    )


def Project_reports(requesdst, project):
    """
    sample main page
    """
    template_name = "result_display/allreports_table.html"

    all_reports = FinalReport.objects.filter(run__project__name=project)

    return render(
        requesdst,
        template_name,
        {"all_reports": all_reports, "project": project},
    )


def display_fastqc_report(requesdst, name: str, report_source: str):
    """display input fastqc report"""

    template_name = "fastqc_html/fastqc_report.html"

    sample = Sample.objects.get(name=name)

    report = QC_REPORT.objects.get(sample=sample, report_source=report_source)

    return render(
        requesdst,
        template_name,
        {"fastqc_input": report, "sample_main": True},
    )


def Sample_detail(requesdst, project="", sample="", name=""):
    """
    home page
    """
    template_name = "result_display/sample_detail.html"
    project_main = Projects.objects.get(name=project)
    sample_main = Sample.objects.get(name=sample, project__name=project)
    #
    run_main = RunMain.objects.get(project__name=project, sample=sample_main, name=name)
    #
    run_detail = RunDetail.objects.get(sample=sample_main, run=run_main)
    #
    run_assembly = RunAssembly.objects.get(sample=sample_main, run=run_main)
    #
    contig_classification = ContigClassification.objects.get(
        sample=sample_main, run=run_main
    )
    #
    read_classification = ReadClassification.objects.get(
        sample=sample_main, run=run_main
    )
    #
    run_remap = RunRemapMain.objects.get(sample=sample_main, run=run_main)
    #
    reference_remap_main = ReferenceMap_Main.objects.filter(
        sample=sample_main, run=run_main
    )
    #
    final_report = FinalReport.objects.filter(sample=sample_main, run=run_main)
    #

    context = {
        "project": project,
        "run_name": name,
        "sample": sample,
        "run_main": run_main,
        "run_detail": run_detail,
        "assembly": run_assembly,
        "contig_classification": contig_classification,
        "read_classification": read_classification,
        "run_remap": run_remap,
        "reference_remap_main": reference_remap_main,
        "final_report": final_report,
    }

    return render(
        requesdst,
        template_name,
        context,
    )


def IGV_display(requestdst):
    """display python plotly app"""
    template_name = "result_display/IGV.html"

    if requestdst.method == "POST":
        form = IGVform(requestdst.POST)
        if form.is_valid():
            print(form.cleaned_data)
            project_name = form.cleaned_data.get("project_name")
            sample_name = form.cleaned_data.get("sample_name")
            run_name = form.cleaned_data.get("run_name")
            reference = form.cleaned_data.get("reference")
            unique_id = form.cleaned_data.get("unique_id")

            data = {"is_ok": False}
            print(project_name)
            print(Sample.objects.filter(name=sample_name, project__name=project_name))
            sample = Sample.objects.get(project__name=project_name, name=sample_name)
            run = RunMain.objects.get(name=run_name, sample=sample)
            ref_map = ReferenceMap_Main.objects.get(
                reference=reference, sample=sample, run=run
            )
            final_report = FinalReport.objects.get(
                sample=sample, run=run, unique_id=unique_id
            )

            path_name_bam = ref_map.bam_file_path
            path_name_bai = ref_map.bai_file_path
            path_name_reference = ref_map.fasta_file_path
            path_name_reference_index = ref_map.fai_file_path
            reference_name = final_report.reference_contig_str

            data["is_ok"] = True

            data["path_reference"] = path_name_reference
            data["path_reference_index"] = path_name_reference_index
            data["path_bam"] = path_name_bam
            data["path_bai"] = path_name_bai

            data["reference_name"] = sample_name
            data["sample_name"] = final_report.reference_contig_str

            #### other files
            data["bam_file_id"] = mark_safe(
                '<strong>Bam file:</strong> <a href="{}" filename="{}">{}</a>'.format(
                    "download_file",
                    os.path.basename(path_name_bam),
                    os.path.basename(path_name_bam),
                )
            )
            data["bai_file_id"] = mark_safe(
                '<strong>Bai file:</strong> <a href="{}" filename="{}">{}</a>'.format(
                    path_name_bai,
                    os.path.basename(path_name_bai),
                    os.path.basename(path_name_bai),
                )
            )
            data["reference_id"] = mark_safe(
                '<strong>Reference:</strong> <a href="{}" filename="{}">{}</a>'.format(
                    path_name_reference,
                    os.path.basename(path_name_reference),
                    os.path.basename(path_name_reference),
                )
            )
            data["reference_index_id"] = mark_safe(
                '<strong>Ref. index:</strong> <a href="{}" filename="{}">{}</a>'.format(
                    path_name_reference_index,
                    os.path.basename(path_name_reference_index),
                    os.path.basename(path_name_reference_index),
                )
            )

            return render(
                requestdst,
                template_name,
                {
                    "sample_name": final_report.reference_contig_str,
                    "run_name": run_name,
                    "reference": reference,
                    "data": data,
                },
            )
    else:
        form = IGVform()


def download_file_igv(requestdst):
    """download fasta file"""
    print(requestdst.method)
    if requestdst.method == "POST":
        form = download_form(requestdst.POST)

        if form.is_valid():
            filepath = form.cleaned_data.get("file_path")

            BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
            DLDIR = os.path.join(BASE_DIR, STATICFILES_DIRS[0], "igv_files")
            filepath = os.path.join(DLDIR, filepath)

            if not os.path.exists(filepath):
                return HttpResponseNotFound(f"file {filepath} not found")

            path = open(filepath, "rb")
            # Set the mime type
            mime_type, _ = mimetypes.guess_type(filepath)
            # Set the return value of the HttpResponse
            response = HttpResponse(path, content_type=mime_type)
            # Set the HTTP header for sending to browser
            response[
                "Content-Disposition"
            ] = "attachment; filename=%s" % os.path.basename(filepath)
            # Return the response value
            return response


def download_file(requestdst):
    """download fasta file"""
    if requestdst.method == "POST":
        form = download_form(requestdst.POST)

        if form.is_valid():
            filepath = form.cleaned_data.get("file_path")

            BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

            print(BASE_DIR)
            print(filepath)

            filepath = BASE_DIR + filepath

            if "//" in filepath:
                filepath = "/" + filepath.split("//")[1]

            if not os.path.isfile(filepath):
                return HttpResponseNotFound(f"file {filepath} not found")

            path = open(filepath, "rb")
            # Set the mime type
            mime_type, _ = mimetypes.guess_type(filepath)
            # Set the return value of the HttpResponse
            response = HttpResponse(path, content_type=mime_type)
            # Set the HTTP header for sending to browser
            response[
                "Content-Disposition"
            ] = "attachment; filename=%s" % os.path.basename(filepath)
            # Return the response value
            return response


def Scaffold_Remap(requesdst, project="", sample="", run="", reference=""):
    """
    home page
    """
    template_name = "result_display/scaffold_remap.html"
    ##

    sample_main = Sample.objects.get(project__name=project, name=sample)
    #
    run_main = RunMain.objects.get(sample=sample_main, name=run)

    try:
        ref_main = ReferenceMap_Main.objects.get(
            reference=reference, sample=sample_main, run=run_main
        )
        map_db = ReferenceContigs.objects.filter(
            reference=ref_main,
            run=run_main,
        )
    except ReferenceMap_Main.DoesNotExist:

        return Http404("Sample not found")

    return render(
        requesdst,
        template_name,
        {
            "table": ContigTable(map_db),
        },
    )
