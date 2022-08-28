import os
import shlex
import subprocess
from webbrowser import BackgroundBrowser

from background_task import background
from django.http import HttpResponseRedirect, JsonResponse
from django.shortcuts import render
from django.urls import reverse
from django.views.generic import ListView
from django.views.generic.edit import FormView
from django_tables2 import RequestConfig
from result_display.constants_settings import ConstantsSettings
from result_display.models import (
    ContigClassification,
    FinalReport,
    Projects,
    ReadClassification,
    ReferenceMap_Main,
    RunAssembly,
    RunDetail,
    RunMain,
    RunRemapMain,
    Sample,
    SampleQC,
)
from result_display.tables import SampleQCTable

from product.file_management import Ephemeral_Project_Manager
from product.forms import UploadFileForm
from product.models import Fastq_Input, Processed, Submitted
from product.tables import RunMainTable, SampleTable
from product.utils import Run_Main_from_Fastq_Input


class Upload_file(FormView):
    template_name = "product/file_upload.html"
    success_url = "sample_view/"
    form_class = UploadFileForm

    manager = Ephemeral_Project_Manager()

    def post(self, request, *args, **kwargs):
        form_class = self.get_form_class()
        form = self.get_form(form_class)

        if form.is_valid():

            file_r1 = form.cleaned_data["file_r1"]
            file_r2 = form.cleaned_data["file_r2"]
            technology = form.cleaned_data["technology"]

            # Save the file to the media folder
            file_r1 = self.manager.handle_uploaded_file(file_r1)
            file_r2 = self.manager.handle_uploaded_file(file_r2) if file_r2 else None
            # Save the file to the database

            new_fastq_input = Fastq_Input(
                technology=technology,
                file_r1=file_r1,
                file_r2=file_r2,
                name=request.user.username,
                project_name=form.cleaned_data["project_name"],
            )

            new_fastq_input.save()
            # Run the main script
            # self.manager.submit_job(new_fastq_input)
            request.session["submit_index"] = new_fastq_input.pk

            url = reverse(
                "televir_submit_job",
                kwargs={
                    "project_name": new_fastq_input.project_name,
                    "pk": new_fastq_input.pk,
                },
            )

            return HttpResponseRedirect(url)
        else:
            return self.form_invalid(form)


@background(schedule=1)
def add_deployment_task(pk):
    # fastq_input = Fastq_Input.objects.get(pk=pk)
    run_task = Run_Main_from_Fastq_Input(pk)
    run_task.Submit()


def process_tasks():
    process_tasks_cmd = "python manage.py process_tasks --duration 60"
    process_tasks_args = shlex.split(process_tasks_cmd)
    process_tasks_subprocess = subprocess.Popen(process_tasks_args)


from background_task.models import Task


def simplify_name(name):
    return (
        name.replace("_", "_")
        .replace("-", "_")
        .replace(" ", "_")
        .replace(".", "_")
        .lower()
    )


def submit_view(request, project_name, pk):
    """
    home page
    """

    template_name = "product/sample_main.html"
    fastq_input = Fastq_Input.objects.get(pk=pk)

    sample_name = os.path.basename(fastq_input.file_r1.path)
    sample_name = simplify_name(sample_name)

    tasks = Task.objects.filter(verbose_name=f"Run_deployment_task_{pk}")

    input_submitted = Submitted.objects.filter(fastq_input=fastq_input)
    if len(input_submitted) == 0:

        if len(tasks) == 0:
            print("No input submitted")
            print("pk: " + str(pk))
            add_deployment_task(pk, verbose_name=f"Run_deployment_task_{pk}")
            process_tasks()
            return render(
                request,
                template_name,
                {
                    "waiting": True,
                    "error": False,
                    "project_name": project_name,
                },
            )
        else:
            print("Input submitted")
            return render(
                request,
                template_name,
                {
                    "waiting": True,
                    "error": False,
                    "project_name": project_name,
                },
            )

    input_processed = Processed.objects.filter(fastq_input=fastq_input)
    if len(input_processed) == 0:

        if len(tasks) == 0:
            print("No input processed and task is running")

            return render(
                request,
                template_name,
                {
                    "waiting": False,
                    "error": True,
                    "project_name": project_name,
                },
            )

        else:
            print("Input running")
            return render(
                request,
                template_name,
                {
                    "waiting": True,
                    "error": False,
                    "project_name": project_name,
                },
            )

    project = Projects.objects.get(name=project_name)
    sample = Sample.objects.get(project=project, name=sample_name)

    try:
        runs = RunMain.objects.filter(sample=sample, project=project)
    except RunMain.DoesNotExist:
        runs = None

    try:
        sampleqc = SampleQC.objects.filter(sample=sample)

    except SampleQC.DoesNotExist:
        sampleqc = None

    sampleqc_table = SampleQCTable(sampleqc)
    runs = RunMainTable(runs)
    RequestConfig(
        request, paginate={"per_page": ConstantsSettings.PAGINATE_NUMBER}
    ).configure(runs)

    return render(
        request,
        template_name,
        {
            "waiting": False,
            "error": False,
            "runs": runs,
            "qc_table": sampleqc_table,
            "sampleqc": [list(sampleqc.values())][0][0],
            "name": sample_name,
            "project_main": True,
            "project_name": project_name,
        },
    )


class Project_page(ListView):
    """Project page"""

    template_name = "product/projects_view.html"

    def get_queryset(self):
        return None

    def get_context_data(self, **kwargs):
        context = super(Project_page, self).get_context_data(**kwargs)

        projects = Projects.objects.filter(project_type=Projects.EXTERNAL)
        context["projects"] = projects

        return context


def ProjectView(request, project_name):
    """
    home page
    """

    template_name = "product/project_page.html"

    samples = Sample.objects.filter(project__name=project_name)

    samples = SampleTable(samples)
    context = {}
    context["samples"] = samples
    context["project_name"] = project_name

    return render(
        request,
        template_name,
        context,
    )


def Sample_Main(requesdst, project_name, sample_name):
    """Sample page"""
    template_name = "product/sample_main.html"

    project = Projects.objects.get(name=project_name, project_type=Projects.EXTERNAL)
    sample = Sample.objects.get(project=project, name=sample_name)

    try:
        runs = RunMain.objects.filter(sample=sample, project=project)
    except RunMain.DoesNotExist:
        runs = None

    try:
        sampleqc = SampleQC.objects.filter(sample=sample)

    except SampleQC.DoesNotExist:
        sampleqc = None

    sampleqc_table = SampleQCTable(sampleqc)
    runs = RunMainTable(runs)
    RequestConfig(
        requesdst, paginate={"per_page": ConstantsSettings.PAGINATE_NUMBER}
    ).configure(runs)

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


def Sample_detail(requesdst, project="", sample="", name=""):
    """
    home page
    """
    template_name = "product/sample_detail.html"
    project_main = Projects.objects.get(name=project)
    sample_main = Sample.objects.get(name=sample, project__name=project)
    #
    run_main = RunMain.objects.get(project__name=project, sample=sample_main, name=name)
    #
    run_detail = RunDetail.objects.get(sample=sample_main, run=run_main)
    #
    run_assembly = RunAssembly.objects.get(sample=sample_main, run=run_main)
    #
    print(run_assembly)
    run_remap = RunRemapMain.objects.get(sample=sample_main, run=run_main)
    #
    read_classification = ReadClassification.objects.get(
        sample=sample_main, run=run_main
    )
    #
    final_report = FinalReport.objects.filter(sample=sample_main, run=run_main)
    #
    contig_classification = ContigClassification.objects.get(
        sample=sample_main, run=run_main
    )
    #

    reference_remap_main = ReferenceMap_Main.objects.filter(
        sample=sample_main, run=run_main
    )
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
