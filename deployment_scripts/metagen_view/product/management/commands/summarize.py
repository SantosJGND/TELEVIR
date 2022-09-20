import os
from datetime import date

import pandas as pd
from django.contrib.auth.models import User
from django.core.management.base import BaseCommand
from product.utils import Run_Main_from_Fastq_Input
from result_display.models import FinalReport, Projects, RunDetail, RunMain, Sample


def read_parameters(row):
    try:
        params = pd.read_csv(row.params_file_path, sep="\t")
        params["run"] = row.run_name
        params["run_id"] = row.run_id
        params["sample"] = row.sample_name
        params["sample_id"] = row.sample_id
        params["project"] = row.project_name
        params["project_id"] = row.project_id
        return params

    except FileNotFoundError:
        return None


def collect_parameters(queryset_df):
    all_parameters = []
    unique_runs_df = queryset_df.drop_duplicates(subset=["run_id"])
    for index, row in unique_runs_df.iterrows():
        params = read_parameters(row)
        if params is not None:
            all_parameters.append(params)
    return pd.concat(all_parameters)


class Command(BaseCommand):
    help = "deploy run"

    def add_arguments(self, parser):
        parser.add_argument(
            "-u",
            "--user",
            type=str,
            help="user to get info from",
        )
        parser.add_argument(
            "-p",
            "--project",
            default="all",
            type=str,
            help="project to target",
        )

    def handle(self, *args, **options):
        ###
        #

        user = User.objects.get(username=options.get("user"))

        if options.get("project") == "all":

            projects = Projects.objects.filter(created_by=user)

        else:
            projects = Projects.objects.get(
                name=options.get("project"), created_by=user
            )

        samples = Sample.objects.filter(project__in=projects)
        mainruns = RunMain.objects.filter(sample__in=samples)
        # rundetails = pd.DataFrame.from_records(
        #    RunDetail.objects.filter(run_main__in=mainruns).values()
        # )
        reports = FinalReport.objects.filter(run__in=mainruns)

        ####
        projects = pd.DataFrame(projects.values()).rename(
            {"name": "project_name"}, axis=1
        )
        samples = pd.DataFrame(samples.values()).rename({"name": "sample_name"}, axis=1)
        mainruns = pd.DataFrame(mainruns.values()).rename({"name": "run_name"}, axis=1)
        # rundetails= pd.DataFrame(rundetails.values())
        reports = pd.DataFrame(reports.values()).rename({"id": "report_id"}, axis=1)

        ####
        reports_df = (
            pd.merge(mainruns, reports, left_on="id", right_on="run_id")
            .drop("id", axis=1)
            .drop("sample_id_y", axis=1)
            .rename({"sample_id_x": "sample_id"}, axis=1)
        )

        reports_df = (
            pd.merge(samples, reports_df, left_on="id", right_on="sample_id")
            .drop("id", axis=1)
            .drop("project_id_y", axis=1)
            .rename({"project_id_x": "project_id"}, axis=1)
        )
        reports_df = pd.merge(
            projects, reports_df, left_on="id", right_on="project_id"
        ).drop("id", axis=1)

        ####
        all_parameters = collect_parameters(reports_df)

        ####
        reports_df.to_csv("all_reports.tsv", index=False, sep="\t")

        all_parameters.to_csv("all_parameters.tsv", index=False, sep="\t")

        # all_parameters=
