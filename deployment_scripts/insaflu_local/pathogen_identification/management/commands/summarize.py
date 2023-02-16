
import pandas as pd
from django.contrib.auth.models import User
from django.core.management.base import BaseCommand
from pathogen_identification.models import (FinalReport, ParameterSet,
                                            PIProject_Sample, Projects,
                                            RunMain)
from pathogen_identification.utilities.utilities_pipeline import Utils_Manager


def read_parameters(row):
    try:
        params = pd.read_csv(row.params_file_path, sep="\t")
        params["run"] = row.run_name
        params["run_id"] = row.run_id
        params["sample"] = row.sample_name_x
        params["sample_id"] = row.sample_id
        params["project"] = row.project_name
        params["project_id"] = row.project_id
        return params

    except FileNotFoundError:
        return None


def collect_parameters(queryset_df: pd.DataFrame) -> pd.DataFrame:
    all_parameters = []
    unique_runs_df = queryset_df.drop_duplicates(subset=["run_id"])
    for index, row in unique_runs_df.iterrows():
        params = read_parameters(row)
        if params is not None:
            all_parameters.append(params)
    all_parameters = pd.concat(all_parameters).rename(
        columns={"sample": "sample_name"})
    return all_parameters


def collect_parameters_project(project: Projects):
    utils = Utils_Manager()
    technology = project.technology
    user = project.owner
    local_tree = utils.generate_project_tree(technology, project, user)
    local_paths = local_tree.get_all_graph_paths_explicit()

    tree_makeup = local_tree.makeup

    pipeline_tree = utils.generate_software_tree(technology, tree_makeup)
    all_paths = utils.get_all_technology_pipelines(technology, tree_makeup)
    pipeline_tree_index = utils.get_software_tree_index(
        technology, tree_makeup)

    parameterset = ParameterSet.objects.filter(
        project=project, status=ParameterSet.STATUS_FINISHED)

    project_params = []
    print(project.name)
    print(len(parameterset))
    print(all_paths.keys())
    for ps in parameterset:
        print(ps.leaf.pk)
        params = all_paths.get(ps.leaf.index, None)
        print(params)

        if params is None:
            continue
        run = RunMain.objects.get(parameter_set=ps)

        params["run"] = run.name
        params["run_id"] = run.pk
        params["sample"] = ps.sample.name
        params["sample_id"] = ps.sample.pk
        params["project"] = ps.project.name
        params["project_id"] = ps.project.pk

        project_params.append(params)
    if len(project_params):
        project_params = pd.concat(project_params)

    return project_params


def projects_params_summary(project_query):
    all_params = []

    for project in project_query:
        project_params = collect_parameters_project(project)
        if len(project_params):
            all_params.append(project_params)

    if len(all_params):
        all_params = pd.concat(all_params)
    return all_params


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

        parser.add_argument(
            "--tech",
            default="nanopore",
            type=str,
            help="project to target",
        )

    def handle(self, *args, **options):
        ###
        #

        user = User.objects.get(username=options.get("user"))

        if options.get("project") == "all":

            projects_query = Projects.objects.filter(owner=user)

        else:
            projects_query = Projects.objects.filter(
                name=options.get("project"), owner=user
            )

        samples = PIProject_Sample.objects.filter(project__in=projects_query)
        mainruns = RunMain.objects.filter(sample__in=samples)

        reports = FinalReport.objects.filter(run__in=mainruns)

        ####
        projects = pd.DataFrame(projects_query.values()).rename(
            {"name": "project_name"}, axis=1
        )
        samples = pd.DataFrame(samples.values()).rename(
            {"name": "sample_name"}, axis=1)
        mainruns = pd.DataFrame(mainruns.values()).rename(
            {"name": "run_name"}, axis=1)
        # rundetails= pd.DataFrame(rundetails.values())
        reports = pd.DataFrame(reports.values()).rename(
            {"id": "report_id"}, axis=1)

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

        all_parameters = projects_params_summary(projects_query)

        ####
        reports_df.to_csv("all_reports.tsv", index=False, sep="\t")

        all_parameters.to_csv("all_parameters.tsv", index=False, sep="\t")

        #
