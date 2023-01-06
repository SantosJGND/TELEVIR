"""
Created on Jan 7, 2018

@author: joao santos
"""

import pathogen_detection.views as result_views
from django.urls import path, re_path
from pathogen_detection import ajax_views
from product import views

urlpatterns = [
    path("", views.entry_page, name="entry_page"),
    path(
        "sample_upload/",
        views.Upload_file.as_view(),
        name="televir_upload_data",
    ),
    path("all_projects/", views.Project_page.as_view(), name="televir_projects_main"),
    # path("login", views.LoginView.as_view(), name="televir_login"),
    path(
        "submit_job/<slug:project_name>/<slug:pk>",
        views.submit_view,
        name="televir_submit_job",
    ),
    path("<slug:project_name>", views.ProjectView, name="televir_project_samples"),
    path(
        "project_<project>/all_reports",
        views.Project_reports,
        name="all_project_reports",
    ),
    path(
        "projects/<slug:project_name>/<slug:sample_name>/",
        views.Sample_Main,
        name="televir_sample_main",
    ),
    path(
        "project_<project>/<sample_name>/all_reports",
        views.Sample_reports,
        name="all_sample_reports",
    ),
    path(
        "project_<project>/all_QC",
        views.sample_QCall,
        name="all_project_qc",
    ),
    path(
        "project_<project>/sample_<name>/<report_source>_fastqc_report",
        views.display_fastqc_report,
        name="display_fastqc_report",
    ),
    ##
    path(
        "<slug:project>/<slug:sample>/<slug:name>",
        views.Sample_detail,
        name="televir_sample_detail",
    ),
    path(
        "<slug:project>/sample_<slug:sample>/<slug:run>/<slug:reference>",
        views.Scaffold_Remap,
        name="scaffold_remap",
    ),
]
