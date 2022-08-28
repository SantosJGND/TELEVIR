"""
Created on Jan 7, 2018

@author: joao santos
"""

import result_display.views as result_views
from django.urls import path, re_path
from result_display import ajax_views

from product import views

urlpatterns = [
    path(
        "",
        views.Upload_file.as_view(),
        name="televir_upload_data",
    ),
    path(
        "submit_job/<slug:project_name>/<slug:pk>",
        views.submit_view,
        name="televir_submit_job",
    ),
    path("all_projects", views.Project_page.as_view(), name="televir_projects_main"),
    path("<slug:project_name>", views.ProjectView, name="televir_project_samples"),
    path(
        "projects/<slug:project_name>/<slug:sample_name>/",
        views.Sample_Main,
        name="televir_sample_main",
    ),
    ##
    path(
        "<slug:project>/<slug:sample>/<slug:name>",
        views.Sample_detail,
        name="televir_sample_detail",
    ),
]
