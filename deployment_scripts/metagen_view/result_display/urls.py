"""
Created on Jan 7, 2018

@author: joao santos
"""

from django.urls import path

from result_display import ajax_views, views

urlpatterns = [
    path("", views.Project_page.as_view(), name="projects_main"),
    # path("login", views.LoginView.as_view(), name="inhouse_login"),
    path(
        "igv_display",
        views.IGV_display,
        name="igv_browser",
    ),  ## get values for IGV
    path(
        "show_igv_<slug:sample_name>/<slug:run_name>/<slug:reference>",
        ajax_views.show_igv,
        name="show_igv",
    ),  ## get values for IGV
    path("download_file", views.download_file, name="download_file"),  ##
    path("download_file_igv", views.download_file_igv, name="download_file_igv"),
    path("<slug:project_name>", views.MainPage, name="project_samples"),
    path(
        "<slug:project_name>/sample_<slug:sample_name>",
        views.Sample_main,
        name="sample_main",
    ),
    path(
        "<slug:project>/sample_<slug:sample>/<slug:name>",
        views.Sample_detail,
        name="sample_detail",
    ),
    ##
]
