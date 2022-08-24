import os
from typing import Type

from result_display.models import Projects


def Update_project(project_directory_path, user: str = "admin", submit_index: int = 1):
    """Updates the project"""
    project_directory_path = os.path.dirname(project_directory_path)
    project_name = os.path.basename(project_directory_path)
    project_name_simple = project_name.replace(".", "_").replace(":", "_")

    try:
        project = Projects.objects.get(name=project_name)
    except Projects.DoesNotExist:
        project = Projects(
            name=project_name,
            full_path=project_directory_path,
            submit_index=submit_index,
        )
        project.save()
