import os
from typing import Type

from django.contrib.auth.models import User
from result_display.models import Projects


def Update_project(project_directory_path, user: str = "admin", submit_index: int = 1):
    """Updates the project"""
    project_directory_path = os.path.dirname(project_directory_path)
    project_name = os.path.basename(project_directory_path)
    project_name_simple = project_name.replace(".", "_").replace(":", "_")

    try:
        user = User.objects.get(username=user)
    except User.DoesNotExist:
        user = User.objects.get(username="admin")

    try:
        project = Projects.objects.get(
            name=project_name,
            created_by=user,
        )

    except Projects.DoesNotExist:
        project = Projects(
            name=project_name,
            full_path=project_directory_path,
            project_type=Projects.EXTERNAL,
            submit_index=submit_index,
            created_by=user,
        )
        project.save()
