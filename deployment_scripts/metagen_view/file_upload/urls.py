from django.conf import settings
from django.conf.urls.static import static
from django.urls import path

from file_upload import views as uploader_views

urlpatterns = [
    path("", uploader_views.UploadView.as_view(), name="fileupload"),
] + static(settings.MEDIA_URL, document_root=settings.MEDIA_ROOT)
