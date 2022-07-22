"""metagen_view URL Configuration

The `urlpatterns` list routes URLs to views. For more information please see:
    https://docs.djangoproject.com/en/4.0/topics/http/urls/
Examples:
Function views
    1. Add an import:  from my_app import views
    2. Add a URL to urlpatterns:  path('', views.home, name='home')
Class-based views
    1. Add an import:  from other_app.views import Home
    2. Add a URL to urlpatterns:  path('', Home.as_view(), name='home')
Including another URLconf
    1. Import the include() function: from django.urls import include, path
    2. Add a URL to urlpatterns:  path('blog/', include('blog.urls'))
"""
from django.conf import settings
from django.conf.urls import include
from django.conf.urls.static import static
from django.contrib import admin
from django.urls import path, re_path
from result_display import igv_app

urlpatterns = [
    re_path("^admin/", admin.site.urls),
    re_path("^", include("result_display.urls")),
    path("file_upload/", include("file_upload.urls")),
    path("django_plotly_dash/", include("django_plotly_dash.urls")),
] + static(settings.STATIC_URL, document_root=settings.STATIC_ROOT)
