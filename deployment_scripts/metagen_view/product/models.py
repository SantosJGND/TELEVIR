from django.contrib.auth.models import User
from django.db import models
from django.urls import reverse

# Create your models here.


class Fastq_Input(models.Model):
    CATEGORIES = (("ILL", "Illumina"), ("ONT", "nanopore"))
    technology = models.CharField(max_length=5, choices=CATEGORIES)
    file_r1 = models.FileField(upload_to="files/")
    file_r2 = models.FileField(upload_to="files/", blank=True, null=True)

    name = models.CharField(max_length=254)
    project_name = models.CharField(max_length=254)
    date_created = models.DateTimeField(auto_now_add=True)
    date_modified = models.DateTimeField(auto_now=True)
    submitted_by = models.ForeignKey(
        User,
        on_delete=models.CASCADE,
        related_name="submitted_by",
        null=True,
        blank=True,
    )

    def __str__(self):
        return self.name


class Submitted(models.Model):
    fastq_input = models.ForeignKey(
        Fastq_Input, on_delete=models.CASCADE, related_name="submitted_fastq_input"
    )
    date_submitted = models.DateTimeField(auto_now_add=True)

    def __str__(self):
        return self.fastq_input.name


class Processed(models.Model):

    fastq_input = models.ForeignKey(Fastq_Input, on_delete=models.CASCADE)
    date_time = models.DateTimeField(auto_now_add=True)

    def __str__(self):
        return self.fastq_input.name
