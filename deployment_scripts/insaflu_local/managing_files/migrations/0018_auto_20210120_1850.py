# -*- coding: utf-8 -*-
# Generated by Django 1.11.18 on 2021-01-20 18:50
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('managing_files', '0017_sample_type_of_fastq'),
    ]

    operations = [
        migrations.AddField(
            model_name='project',
            name='last_change_date',
            field=models.DateTimeField(blank=True, null=True, verbose_name='Last change date'),
        ),
        migrations.AlterField(
            model_name='project',
            name='creation_date',
            field=models.DateTimeField(auto_now_add=True, verbose_name='Uploaded date'),
        ),
    ]