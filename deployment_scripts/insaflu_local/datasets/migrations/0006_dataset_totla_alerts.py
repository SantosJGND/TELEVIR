# -*- coding: utf-8 -*-
# Generated by Django 1.11.18 on 2022-06-08 09:29
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('datasets', '0005_datasetconsensus_is_project_sample_finished'),
    ]

    operations = [
        migrations.AddField(
            model_name='dataset',
            name='totla_alerts',
            field=models.IntegerField(default=0),
        ),
    ]