# -*- coding: utf-8 -*-
# Generated by Django 1.11.18 on 2021-10-17 12:25
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('managing_files', '0022_project_number_passed_sequences'),
    ]

    operations = [
        migrations.AlterField(
            model_name='software',
            name='version',
            field=models.CharField(max_length=200),
        ),
    ]