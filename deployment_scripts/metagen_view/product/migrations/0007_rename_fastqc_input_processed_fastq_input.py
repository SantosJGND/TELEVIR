# Generated by Django 3.2.14 on 2022-08-24 11:47

from django.db import migrations


class Migration(migrations.Migration):

    dependencies = [
        ('product', '0006_processed_submitted'),
    ]

    operations = [
        migrations.RenameField(
            model_name='processed',
            old_name='fastqc_input',
            new_name='fastq_input',
        ),
    ]