# Generated by Django 3.2.14 on 2022-08-23 12:29

from django.db import migrations


class Migration(migrations.Migration):

    dependencies = [
        ('product', '0002_alter_fastq_input_file_r2'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='fastq_input',
            name='email',
        ),
    ]