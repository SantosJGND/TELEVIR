# Generated by Django 3.2.14 on 2022-08-05 00:49

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('result_display', '0017_runremapmain_remap_plan'),
    ]

    operations = [
        migrations.AlterField(
            model_name='rundetail',
            name='input',
            field=models.CharField(blank=True, max_length=200, null=True),
        ),
        migrations.AlterField(
            model_name='rundetail',
            name='merged_files',
            field=models.CharField(blank=True, max_length=200, null=True),
        ),
        migrations.AlterField(
            model_name='rundetail',
            name='processed',
            field=models.CharField(blank=True, max_length=200, null=True),
        ),
        migrations.AlterField(
            model_name='rundetail',
            name='processing_final',
            field=models.CharField(blank=True, max_length=200, null=True),
        ),
        migrations.AlterField(
            model_name='rundetail',
            name='sift_removed_pprc',
            field=models.CharField(blank=True, max_length=200, null=True),
        ),
    ]
