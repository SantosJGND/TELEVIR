# Generated by Django 3.2.14 on 2022-09-15 00:36

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('result_display', '0025_auto_20220828_1106'),
    ]

    operations = [
        migrations.AlterField(
            model_name='finalreport',
            name='covplot',
            field=models.CharField(blank=True, max_length=250, null=True),
        ),
        migrations.AlterField(
            model_name='finalreport',
            name='refa_dotplot',
            field=models.CharField(blank=True, max_length=250, null=True),
        ),
    ]
