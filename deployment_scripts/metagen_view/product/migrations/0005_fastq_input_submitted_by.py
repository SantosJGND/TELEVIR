# Generated by Django 3.2.14 on 2022-08-23 15:15

from django.conf import settings
from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    dependencies = [
        migrations.swappable_dependency(settings.AUTH_USER_MODEL),
        ('product', '0004_remove_fastq_input_description'),
    ]

    operations = [
        migrations.AddField(
            model_name='fastq_input',
            name='submitted_by',
            field=models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.CASCADE, related_name='submitted_by', to=settings.AUTH_USER_MODEL),
        ),
    ]
