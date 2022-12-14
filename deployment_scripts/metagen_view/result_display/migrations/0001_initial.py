# Generated by Django 3.2.14 on 2022-07-18 16:31

from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    initial = True

    dependencies = [
    ]

    operations = [
        migrations.CreateModel(
            name='Projects',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('name', models.CharField(blank=True, db_index=True, max_length=29, null=True)),
                ('full_path', models.CharField(blank=True, db_index=True, max_length=200, null=True)),
            ],
            options={
                'ordering': ['name'],
            },
        ),
        migrations.CreateModel(
            name='RunMain',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('suprun', models.CharField(blank=True, max_length=100, null=True)),
                ('name', models.CharField(blank=True, db_index=True, max_length=100, null=True)),
                ('params_file_path', models.CharField(blank=True, max_length=100, null=True)),
                ('enrichment', models.CharField(blank=True, max_length=10, null=True)),
                ('enrichment_performed', models.BooleanField(blank=True, null=True)),
                ('host_depletion', models.CharField(blank=True, max_length=10, null=True)),
                ('host_depletion_performed', models.BooleanField(blank=True, null=True)),
                ('processed_reads_r1', models.CharField(blank=True, max_length=200, null=True)),
                ('processed_reads_r2', models.CharField(blank=True, max_length=200, null=True)),
                ('reads_after_processing', models.CharField(blank=True, max_length=100, null=True)),
                ('reads_proc_percent', models.CharField(blank=True, max_length=100, null=True)),
                ('assembly_performed', models.CharField(blank=True, max_length=10, null=True)),
                ('assembly_method', models.CharField(blank=True, max_length=10, null=True)),
                ('assembly_max', models.CharField(blank=True, max_length=100, null=True)),
                ('read_classification', models.CharField(blank=True, max_length=10, null=True)),
                ('contig_classification', models.CharField(blank=True, max_length=15, null=True)),
                ('remap', models.CharField(blank=True, max_length=10, null=True)),
                ('finished', models.CharField(blank=True, max_length=10, null=True)),
                ('runtime', models.CharField(blank=True, max_length=100, null=True)),
                ('report', models.CharField(blank=True, max_length=100, null=True)),
                ('project', models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.CASCADE, to='result_display.projects')),
            ],
            options={
                'ordering': ['name'],
            },
        ),
        migrations.CreateModel(
            name='Sample',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('name', models.CharField(blank=True, db_index=True, max_length=100, null=True)),
                ('name_extended', models.CharField(blank=True, db_index=True, max_length=100, null=True)),
                ('type', models.CharField(blank=True, max_length=10, null=True)),
                ('combinations', models.IntegerField(blank=True, null=True)),
                ('input', models.TextField(blank=True, null=True)),
                ('technology', models.CharField(blank=True, max_length=100, null=True)),
                ('report', models.CharField(blank=True, max_length=100, null=True)),
                ('project', models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.CASCADE, to='result_display.projects')),
            ],
            options={
                'ordering': ['name'],
            },
        ),
        migrations.CreateModel(
            name='SampleQC',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('software', models.CharField(blank=True, max_length=100, null=True)),
                ('qc_type', models.CharField(blank=True, max_length=100, null=True)),
                ('encoding', models.CharField(blank=True, max_length=100, null=True)),
                ('input_reads', models.CharField(blank=True, max_length=100, null=True)),
                ('processed_reads', models.CharField(blank=True, max_length=100, null=True)),
                ('percent_passed', models.FloatField(blank=True, null=True)),
                ('sequence_length', models.CharField(blank=True, max_length=100, null=True)),
                ('percent_gc', models.FloatField(blank=True, null=True)),
                ('input_fastqc_report', models.FileField(blank=True, null=True, upload_to='input_fastqc_report')),
                ('processed_fastqc_report', models.FileField(blank=True, null=True, upload_to='processed_fastqc_report')),
                ('sample', models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.CASCADE, to='result_display.sample')),
            ],
            options={
                'ordering': ['sample'],
            },
        ),
        migrations.CreateModel(
            name='RunRemapMain',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('performed', models.BooleanField(default=False)),
                ('method', models.CharField(blank=True, max_length=10, null=True)),
                ('merged_log', models.CharField(blank=True, max_length=200, null=True)),
                ('remap_plan', models.CharField(blank=True, max_length=200, null=True)),
                ('found_total', models.IntegerField(blank=True, null=True)),
                ('coverage_minimum', models.IntegerField(blank=True, null=True)),
                ('coverage_maximum', models.IntegerField(blank=True, null=True)),
                ('success', models.IntegerField(blank=True, null=True)),
                ('run', models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.CASCADE, to='result_display.runmain')),
                ('sample', models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.CASCADE, to='result_display.sample')),
            ],
            options={
                'ordering': ['method'],
            },
        ),
        migrations.AddField(
            model_name='runmain',
            name='sample',
            field=models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.CASCADE, to='result_display.sample'),
        ),
        migrations.CreateModel(
            name='RunIndex',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('name', models.CharField(blank=True, db_index=True, max_length=100, null=True)),
                ('project', models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.CASCADE, to='result_display.projects')),
                ('sample', models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.CASCADE, to='result_display.sample')),
            ],
        ),
        migrations.CreateModel(
            name='RunDetail',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('name', models.CharField(blank=True, db_index=True, max_length=100, null=True)),
                ('max_depth', models.FloatField(blank=True, null=True)),
                ('max_depthR', models.FloatField(blank=True, null=True)),
                ('max_gaps', models.IntegerField(blank=True, null=True)),
                ('max_prop', models.FloatField(blank=True, null=True)),
                ('max_mapped', models.IntegerField(blank=True, null=True)),
                ('input', models.CharField(blank=True, max_length=100, null=True)),
                ('processed', models.CharField(blank=True, max_length=100, null=True)),
                ('processed_percent', models.FloatField(blank=True, null=True)),
                ('sift_preproc', models.BooleanField(blank=True, null=True)),
                ('sift_remap', models.BooleanField(blank=True, null=True)),
                ('sift_removed_pprc', models.CharField(blank=True, max_length=100, null=True)),
                ('processing_final', models.CharField(blank=True, max_length=100, null=True)),
                ('processing_final_percent', models.FloatField(blank=True, null=True)),
                ('merged', models.BooleanField(blank=True, null=True)),
                ('merged_number', models.IntegerField(blank=True, null=True)),
                ('merged_files', models.CharField(blank=True, max_length=100, null=True)),
                ('run', models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.CASCADE, to='result_display.runmain')),
                ('sample', models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.CASCADE, to='result_display.sample')),
            ],
            options={
                'ordering': ['name'],
            },
        ),
        migrations.CreateModel(
            name='RunAssembly',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('performed', models.BooleanField(default=False)),
                ('assembly_contigs', models.CharField(blank=True, max_length=200, null=True)),
                ('method', models.CharField(blank=True, max_length=10, null=True)),
                ('contig_number', models.IntegerField(blank=True, null=True)),
                ('contig_max', models.CharField(blank=True, max_length=100, null=True)),
                ('contig_min', models.CharField(blank=True, max_length=100, null=True)),
                ('contig_mean', models.CharField(blank=True, max_length=100, null=True)),
                ('contig_trim', models.CharField(blank=True, max_length=100, null=True)),
                ('run', models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.CASCADE, to='result_display.runmain')),
                ('sample', models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.CASCADE, to='result_display.sample')),
            ],
            options={
                'ordering': ['method'],
            },
        ),
        migrations.CreateModel(
            name='ReferenceMap_Main',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('reference', models.CharField(blank=True, db_index=True, max_length=100, null=True)),
                ('taxid', models.CharField(blank=True, max_length=100, null=True)),
                ('reference_contig_str', models.CharField(blank=True, max_length=100, null=True)),
                ('report', models.CharField(blank=True, max_length=100, null=True)),
                ('plotly_dotplot', models.TextField(blank=True, null=True)),
                ('bam_file_path', models.CharField(blank=True, max_length=100, null=True)),
                ('bai_file_path', models.CharField(blank=True, max_length=100, null=True)),
                ('fasta_file_path', models.CharField(blank=True, max_length=100, null=True)),
                ('fai_file_path', models.CharField(blank=True, max_length=100, null=True)),
                ('run', models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.CASCADE, to='result_display.runmain')),
                ('sample', models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.CASCADE, to='result_display.sample')),
            ],
            options={
                'ordering': ['reference'],
            },
        ),
        migrations.CreateModel(
            name='ReferenceContigs',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('contig', models.CharField(blank=True, max_length=100, null=True)),
                ('depth', models.CharField(blank=True, max_length=100, null=True)),
                ('depthr', models.CharField(blank=True, max_length=100, null=True)),
                ('coverage', models.CharField(blank=True, max_length=100, null=True)),
                ('reference', models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.CASCADE, to='result_display.referencemap_main')),
                ('run', models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.CASCADE, to='result_display.runmain')),
            ],
            options={
                'ordering': ['reference'],
            },
        ),
        migrations.CreateModel(
            name='ReadClassification',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('performed', models.BooleanField(default=False)),
                ('method', models.CharField(blank=True, max_length=10, null=True)),
                ('read_classification_report', models.CharField(blank=True, max_length=200, null=True)),
                ('classification_number', models.IntegerField(blank=True, null=True)),
                ('classification_minhit', models.IntegerField(blank=True, null=True)),
                ('success', models.BooleanField(default=False)),
                ('run', models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.CASCADE, to='result_display.runmain')),
                ('sample', models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.CASCADE, to='result_display.sample')),
            ],
            options={
                'ordering': ['method'],
            },
        ),
        migrations.CreateModel(
            name='QC_REPORT',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('report_source', models.CharField(blank=True, max_length=100, null=True)),
                ('QC_report', models.CharField(blank=True, max_length=100, null=True)),
                ('sample', models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.CASCADE, to='result_display.sample')),
            ],
            options={
                'ordering': ['sample'],
            },
        ),
        migrations.CreateModel(
            name='FinalReport',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('reference', models.CharField(blank=True, db_index=True, max_length=100, null=True)),
                ('unique_id', models.CharField(blank=True, max_length=100, null=True)),
                ('reference_length', models.IntegerField(blank=True, null=True)),
                ('taxid', models.CharField(blank=True, max_length=100, null=True)),
                ('simple_id', models.CharField(blank=True, max_length=100, null=True)),
                ('description', models.CharField(blank=True, max_length=100, null=True)),
                ('ref_db', models.CharField(blank=True, max_length=100, null=True)),
                ('reference_contig_str', models.CharField(blank=True, max_length=100, null=True)),
                ('accid', models.CharField(blank=True, max_length=100, null=True)),
                ('coverage', models.FloatField(blank=True, null=True)),
                ('depth', models.FloatField(blank=True, null=True)),
                ('depthR', models.FloatField(blank=True, null=True)),
                ('mapped_reads', models.IntegerField(blank=True, null=True)),
                ('ref_proportion', models.FloatField(blank=True, null=True)),
                ('mapped_proportion', models.FloatField(blank=True, null=True)),
                ('ngaps', models.IntegerField(blank=True, null=True)),
                ('success', models.CharField(blank=True, max_length=100, null=True)),
                ('refa_dotplot', models.TextField(blank=True, null=True)),
                ('refa_dotplot_exists', models.BooleanField(default=False)),
                ('covplot', models.TextField(blank=True, null=True)),
                ('covplot_exists', models.BooleanField(default=False)),
                ('bam_path', models.CharField(blank=True, max_length=100, null=True)),
                ('bai_path', models.CharField(blank=True, max_length=100, null=True)),
                ('reference_path', models.CharField(blank=True, max_length=100, null=True)),
                ('reference_index_path', models.CharField(blank=True, max_length=100, null=True)),
                ('reference_assembly_paf', models.CharField(blank=True, max_length=100, null=True)),
                ('mapped_scaffolds_path', models.CharField(blank=True, max_length=100, null=True)),
                ('mapped_scaffolds_index_path', models.CharField(blank=True, max_length=100, null=True)),
                ('run', models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.CASCADE, to='result_display.runmain')),
                ('sample', models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.CASCADE, to='result_display.sample')),
            ],
        ),
        migrations.CreateModel(
            name='ContigClassification',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('performed', models.BooleanField(default=False)),
                ('method', models.CharField(blank=True, max_length=10, null=True)),
                ('contig_classification_report', models.CharField(blank=True, max_length=200, null=True)),
                ('classification_number', models.IntegerField(blank=True, null=True)),
                ('classification_minhit', models.IntegerField(blank=True, null=True)),
                ('success', models.BooleanField(default=False)),
                ('run', models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.CASCADE, to='result_display.runmain')),
                ('sample', models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.CASCADE, to='result_display.sample')),
            ],
            options={
                'ordering': ['method'],
            },
        ),
    ]
