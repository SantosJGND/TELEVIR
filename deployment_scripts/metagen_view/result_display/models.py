import codecs
import os
from turtle import back

from django.contrib.auth.models import User
from django.db import models
from django.utils.safestring import mark_safe

# Create your models here.


class Projects(models.Model):

    INHOUSE = 1
    EXTERNAL = 2

    name = models.CharField(max_length=100, blank=True, db_index=True, null=True)
    full_path = models.CharField(max_length=200, blank=True, db_index=True, null=True)
    created_by = models.ForeignKey(
        User, on_delete=models.CASCADE, blank=True, null=True
    )

    project_type = models.IntegerField(
        choices=((INHOUSE, "In-House"), (EXTERNAL, "External")), default=INHOUSE
    )
    date_created = models.DateTimeField(auto_now_add=True)
    date_modified = models.DateTimeField(auto_now_add=True)

    submit_index = models.IntegerField(default=0)

    class Meta:
        ordering = ["name"]

    def __str__(self):
        return self.name


class Sample(models.Model):
    """
    Main sample information. Connects to the RunMain and QC models.
    """

    project = models.ForeignKey(
        Projects,
        null=True,
        blank=True,
        on_delete=models.CASCADE,
    )

    name = models.CharField(
        max_length=100, db_index=True, blank=True, null=True
    )  # Create your models here. # Name of the sample
    name_extended = models.CharField(
        max_length=100, db_index=True, blank=True, null=True
    )  ## extra name to show in the settings HTML table

    type = models.CharField(
        max_length=10, blank=True, null=True
    )  # sample type: SE or PE
    combinations = models.IntegerField(blank=True, null=True)  # number of combinations
    input = models.TextField(blank=True, null=True)  # input files
    technology = models.CharField(
        max_length=100,
        name="technology",
        blank=True,
        null=True,
    )  # encoding

    report = models.CharField(max_length=100, blank=True, null=True)  # report file

    class Meta:
        ordering = [
            "name",
        ]

    def __str__(self):
        return self.name


class SampleQC(models.Model):
    """Results of sample quality control."""

    # project = models.ForeignKey(
    #    Projects,
    #    null=True,
    #    blank=True,
    #    on_delete=models.CASCADE,
    # )

    sample = models.ForeignKey(
        Sample, blank=True, null=True, on_delete=models.CASCADE
    )  ## sample
    software = models.CharField(max_length=100, blank=True, null=True)  # software used

    qc_type = models.CharField(
        max_length=100, name="qc_type", blank=True, null=True
    )  # qc type
    encoding = models.CharField(
        max_length=100, name="encoding", blank=True, null=True
    )  # encoding
    input_reads = models.CharField(
        max_length=100, name="input_reads", blank=True, null=True
    )  # input reads

    processed_reads = models.CharField(
        max_length=100, name="processed_reads", blank=True, null=True
    )  # processed reads

    percent_passed = models.FloatField(
        name="percent_passed", blank=True, null=True
    )  # percent passed

    sequence_length = models.CharField(
        max_length=100, blank=True, null=True
    )  # Read length distribution after processing
    percent_gc = models.FloatField(
        name="percent_gc", blank=True, null=True
    )  # percent GC content after filtering.
    # report = tables.LinkColumn("report", text="Report", args=["pk"])
    input_fastqc_report = models.FileField(
        upload_to="input_fastqc_report", blank=True, null=True
    )  # input fastqc report

    processed_fastqc_report = models.FileField(
        upload_to="processed_fastqc_report", blank=True, null=True
    )  # processed fastqc report

    def render_input_fastqc(self):
        html_path = self.input_fastqc_report.path
        print(html_path)

        if os.path.exists(html_path):
            print("file exists")
            html_string = codecs.open(html_path, "r").read()
            return mark_safe(html_string)
        return None

    def render_processed_fastqc(self):
        html_path = self.processed_fastqc_report.path
        if os.path.exists(html_path):
            html_string = codecs.open(html_path, "r").read()
            return mark_safe(html_string)
        return None

    class Meta:
        ordering = [
            "sample",
        ]

    def __str__(self):
        return self.name


class QC_REPORT(models.Model):

    RAW = "input"
    PROCESSED = "processed"

    sample = models.ForeignKey(
        Sample, blank=True, null=True, on_delete=models.CASCADE
    )  ## sample

    report_source = models.CharField(max_length=200, blank=True, null=True)  # qc type

    QC_report = models.CharField(max_length=250, blank=True, null=True)  # qc type

    class Meta:
        ordering = [
            "sample",
        ]

    def __str__(self):
        return self.QC_report


class RunIndex(models.Model):
    project = models.ForeignKey(
        Projects,
        null=True,
        blank=True,
        on_delete=models.CASCADE,
    )
    sample = models.ForeignKey(Sample, blank=True, null=True, on_delete=models.CASCADE)
    name = models.CharField(
        max_length=100, db_index=True, blank=True, null=True
    )  # Create your models here.


class RunMain(models.Model):

    project = models.ForeignKey(
        Projects,
        null=True,
        blank=True,
        on_delete=models.CASCADE,
    )

    suprun = models.CharField(max_length=100, blank=True, null=True)
    name = models.CharField(
        max_length=100, db_index=True, blank=True, null=True
    )  # Create your models here.

    sample = models.ForeignKey(Sample, blank=True, null=True, on_delete=models.CASCADE)

    params_file_path = models.CharField(max_length=250, blank=True, null=True)

    processed_reads_r1 = models.CharField(
        max_length=200, blank=True, null=True
    )  # processed reads
    processed_reads_r2 = models.CharField(
        max_length=200, blank=True, null=True
    )  # processed reads

    enrichment = models.CharField(
        max_length=20, blank=True, null=True
    )  # enrichment method if any
    enrichment_performed = models.BooleanField(
        blank=True, null=True
    )  # enrichment performed
    enrichment_args = models.CharField(
        max_length=50, blank=True, null=True
    )  # enrichment args

    enrichment_db = models.CharField(
        max_length=200, blank=True, null=True
    )  # enrichment db if any

    host_depletion = models.CharField(
        max_length=20, blank=True, null=True
    )  # host depletion method if any
    host_depletion_performed = models.BooleanField(
        blank=True, null=True
    )  # host depletion performed

    host_depletion_args = models.CharField(
        max_length=50, blank=True, null=True
    )  # enrichment args

    host_depletion_db = models.CharField(
        max_length=200, blank=True, null=True
    )  # enrichment db if any

    reads_after_processing = models.CharField(
        max_length=100, blank=True, null=True
    )  # reads after processing
    reads_proc_percent = models.CharField(
        max_length=100, blank=True, null=True
    )  # percent of reads after processing

    assembly_performed = models.CharField(
        max_length=10, blank=True, null=True
    )  # assembly method if any
    assembly_method = models.CharField(
        max_length=50, blank=True, null=True
    )  # assembly method if any

    assembly_max = models.CharField(
        max_length=100, blank=True, null=True
    )  # max length of contig.
    read_classification = models.CharField(
        max_length=50, blank=True, null=True
    )  # read classification method if any
    read_classification_performed = models.BooleanField(
        blank=True, null=True
    )  # read classification performed

    contig_classification = models.CharField(max_length=50, blank=True, null=True)
    contig_classification_performed = models.BooleanField(blank=True, null=True)

    remap = models.CharField(
        max_length=50, blank=True, null=True
    )  # remap method if any
    remap_performed = models.BooleanField(blank=True, null=True)
    remap_args = models.CharField(max_length=50, blank=True, null=True)

    finished = models.CharField(max_length=10, blank=True, null=True)  # SE or PE
    runtime = models.CharField(max_length=100, blank=True, null=True)

    report = models.CharField(max_length=200, blank=True, null=True)

    static_dir = models.CharField(max_length=250, blank=True, null=True)

    class Meta:

        ordering = [
            "name",
        ]

    def __str__(self):
        return self.name


class RunDetail(models.Model):

    run = models.ForeignKey(RunMain, blank=True, null=True, on_delete=models.CASCADE)

    name = models.CharField(
        max_length=100, db_index=True, blank=True, null=True
    )  # Create your models here.
    sample = models.ForeignKey(Sample, blank=True, null=True, on_delete=models.CASCADE)
    run = models.ForeignKey(RunMain, blank=True, null=True, on_delete=models.CASCADE)

    max_depth = models.FloatField(blank=True, null=True)
    max_depthR = models.FloatField(blank=True, null=True)
    max_gaps = models.IntegerField(blank=True, null=True)
    max_prop = models.FloatField(blank=True, null=True)
    max_mapped = models.IntegerField(blank=True, null=True)

    input = models.CharField(max_length=300, blank=True, null=True)
    processed = models.CharField(max_length=300, blank=True, null=True)
    processed_percent = models.FloatField(blank=True, null=True)
    sift_preproc = models.BooleanField(blank=True, null=True)
    sift_remap = models.BooleanField(blank=True, null=True)

    sift_removed_pprc = models.CharField(max_length=300, blank=True, null=True)
    processing_final = models.CharField(max_length=300, blank=True, null=True)
    processing_final_percent = models.FloatField(blank=True, null=True)
    merged = models.BooleanField(blank=True, null=True)
    merged_number = models.IntegerField(blank=True, null=True)
    merged_files = models.CharField(max_length=300, blank=True, null=True)

    class Meta:
        ordering = [
            "name",
        ]

    def __str__(self):
        return self.name


class RunAssembly(models.Model):

    run = models.ForeignKey(RunMain, blank=True, null=True, on_delete=models.CASCADE)
    sample = models.ForeignKey(Sample, blank=True, null=True, on_delete=models.CASCADE)

    performed = models.BooleanField(default=False)

    assembly_contigs = models.CharField(max_length=200, blank=True, null=True)

    method = models.CharField(
        max_length=20, blank=True, null=True
    )  # assembly method if any

    args = models.CharField(max_length=50, blank=True, null=True)  # assembly args

    contig_number = models.IntegerField(blank=True, null=True)

    contig_max = models.CharField(
        max_length=100, blank=True, null=True
    )  # max length of contig.
    contig_min = models.CharField(
        max_length=100, blank=True, null=True
    )  # min length of contig.
    contig_mean = models.CharField(
        max_length=100, blank=True, null=True
    )  # mean length of contig.

    contig_trim = models.CharField(
        max_length=100, blank=True, null=True
    )  # contig min len filter

    class Meta:
        ordering = [
            "method",
        ]

    def __str__(self):
        return self.method


class ReadClassification(models.Model):

    run = models.ForeignKey(RunMain, blank=True, null=True, on_delete=models.CASCADE)
    sample = models.ForeignKey(Sample, blank=True, null=True, on_delete=models.CASCADE)

    performed = models.BooleanField(default=False)

    method = models.CharField(
        max_length=20, blank=True, null=True
    )  # read classification method if any

    args = models.CharField(
        max_length=250, blank=True, null=True
    )  # read classification args
    db = models.CharField(
        max_length=250, blank=True, null=True
    )  # read classification db if any

    read_classification_report = models.CharField(
        max_length=250, blank=True, null=True
    )  # read classification report

    classification_number = models.IntegerField(blank=True, null=True)
    classification_minhit = models.IntegerField(blank=True, null=True)

    success = models.BooleanField(default=False)

    class Meta:
        ordering = [
            "method",
        ]

    def __str__(self):
        return self.method


class ContigClassification(models.Model):

    run = models.ForeignKey(RunMain, blank=True, null=True, on_delete=models.CASCADE)
    sample = models.ForeignKey(Sample, blank=True, null=True, on_delete=models.CASCADE)

    performed = models.BooleanField(default=False)

    method = models.CharField(
        max_length=250, blank=True, null=True
    )  # contig classification method if any

    args = models.CharField(
        max_length=250, blank=True, null=True
    )  # read classification args
    db = models.CharField(
        max_length=250, blank=True, null=True
    )  # read classification db if any

    contig_classification_report = models.CharField(
        max_length=300, blank=True, null=True
    )  # contig classification report

    classification_number = models.IntegerField(blank=True, null=True)
    classification_minhit = models.IntegerField(blank=True, null=True)

    success = models.BooleanField(default=False)

    class Meta:
        ordering = [
            "method",
        ]

    def __str__(self):
        return self.method


class RunRemapMain(models.Model):

    run = models.ForeignKey(RunMain, blank=True, null=True, on_delete=models.CASCADE)
    sample = models.ForeignKey(Sample, blank=True, null=True, on_delete=models.CASCADE)
    merged_log = models.CharField(max_length=200, blank=True, null=True)
    performed = models.BooleanField(default=False)

    method = models.CharField(
        max_length=250, blank=True, null=True
    )  # remap method if any

    found_total = models.IntegerField(blank=True, null=True)

    coverage_minimum = models.IntegerField(blank=True, null=True)
    coverage_maximum = models.IntegerField(blank=True, null=True)

    success = models.IntegerField(blank=True, null=True)
    remap_plan = models.CharField(max_length=200, blank=True, null=True)

    class Meta:
        ordering = [
            "method",
        ]

    def __str__(self):
        return self.method


class ReferenceMap_Main(models.Model):

    reference = models.CharField(
        max_length=100, db_index=True, blank=True, null=True
    )  # Create your models here.

    run = models.ForeignKey(RunMain, blank=True, null=True, on_delete=models.CASCADE)

    sample = models.ForeignKey(Sample, blank=True, null=True, on_delete=models.CASCADE)
    taxid = models.CharField(max_length=20, blank=True, null=True)
    reference_contig_str = models.CharField(max_length=100, blank=True, null=True)

    report = models.CharField(max_length=200, blank=True, null=True)

    plotly_dotplot = models.TextField(blank=True, null=True)

    bam_file_path = models.CharField(max_length=250, blank=True, null=True)
    bai_file_path = models.CharField(max_length=250, blank=True, null=True)
    fasta_file_path = models.CharField(max_length=250, blank=True, null=True)
    fai_file_path = models.CharField(max_length=250, blank=True, null=True)

    class Meta:

        ordering = [
            "reference",
        ]

    def __str__(self):
        return self.reference


class FinalReport(models.Model):

    run = models.ForeignKey(RunMain, blank=True, null=True, on_delete=models.CASCADE)
    sample = models.ForeignKey(Sample, blank=True, null=True, on_delete=models.CASCADE)
    reference = models.CharField(
        max_length=100, db_index=True, blank=True, null=True
    )  # Create your models here.

    unique_id = models.CharField(max_length=20, blank=True, null=True)

    reference_length = models.IntegerField(blank=True, null=True)

    taxid = models.CharField(max_length=20, blank=True, null=True)
    simple_id = models.CharField(max_length=20, blank=True, null=True)
    description = models.CharField(max_length=150, blank=True, null=True)

    ref_db = models.CharField(max_length=400, blank=True, null=True)
    reference_contig_str = models.CharField(max_length=100, blank=True, null=True)

    accid = models.CharField(max_length=20, blank=True, null=True)
    coverage = models.FloatField(blank=True, null=True)
    depth = models.FloatField(blank=True, null=True)
    depthR = models.FloatField(blank=True, null=True)
    mapped_reads = models.IntegerField(blank=True, null=True)
    ref_proportion = models.FloatField(blank=True, null=True)
    mapped_proportion = models.FloatField(blank=True, null=True)
    ngaps = models.IntegerField(blank=True, null=True)
    mapping_success = models.CharField(max_length=20, blank=True, null=True)
    classification_success = models.CharField(max_length=20, blank=True, null=True)
    windows_covered = models.CharField(max_length=20, blank=True, null=True)

    refa_dotplot = models.CharField(max_length=350, blank=True, null=True)
    refa_dotplot_exists = models.BooleanField(default=False)
    covplot = models.CharField(max_length=350, blank=True, null=True)
    covplot_exists = models.BooleanField(default=False)
    bam_path = models.CharField(max_length=350, blank=True, null=True)
    bai_path = models.CharField(max_length=350, blank=True, null=True)
    reference_path = models.CharField(max_length=350, blank=True, null=True)
    reference_index_path = models.CharField(max_length=350, blank=True, null=True)
    reference_assembly_paf = models.CharField(max_length=350, blank=True, null=True)
    mapped_scaffolds_path = models.CharField(max_length=350, blank=True, null=True)
    mapped_scaffolds_index_path = models.CharField(
        max_length=350, blank=True, null=True
    )


class ReferenceContigs(models.Model):

    reference = models.ForeignKey(
        ReferenceMap_Main, blank=True, null=True, on_delete=models.CASCADE
    )
    run = models.ForeignKey(RunMain, blank=True, null=True, on_delete=models.CASCADE)
    contig = models.CharField(max_length=100, blank=True, null=True)
    # length = models.CharField(max_length=100, blank=True, null=True)
    depth = models.CharField(max_length=100, blank=True, null=True)
    depthr = models.CharField(max_length=100, blank=True, null=True)
    coverage = models.CharField(max_length=100, blank=True, null=True)

    class Meta:

        ordering = [
            "reference",
        ]

    def __str__(self):
        return self.contig
