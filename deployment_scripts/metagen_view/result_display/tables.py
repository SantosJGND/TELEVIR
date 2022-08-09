from typing import DefaultDict

import django_tables2 as tables

from result_display.models import Projects, ReferenceContigs, RunMain, Sample, SampleQC


class SampleTable(tables.Table):
    class Meta:
        model = Sample
        # attrs = {
        #    "class": "semantic-ui-table",
        # }
        # template_name = "django_tables2/bootstrap4.html"
        attrs = {"class": "paleblue"}
        fields = (
            "name",
            "combinations",
            "report",
            "input",
            "technology",
            "type",
        )

    def render_combinations(self, record):
        return RunMain.objects.filter(
            sample__name=record.name, project__name=record.project
        ).count()

    report = tables.LinkColumn(
        "sample_main", text="Report", args=[tables.A("project__name"), tables.A("name")]
    )


class SampleQCTable(tables.Table):
    class Meta:
        model = SampleQC
        attrs = {
            "class": "paleblue",
        }
        fields = (
            "encoding",
            "input_reads",
            "processed_reads",
            "sequence_length",
            "percent_gc",
        )


class ContigTable(tables.Table):
    class Meta:
        model = ReferenceContigs
        attrs = {
            "class": "paleblue",
        }
        fields = ("contig", "depth", "depthr", "coverage")


class RunMainTable(tables.Table):
    class Meta:
        model = RunMain
        attrs = {
            "class": "paleblue",
        }
        fields = (
            "name",
            "enrichment",
            "host_depletion",
            "assembly_method",
            "read_classification",
            "contig_classification",
            "finished",
            "runtime",
        )

    report = tables.LinkColumn(
        "sample_detail",
        text="Details",
        args=[tables.A("project"), tables.A("sample"), tables.A("name")],
    )

    def render_runtime(self, record):
        return float(record.runtime.split()[0])
