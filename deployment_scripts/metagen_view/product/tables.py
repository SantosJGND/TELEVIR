from typing import DefaultDict

import django_tables2 as tables
from result_display.models import RunMain, Sample


class SampleTable(tables.Table):
    class Meta:
        model = Sample

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
        "televir_sample_main",
        text="detail",
        args=[tables.A("project__name"), tables.A("name")],
    )


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
        "televir_sample_detail",
        text="Details",
        args=[tables.A("project"), tables.A("sample"), tables.A("name")],
    )

    def render_runtime(self, record):
        return float(record.runtime.split()[0])
