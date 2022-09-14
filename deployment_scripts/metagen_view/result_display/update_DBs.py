import os
from typing import Type

from django.contrib.auth.models import User
from django.core.files import File

from result_display.models import (
    QC_REPORT,
    ContigClassification,
    FinalReport,
    Projects,
    ReadClassification,
    ReferenceContigs,
    ReferenceMap_Main,
    RunAssembly,
    RunDetail,
    RunMain,
    RunRemapMain,
    Sample,
    SampleQC,
)
from result_display.samples import Reference_Map, RunMetaData, SampleMetaData


####################################################################################################################
####################################################################################################################
def Update_project(
    project_directory_path,
    ser: str = "admin",
):
    """Updates the project"""
    project_directory_path = os.path.dirname(project_directory_path)
    project_name = os.path.basename(project_directory_path)
    project_name_simple = project_name.replace(".", "_").replace(":", "_")
    try:
        user = User.objects.get(username=ser)
    except User.DoesNotExist:
        user = User.objects.get(username="admin")

    try:
        project = Projects.objects.get(
            name=project_name, project_type=Projects.EXTERNAL
        )
    except Projects.DoesNotExist:
        project = Projects(
            name=project_name,
            full_path=project_directory_path,
            project_type=Projects.INHOUSE,
            created_by=user,
        )
        project.save()


def Update_Sample(sample_class: Type[SampleMetaData]):
    """
    Update Sample class.

    :param sample_class:
    :return: None
    """
    if os.path.isdir(sample_class.sample_dir):

        try:
            Sample.objects.get(
                name_extended=sample_class.sample_name,
                project__name=sample_class.project_name,
            )
        except Sample.DoesNotExist:

            Update_sample(sample_class)

        sample = Sample.objects.get(
            name_extended=sample_class.sample_name,
            project__name=sample_class.project_name,
        )

        try:
            SampleQC.objects.get(sample=sample)
        except SampleQC.DoesNotExist:
            Update_sample_qc(sample_class)


def Update_sample(sample_class: Type[SampleMetaData]):
    """update sample_class.
    :param sample_class:
    :return: None
    """
    try:
        sample = Sample.objects.get(
            name_extended=sample_class.sample_name,
            project__name=sample_class.project_name,
        )
    except Sample.DoesNotExist:
        #
        project = Projects.objects.get(name=sample_class.project_name)

        sample = Sample(
            project=project,
            name_extended=sample_class.sample_name,
            name=os.path.splitext(sample_class.sample_name)[0],
            technology=sample_class.sample_technology,
            type=sample_class.sample_type,  # SE or PE
            combinations=len(sample_class.runs),
            input=sample_class.sample_input,  # input files,
            report="report",
        )
        sample.save()


def Update_sample_qc(sample_class: Type[SampleMetaData]):
    """update sample_class.qc_data.
    :param sample_class:
    :return: None
    """

    sample = Sample.objects.get(
        project__name=sample_class.project_name, name_extended=sample_class.sample_name
    )

    reads_before_processing = sample_class.qcdata["input"].loc["Total_Sequences"][
        "value"
    ]
    reads_after_processing = sample_class.qcdata["processed"].loc["Total_Sequences"][
        "value"
    ]

    percent_passed = (int(reads_after_processing) / int(reads_before_processing)) * 100

    input_report = open(sample_class.input_fastqc_report, "r")
    processed_report = open(sample_class.processed_fastqc_report, "r")

    try:

        sampleqc = SampleQC(
            sample=sample,
            software=sample_class.qc_soft,
            encoding=sample_class.sample_technology,
            input_reads=f"{int(reads_before_processing):,}",
            processed_reads=f"{int(reads_after_processing):,}",
            percent_passed=round(percent_passed, 2),
            sequence_length=sample_class.qcdata["processed"].loc["Sequence length"][
                "value"
            ],
            percent_gc=sample_class.qcdata["processed"].loc["%GC"]["value"],
            input_fastqc_report=File(
                input_report, name=os.path.basename(input_report.name)
            ),
            processed_fastqc_report=File(
                processed_report, name=os.path.basename(processed_report.name)
            ),
        )

        sampleqc.save()

    except:
        print("failed to input")
    finally:
        input_report.close()
        processed_report.close()


def Update_QC_report(sample_class: Type[SampleMetaData]):
    """
    Update QC data for sample_class.

    :param sample_class:
    :return: None
    """
    sample = Sample.objects.get(
        project__name=sample_class.project_name, name_extended=sample_class.sample_name
    )

    try:
        qc_report = QC_REPORT.objects.get(sample=sample, report_source=QC_REPORT.RAW)
    except QC_REPORT.DoesNotExist:
        qc_report = QC_REPORT(
            sample=sample,
            report_source=QC_REPORT.RAW,
            QC_report=sample_class.input_fastqc_report,
        )
        qc_report.save()

    try:
        qc_report = QC_REPORT.objects.get(
            sample=sample, report_source=QC_REPORT.PROCESSED
        )
    except QC_REPORT.DoesNotExist:
        qc_report = QC_REPORT(
            sample=sample,
            report_source=QC_REPORT.PROCESSED,
            QC_report=sample_class.processed_fastqc_report,
        )
        qc_report.save()


def Update_Sample_Runs(sample_class: Type[RunMetaData]):
    """get run data
    Update ALL run TABLES:
    - RunMain,
    - RunDetail,
    - RunAssembly,
    - ReadClassification,
    - ContigClassification,
    - RunRemapMain,
    - ReferenceMap_Main
    - ReferenceContigs
    - FinalReport,

    :param sample_class:
    :return: run_data
    """

    for run in sample_class.runs:

        Update_RunMain(run)
        Update_Sample_Runs_DB(run)
        Update_RefMap_DB(run)


def Update_RefMap_DB(run_class: Type[RunMetaData]):
    """
    Update Remap TABLES with info on this run.

    :param run_class:
    :return: run_data
    """

    for ref_map in run_class.references:

        Update_ReferenceMap(
            ref_map, run_class, sample_name=run_class.sample, name=run_class.name
        )


def Update_ReferenceMap(
    ref_map: Type[Reference_Map],
    run_class: Type[RunMetaData],
    sample_name="",
    name="",
):
    """
    Updates the reference map data to TABLES.
    - ReferenceMap_Main,
    - ReferenceContigs
    """
    sample = Sample.objects.get(
        name_extended=sample_name, project__name=run_class.project_name
    )
    run = RunMain.objects.get(
        project__name=run_class.project_name,
        suprun=run_class.suprun.name,
        name=name,
        sample=sample,
    )

    run_context = run_class.get_context()

    try:
        map_db = ReferenceMap_Main.objects.get(
            reference=ref_map.reference, sample=sample, run=run
        )
    except ReferenceMap_Main.DoesNotExist:
        map_db = ReferenceMap_Main(
            reference=ref_map.reference,
            sample=sample,
            run=run,
            taxid=ref_map.taxid,
            report=run_class.report,
            plotly_dotplot=run_context["remap"]["plotly_dotplots"],
            bam_file_path=ref_map.bam_path,
            bai_file_path=ref_map.bam_index_path,
            fasta_file_path=ref_map.reference_path,
            fai_file_path=ref_map.reference_index_path,
        )
        map_db.save()

    remap_stats = ref_map.remap_stats.set_index("ID")

    for seqid, row in remap_stats.iterrows():
        try:
            map_db_seq = ReferenceContigs.objects.get(
                reference=map_db, run=run, contig=seqid
            )
        except ReferenceContigs.DoesNotExist:
            map_db_seq = ReferenceContigs(
                contig=seqid,
                reference=map_db,
                run=run,
                depth=row["Hdepth"],
                depthr=row["HdepthR"],
                coverage=row["coverage"],
            )
            map_db_seq.save()

        map_db_seq.report = run_class.report
        map_db_seq.save()


def Update_RunMain(run_class: Type[RunMetaData]):
    """update run data for run_class. Update run_class.run_data.

    :param run_class:
    :return: None
    """
    project = Projects.objects.get(name=run_class.project_name)
    sample = Sample.objects.get(
        name_extended=run_class.sample, project__name=run_class.project_name
    )

    reads_after_processing = run_class.suprun.read_summary.loc["CLEAN"][0]
    reads_proc_percent = (
        reads_after_processing / run_class.suprun.read_summary.loc["INPUT"][0]
    ) * 100

    enrichment = run_class.suprun.conf["ENRICH"] == "true"
    host_depletion = run_class.suprun.conf["DEPLETE"] == "true"

    host_depletion_method = "None"
    enrichment_method = "None"
    if enrichment:
        enrichment_method = run_class.suprun.conf["HD"]

    if host_depletion:
        host_depletion_method = run_class.suprun.conf["HD"]

    try:
        runmain = RunMain.objects.get(
            project__name=run_class.project_name,
            suprun=run_class.suprun.name,
            sample=sample,
            name=run_class.name,
        )
    except RunMain.DoesNotExist:

        runmain = RunMain(
            suprun=run_class.suprun.name,
            project=project,
            sample=sample,
            name=run_class.name,
            params_file_path=run_class.params_file_path,
            processed_reads_r1=run_class.processed_reads_r1,
            processed_reads_r2=run_class.processed_reads_r2,
            assembly_performed=run_class.suprun.conf["ASSEMBLE"],
            assembly_method=run_class.suprun.conf["ASSEMBLY_SOFT"],
            reads_after_processing=f"{reads_after_processing:,}",
            reads_proc_percent=round(reads_proc_percent, 2),
            host_depletion=host_depletion_method,
            host_depletion_performed=host_depletion,
            enrichment_performed=enrichment,
            enrichment=enrichment_method,
            assembly_max=f"{run_class.suprun.assembly_summary['length'].max():,}",
            remap=run_class.conf["REMAP"],
            read_classification=run_class.conf["CLASSM"],
            contig_classification=run_class.suprun.conf["ASSEMBLE_CLASS"],
            runtime=f"{run_class.runtime / 60:.2f} m",
            finished=str(run_class.finished),
            report="report",
        )

        runmain.save()


def Update_Sample_Runs_DB(run_class: Type[RunMetaData]):
    """
    Update ALL run TABLES for one run_class.:
    - RunMain,
    - RunDetail,
    - RunAssembly,
    - ReadClassification,
    - ContigClassification,
    - RunRemapMain,
    - ReferenceMap_Main
    - ReferenceContigs
    - FinalReport,

    :param run_class:
    :return: run_data
    """
    sample = Sample.objects.get(
        project__name=run_class.project_name, name_extended=run_class.sample
    )

    run_context = run_class.context

    try:
        runmain = RunMain.objects.get(
            project__name=run_class.project_name,
            suprun=run_class.suprun.name,
            sample=sample,
            name=run_class.name,
        )
    except RunMain.DoesNotExist:
        return

    try:
        run_detail = RunDetail.objects.get(run=runmain, sample=sample)
    except RunDetail.DoesNotExist:

        run_detail = RunDetail(
            run=runmain,
            sample=sample,
            max_depth=run_class.run_detail_report.max_depth,  #
            max_depthR=run_class.run_detail_report.max_depthR,  #
            max_gaps=run_class.run_detail_report.max_gaps,  #
            max_prop=run_class.run_detail_report.max_prop,  #
            max_mapped=run_class.run_detail_report.max_mapped,  #
            input=run_class.run_detail_report.input,  #
            processed=run_class.run_detail_report.processed,  #
            processed_percent=run_class.run_detail_report.processed_percent,  #
            sift_preproc=run_class.run_detail_report.sift_preproc,  #
            sift_remap=run_class.run_detail_report.sift_remap,  #
            sift_removed_pprc=run_class.run_detail_report.sift_removed_pprc,
            processing_final=run_class.run_detail_report.processing_final,  #
            processing_final_percent=run_class.run_detail_report.processing_final_percent,  #
            merged=run_class.run_detail_report.merged,  #
            merged_number=run_class.run_detail_report.merged_number,  #
            merged_files=run_class.run_detail_report.merged_files,  #
        )

        run_detail.save()

    try:
        run_assembly = RunAssembly.objects.get(run=runmain, sample=sample)
    except RunAssembly.DoesNotExist:

        run_assembly = RunAssembly(
            run=runmain,
            sample=sample,
            performed=run_class.assembly_report.performed,
            method=run_class.assembly_report.assembly_soft,
            contig_number=run_class.assembly_report.assembly_number,
            contig_max=run_class.assembly_report.assembly_max,
            contig_min=run_class.assembly_report.assembly_min,
            contig_mean=run_class.assembly_report.assembly_mean,
            contig_trim=run_class.assembly_report.assembly_trim,
            assembly_contigs=run_class.assembly_contigs,
        )
        run_assembly.save()

    try:
        read_classification = ReadClassification.objects.get(run=runmain, sample=sample)
    except ReadClassification.DoesNotExist:

        read_classification = ReadClassification(
            run=runmain,
            sample=sample,
            read_classification_report=run_class.read_classification_path,
            performed=run_class.read_classification_results.performed,
            method=run_class.read_classification_results.method,
            classification_number=run_class.read_classification_results.classification_number,
            classification_minhit=run_class.read_classification_results.classification_minhit,
            success=run_class.read_classification_results.success,
        )
        read_classification.save()

    try:
        contig_classification = ContigClassification.objects.get(
            run=runmain, sample=sample
        )
    except ContigClassification.DoesNotExist:

        contig_classification = ContigClassification(
            run=runmain,
            sample=sample,
            contig_classification_report=run_class.assembly_classification_path,
            performed=run_context["class_contigs"]["performed"],
            method=run_context["class_contigs"]["assembly_class_soft"],
            classification_number=run_context["class_contigs"]["assembly_class_number"],
            classification_minhit=run_context["class_contigs"]["assembly_min_hit"],
            success=run_context["class_contigs"]["success"],
        )
        contig_classification.save()

    try:
        remap_main = RunRemapMain.objects.get(run=runmain, sample=sample)
    except RunRemapMain.DoesNotExist:
        remap_main = RunRemapMain(
            run=runmain,
            sample=sample,
            merged_log=run_class.merged_log_path,
            remap_plan=run_class.remap_plan_path,
            performed=run_context["remap"]["performed"],
            method=run_context["remap"]["remap_soft"],
            found_total=run_context["remap"]["found_total"],
            coverage_maximum=run_context["remap"]["cov_max"],
            coverage_minimum=run_context["remap"]["cov_min"],
            success=run_context["remap"]["success"],
        )
        remap_main.save()

    for i, row in run_context["final_data"].iterrows():
        try:
            report_row = FinalReport.objects.get(
                run=runmain,
                sample=sample,
                unique_id=row["unique_id"],
            )
        except FinalReport.DoesNotExist:
            report_row = FinalReport(
                run=runmain,
                sample=sample,
                unique_id=row["unique_id"],
                reference_length=row["contig_length"],
                taxid=row["taxid"],
                accid=row["ID"],
                reference_contig_str=row["contig_string"],
                simple_id=row["simple_id"],
                description=row["description"],
                ref_db=row["refdb"],
                coverage=row["coverage"],
                depth=row["Hdepth"],
                depthR=row["HdepthR"],
                ngaps=row["ngaps"],
                mapped_reads=row["mapped"],
                ref_proportion=row["ref_prop"],
                mapped_proportion=row["mapped_prop"],
                success=row["success"],
                refa_dotplot_exists=row["refa_dotplot_exists"],
                covplot_exists=row["covplot_exists"],
                refa_dotplot=row["refa_dotplot_path"],
                covplot=row["covplot"],
                bam_path=row["bam_path"],
                bai_path=row["bai_path"],
                reference_path=row["reference_path"],
                reference_index_path=row["reference_index_path"],
                reference_assembly_paf=row["reference_assembly_paf"],
                mapped_scaffolds_path=row["mapped_scaffolds_path"],
                mapped_scaffolds_index_path=row["mapped_scaffolds_index_path"],
            )

            report_row.save()
