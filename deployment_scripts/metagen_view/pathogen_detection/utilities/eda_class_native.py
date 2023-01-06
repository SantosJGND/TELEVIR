import itertools as it
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


class Validator:

    validation: pd.DataFrame

    def __init__(self, filepath):

        self.validation_set = self.load_validation(filepath)

    @staticmethod
    def check_content(df):
        checksum = 0
        if "taxid" not in df.columns:
            df["taxid"] = np.nan
            checksum += 1

        if "accid" not in df.columns:
            df["accid"] = np.nan
            checksum += 1

        if "description" not in df.columns:
            df["description"] = np.nan
            checksum += 1

        if checksum == 3:
            print(
                "all columns absent. provide at least one validator: taxid, accid or description."
            )

        return df

    def load_validation(self, file_path):

        df = pd.read_csv(file_path, sep="\t")

        df = self.check_content(df)

        df = df.dropna()
        df["taxid"] = df.taxid.apply(lambda x: [y for y in x.split(";")])
        df["accid"] = df.accid.apply(lambda x: [y for y in x.split(";")])

        df.set_index("sample_name", inplace=True)

        return df

    def assess(self, x):

        if x.taxid in self.validation_set.loc[x.sample_name].taxid:
            return True
        if x.accid in self.validation_set.loc[x.sample_name].accid:
            return True

        if self.validation_set.loc[x.sample_name].description[0] in x.description:
            return True

        return False


class run_eda:

    data_total: pd.DataFrame
    softs: pd.DataFrame
    source_total: dict
    sources: list

    def __init__(self, validator: Validator, report_file: str, params_file: str):
        self.validator = validator
        self.report_file = report_file
        self.params_file = params_file
        self.dir = os.path.dirname(report_file)

    def run_input(self) -> None:
        """
        return report data sets
        """

        self.data_total = pd.read_csv(self.report_file, sep="\t")
        self.softs = pd.read_csv(self.params_file, sep="\t")

    def remove_samples_unvalidated(self):
        absent_samples_idx = ~self.data_total.sample_name.isin(
            self.validator.validation_set.index
        )
        absent_samples = self.data_total.sample_name[absent_samples_idx].unique()
        absent_runs = self.data_total.run_id[absent_samples_idx].unique()

        self.data_total = self.data_total[
            ~self.data_total.sample_name.isin(absent_samples)
        ]

        self.softs = self.softs[~self.softs.run_id.isin(absent_runs)]

    def complement_input(self):

        self.data_total["runid"] = self.data_total.apply(
            lambda x: f"{x.project_id}_{x.sample_id}_{x.run_id}", axis=1
        )

        self.softs["runid"] = self.softs.apply(
            lambda x: f"{x.project_id}_{x.sample_id}_{x.run_id}", axis=1
        )

        self.data_total["source"] = self.data_total.apply(
            lambda x: f"{x.project_id}_{x.sample_id}", axis=1
        )

        self.softs["source"] = self.softs.apply(
            lambda x: f"{x.project_id}_{x.sample_id}", axis=1
        )

    def split(self):

        self.source_total = {
            source: self.data_total[self.data_total.source == source].reset_index(
                drop=True
            )
            for source in self.data_total.source.unique()
        }

        for source in self.data_total.source.unique():
            self.source_total[source] = self.data_total[
                self.data_total.source == source
            ].reset_index(drop=True)

        self.sources = list(self.source_total.keys())

        self.args_dict = {
            runid: self.softs[self.softs.runid == runid].reset_index(drop=True)
            for runid in self.softs.runid.unique()
        }

        return self

    def run_summary(self, stotal):
        run_summaries = []

        for run in stotal.runid.unique():
            source = stotal.source.unique()[0]
            rda = stotal[stotal.runid == run].reset_index(drop=True)
            found_success = np.nansum(rda.success == True) > 0
            runtime = max(rda.runtime)
            percent_over = 0
            Hdepth = 0
            HdepthR = 0
            mapped_reads = 0
            ref_proportion = 0
            mapped_proportion = 0
            ngaps = 0
            assh = False
            dsuc = rda[rda.success == True]
            focus = 0

            if len(dsuc):
                focus = (
                    dsuc[dsuc.classification_success == True].shape[0] / rda.shape[0]
                )
                percent_over = np.mean(dsuc["coverage"])
                Hdepth = np.mean(dsuc["depth"])
                HdepthR = np.mean(dsuc["depthR"])
                mapped_reads = np.mean(dsuc["mapped_reads"])
                ref_proportion = np.mean(dsuc["ref_proportion"])
                mapped_proportion = np.mean(dsuc["mapped_proportion"])
                ngaps = np.mean(dsuc["ngaps"])
                assh = sum(dsuc.classification_success.str.contains("contigs")) > 0
                readh = sum(dsuc.classification_success.str.contains("reads")) > 0

            run_summaries.append(
                [
                    rda.sample_name.unique()[0],
                    source,
                    run,
                    found_success,
                    runtime,
                    percent_over,
                    Hdepth,
                    HdepthR,
                    mapped_reads,
                    ref_proportion,
                    mapped_proportion,
                    ngaps,
                    readh,
                    assh,
                    focus,
                ]
            )

        run_summaries = pd.DataFrame(
            run_summaries,
            columns=[
                "sample_name",
                "source",
                "runid",
                "success",
                "runtime",
                "coverage",
                "depth",
                "depthR",
                "mapped_reads",
                "ref_proportion",
                "mapped_proportion",
                "ngaps",
                "rhelp",
                "ahelp",
                "focus",
            ],
        )
        return run_summaries

    def summarize_runs(self):
        if os.path.isfile(f"run_summaries.tsv"):
            pass
            self.summaries = pd.read_csv(f"run_summaries.tsv", sep="\t")

        else:
            # self.split()
            self.summaries = {}
            self.descriptions = {}
            print(self.sources)
            for file in self.sources:

                stotal = self.source_total[file]

                stotal["success"] = stotal.apply(self.validator.assess, axis=1)
                run_summaries = self.run_summary(stotal)
                self.summaries[file] = run_summaries

            self.summaries = pd.concat(list(self.summaries.values()), axis=0)
            self.summaries.to_csv(f"run_summaries.tsv", sep="\t")

    def summarize_samples(self):
        if os.path.isfile(f"run_summaries.tsv"):
            pass
            self.summaries = pd.read_csv(f"run_summaries.tsv", sep="\t")

        else:
            # self.split()
            self.summaries = {}
            self.descriptions = {}
            print(self.sources)
            for file in self.sources:

                stotal = self.source_total[file]

                stotal["success"] = stotal.apply(self.validator.assess, axis=1)
                run_summaries = self.run_summary(stotal)
                self.summaries[file] = run_summaries

            self.summaries = pd.concat(list(self.summaries.values()), axis=0)
            self.summaries.to_csv(f"run_summaries.tsv", sep="\t")

    def describe(self, file=""):
        self.descriptions = {}

        for file in self.sources:
            stotal = self.source_total[file]
            sample = stotal.sample_name.unique()[0]
            infection = self.validator.validation_set.loc[sample][0]
            run_summaries = self.summaries[self.summaries.sample_name == sample]
            ssums = run_summaries[run_summaries.success == True]
            ###
            psc_runs = sum(run_summaries.success == True) / run_summaries.shape[0]
            psc_runs = round(psc_runs * 100, 3)
            report_runs = f"{run_summaries.shape[0]}"  # f"total number of runs: {run_summaries.shape[0]}"
            report_success = f"{psc_runs} %,  ({ssums.shape[0]})"  # f"N runs that find {infection} : {ssums.shape[0]},  ({psc_runs} %)"
            ##
            scruns = stotal[stotal.description.str.contains(infection, na=False)]
            species_unique = "\n".join(
                [x + ";" for x in scruns.description.unique()]
            )  # "species found with this description: " + ", ".join(list(scruns.description.unique()))
            ##
            allcover = (ssums.success == True) & (ssums.ahelp == True)
            allcover = sum(allcover)
            if ssums.shape[0]:
                succass = round(100 * allcover / ssums.shape[0], 3)
            else:
                succass = 0

            succtotal = round(100 * allcover / run_summaries.shape[0], 3)

            report1 = f"{succass} %"  # f"percent of successfull runs where assembly is also successful : {succass} %"
            report2 = f"{succtotal} %"  # f"percent of runs that find {infection} with both reads and assembly : {succtotal} %"

            self.descriptions[file] = [
                infection,
                report_runs,
                report_success,
                species_unique,
                report1,
                report2,
            ]

    def pretty_print(self, file):
        text = self.descriptions[file]
        return text

    def describe_print(self, file=""):
        text = [
            [
                "infection",
                "Number of runs",
                "success",
                "species",
                "% success assembled",
                "% success + assembly",
            ]
        ]

        if file:
            header = ["summary \ source"] + [file]
            text.append(self.pretty_print(file))
        else:
            text.extend([self.pretty_print(file) for file in self.sources])
            header = ["summary \ source"] + self.sources

        text = pd.DataFrame(text).T
        text.columns = header
        return text

    def eda_plots(self, file, cols=["Hdepth%", "coverage", "runtime"]):

        run_summaries = self.summaries[file]
        ssums = run_summaries[run_summaries.success == True]

        for col in cols:
            sns.histplot(data=ssums, x=col).set_title(file)
            plt.show()

    def combine_data(self):
        if os.path.isfile(f"combined_reports.tsv"):
            combdat = pd.read_csv("combined_reports.tsv", sep="\t")

        else:
            total_softs = self.softs
            total_reports = self.summaries

            # for s, g in total_reports.items():
            #    g["source"] = s

            combdat = (
                pd.merge(total_softs, total_reports, on="runid")
                .rename({"source_x": "source"}, axis=1)
                .drop("source_y", axis=1)
            )

            combdat["complete"] = (combdat["success"] == True) & (
                combdat["ahelp"] == True
            ) & combdat["rhelp"] == True

            combdat.to_csv("combined_reports.tsv", sep="\t", header=True, index=False)

        self.combdat = combdat

    def combine_data_full(self):
        if os.path.isfile(f"reports.comb_full.tsv"):
            combdat_full = pd.read_csv(f"reports.comb_full.tsv", sep="\t")

        else:
            args_all = []
            for x, g in self.args_dict.items():
                newg = g.copy()
                newg["id"] = [x]
                newg.reset_index(inplace=True, drop=True)

                if len(args_all):
                    args_all = pd.concat([args_all, newg], axis=0)
                    args_all.reset_index(inplace=True, drop=True)
                else:
                    args_all = newg

            combdat_full = pd.merge(args_all, self.summaries, on="runid")
            combdat_full["complete"] = (
                (combdat_full["success"] == True)
                & (combdat_full["ahelp"] == True)
                & (combdat_full["rhelp"] == True)
            )

            combdat_full.to_csv(
                f"reports.comb_full.tsv", sep="\t", header=True, index=False
            )

        self.combdat_full = combdat_full

    def process(self, df, filter=False):
        dt = []
        df["runtime"] = df["runtime"].apply(lambda x: x.split(" ")[0])
        df["runtime"] = df["runtime"].astype(float)
        print(df.columns)

        for source in df.source.unique():
            sub = df[df.source == source].reset_index(drop=True)

            sub = sub.replace("NA", np.NaN)
            for cols in ["runtime", "coverage", "depth", "depthR"]:
                # print(sub[cols])
                if np.nanstd(sub[cols]) == 0:
                    sub[cols] = np.NaN
                else:
                    sub[cols] = (sub[cols] - np.nanmean(sub[cols])) / np.nanstd(
                        sub[cols]
                    )

                if filter:
                    sub = sub[(sub[cols] < 4) & (sub[cols] > -4)]

            dt.append(sub)

        dt = pd.concat(dt, axis=0)
        dt_cols = list(dt.columns)

        dt.columns = dt_cols
        return dt

    def runid_summary(self, combdat_full_process):
        run_assess = []
        for rid in combdat_full_process.source.unique():
            rclub = combdat_full_process[
                combdat_full_process.source == rid
            ].reset_index(drop=True)

            csuc = sum(rclub.success) / rclub.shape[0]
            cass = sum(rclub.ahelp) / rclub.shape[0]
            ccov = np.median(rclub["coverage"])
            cdepth = np.median(rclub["depth"])
            cdepthr = np.median(rclub["depthR"])
            cruntime = np.median(rclub[~rclub.runtime.isnull()]["runtime"])
            crsucc = sum(~rclub.runtime.isnull()) / rclub.shape[0]

            runsum = pd.DataFrame(
                [
                    [
                        rid,
                        csuc,
                        cass,
                        ccov,
                        cdepth,
                        cdepthr,
                        cruntime,
                        crsucc,
                        sum(rclub.complete) / rclub.shape[0],
                    ]
                ],
                columns=[
                    "source",
                    "success",
                    "ahelp",
                    "coverage",
                    "depth",
                    "depthR",
                    "runtime",
                    "finished",
                    "complete",
                ],
            )

            run_assess.append(runsum)

        run_assess = pd.concat(run_assess).reset_index(drop=True)
        return run_assess

    def interactive_1d(self, cols=["Hdepth%", "%>2", "runtime"]):
        @interact(file=self.sources)
        def eda_plots_interact(file):
            f, axes = plt.subplots(1, 3)
            f.set_size_inches(15, 7)

            run_summaries = self.summaries[self.summaries.source == file]
            ssums = run_summaries[run_summaries.success == True]

            for ix, col in enumerate(cols):
                sns.histplot(data=ssums, x=col, ax=axes[ix]).set_title(col)

    def interactive_2d(self, cols=["Hdepth%", "%>2", "runtime"], kind="hist"):
        @interact(file=self.sources)
        def eda_plots_interact(file):
            run_summaries = self.summaries[file]
            ssums = run_summaries[run_summaries.success == True]
            f, axes = plt.subplots(1, 3)
            f.set_size_inches(15, 7)

            for ix, comb in enumerate(it.combinations(cols, 2)):
                plt.hist2d(
                    ssums[comb[0]], ssums[comb[1]], bins=(15, 15), cmap=plt.cm.BuPu
                )
                axes[ix].set(xlabel=comb[0], ylabel=comb[1])

            plt.show()


class plot_interact:
    def single_data_plot(scomb, module, rx):

        # scomb=scomb.sort_values(rx, ascending= False)
        source = scomb.source.unique()[0]

        sns.set(rc={"figure.figsize": (15, 8)})

        sns.boxplot(data=scomb, x=module, y=rx).set_title(source)

        plt.show()

    def explain(self, dataset, Rx, modules):
        dropdown_source = widgets.Dropdown(
            options=sorted(dataset.source.unique()),
            value=sorted(dataset.source.unique())[0],
            description="source:",
        )
        dropdown_variable = widgets.Dropdown(
            options=Rx, value=Rx[0], description="Var:"
        )
        dropdown_module = widgets.Dropdown(
            options=modules, value=modules[0], description="Module:"
        )

        def dropdown_source_eventhandler(change):
            """
            Eventhandler for the state dropdown widget
            """
            # display(input_widgets)
            source_choice = change.new
            # output_by_state(source_choice, dropdown_variable.value, dropdown_module.value)
            # IPython.display.clear_output(wait=True)

        def dropdown_variable_eventhandler(change):
            """
            Eventhandler for the question dropdown widget
            """
            # display(input_widgets)
            variable_choice = change.new
            # output_by_state(dropdown_source.value, variable_choice, dropdown_module.value)
            # IPython.display.clear_output(wait=True)

        def dropdown_module_eventhandler(change):
            """
            Event handler for the stratification dropdown widget
            """
            # display(input_widgets)
            module_choice = change.new
            # output_by_state(dropdown_source.value, dropdown_variable.value, module_choice)
            #

        dropdown_source.observe(dropdown_source_eventhandler, names="value")
        dropdown_variable.observe(dropdown_variable_eventhandler, names="value")
        dropdown_module.observe(dropdown_module_eventhandler, names="value")

        @interact(
            source=dropdown_source, variable=dropdown_variable, module=dropdown_module
        )
        def output_by_state(source, variable, module):
            """
            Takes in a state value, the specific question from the for loop and the specified stratification category.
            This function is called by the dropdown handlers below to pull the data based on user-input.
            """
            # IPython.display.clear_output(wait=True)

            output_data = dataset[dataset.source == source]
            output_data = single_data_plot(output_data, module, variable)


def single_data_plot(scomb, module, rx):

    # scomb=scomb.sort_values(rx, ascending= False)
    source = scomb.source.unique()[0]

    sns.set(rc={"figure.figsize": (15, 8)})

    sns.boxplot(data=scomb, x=module, y=rx).set_title(source)

    plt.show()


def dropdown_menu_eda(dataset, Rx, modules):

    output = widgets.Output()
    dropdown_source = widgets.Dropdown(
        options=sorted(dataset.source.unique()),
        value=sorted(dataset.source.unique())[0],
        description="source:",
    )
    dropdown_variable = widgets.Dropdown(options=Rx, value=Rx[0], description="Var:")
    dropdown_module = widgets.Dropdown(
        options=modules, value=modules[0], description="Module:"
    )

    def output_by_state(source, variable, module):
        """
        Takes in a state value, the specific question from the for loop and the specified stratification category.
        This function is called by the dropdown handlers below to pull the data based on user-input.
        """

        output_data = dataset[dataset.source == source]

        output_data = single_data_plot(output_data, module, variable)

        with output:
            display(output_data)

    def dropdown_source_eventhandler(change):
        """
        Eventhandler for the state dropdown widget
        """
        display(input_widgets)
        source_choice = change.new
        output_by_state(source_choice, dropdown_variable.value, dropdown_module.value)
        IPython.display.clear_output(wait=True)

    def dropdown_variable_eventhandler(change):
        """
        Eventhandler for the question dropdown widget
        """
        display(input_widgets)
        variable_choice = change.new
        output_by_state(dropdown_source.value, variable_choice, dropdown_module.value)
        IPython.display.clear_output(wait=True)

    def dropdown_module_eventhandler(change):
        """
        Event handler for the stratification dropdown widget
        """
        display(input_widgets)
        module_choice = change.new
        output_by_state(dropdown_source.value, dropdown_variable.value, module_choice)
        IPython.display.clear_output(wait=True)

    dropdown_source.observe(dropdown_source_eventhandler, names="value")
    dropdown_variable.observe(dropdown_variable_eventhandler, names="value")
    dropdown_module.observe(dropdown_module_eventhandler, names="value")

    input_widgets = widgets.HBox([dropdown_source, dropdown_variable, dropdown_module])

    display(input_widgets)
    output_by_state(sorted(dataset.source.unique())[0], Rx[0], modules[0])

    IPython.display.clear_output(wait=True)
