import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

import seaborn as sns
import os

from ipywidgets import interact, interactive, fixed, interact_manual
import ipywidgets as widgets

class run_eda():
    
    def __init__(self, validation):
        self.validation= validation
    
    
    def run_input(self, rundir):
        """
        return report data sets
        """
        dirs= next(os.walk(rundir), (None, None, []))[1]
        files= {}

        data= []

        for dir in dirs:
            #print(dir)
            files[dir]= {}
            outd= f"{rundir}{dir}/output/"
            outdirs= os.listdir(outd)

            logs= [f"{rundir}{dir}/run_{x}/" for x in outdirs]
            #outdirs= [f"{outd}/{x}/" for x in outdirs]

            for suprun in outdirs:

                odir= f"{outd}/{suprun}/"
                reports= os.listdir(odir)
                if reports:

                    for rep in reports:

                        rep= pd.read_csv(f"{odir}/{rep}", sep= "\t")
                        ocols= list(rep.columns)
                        low_run= rep.suffix.unique()[0]
                        rep["run"]= [f"{dir}_{suprun}_{low_run}"] * rep.shape[0]
                        rep["source"]= [dir] * rep.shape[0]

                        rep= rep[["source", "run"] + ocols]

                        if len(data): 
                            data= pd.concat((data, rep), axis= 0)
                        else:
                            data= rep
        
        self.files= files
        self.odat= data
        self.data_total= data[data.ID == "total"]
        
    
    def split(self):
        self.source_total= {}
        
        for source in self.data_total.source.unique():
            
            self.source_total[source]= self.data_total[self.data_total.source == source].reset_index(drop= True)
        
        self.sources= list(self.source_total.keys())
        
        return self
    
    def run_summary(self, stotal):
        run_summaries= []
        for run in stotal.run.unique():
            rda= stotal[stotal.run == run].reset_index(drop= True)
            found_success= sum(rda.success) > 0
            runtime= max(rda.time)
            percent_over= 0
            Hdepth= 0
            assh= False
            dsuc= rda[rda.success == True]

            if len(dsuc):
                percent_over= np.mean(rda["%>2"])
                Hdepth= np.mean(rda["Hdepth%"])
                assh= sum(rda.aclass == True) > 0

            run_summaries.append([run, found_success, runtime, percent_over, Hdepth, assh])

        run_summaries= pd.DataFrame(run_summaries, 
                               columns= ["run", "success", "runtime", "%>2", "Hdepth%", "ass_help"]
                                   )
        return run_summaries
    
    
    def summarize(self):
        self.summaries= {}
        self.descriptions= {}
        for file in self.sources:
            
            infection=self.validation[file][0]

            stotal= self.source_total[file]
            stotal["success"]= stotal.description.str.contains(infection)
            run_summaries= self.run_summary(stotal)
            self.summaries[file]= run_summaries
    
    def describe(self, file= ""):
        self.descriptions= {}
        
        for file in self.sources:
            stotal= self.source_total[file]
            infection=self.validation[file][0]
            run_summaries= self.summaries[file]
            ssums= run_summaries[run_summaries.success == True]
            ###
            psc_runs= sum(run_summaries.success) / run_summaries.shape[0]
            psc_runs= round(psc_runs * 100, 3)
            report_runs= f"total number of runs: {run_summaries.shape[0]}"
            report_success= f"N runs that find {infection} : {ssums.shape[0]},  ({psc_runs} %)"
            ##
            scruns= stotal[stotal.success == True]
            species_unique= "species found with this description: " + ", ".join(list(scruns.description.unique()))
            ##
            allcover= (ssums.success == True) & (ssums.ass_help == True)
            allcover= sum(allcover)
            succass= round(100 * allcover / ssums.shape[0], 3)
            succtotal= round(100 * allcover / run_summaries.shape[0], 3)

            report1= f"percent of successfull runs where assembly is also successful : {succass} %"
            report2= f"percent of runs that find {infection} with both reads and assembly : {succtotal} %"
            
            self.descriptions[file]= [report_runs, report_success, species_unique, report1, report2]
    
    def pretty_print(self, file): 
        print(f"source: {file}, infection: f{self.validation[file]}")
        print("\n".join(self.descriptions[file]) + "\n")
        
    def describe_print(self, file= ""):
        
        if file:
            self.pretty_print(file)
        else:
            for file in self.sources:
                self.pretty_print(file)
    
    def eda_plots(self, file, cols= ["Hdepth%", "%>2", "runtime"]):
        
        run_summaries= self.summaries[file]
        ssums= run_summaries[run_summaries.success == True]
        
        for col in cols: 
            sns.histplot(data= ssums, x= col).set_title(file)
            plt.show()
    
    def interactive(self, cols= ["Hdepth%", "%>2", "runtime"]):
        @interact (file= self.sources)
        def eda_plots_interact(file):

            run_summaries= self.summaries[file]
            ssums= run_summaries[run_summaries.success == True]

            for col in cols: 
                sns.histplot(data= ssums, x= col).set_title(file)
                plt.show()

