import copy
import itertools as it
import os
import subprocess
import time
from threading import Thread

import numpy as np
import pandas as pd


class metaclass_run():
    """
    class to setup metagenomics classification run.
    """
    def __init__(self, id= "1", rdir= "", child= ""):

        if len(rdir):
            self.parse_config(rdir, child= child)

        else:
            self.id = id
            self.actions = {
                "CLEAN": False,
                "QCONTROL": True,
                "ASSEMBLE": True,
                "DEPLETE": True,
                "CLASSIFY": False,
                "REMAP": False
            }
            self.r1 = ""
            self.r2 = ""
            self.paired = False

    def parse_config(self, rdir, child= ""):
        """
        function to return class instance by parsong config in given directory.
        :param rdir: previously spawned directory with config file.
        :param child: ID of child config file present in rdir.
        :return: self
        """
        if len(child):
            confgf= rdir + "confch_{}.sh".format(child)
            self.configf= confgf
            self.main= rdir + "main_{}.sh".format(child)
        else:
            confgf= rdir + "config.sh"
            self.configf= confgf
            self.main= rdir + "main.sh"
        settings= {}

        with open(confgf, "r") as fp:
            ln= fp.readline().strip()
            while ln:
                if ln[0] == "#":
                    ln= fp.readline().strip()
                    continue
                ln=ln.split("=")
                settings[ln[0]]= ln[1].strip('"')

                ln = fp.readline().strip()

        self.actions= {}
        self.config= {}
        self.params= {}

        for dr in ACTIONS.keys():
            if dr in settings.keys(): self.actions[dr]= settings[dr] == "true"

        for dr in SOFTWARE.keys():
            if dr in settings.keys(): self.config[dr]= [settings[dr]]

        pooled_params= {**self.ARGS_ENRICH, **self.ARGS_ASS, **self.ARGS_CLASS, **self.ARGS_REMAP}
        for dr in pooled_params.keys():
            if dr in settings.keys(): self.params= [settings[dr]]

        self.dir= settings["RDIR"]
        self.id= settings["SUFFIX"]
        self.r1= settings["INPUT"]
        self.config= pd.DataFrame(self.config)
        self.params= pd.DataFrame(self.params)

        if "PAIR" in settings.keys():
            self.r2= settings["PAIR"]
            self.paired= True


    def config_get(self,config_frame):
        """
        function to incorporate SOFTWARE data frame into self.config.
        :param config_frame: data FRAME. columns= modules, single row, elemets = software.
        :return: self with updated config.
        """

        self.config= config_frame

        if "HD" in self.config.columns:
            if len(self.config["HD"][0]) == 0:
                self.actions["DEPLETE"]= False

    def params_get(self,param_frame):
        """
        function to incorporate parameter data frame.
        :param param_frame:
        :return:
        """
        self.params= param_frame

    def write_config(self, config="config.sh"):
        """
        write config file based on actions, dirs, config and params dfs.
        :param config:
        :return:
        """
        #
        self.configf = config
        lines = ["#!/bin/bash"]
        #
        lines.append("## INPUT")
        lines.append('{}="{}"'.format("INPUT", self.r1))
        if self.r2:lines.append('{}="{}"'.format("PAIR", self.r2))
        else: lines.append('{}=``'.format("PAIR", self.r2))

        lines.append("##")
        #
        for dr, path in SOURCE.items():
            lines.append('{}="{}"'.format(dr, path))
        lines.append("##")
        #
        for dr, path in DIRS.items():

            lines.append('{}="{}"'.format(dr, self.dir + path))
        lines.append("##")
        #
        for dr, g in self.actions.items():
            lines.append('{}="{}"'.format(dr, str(g).lower()))
        lines.append("##")
        #
        for dr in self.config.columns:
            lines.append('{}="{}"'.format(dr, self.config[dr][0]))
        lines.append("##")
        #
        self.params= self.params.loc[:,~self.params.columns.duplicated(keep= "first")]
        for dr in self.params.columns:
            lines.append('{}="{}"'.format(dr, self.params[dr][0]))
        lines.append("##")
        #
        lines.append('SUFFIX="{}"'.format(self.id))
        lines.append('RDIR="{}"'.format(self.dir))

        ###
        with open(config, "w") as fn:
            fn.write("\n".join(lines))

    def prep_env(self, rdir= ""):
        """
        from main directory bearing scripts, params.py and main.sh, create metagenome run directory
        :return:
        """
        #
        if not rdir:
            rdir= os.getcwd()

        self.dir = rdir + "run_{}/".format(self.id)
        os.system("mkdir -p " + self.dir)
        os.system("cp main.sh {}".format(self.dir))
        os.system("cp -r {} {}".format("scripts/", self.dir))
        for dir in DIRS.values():
            os.system("mkdir -p " + self.dir + dir)

    def spawn(self, fofn= "", config="config.sh", source="main.sh", sink="main.sh"):
        """
        write config file, copy and modify instance main.sh file.
        :param fofn: file of fastas, one per line, max two.
        :param config: name of config file
        :param source: name of main*.sh file to copy and modify.
        :param sink: destinaion main*.sh file. to become instance main.
        :return:
        """
        #
        if len(fofn):
            self.read_fofn(fofn)

        self.write_config(config=config)
        self.configf = config
        #
        if source != sink:
            os.system("""cp {} {}""".format(source, sink))
        #
        os.system("sed -i 's#.*CONF_FILE=.*#CONF_FILE={}#g' {}".format(config, sink))
        self.main = sink

    def read_fofn(self, fofn):
        """
        read file of files.
        :param fofn: file path
        :return:
        """
        d = 0
        with open(fofn, "r") as fn:
            ln = fn.readline().strip()
            while ln:
                if d == 0:
                    self.r1 = ln
                if d == 1:
                    self.r2 = ln
                    self.paired = True
                d += 1
                ln = fn.readline().strip()

    def run_main(self):
        """
        deploy meta classsification run
        :return:
        """
        #
        os.system("{}".format(self.main))

    def clean(self, delete= True, outf= ""):

        trash= [
            "remap/", "reads/",  "scripts/", "host_depletion/", "classification/"
        ]

        if outf:
            os.makedirs(outf + self.id, exist_ok= True)
            subprocess.Popen("mv {}*report.tsv {}".format(self.dir, outf + self.id + "/"),shell= True)

        subprocess.Popen("mv {}*tsv {}".format(self.dir,self.dir + DIRS["OUTD"]), shell= True)
        subprocess.Popen("mv {}*sh {}".format(self.dir,self.dir + DIRS["OUTD"]), shell= True)

        if delete:
            for dr in trash: subprocess.run(["rm", "-r", self.dir + dr])


class meta_orchestra():
    def __init__(self, sup= 1, down= 1, smax= 3, odir= ""):
        self.sup= sup
        self.down= down
        self.processes= {}
        self.smax= smax
        self.odir= odir


    def clean(self, delete: bool= True):
        for sid, sac in self.projects.items(): sac.clean(delete= delete, outf= self.outd)

    def sample_main(self, sample=1, cols=["QC", "HD", "ASSEMBLY_SOFT", "ASSEMBLE_CLASS"]):
        """
        sample module / software combinations from dictionaries in params.py. random.
        :param sample: how many combinations to sample. corresponds to number of metclass run directories to be created.
        :param cols: keys to sample combinations of in SOFTWARE dict in params.py
        :return: data frame.
        """
        venue= [self.SOFTWARE[x] for x in cols]
        venues= list(it.product(*venue))
        #
        if sample > 0 and sample < len(venues):

            vex = np.random.choice(list(range(len(venues))), sample, replace= False)
            venues= [venues[x] for x in vex]

        venues = pd.DataFrame(venues)
        venues.columns = cols
        #
        return venues

    def params_extract(self, show, pooled_params, sample=1):
        """
        takes list of software, which might or not have entries in the argument dictionaries.
        """
        new_features = []
        nvens = []
        #
        for soft in show:
            if soft in pooled_params.keys():
                for c, g in pooled_params[soft].items():
                    new_features.append(c)
                    nvens.append(g)
        #
        nvens = list(it.product(*nvens))
        if sample:
            vex = np.random.choice(list(range(len(nvens))), sample, replace= False)
            nvens= [nvens[x] for x in vex]
        nvens = pd.DataFrame(nvens).reset_index(drop=True)
        nvens.columns = new_features
        return nvens

    def record_runs(self, outf= "summary.tsv"):
        """
        record instance software and parameter informationo and print to common csv.
        :param outf: common csv file path.
        :return: None
        """
        simple= []
        for sd,tb in self.processes.items():
            tb.to_csv(sd + ".args.tsv", sep= "\t", header= True, index= False)
            simple.append(tb[["DIR","ID","rclass_report","aclass_report"]])

        simple= pd.concat(simple, axis= 0)
        header= not os.path.isfile(outf)
        simple.to_csv(outf, header= header, sep= "\t", index= False, mode= "a")

    def run_proc(self, srun, sac, hdconf, paramCombs):
        """
        prepare and deploy instance of metaclassification run with only post assembly streps.

        :param sac: metaclassification run class instance.
        :param hdconf: data frame of software combinations. 1 row.
        :param srun: row index in parameter df.
        :return: None
        """

        child = copy.deepcopy(sac)
        child.id = "child_{}".format("-".join(np.array(srun,dtype= str)))

        child.actions["QCONTROL"]= False
        child.actions["DEPLETE"] = False
        child.actions["ASSEMBLE"] = False
        child.actions["CLASSIFY"]= True
        child.actions["REMAP"]= True
        child.actions["CLEAN"]= True

        if len(paramCombs[srun[0]]) == 0:
            return self

        child.config_get(hdconf.loc[[srun[0]]].reset_index(drop=True))
        child.params_get(paramCombs[srun[0]].loc[[srun[1]]].reset_index(drop=True))
        #
        child.spawn(fofn= child.dir + DIRS["LOGD"] + "{}_latest.fofn".format(sac.id),
                        config=child.dir + "confch_{}.sh".format(child.id),
                    source=child.main, sink= child.dir + "main_{}.sh".format(child.id))

        cstart= time.perf_counter()
        child.run_main()
        cend= time.perf_counter()
        child.exect= cend - cstart
        #
        ### gather report / output, mix parent and child params.
        child_params= pd.concat([sac.config, child.config, sac.params, child.params], axis= 1)
        child_params["rclass_report"]= DIRS["CLASSD"] + "reads/{}.report.tsv".format(child.id)

        if sac.actions["ASSEMBLE"]== True:
            child_params["aclass_report"]= DIRS["CLASSD"] + "assembly/{}.report.tsv".format(child.id)
        else: child_params["aclass_report"]= "NA"

        child_params["ID"]= child.id
        child_params["DIR"]= child.dir
        child_params["TOD"]= sac.exect + child.exect

        self.processes[child.dir + child.id]= child_params

    def low_run(self,prj):
        """
        deploy children of current metaclass run instance, controls parallel thread deployment tto not exceed self.smax.
        :param prj: metaclass_run instance.
        :return: None
        """
        hdconf = self.sample_main(sample=self.down, cols=["CLASSM", "REMAP_SOFT"])
        params2 = {**self.ARGS_CLASS, **self.ARGS_REMAP}
        paramCombs = [self.params_extract(hdconf.iloc[idx], params2, sample=0) for idx in range(hdconf.shape[0])]


        possibilities= [[(x,y) for y in range(paramCombs[x].shape[0])] for x in range(hdconf.shape[0])]
        suprelay= list(it.chain(*possibilities))

        print("low run possibilities: {}".format(len(suprelay)))

        if self.down > 0 and self.down < len(suprelay):
            vex= np.random.choice(list(range(len(suprelay))),self.down, replace= False)
            suprelay= [suprelay[x] for x in vex]
        #
        rlow= len(suprelay)
        #
        if rlow > self.smax:
            pack = list(np.arange(0, rlow, self.smax)) + [rlow]
            lset = [list(range(pack[x], pack[x + 1])) for x in range(len(pack) - 1)]
            #
            for a_args in lset:
                threads = [Thread(target=self.run_proc, args=(suprelay[x], prj, hdconf.copy(), paramCombs.copy())) for x in a_args]
                for th in threads: th.start()
                for th in threads: th.join()

        else:
            a_args = list(range(self.down))
            threads = [Thread(target=self.run_proc, args=(suprelay[x], prj, hdconf.copy(), paramCombs.copy())) for x in a_args]
            for th in threads: th.start()
            for th in threads: th.join()


    def low_deploy(self, sac= ""):
        """
        deploy post assembly steps across metaclass runs.
        :param sac: metaclass_run instance.
        :return: self
        """
        if sac:
            sel= [x for x,g in self.projects.items() if g.dir == sac][0]
            self.low_run(sac)

        else:
            for idr, prj in self.projects.items():
                self.low_run(prj)

        return self

    def param_prepare(self, technology, filters= {}):

        try:
            param_module = __import__("metaparams.params_%s" % technology, fromlist=["metaparams"])
        except ImportError:
            # Display error message
            print("couldnt find params")

        self.QC= param_module.ARGS_QC
        self.HD= param_module.ARGS_ENRICH
        self.ASSM= param_module.ARGS_ASS
        self.CLASSM= param_module.ARGS_CLASS
        self.REMAP= param_module.ARGS_REMAP
        self.SOFTWARE= param_module.SOFTWARE

        key2arg= {
            "QC": self.QC,
            "HD": self.HD,
            "ASSM": self.ASSM
            "ASSM_CLASS": self.ARGS_
            "CLASSM": self.CLASSM,
            "REMAP": self.CLASSM
        }

        for rset, sof in filters.items():
            self. SOFTWARE[rset]= [sof]



    def sup_deploy(self, fofn):
        """
        deploy metaclass runs.
        :return: self
        """
        conf= self.sample_main(sample= 0)
        params1 = {**self.ARGS_ENRICH, **self.ARGS_ASS}
        paramCombs = [self.params_extract(conf.iloc[idx], params1, sample=0) for idx in range(conf.shape[0])]

        possibilities= [[(x,y) for y in range(paramCombs[x].shape[0])] for x in range(conf.shape[0])]
        suprelay= list(it.chain(*possibilities))
        print("sup run possibilities: {}".format(len(suprelay)))

        if self.sup > 0 and self.sup < len(suprelay):
            vex= np.random.choice(list(range(len(suprelay))),self.sup, replace= False)
            suprelay= [suprelay[x] for x in vex]
        #
        rsup= len(suprelay)

        self.rdir= self.odir + os.path.basename(fofn) + "/"
        os.makedirs(self.rdir, exist_ok= True)
        self.outd= self.rdir + "output/"
        os.makedirs(self.outd, exist_ok= True)

        self.projects= {}

        if rsup > self.smax:
            pack= list(np.arange(0,rsup, self.smax)) + [rsup]
            lset= [list(range(pack[x],pack[x+1])) for x in range(len(pack)-1)]

            for a_args in lset:
                threads = [Thread(target=self.sup_run, args=(suprelay[x], fofn, conf.copy(), paramCombs.copy())) for x in a_args]
                for th in threads: th.start()
                for th in threads: th.join()

        else:
            a_args= list(range(rsup))
            threads= [Thread(target=self.sup_run,args= (suprelay[x], fofn,conf.copy(), paramCombs.copy())) for x in a_args]
            for th in threads: th.start()
            for th in threads: th.join()

        return self


    def sup_run(self, idx, fofn, conf, paramcomb):
        """
        prepare metaclass_run for deployment. create directories, config and main files, set up software and modules.
        :param idx: row index in software db.
        :param conf: software combination df. one line.
        :return: self
        """
        nrun = metaclass_run(id= "-".join(np.array(idx, dtype= str)))
        #
        nrun.config_get(conf.loc[[idx[0]]].reset_index(drop=True))
        #
        nrun.params_get(paramcomb[idx[0]].loc[[idx[1]]].reset_index(drop=True))
        #
        nrun.prep_env(rdir= self.rdir)

        nrun.spawn(fofn= fofn,
                   config= nrun.dir + "config.sh",
                   sink= nrun.dir + "main.sh")  ### fofn currently flobal variable.

        supstart= time.perf_counter()
        nrun.run_main()
        supend= time.perf_counter()
        nrun.exect= supend - supstart

        self.projects[nrun.id]= nrun

        return self
