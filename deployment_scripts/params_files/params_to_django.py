import pandas as pd
import params_illumina as pai
import params_nanopore as pan

Technology= pd.DataFrame([["illumina"], ["nanopore"]],
                         columns= ["name"])
Module= pd.DataFrame(list(pai.SOFTWARE.keys()), columns= ["name"])

Software= pd.DataFrame(columns= ["tech", "module", "name"])
Arg= pd.DataFrame(columns= ["tech", "module", "soft", "name"])
Flag= pd.DataFrame(columns= ["tech", "module", "soft", "arg", "name"])

source_dict= {
    "illumina": pai,
    "nanopore": pan
}

if __name__ == "__main__":
    for tech, source in source_dict.items():
        soft_mdict = {
            "QC": source.ARGS_HD,
            "HD": source.ARGS_HD,
            "ASSEMBLY_SOFT": source.ARGS_ASS,
            "ASSEMBLE_CLASS": source.ARGS_CLASS,
            "CLASSM": source.ARGS_CLASS,
            "REMAP_SOFT": source.ARGS_REMAP
        }

        for mod, softl in source.SOFTWARE.items():
            for sf in softl:
                nrow= pd.DataFrame({
                    "tech": [tech],
                    "module": [mod],
                    "name": [sf]
                })
                Software= pd.concat([Software,nrow],ignore_index=True, axis= 0)

        for mod, softl in source.SOFTWARE.items():
            argdict = soft_mdict[mod]
            for soft in softl:
                if soft in argdict.keys():
                    for arg, flags in argdict[soft].items():
                        nrow= pd.DataFrame({
                            "tech": [tech],
                            "module": [mod],
                            "soft": [soft],
                            "name": [arg]
                        })
                        Arg= pd.concat([Arg,nrow], ignore_index= True, axis= 0)

                        for fl in flags:
                            nrow= pd.DataFrame({
                                "tech": [tech],
                                "module": [mod],
                                "soft": [soft],
                                "arg": [arg],
                                "name": [fl]
                            })
                        Flag= pd.concat([Flag,nrow], ignore_index= True, axis= 0)

                else:
                    print(soft)


    print(Module)
