#!/bin/bash

CONF_FILE=config.sh
source $CONF_FILE
source $CONDA"etc/profile.d/conda.sh"

export CONDA
export ENVSDIR
export CONF_FILE

##
SECONDS=0
TIMESTAMP=`date +"%Y-%m-%d %T"`
WDH=""
FILTD=$WDH$FILTD
ASSD=$WDH$ASOUT$ASSEMBLY_SOFT"/"
export WDH

# PARAMS
OFQ=""

print_usage() {
    printf "Usage: ..."
}

## INPUT ##
if [ -z $INPUT ]; then echo "missing input"; exit 1; fi

IFILE=`basename $INPUT`
INAME=${IFILE%.f*}
CLEAN_R1=$WDH$CLEAND$INAME.fastq.gz
FILT_R1=$FILTD$INAME.$HD.filt.fq.gz
C1=$INPUT
echo $C1 >  $LOGD$SUFFIX"_latest.fofn"

export INAME

if [ ! -z $PAIR ]; then
    PNAME=`basename $PAIR`
    PNAME=${PNAME%.f*}
    CLEAN_R2=$WDH$CLEAND$PNAME.fastq.gz
    FILT_R2=$FILTD$PNAME.$HD.filt.fq.gz
    C2=$PAIR
    echo $C2 >> $LOGD$SUFFIX"_latest.fofn"
    export PNAME
fi

###############################################
### PREPROCESS
echo $INPUT $PAIR

if $QCONTROL; then
    echo "QUALITY CONTROL"
    if [ ! -f $CLEAN_R1 ]; then
        
        echo $PAIR
        $BIN"preprocess/process_data.sh" -a $INPUT -b $PAIR #2>> $LOGD$SUFFIX.preproc.log
        
    fi
    
    CN=`echo $(zcat ${C1} | wc -l)/4|bc`
    
    echo -e "INPUT\t$CN" >> $LOGD"reads_latest.stats"
    
    C1=$CLEAN_R1
    C2=$CLEAN_R2
    echo $C1 >  $LOGD$SUFFIX"_latest.fofn"
    echo $C2 >> $LOGD$SUFFIX"_latest.fofn"
    
    CN=`echo $(zcat ${C1} | wc -l)/4|bc`
    echo -e "QC\t$CN" >> $LOGD"reads_latest.stats"
    
fi

################################################
##### HOST DEPLETION

if $ENRICH; then
    echo "ENRICHMENT"
    
    $BIN"classification/class_selector.sh" -a $C1 \
    -d $HDOUT -o $FILTD -m $HD -p f -s $SUFFIX -b $C2 #2>> $LOGD$SUFFIX.depletion.log
    
    C1=$FILT_R1
    C2=$FILT_R2
    echo $C1 >  $LOGD$SUFFIX"_latest.fofn"
    echo $C2 >> $LOGD$SUFFIX"_latest.fofn"
    
    CN=`echo $(zcat ${C1} | wc -l)/4|bc`
    echo -e "ENRICHMENT\t$CN" >> $LOGD"reads_latest.stats"
    
fi

if $DEPLETE && [ ! -f $REFERENCE ]; then
    echo "DEPLETION"
    
    conda activate $ENVSDIR"remap/remap"
    
    if [ ! -f $REFERENCE ]; then
        bwa mem -t 8 $REFERENCE $C1 $C2 | samtools view -Sb -f 4 - > samtools fastq -1 $FILTD$INAME.$HD.filt.fq -2 $FILTD$PNAME.$HD.filt.fq -0 /dev/null -s /dev/null -
        gzip $FILTD$INAME.$HD.filt.fq
        gzip $FILTD$PNAME.$HD.filt.fq
        
        C1=$FILT_R1
        C2=$FILT_R2
        echo $C1 >  $LOGD$SUFFIX"_latest.fofn"
        echo $C2 >> $LOGD$SUFFIX"_latest.fofn"
        
    else
        
        bwa mem -t 8 $REFERENCE $C1 | samtools view -Sb -f 4 - > samtools fastq -1 $FILTD$INAME.$HD.filt.fq -0 /dev/null -s /dev/null -
        gzip $FILTD$INAME.$HD.filt.fq
        
        C1=$FILT_R1
        echo $C1 >  $LOGD$SUFFIX"_latest.fofn"
        
    fi
    
    CN=`echo $(zcat ${C1} | wc -l)/4|bc`
    echo -e "DEPLETION\t$CN" >> $LOGD"read_latest.stats"
    
    
fi


if [ `zcat $C1 | wc -l` -gt 0 ]; then
    
    echo $C1 >  $LOGD$SUFFIX"_latest.fofn"
    
    if [ ! -z $PAIR ]; then
        echo "CLEANING"
        
        echo $C2 >> $LOGD$SUFFIX"_latest.fofn"
        
        $BIN"preprocess/clean_unpaired.sh" -a $C1 -b $C2
        
    fi
    
fi



CN=`echo $(zcat ${C1} | wc -l)/4|bc`
echo -e "CLEAN\t$CN" >> $LOGD"reads_latest.stats"

#################################################
##### ASSEMBLY
if $ASSEMBLE; then
    echo "ASSEMBLY"
    touch $LOGD"assembly.fofn"
    
    if [ $ASSEMBLY_SOFT == "spades" ]; then
        
        $BIN"assembly/spades.sh" -a $C1 -d $ASSD -o $ASSD -b $C2 2>> $LOGD$SUFFIX.assembly.log
        
        elif [ $ASSEMBLY_SOFT == "velvet" ]; then
        
        $BIN"assembly/velvet.sh" -a $C1 -d $ASSD -o $ASSD -b $C2 2>> $LOGD$SUFFIX.assembly.log
        
    fi
    
    ##################################################
    ### filter assembly
    ASSEMBLY=$WDH$ASOUT$ASSEMBLY_SOFT"/"$INAME"/"$INAME.scaffolds.fasta.gz
    ASSEMBLY_FILTERED=$WDH$ASOUT$ASSEMBLY_SOFT"/"$INAME.filtered.scaffolds.fasta
    
    conda activate $ENVSDIR"Pyenv/pyenv"
    
    python $BIN"assembly/filter_fasta.py" \
    --input $ASSEMBLY \
    --output $ASSEMBLY_FILTERED \
    --length $ASSEMBLY_LTRIM \
    --dir $LOGD
    
    gzip $ASSEMBLY_FILTERED
    ASSEMBLY_FILTERED=$ASSEMBLY_FILTERED.gz
    
    conda deactivate
    
    ###################################################
    ##### CLASSIFICATION
    
    if [ `zgrep "^>" $ASSEMBLY_FILTERED | wc -l` -gt 0 ]; then
        echo "ASSEMBLY CLASSIFICATION"
        echo $CLASSD"assembly/"
        
        $BIN"classification/class_selector.sh" -a $ASSEMBLY \
        -d $CLASSD"assembly/" \
        -o $CLASSD"assembly/" \
        -m $ASSEMBLE_CLASS \
        -p r -s "assembly" 2>> $LOGD$SUFFIX.aclass.log
        
        echo $ASSEMBLY_FILTERED > $LOGD"assembly.fofn"
        
        if $VIRSORT; then
            echo "VIRSORT2"
            
            $BIN"remap/virsorter2.sh" -a $ASSEMBLY_FILTERED
            
        fi
    else
        echo $ASSEMBLY_FILTERED > $LOGD"assembly_missing.fofn"
        
    fi
    
    
fi
####################### END OF ASSSEMBLY
if $CLASSIFY; then
    
    $BIN"classification/class_selector.sh" -a $C1 -d $CLASSD"reads/" \
    -m $CLASSM -p r -s $SUFFIX -b $C2 2>> $LOGD$SUFFIX.rclass.log
    
fi

##### MAPPING

if $REMAP; then
    mkdir -p $LOGD$SUFFIX"/"
    
    conda activate $ENVSDIR"Pyenv/pyenv"
    
    python $BIN"remap/merge_results.py" \
    --classd $CLASSD \
    --suffix $SUFFIX \
    --taxd $METAD \
    --logd $LOGD$SUFFIX"/" \
    --max_depth $REF_KEEP
    
    conda deactivate
    
    #####
    
    PATTERN=$CLASSD$SUFFIX*.targets
    
    if ls $PATTERN 1> /dev/null 2>&1; then
        for REF in $PATTERN; do
            REFN=`basename $REF`
            REFN=${REFN%.targets}
            REMAP_REF=`cat $CLASSD$REFN.file`
            mkdir -p $LOGD$SUFFIX"/"$REFN"/"
            echo $acc
            echo $REF
            
            conda activate $ENVSDIR"remap/remap"
            >$REMD$REFN.fasta
            for acc in `cat $REF`; do
                samtools faidx $REF_FASTA$REMAP_REF $acc >> $REMD$REFN.fasta
            done
            conda deactivate
            
            echo "ref here: " $REMD$REFN.fasta
            
            if [ $REMAP_SOFT == "snippy" ]; then
                $BIN"remap/snippy.sh" -a $C1 -r $REMD$REFN.fasta -d $REMD -s "ref" -b $C2 2>> $LOGD$SUFFIX".remap.log"
                
                elif [ $REMAP_SOFT == "rematch" ]; then
                
                $BIN"remap/rematch.sh" -a $C1 -r $REMD$REFN.fasta -d $REMD -b $C2 2>> $LOGD$SUFFIX".remap.log"
                
                elif [ $REMAP_SOFT == "bowtie" ]; then
                
                $BIN"remap/bowtie.sh" -a $C1 -r $REMD$REFN.fasta -d $REMD -b $C2  2>> $LOGD$SUFFIX".remap.log"
                
                elif [ $REMAP_SOFT == "minimap-rem" ]; then
                
                $BIN"remap/minimap2.sh" -a $C1 -r $REMD$REFN.fasta -d $REMD -b $C2 2>> $LOGD$SUFFIX".remap.log"
                
            fi
            
            if [ -s $LOGD"assembly.fofn" ]; then
                ASSEMBLY=`cat $LOGD"assembly.fofn"`
                gunzip $ASSEMBLY
                ASSEMBLY=${ASSEMBLY%.gz}
                
                $BIN"remap/hd_minimapASM.sh" -a $ASSEMBLY -f \
                -q $REMD$REFN.fasta \
                -s $REFN \
                -d $REMD$SUFFIX"/"$REFN"/" 2>> $LOGD$SUFFIX".aremap.log"
                
                conda activate $ENVSDIR"remap/remap"
                bgzip $ASSEMBLY
                conda deactivate
                
                if [ -s $REMD$SUFFIX"/"$REFN"/minimap2/"$REFN.paf ]; then
                    
                    conda activate $ENVSDIR"remap/Renv"
                    Rscript --vanilla $BIN"remap/pafr.R" $REMD$SUFFIX"/"$REFN"/minimap2/"$REFN.paf \
                    $REMD$SUFFIX"/"$REFN"/minimap2/"
                    conda deactivate
                    
                    #conda activate $ENVSDIR"Pyenv/pyenv"
                    #python $BIN"remap/dotplot.py" \
                    #$REMD$SUFFIX"/"$REFN"/minimap2/"$REFN.paf \
                    #$REMD$SUFFIX"/"$REFN"/minimap2/"
                    #conda deactivate
                    
                fi
                
                mapped_scaffolds=$REMD$SUFFIX"/"$REFN"/minimap2/"$REFN.filt.fa
                
                if [ `zgrep "^>" $mapped_scaffolds.gz | wc -l` -gt 0 ]; then
                    
                    R1_MAP=$C1
                    if [ -f $REMD$SUFFIX"/"$REFN"/snps.R1.kept.fq.gz" ]; then
                        echo "REFERENCE MAPPED READS FOUND, MAPPING TO ASSEMBLY"
                        echo "R1_MAP" $REMD$SUFFIX"/"$REFN"/snps.R1.kept.fq.gz"
                        R1_MAP=$REMD$SUFFIX"/"$REFN"/snps.R1.kept.fq.gz"
                    fi
                    
                    R2_MAP=$C2
                    if [ -f $REMD$SUFFIX"/"$REFN"/snps.R2.kept.fq.gz" ]; then
                        R2_MAP=$REMD$SUFFIX"/"$REFN"/snps.R2.kept.fq.gz"
                    fi
                    
                    gunzip $mapped_scaffolds.gz
                    mv $mapped_scaffolds $REMD$SUFFIX"/"$REFN"/assembly".fasta
                    
                    echo "ATTENTION"
                    echo $mapped_scaffolds $REMD$SUFFIX"/"$REFN"/assembly".fasta
                    
                    $BIN"remap/snippy.sh" -a $R1_MAP \
                    -r $REMD$SUFFIX"/"$REFN"/assembly".fasta \
                    -d $REMD$SUFFIX"/"$REFN"/" -s "assembly" \
                    -b $R2_MAP 2>> $LOGD$SUFFIX".rremap.log"
                    
                    if [ -s $REMD$SUFFIX"/"$REFN"/"$SUFFIX"/assembly/snps.stats" ]; then
                        mv $REMD$SUFFIX"/"$REFN"/"$SUFFIX"/assembly/snps.stats" $LOGD$SUFFIX"/"$REFN"/scaffold_remap.stats"
                    fi
                    conda activate $ENVSDIR"remap/remap"
                    bgzip $REMD$SUFFIX"/"$REFN"/assembly".fasta
                    conda deactivate
                    
                fi
                
            fi
            
            #rm $REMD$REFN.fasta
            #rm $REMD$REFN.fasta.fai
            
        done
        
        
        conda activate $ENVSDIR"Pyenv/pyenv"
        python $BIN"remap/merge_reports.py" \
        --suffix $SUFFIX \
        --classd $CLASSD \
        --remd $REMD \
        --logdir $LOGD \
        --taxd $METAD \
        --output $RDIR$SUFFIX.final.report.tsv \
        --seconds $SECONDS
        
        conda deactivate
        
    fi
    
fi

##### SUMMARY & REPORT PREP

####### CLEAN

if $CLEAN ; then
    
    rm -r $CLASSD$SUFFIX* $CLASSD"reads/"$CLASSM"/" $CLASSD"reads/"$CLASSM"/"
    
fi


