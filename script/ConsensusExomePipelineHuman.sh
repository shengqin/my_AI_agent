#!/bin/bash

#Wrapper for consensus variant calling from input paired tumor/normal sorted deduped bam files
#required arguments 1-sample list input file 2-output directory #sample input file should be in the format:
#sample_name-T1	fullpath/sample_tumor_sorted_deduped.bam	fullpath/sample_normal_sorted_deduped.bam
#each sample name should be unique, use "-T2" etc. to designate multiple samples from the same patient
#prior to running - modify the Settings below

##########     Settings - please modify this section as desired    #################
CNV_kit_switch="off" #"on" will run CNVkit for copynumber
CNV_reference_switch="on" #on will generate pooled reference for all samples in directory
CNV_target_antitarget="on" #will create target/antitarget coverage
main_pipeline="on" #turn off and on the main consensus and annotation portion of the srcipt, useful for when running filtering/plotting at the end
filter_and_plot="off" #if running several samples in parallel into the same results folder, wait until they are done to run this
Absolute_switch="off" #run this after monte_carlo monitoring
indel_switch="on" #for use with nondeduped bams to generate a list of indels.

export scripts_dir="/oak/stanford/groups/emoding/scripts/exome_pipeline_scripts/wespipeline/scripts" #Important first directory specification to show where required scripts are
export index_dir="/oak/stanford/groups/emoding/scripts/exome_pipeline_scripts/indices"
export dependencies_dir="/oak/stanford/groups/emoding/scripts/exome_pipeline_scripts/dependencies" #Specify where indices and dependencies are
export selector_bed="$index_dir/xGen_Exome_Research_Panel_v2_targets.bed"  #enter full path of the capture selector bed; default is MedExome
export selector_padded_bed="$index_dir/xGen_Exome_Research_Panel_v2_probes.bed" #enter full path of the capture selector bed; default is MedExome
export cnv_selector="$index_dir/MedExome.gene_label.bed" #modified capture .bed file with only gene name in 4th column; needed for better visualization on cnv diagram- if unavailable can use the selector bed file

export bam_files=$3 #input from command line specifying bam file location
echo $bam_files
export nondeduped_bam_files=$4 #input from command line specifying nondeduped bam file location
echo $nondeduped_bam_files
export results_folder=$5 #Where should the results go?


#######   Should not require modification past this line    ##########################################

#required scripts
	#MSK_process_variants.sh tag
	#MSK_process_annovar_output.sh annovar.multianno.txt annovar.coding.file tag
	#MSK_run_varscan.sh	
	#merge_VCFs.R
	#MSK_run_stats.sh
	#plot_coverage.R 
	#MSK_polysolver.sh

#Dependenciesexport 
export annovar_dir="$dependencies_dir/annovar_2015-03-22" #dir to annovar scripts
export samtools="$dependencies_dir/samtools-0.1.19/samtools" #samtools command
export varscan2_4_1_dir="$dependencies_dir/varscan-master"
export bam_readcount="$dependencies_dir/bam-readcount/bin/bam-readcount"
export strelka_dir="$dependencies_dir/strelka_workflow-1.0.14"
export GATK_dir="$dependencies_dir/gatk-4.1.9.0/"
export strelka_exome_config="$strelka_dir/etc/strelka_config_bwa_exome.ini" #this file may need to be modified for non-exome variant calling
export mutect_jar="$dependencies_dir/mutect/mutect-1.1.7.jar"
export cnvkit_script="$dependencies_dir/cnvkit-0.9.6/cnvkit.py"
export plot_coverage="$scripts_dir/plot_coverage.R"
export PATH=$PATH:$scripts_dir/dependencies/ensembl-vep/:$dependencies_dir/bin/bin/ #Add some tools to the PATH
export vep="$dependencies_dir/ensembl-vep/vep"
export convert2bed_script="$dependencies_dir/bin/convert2bed"
export bam2freq_script="$scripts_dir/bam-snvfreq-single.py"
export annotate_vars_script="$scripts_dir/annotate_var_list_consensus_ajay.py"
export filter_and_plot_script="$scripts_dir/ExomePlottingHuman.R"
export run_absolute="$scripts_dir/RunAbsolute.R"

#indices/references
export reference_fasta="$index_dir/hg19.fa" 
export db="$annovar_dir/humandb" #dir to annovar human db files
export cnv_access_file="$dependencies_dir/CNVkit/cnvkit-master/data/access-10kb.hg19.bed" #for cnvkit
export repeatmasker="$index_dir/hg19.masked.fixed.bed" #Repeatmasker reference
export genomad="$index_dir/af.only.genomad.hg19.chr_v2.bed" #POP AF reference
export annotate_genomad_script="$scripts_dir/annotate_genomAD.R" #for annotating POPAF
export refFlat="$index_dir/hg19refFlat.txt" #for CNVkit

#echo $PATH

if [ "$#" -ne "5" ]; then
        echo "Incorrect number of arguments.  Please run as"
        echo "placeholder"
        input_mode="multi"
        exit
else
        export input_file=$1
        export starting_dir=$(echo $2 | sed 's/\/$//') #removes last slash if present in directory name
        if [ $2 == "." ]; then
                starting_dir=$PWD
        fi
        export working_dir="$starting_dir/exome_pipeline"
fi


#functions
create_dir_structure() { #make directories; Argument 1 is space delimited list of samples
mkdir -p "$working_dir"
mkdir -p "$working_dir/bam_files"
mkdir -p "$working_dir/variants"
#export sample_table="$working_dir/sample_table.temp"
#awk '{print $1,$1,$2,$3}' $input_file | awk '{gsub(/-T.*/,"",$2)}1' > $working_dir/sample_table.temp #sample list file for downstream scripts
export tumor_samples=$(awk '{print $1}' $input_file)
export normal_samples=$(awk '{print $1}' $input_file | sed 's/-T.*//')

for sample in $tumor_samples
        do
        mkdir -p "$working_dir/$sample"
done
export bam_dir="$working_dir/bam_files"
export HLA_dir="$working_dir/HLA_files"
export run_stat_dir="$working_dir/run_stats"
export variant_dir="$working_dir/variants"
export cnv_dir="$working_dir/CNVkit"
}

echo $tumor_samples

########## added this as a test ###########
ml load biology bedtools


###### Real Functions Begin Here ######

call_varscan_variants() {
    for tag in $tumor_samples
do
	$scripts_dir/MSK_run_varscan.sh $varscan2_4_1_dir $reference_fasta $bam_readcount $working_dir $tag $bam_files $scripts_dir
	mkdir -p $results_folder/somatic_and_germline
	cp $tag.vcf.snp.Germline.hc $results_folder/somatic_and_germline
done
}

run_mutect_caller() {
mutect_jar=$1
reference_fasta=$2
working_dir=$3
tag=$4

tumorbam=$bam_files/$tag/${tag}-T*.bam
normalbam=$bam_files/$tag/${tag}-N*.bam

echo $reference_fasta
echo $tumorbam
echo $normalbam
echo $6

cd $working_dir/$tag
$GATK_dir/gatk Mutect2 -R $reference_fasta -I $tumorbam -I $normalbam -O $tag.mutect.raw.out -normal normal -L $selector_padded_bed

mkdir -p $results_folder/somatic_and_germline

$GATK_dir/gatk FilterMutectCalls -V $tag.mutect.raw.out -O $tag.mutect.pass.raw.out -R $reference_fasta
cp $tag.mutect.pass.raw.out $results_folder/somatic_and_germline

bcftools view -f PASS $tag.mutect.pass.raw.out > $tag.passed.vcf
bcftools view --no-header $tag.passed.vcf > $tag.mutect.passed.vcf

awk -F'[\t:]' '{print 100*$25}' $tag.mutect.passed.vcf > $tag.mutect.vaf.temp

paste $tag.mutect.passed.vcf $tag.mutect.vaf.temp > $tag.mutect.vaf.pre.temp

awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6,$10,$11,$NF, "mutect"}' $tag.mutect.vaf.pre.temp > $tag.mutect.convert2annovar.input

cd $working_dir
}
export -f run_mutect_caller

call_mutect_variants() {
for i in $tumor_samples
do
    run_mutect_caller $mutect_jar $reference_fasta $working_dir $i
done
}

run_strelka(){
strelka_dir=$1
strelka_exome_config=$2
reference_fasta=$3
working_dir=$4
tag=$5

tumorbam=$6
normalbam=$7

echo $tumorbam
echo $normalbam

rm -r $working_dir/$tag/strelka_workdir
mkdir $working_dir/$tag/strelka_workdir
cd $working_dir/$tag/strelka_workdir

echo $2

cp $strelka_exome_config config.ini
$strelka_dir/bin/configureStrelkaWorkflow.pl --normal=$normalbam --tumor=$tumorbam --ref=$reference_fasta --config=$strelka_exome_config --output-dir=$working_dir/$tag/strelka_workdir/results 
cd results
make -j 10
cd $working_dir/$tag
cat $working_dir/$tag/strelka_workdir/results/results/passed.somatic.snvs.vcf $working_dir/$tag/strelka_workdir/results/results/passed.somatic.indels.vcf | grep -v "#" > $working_dir/$tag/$tag.strelka.raw.vcf
export strelka_snv_count=$(grep -v "#" $working_dir/$tag/strelka_workdir/results/results/passed.somatic.snvs.vcf | wc -l)
}
export -f run_strelka

call_strelka_variants() {
echo $tumor_samples
for tag in $tumor_samples
do
echo $tag
run_strelka $strelka_dir $strelka_exome_config $reference_fasta $working_dir $tag $bam_files/$tag/${tag}-T*.bam $bam_files/$tag/${tag}-N*.bam
done
}

displaytime (){
  local T=$1
  local D=$((T/60/60/24))
  local H=$((T/60/60%24))
  local M=$((T/60%60))
  local S=$((T%60))
  (( $D > 0 )) && printf '%d days ' $D
  (( $H > 0 )) && printf '%d hours ' $H
  (( $M > 0 )) && printf '%d minutes ' $M
  (( $D > 0 || $H > 0 || $M > 0 )) && printf 'and '
  printf '%d seconds\n' $S
}

create_log() {
#rm $working_dir/log.txt
start=$(date)
echo -e "Start time $start" > $working_dir/$logfile
#echo "Demultiplexing duration ="<(displaytime $demux_duration) >> log.txt
#echo "Mapping duration ="<(displaytime $mapping_duration) >> log.txt
#echo "Deduping duration ="<(displaytime $rm_dup_duration) >> log.txt
#echo "Realign duration ="<(displaytime $realign_duration) >> log.txt
#echo "Recalibrate duration ="<(displaytime $recalibrate_bases_duration) >> log.txt
#echo "Stats calculation duration ="<(displaytime $runstats_duration) >> log.txt
}

consensus_variants() {
for i in $tumor_samples
do
	$scripts_dir/process_variants.sh $i $working_dir $annovar_dir $variant_dir $scripts_dir
	sed -i.bak $'s/\t/    /g' $working_dir/$tag/$tag.merged.avinput #replace merged input delimiter with spaces for compatibility later on
done
}

run_annovar() {
tag=$1
avinputfile="$working_dir/$tag/$tag.merged.avinput"
consensusavinputfile="$working_dir/$tag/$tag.consensus.avinput"
germlineavinputfile="$working_dir/$tag/$tag.germline.avinput"
protocol=refGene,snp138,cosmic70,clinvar_20140929,1000g2014oct_all,esp6500si_all,exac02
operation=g,f,f,f,f,f,f
perl $annovar_dir/table_annovar.pl $avinputfile $db --buildver hg19 -protocol $protocol -operation $operation -nastring "." --otherinfo

perl $annovar_dir/table_annovar.pl $consensusavinputfile $db --buildver hg19 -protocol $protocol -operation $operation -nastring "." --otherinfo

perl $annovar_dir/table_annovar.pl $germlineavinputfile $db --buildver hg19 -protocol $protocol -operation $operation -nastring "." --otherinfo

perl $annovar_dir/coding_change.pl --includesnp $consensusavinputfile.refGene.exonic_variant_function $db/hg19_refGene.txt $db/hg19_refGeneMrna.fa > $working_dir/$tag/$tag.coding\
change

perl $annovar_dir/coding_change.pl --includesnp $avinputfile.refGene.exonic_variant_function $db/hg19_refGene.txt $db/hg19_refGeneMrna.fa > $working_dir/$tag/$tag.merged.codingch\
ange

perl $annovar_dir/coding_change.pl --includesnp $germlineavinputfile.refGene.exonic_variant_function $db/hg19_refGene.txt $db/hg19_refGeneMrna.fa > $working_dir/$tag/$tag.germlin\
e.codingchange

#cp $working_dir/$tag/$tag.*.avinput.hg19_multianno.txt $variant_dir/
}
export -f run_annovar

annotate_annovar() {
parallel --jobs 4 run_annovar ::: $tumor_samples
}

cnv_analysis() {
#create sample file with switch - for single vs multiplex; creates sample file rather than just variable for future scripts
mkdir -p "$starting_dir/exome_pipeline/CNVkit"
mkdir -p "$results_folder/cnv"
cd $cnv_dir

 pooled_reference=$scripts_dir/dependencies/Reference.cnn
 tumorbam=$bam_files/*/*-T*.bam
 normalbam=$bam_files/*/*-N*.bam
 allbams=$bam_files/*/*.bam

 $cnvkit_script autobin $allbams -t $selector_padded_bed -g $cnv_access_file --annotate $refFlat --short-names #AUTOBIN

 autobins=${selector_padded_bed##*/}

            if [ $CNV_target_antitarget == "on" ] #turn on for first time run
                then
                for i in $allbams
                do
                    samp=${i##*/}
                    $cnvkit_script coverage $i ${autobins%.*}.target.bed -o ${samp%%.*}.targetcoverage.cnn
                    $cnvkit_script coverage $i ${autobins%.*}.antitarget.bed -o ${samp%%.*}.antitargetcoverage.cnn
                done
                fi

            if [ $CNV_reference_switch == "on" ] #turn on for first time run
                then
                $cnvkit_script reference *-N*coverage.cnn -f $reference_fasta -o $scripts_dir/dependencies/Reference.cnn
                fi

        pooled_reference=$scripts_dir/dependencies/Reference.cnn

        for i in $allbams
        do
                samp=${i##*/}
                echo $i
                $cnvkit_script fix ${samp%%.*}.targetcoverage.cnn ${samp%%.*}.antitargetcoverage.cnn $pooled_reference -o ${samp%%.*}.cnr
                $cnvkit_script segment ${samp%%.*}.cnr -o ${samp%%.*}.cns
                $cnvkit_script scatter ${samp%%.*}.cnr -s ${samp%%.*}.cns -o ${samp%%.*}.scatter.pdf
                $cnvkit_script diagram -s ${samp%%.*}.cns ${samp%%.*}.cnr
                $cnvkit_script genemetrics ${samp%%.*}.cnr > ${samp%%.*}.genemetrics_no_cns.txt
                $cnvkit_script genemetrics ${samp%%.*}.cnr -s ${samp%%.*}.cns > ${samp%%.*}.genemetrics.txt
                $cnvkit_script call ${samp%%.*}.cns -m threshold -t=-1.1,-0.4,0.3,0.7 -o ${samp%%.*}.call.cns
                $cnvkit_script breaks ${samp%%.*}.cnr ${samp%%.*}.cns > ${samp%%.*}.gene.breaks.txt
                $cnvkit_script scatter -s ${samp%%.*}.cns -c chr1 -g Phlda3 -o ${samp%%.*}.chr1.scatter.pdf

		cp ${samp%%.*}.call.cns $results_folder/cnv
                cp ${samp%%.*}.scatter.pdf $results_folder/cnv
                cp ${samp%%.*}.cns $results_folder/cnv
                cp ${samp%%.*}.cnr $results_folder/cnv
                cp ${samp%%.*}.genemetrics_no_cns.txt $results_folder/cnv
                cp ${samp%%.*}.genemetrics.txt $results_folder/cnv
                cp ${samp%%.*}.gene.breaks.txt $results_folder/cnv
                cp ${samp%%.*}.chr1.scatter $results_folder/cnv
        done

        $cnvkit_script heatmap *.cnr -d -o All.pdf
        $cnvkit_script heatmap *.cnr -d -c chr1 -o All.chr1.pdf

}

annotate_depth() {
mkdir -p $results_folder/variants
for tag in $tumor_samples
do
    
    cd $working_dir/$tag

    if [ $indel_mode == "on" ]
    then
	#Add in indels from annovar
	echo "Exporting annotated INDELS from nondeduped BAMS"
	Rscript $scripts_dir/ExportIndels.R $tag.consensus.avinput.hg19_multianno.txt
	awk -v OFS="\t" '{print $0,"0","30","2","30"}' $tag.consensus.avinput.hg19_multianno.txt.indels > $tag.consensus.avinput.hg19_multianno.txt.indels.anno
	cp $tag.consensus.avinput.hg19_multianno.txt.indels.anno $working_dir/variants/$tag.consensus.avinput.hg19_multianno.annoGermline.annoTumor.txt
	echo "Done annotating INDELS"
    fi

    if [ $indel_mode != "on" ]
    then
    cp $tag.consensus.avinput.hg19_multianno.txt $variant_dir/
    echo "Generating BED file"
    cat $tag.strelka.raw.vcf | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k2,2n"}' > $tag.strelka.sorted.vcf
    convert2bed -i vcf < $tag.strelka.sorted.vcf > $tag.strelka.bed
    convert2bed -i vcf < $tag.mutect.passed.vcf > $tag.mutect.bed
    convert2bed -i vcf < $tag.varscan.filtered.raw.vcf > $tag.varscan.bed
    bedops --everything $tag.{mutect,strelka,varscan}.bed > $tag.consensus.bed

    echo "Done generating BED file"
    tumorbam=$bam_files/$tag/${tag}-T*.bam
    normalbam=$bam_files/$tag/${tag}-N*.bam

    echo "Starting FREQ Generation"

    #Generate Tumor Freq
    python $bam2freq_script $tumorbam $reference_fasta $tag.consensus.bed 
    #Generate Normal Freq
    python $bam2freq_script $normalbam $reference_fasta $tag.consensus.bed
    
    echo "done generating freqs"


    #Copying and renaming FREQ files because the annotate_vars script expects it in this format
    cp $bam_files/$tag/${tag}-T*freq.paired.Q30.txt $working_dir/$tag/${tag}-T.freq.paired.Q30.txt
    cp $bam_files/$tag/${tag}-N*freq.paired.Q30.txt $working_dir/$tag/${tag}-N.freq.paired.Q30.txt
    
    echo "annotating read depth"
    #Annotate with Normal Reads
    cd $working_dir/variants/
    python $annotate_vars_script $tag.consensus.avinput.hg19_multianno.txt germline $working_dir/$tag/
    python $annotate_vars_script $tag.consensus.avinput.hg19_multianno.annoGermline.txt tumor $working_dir/$tag/ 
    echo "done annotating read depth"
    fi

    echo "adding repeatmasker field"
    #Add repeatmasker field
    if [ $indel_mode != "on" ] #added this to remove header of the annoGermline.annoTumor.txt file since update to sherlock, only needed for snv mode
    then
    cp $working_dir/variants/$tag.consensus.avinput.hg19_multianno.annoGermline.annoTumor.txt $working_dir/variants/$tag.consensus.avinput.hg19_multianno.annoGermline.annoTumor.snv.txt
    tail -n +2 $working_dir/variants/$tag.consensus.avinput.hg19_multianno.annoGermline.annoTumor.snv.txt > $working_dir/variants/$tag.consensus.avinput.hg19_multianno.annoGermline.annoTumor.txt
    fi 
    bedtools intersect -loj -a $working_dir/variants/$tag.consensus.avinput.hg19_multianno.annoGermline.annoTumor.txt -b $repeatmasker > $working_dir/variants/$tag.consensus.annotated.masked.txt
    echo "done adding repeatmasker"

    #echo "adding genomad POP AF"
    #Add genomad information
    #Rscript $annotate_genomad_script $working_dir/variants/$tag.consensus.annotated.masked.txt $genomad $working_dir/variants/$tag.consensus.final.txt
    #echo "done with all annotations"
    
    #i am doing the genomad part inside the if instead so that if nothing is called in the indel step, it does not just write the old consensus.final.txt file from SNV
    if [ $indel_mode != "on" ]
    then    
    echo "adding genomad POP AF"
    #Add genomad information
    Rscript $annotate_genomad_script $working_dir/variants/$tag.consensus.annotated.masked.txt $genomad $working_dir/variants/$tag.consensus.final.txt
    echo "done with all annotations"
    #Copy file to results folder
    mkdir -p $results_folder
    cp $working_dir/variants/$tag.consensus.final.txt $results_folder/variants/$tag.variant.results
    cp $working_dir/variants/$tag.vcf.venn.png $results_folder/variants/$tag.venn.png
    cp $working_dir/variants/$tag.vcf.scatter.png $results_folder/variants/$tag.scatter.png
    cp $working_dir/variants/$tag.consensus.final.txt $working_dir/variants/$tag.variant.snv.results
    cp $working_dir/variants/$tag.vcf.venn.png $working_dir/variants/$tag.venn.snv.png
    cp $working_dir/variants/$tag.vcf.scatter.png $working_dir/variants/$tag.scatter.snv.png
    fi
    
    if [ $indel_mode == "on" ]
    then
    echo "adding genomad POP AF"
    #Add genomad information
    Rscript $annotate_genomad_script $working_dir/variants/$tag.consensus.annotated.masked.txt $genomad $working_dir/variants/$tag.consensus.indel.final.txt
    echo "done with all annotations"
    #copy file to results folder
    mkdir -p $results_folder
    cp $working_dir/variants/$tag.consensus.indel.final.txt $results_folder/variants/$tag.variant.indels.results
    cp $working_dir/variants/$tag.vcf.venn.png $results_folder/variants/$tag.venn.indels.png
    cp $working_dir/variants/$tag.vcf.scatter.png $results_folder/variants/$tag.scatter.indels.png
    cp $working_dir/variants/$tag.consensus.indel.final.txt $working_dir/variants/$tag.variant.indels.results
    cp $working_dir/variants/$tag.vcf.venn.png $working_dir/variants/$tag.venn.indels.png
    cp $working_dir/variants/$tag.vcf.scatter.png $working_dir/variants/$tag.scatter.indels.png    
fi

done
}

filter_and_plot() {
Rscript $filter_and_plot_script $results_folder/variants
cd $results_folder/variants
mkdir -p $results_folder/variants/sorted_beds/
for i in *filtered.bed
do
sort -V -k1,1 -k2,2 $i > sorted_beds/${i%%.*}.filtered.sorted.bed
done
}

ABSOLUTE() {
mkdir -p $results_folder/ABSOLUTE
Rscript $run_absolute $results_folder/variants/ $results_folder/cnv/ $results_folder/ABSOLUTE/
}

##########################Script starts here########################################
start=$(date)
start_time=$(date +%s)
runid=$(date | awk '{print $2$3$NF$4}')
samplefile=samples.$runid.txt
logfile=log.$runid.txt

create_dir_structure #creates directory structure
echo -e "Start time $start" > $working_dir/$logfile

if [ $main_pipeline == "on" ]
       then
    
       indel_mode="off"
       SECONDS=0
       echo "Starting Varscan variant calling" | tee -a $working_dir/$logfile
       call_varscan_variants  #calls variants using varscan plus realign and fpf filter
       varscan_duration=$(displaytime $SECONDS)
       echo "Varscan duration=$varscan_duration" | tee -a $working_dir/$logfile

       SECONDS=0
       echo "Starting Mutect variant calling" | tee -a $working_dir/$logfile
       call_mutect_variants #calls variants using mutect  
       mutect_duration=$(displaytime $SECONDS)
       echo "Mutect duration=$mutect_duration" | tee -a $working_dir/$logfile

       SECONDS=0
       echo "Starting Strelka variant calling" | tee -a $working_dir/$logfile
       call_strelka_variants #calls variants using strelka
       strelka_duration=$(displaytime $SECONDS)
       echo "Strelka duration=$strelka_duration" | tee -a $working_dir/$logfile

       SECONDS=0
       echo "Calling consensus_variants"
       consensus_variants #add VAFs, reformat variant call outputs.  R script to determine consensus calls, create the annovar input files
       echo "finished consensus_variants"
       consensus_duration=$(displaytime $SECONDS) 
       echo "Consensus calling duration=$consensus_duration" | tee -a $working_dir/$logfile

       SECONDS=0
       annotate_annovar #annotate all calls and consensus calls
       annovar_duration=$(displaytime $SECONDS)
       echo "Annovar duration=$annovar_duration" | tee -a $working_dir/$logfile

       SECONDS=0
       annotate_depth
       echo "Done annotating tumor and germline depths to variant output"
       annotate_depth_dur=$(displaytime $SECONDS)
       echo "annotate depth duration=$annotate_depth_dur" | tee -a $working_dir/$logfile   

       #Begin INDEL portion of main pipeline where bam_files variable is switched to the nondeduped bam file folder, and runfolder variable is modified to indicate nondeduped bams. indel mode switch is also changed.

       if [ $indel_switch == "on" ]
       then
	   
	   bam_files=$nondeduped_bam_files
	   $2=$2_nondeduped
	   indel_mode="on"

       SECONDS=0
       echo "Starting Varscan indel variant calling" | tee -a $working_dir/$logfile
       call_varscan_variants  #calls variants using varscan plus realign and fpf filter
       varscan_duration=$(displaytime $SECONDS)
       echo "Varscan duration=$varscan_duration" | tee -a $working_dir/$logfile

       SECONDS=0
       echo "Starting Mutect indel variant calling" | tee -a $working_dir/$logfile
       call_mutect_variants #calls variants using mutect
       mutect_duration=$(displaytime $SECONDS)
       echo "Mutect duration=$mutect_duration" | tee -a $working_dir/$logfile

       SECONDS=0
       echo "Starting Strelka indel variant calling" | tee -a $working_dir/$logfile
       call_strelka_variants #calls variants using strelka
       strelka_duration=$(displaytime $SECONDS)
       echo "Strelka duration=$strelka_duration" | tee -a $working_dir/$logfile

        SECONDS=0
       echo "Calling consensus_variants"
       consensus_variants #add VAFs, reformat variant call outputs.  R script to determine consensus calls, create the annovar input files
       echo "finished consensus_variants"
       consensus_duration=$(displaytime $SECONDS)
       echo "Consensus calling duration=$consensus_duration" | tee -a $working_dir/$logfile

       SECONDS=0
       annotate_annovar #annotate all calls and consensus calls
       annovar_duration=$(displaytime $SECONDS)
       echo "Annovar duration=$annovar_duration" | tee -a $working_dir/$logfile

       SECONDS=0
       annotate_depth
       echo "Done annotating tumor and germline depths to variant output"
       annotate_depth_dur=$(displaytime $SECONDS)
       echo "annotate depth duration=$annotate_depth_dur" | tee -a $working_dir/$logfile
       fi
fi

SECONDS=0
if [ $filter_and_plot == "on" ]
      then
      filter_and_plot
      echo "Done Filtering and Plotting"
      filtration_duration=$(displaytime $SECONDS)
      echo "filtration duration=$filtration_duration" | tee -a $working_dir/$logfile
fi

SECONDS=0
if [ $CNV_kit_switch == "on" ]
	then
	echo "Detecting copy number variation" | tee -a $working_dir/$logfile
	cnv_analysis #run cnvkit of each tumor against pooled germlines background
	cnv_duration=$(displaytime $SECONDS)
	echo "Copy number calling duration=$cnv_duration" | tee -a $working_dir/$logfile
fi

SECONDS=0
if [ $Absolute_switch == "on" ]
    then
    echo "Running ABSOLUTE" | tee -a $working_dir/$logfile
    ABSOLUTE #generate potential ploidy/purity models
    ABSOLUTE_diration=$(displaytime $SECONDS)
    echo "run ABSOLUTE duration=$ABSOLUTE_duration" | tee -a $working_dir/$logfile
fi

#####End script
stop=$(date)
echo -e "Completed time $stop" >> $working_dir/$logfile
stop_time=$(date +%s)
duration=$(($stop_time-$start_time))
total_duration=$(displaytime $duration)
echo "Total pipeline duration=$total_duration" | tee -a $working_dir/$logfile
exit
