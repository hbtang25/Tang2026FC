[TOC]

# 0 Set environment variables

    db=/opt/database/meta
    soft=/opt/miniforge3
    script=/opt/script
    
    wd=~/koji/WGS
    mkdir -p $wd 
    cd $wd

    conda activate meta_py3.9


# 1 Working directory and initial files

    # Create directories
    mkdir -p seq result temp
    

# 2 FastQC quality assessment
    
    mkdir -p result/01_raw_qc
    fastqc seq/*.gz -o result/01_raw_qc -t 64
    
    # Generate comparative reports for multiple samples
    mkdir -p result/02_raw_multiqc
    multiqc -f -d result/01_raw_qc -o result/02_raw_multiqc



# 3 Sequence quality control and host removal Qaulity control and remove host contamination
    
    # Activate the environment
    conda activate meta_kneaddata
    mkdir -p result/03_kneaddata
    
    parallel -j 2 --xapply \
      "kneaddata -i {1} -i {2} \
      -o result/03_kneaddata -v -t 16 --remove-intermediate-output \
      --trimmomatic ${soft}/envs/meta_kneaddata/share/trimmomatic/ \
      --trimmomatic-options 'ILLUMINACLIP:${soft}/envs/meta_kneaddata/share/trimmomatic/adapters/TruSeq3-PE.fa:2:40:15 SLIDINGWINDOW:4:20 MINLEN:50' \
      --bowtie2-options '--very-sensitive --dovetail' \
      -db ${db}/kneaddata/human_genome/Homo_sapiens \
      -db ${db}/kneaddata/glycine_max/glycine_max" \
      ::: seq/*_1.fq.gz ::: seq/*_2.fq.gz > temp/kneaddata.log 2>&1
      
    # Summarize quality control results x
    mkdir -p result/03_kneaddata
    kneaddata_read_count_table --input result/03_kneaddata --output result/03_kneaddata/kneaddata_sum.txt
    cat result/03_kneaddata/kneaddata_sum.txt 


# 4 Quality reassessment after quality control

    conda activate meta_py3.9
    mkdir -p result/04_trim_qc
    fastqc result/03_kneaddata/*_1_kneaddata_paired_* -o result/04_trim_qc -t 4

    mkdir -p result/05_trim_multiqc
    multiqc -f -d result/04_trim_qc result/03_kneaddata -o result/05_trim_multiqc
    


# 5 Reference-based metagenomic analysis workflow Reference-based Metagenomic Pipeline (HUMAnN2)

## 5.1 Merge quality-controlled files as HUMAnN2 input

    mkdir -p temp/concat
    
    parallel -j 4 --xapply \
      "cat result/03_kneaddata/{1}*_1_kneaddata_paired_?.fastq > temp/concat/{1}.fq" \
      ::: `tail -n+2 result/metadata.tsv | cut -f 1`
    

    # Check the number and size of sample files
    ls -sh temp/concat/*.fq 


## 5.2 Calculate taxonomic and functional profiles using HUMAnN2

    # Taxonomic profiling uses MetaPhlAn2, with bowtie2 alignment against nucleotide sequences;
    conda activate meta_humann2
    mkdir -p temp/humann2
    
    # Run multiple tasks simultaneously if sufficient disk space is available
    parallel -j 2 --xapply \
      "humann2 --input {} \
        --output temp/humann2/ --threads 16" \
    ::: temp/concat/*.fq > temp/humann2/humann2.log 2>&1

    # Copy important files to the humann2 directory
    for i in `tail -n+2 result/metadata.tsv|cut -f1`;do
       cp temp/humann2/${i}_humann2_temp/${i}_metaphlan_bugs_list.tsv temp/humann2/
    done

    # Delete temporary files
    rm temp/concat/* temp/humann2/*_humann2_temp

  
## 5.3 Taxonomic composition table
### 5.3.1 Merge sample results

    mkdir -p result/06_metaphlan2
    merge_metaphlan_tables.py temp/humann2/*_metaphlan_bugs_list.tsv | \
      sed 's/_metaphlan_bugs_list//g' > result/06_metaphlan2/taxonomy.tsv


### 5.3.2 Convert to spf format for stamp analysis

    perl ${script}/metaphlan_to_stamp.pl result/06_metaphlan2/taxonomy.tsv > result/06_metaphlan2/taxonomy.spf
    head -n6 result/06_metaphlan2/taxonomy.spf

### 5.3.3 Generate heatmaps

    metaphlan_hclust_heatmap.py --in result/06_metaphlan2/taxonomy_jiangmen.tsv \
        --out result/06_metaphlan2/heatmap_s_jiangmen.pdf \
        -c Purples --tax_lev s --top 5 --minv 0.1 -s log -x 0.6 -y 0.4

    metaphlan_hclust_heatmap.py --in result/06_metaphlan2/taxonomy.tsv \
        --out result/06_metaphlan2/heatmap_s.pdf \
        -c Purples --tax_lev s --top 3 --minv 0.1 -s log -x 0.6 -y 0.4
    

    cat result/06_metaphlan2/taxonomy.tsv |grep 's__' |grep -v -E '\|t__|unclassified' |sed 's/^.*s__//' \
      |sed "1i `sed -n 1p result/06_metaphlan2/taxonomy.tsv`" \
      |sed "2i `sed -n 2p result/06_metaphlan2/taxonomy.tsv`" >result/06_metaphlan2/taxonomy_s.tsv
    
    
    cat result/06_metaphlan2/taxonomy.tsv |grep 'p__' |grep -v '|c__' |sed 's/^.*p__//' \
      |sed "1i `sed -n 1p result/06_metaphlan2/taxonomy.tsv`" \
      |sed "2i `sed -n 2p result/06_metaphlan2/taxonomy.tsv`" >result/06_metaphlan2/taxonomy_p.tsv

    # LEfSe
    conda activate lefse
   
    # Convert to the internal lefse format
    lefse-format_input.py LEfSe.txt input.in -c 1 -o 1000000
    # Run lefse
    run_lefse.py input.in input.res -l 4
    cat input.res |grep 's__' |grep -v 't__' > input_Species.res
    
    # Generate a taxonomic cladogram highlighting differential taxa
    lefse-plot_cladogram.py input.res cladogram.pdf --format pdf
    # lefse-plot_cladogram.py input_Species.res cladogram_Species.pdf --format pdf
    
    # Generate bar plots of all differential features
    # lefse-plot_res.py input.res res.pdf --subclades -5 --format pdf
    lefse-plot_res.py input_Species.res res_Species.pdf  --subclades -5 --format pdf
  
    


## 5.4 Functional composition analysis

    conda activate meta_humann2
    mkdir -p ${wd}/result/07_humann2 && cd ${wd}/result/07_humann2


### 5.4.1 Merge, normalize, and stratify functional profiles

    humann2_join_tables --input ${wd}/temp/humann2 --file_name pathabundance \
      --output pathabundance.tsv
    
    # Remove redundant information from column names
    sed -i 's/_Abundance//g' pathabundance.tsv
     
    # Preview and summarize results
    grep -v -E 'UNMAPPED|UNINTEGRATED|unclassified' pathabundance.tsv \
      > pathabundance_grep.tsv

    # Normalize to relative abundance (relab) or copies per million (cpm)
    humann2_renorm_table --input pathabundance_grep.tsv \
      --units cpm \
      --output pathabundance_grep_cpm.tsv

    head -n5 pathabundance_grep_cpm.tsv

    # Stratified results: taxon-specific functional composition table (stratified) and overall functional composition table (unstratified)
    humann2_split_stratified_table --output ./ \
      --input pathabundance_grep_cpm.tsv 
    # Split stratified tables created:
    # result/humann2/pathabundance_relab_stratified.tsv
    # result/humann2/pathabundance_relab_unstratified.tsv

    # Overview of annotation results
      Rscript ${script}/mean_humann2.R --input pathabundance_grep_cpm_unstratified.tsv \
        --design ../metadata.tsv \
        --group Location --thre 0 \
        --output pathabundance_grep_cpm_unstratified_mean.txt
      
      tail -n+2 pathabundance_grep_cpm_unstratified_mean.txt |sort -t $'\t' -k2,2gr |sed "1i `sed -n 1p pathabundance_grep_cpm_unstratified_mean.txt`" \
        > pathabundance_grep_cpm_unstratified_mean_sort.txt
      rm pathabundance_grep_cpm_unstratified_mean.txt
  
      
      # 1. Generate an overall annotation plot
      Rscript ${script}/bar_humann2.R \
        -i pathabundance_grep_cpm_unstratified_mean_sort.txt \
        -g All -c pathabundance -n 10 --width 400 --height 200 -o ./



### 5.4.2 Convert to ec annotations
    
    for i in `tail -n+2 ${wd}/result/metadata.tsv|cut -f1`;do
      humann2_regroup_table \
      -i ${wd}/temp/humann2/${i}_genefamilies.tsv \
      -g uniref90_level4ec \
      -o ${wd}/temp/humann2/${i}_ec.tsv
    done
    # Merge results and correct sample names
    humann2_join_tables \
      --input ${wd}/temp/humann2/ \
      --file_name ec \
      --output ec.tsv
    
    
   grep -v -E 'UNMAPPED|UNGROUPED|unclassified' ec.tsv \
      > ec_grep.tsv

    # Normalize to relative abundance (relab) or copies per million (cpm)
    humann2_renorm_table --input ec_grep.tsv \
      --units cpm \
      --output ec_grep_cpm.tsv

    head -n5 ec_grep_cpm.tsv

    humann2_split_stratified_table --output ./ \
      --input ec_grep_cpm.tsv 

    head -n3 ec_grep_cpm_stratified.tsv
    head -n3 ec_grep_cpm_unstratified.tsv
    # Two output tables: stratified (functional composition of each taxon) and unstratified (overall functional composition)

    # Overview of annotation results
    Rscript ${script}/mean_humann2.R --input ec_grep_cpm_unstratified.tsv \
      --design ../metadata.tsv \
      --group Location --thre 0 \
      --output ec_grep_cpm_unstratified_mean.txt
      
    tail -n+2 ec_grep_cpm_unstratified_mean.txt |sort -t $'\t' -k2,2nr |sed "1i `sed -n 1p ec_grep_cpm_unstratified_mean.txt`" \
      > ec_grep_cpm_unstratified_mean_sort.txt
    rm ec_grep_cpm_unstratified_mean.txt

      

    # Contributions to individual enzyme EC numbers
    Rscript ${script}/mean_humann2.R --input ec_grep_cpm_stratified.tsv \
      --design ../metadata.tsv \
      --group Location --thre 0 \
      --output ec_grep_cpm_stratified_mean.txt


    for j in {3..5};do
      for i in `cat Enzyme_EC.txt |tail -n+2 |cut -f 2`;do
        group=`cat ec_grep_cpm_unstratified_mean_sort.txt |head -n1 |cut -f ${j}`
        cat <(cat ec_grep_cpm_unstratified_mean_sort.txt |grep -w "EC:${i}" |cut -f 1,${j}) \
          <(grep "EC:$i"\| ec_grep_cpm_stratified_mean.txt |sed 's/^EC.*s__//' \
            |cut -f 1,${j} |sort -k2,2nr |head -n 5 |sed 's/\t/(/;s/$/)/;s/_/ /'|tr '\n' ',' |sed 's/,/, /g;s/, $/\n/') \
          |tr '\n' '\t' |sed 's/\t$/\n/'   >> Enzyme_Contribution_${group}.txt
      done
    done


## 5.5 LEfSe differential abundance analysis and Cladogram

    conda activate lefse_v1.1
    mkdir -p temp/lefse

    # # Preview output data
    head -n3 result/06_metaphlan2/taxonomy.tsv
    # Extract the sample row, convert it to one sample per line, and rename ID to SampleID
    head -n1 result/06_metaphlan2/taxonomy.tsv|tr '\t' '\n'|sed '1 s/ID/SampleID/' >temp/lefse/sampleid
    head -n3 temp/lefse/sampleid

    awk 'BEGIN{OFS=FS="\t"}NR==FNR{a[$1]=$2}NR>FNR{print a[$1]}' result/metadata.tsv temp/lefse/sampleid|tr '\n' '\t'|sed 's/\t$/\n/' > temp/lefse/groupid
    cat temp/lefse/groupid

    awk 'BEGIN{OFS=FS="\t"}NR==FNR{a[$1]=$3}NR>FNR{print a[$1]}' result/metadata.tsv temp/lefse/sampleid|tr '\n' '\t'|sed 's/\t$/\n/' > temp/lefse/subgroupid
    cat temp/lefse/subgroupid

    # Merge grouping information and data
    cat temp/lefse/groupid temp/lefse/subgroupid <(head -n1 result/06_metaphlan2/taxonomy.tsv) <(tail -n+3 result/06_metaphlan2/taxonomy.tsv) > result/06_metaphlan2/lefse.txt
    head -n5 result/06_metaphlan2/lefse.txt


    # Convert to the internal lefse format
    lefse_format_input.py result/06_metaphlan2/lefse.txt temp/lefse/input.in -c 1 -s 2 -u 3 -o 1000000
    # Run lefse
    lefse_run.py temp/lefse/input.in temp/lefse/input.res -l 2
    cat temp/lefse/input.res |grep 's__*' |grep -v 't__' > temp/lefse/input_genus.res

    # Generate a taxonomic cladogram highlighting differential taxa
    lefse-plot_cladogram.py temp/lefse/input.res result/06_metaphlan2/lefse_cladogram.pdf --format pdf --dpi 600
    # Generate bar plots of all differential features
    lefse_plot_res.py temp/lefse/input.res result/06_metaphlan2/lefse_res.pdf --format pdf --dpi 600
    lefse_plot_res.py temp/lefse/input_genus.res result/06_metaphlan2/input_genus.res.pdf --format pdf

    # Batch-generate bar plots of all differential features
    mkdir -p result/06_metaphlan2/features
    lefse_plot_features.py -f diff --archive none --format pdf \
      temp/lefse/input.in temp/lefse/input_genus.res result/06_metaphlan2/features/
