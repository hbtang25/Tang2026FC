
[TOC]

# Set variables and environment
 conda activate tree
 wd=~/project/tree
 cd ${wd}

# bac
rm -rf track*_bac
# Generate default tree parameters; additional styles can be adjusted manually
cat cfg/global_bac.cfg tree2_label_color_bac.txt > track0_bac
# Merge all annotations; more tracks will be generated subsequently to make the tree more complex
cat track*_bac > graphlan_annotate_bac.txt

# Obtain the node IDs for the final visualization
cut -f 6 -d '.' tree1_backbone_bac.txt > tree1_backbone_bac_id
# Annotate the mean abundance of the results
awk 'BEGIN{OFS=FS="\t"} NR==FNR{a[$1]=$3} NR>FNR {print $1,a[$1]}' filtered_taxonomy_bac.txt tree1_backbone_bac_id > tree1_backbone_bac

# Annotate the results filtered at 0.5% as squares in ring 1; m represents magenta and R represents squares in cfg/ring1.cfg
cat cfg/ring1.cfg <(awk '$2=="Firmicutes"' tree1_backbone_bac | cut -f 1 | sed 's/$/\tring_shape\t1\tR/') \
    <(awk '$2=="Proteobacteria"' tree1_backbone_bac | cut -f 1 | sed 's/$/\tring_shape\t1\tv/') \
    <(awk '$2=="Firmicutes"' tree1_backbone_bac | cut -f 1 | sed 's/$/\tring_color\t1\tm/') \
    <(awk '$2=="Proteobacteria"' tree1_backbone_bac | cut -f 1 | sed 's/$/\tring_color\t1\ty/') > track1_bac


cat track*_bac > graphlan_annotate_bac.txt

# Add the mean abundance across all samples as a heatmap in ring 3
# Essentially, the heatmap is represented by the transparency of the circular bands
awk 'BEGIN{OFS=FS="\t"} NR==FNR{a[$1]=$3} NR>FNR {print $1,a[$1]}' otutab_mean_bac.txt tree1_backbone_bac_id > tree2_group1_mean_bac
awk 'BEGIN{OFS=FS="\t"} NR==FNR{a[$1]=$4} NR>FNR {print $1,a[$1]}' otutab_mean_bac.txt tree1_backbone_bac_id > tree2_group2_mean_bac
awk 'BEGIN{OFS=FS="\t"} NR==FNR{a[$1]=$5} NR>FNR {print $1,a[$1]}' otutab_mean_bac.txt tree1_backbone_bac_id > tree2_group3_mean_bac
awk 'BEGIN{OFS=FS="\t"} NR==FNR{a[$1]=$6} NR>FNR {print $1,a[$1]}' otutab_mean_bac.txt tree1_backbone_bac_id > tree2_group4_mean_bac
awk 'BEGIN{OFS=FS="\t"} NR==FNR{a[$1]=$7} NR>FNR {print $1,a[$1]}' otutab_mean_bac.txt tree1_backbone_bac_id > tree2_group5_mean_bac
awk 'BEGIN{OFS=FS="\t"} NR==FNR{a[$1]=$8} NR>FNR {print $1,a[$1]}' otutab_mean_bac.txt tree1_backbone_bac_id > tree2_group6_mean_bac
awk 'BEGIN{OFS=FS="\t"} NR==FNR{a[$1]=$2} NR>FNR {print $1,a[$1]}' otutab_mean_bac.txt tree1_backbone_bac_id > tree2_group0_mean_bac


awk 'BEGIN{OFS=FS="\t"} NR==FNR{a[$1]=$0} NR>FNR {print a[$1]}' filtered_taxonomy_bac.txt tree2_group0_mean_bac > filtered_taxonomy_bac_assigned.txt


# Annotate the results filtered at 0.5% as squares in ring 1; m represents magenta and R represents squares in cfg/ring1.cfg
cat cfg/bar2.cfg <(sed 's/\t/\tring_height\t2\t/g' tree2_group0_mean_bac) > track2_bac
cat cfg/heat3.cfg <(sed 's/\t/\tring_alpha\t3\t/g' tree2_group1_mean_bac) > track3_bac
cat cfg/heat4.cfg <(sed 's/\t/\tring_alpha\t4\t/g' tree2_group2_mean_bac) > track4_bac
cat cfg/heat5.cfg <(sed 's/\t/\tring_alpha\t5\t/g' tree2_group3_mean_bac) > track5_bac
cat cfg/heat6.cfg <(sed 's/\t/\tring_alpha\t6\t/g' tree2_group4_mean_bac) > track6_bac
cat cfg/heat7.cfg <(sed 's/\t/\tring_alpha\t7\t/g' tree2_group5_mean_bac) > track7_bac
cat cfg/heat8.cfg <(sed 's/\t/\tring_alpha\t8\t/g' tree2_group6_mean_bac) > track8_bac

  
# Plot the tree with rectangles in the first ring to display features with abundance greater than ten million
cat track*_bac > graphlan_annotate_bac.txt
graphlan_annotate.py --annot graphlan_annotate_bac.txt tree1_backbone_bac.txt graphlan_bac.xml
graphlan.py graphlan_bac.xml graphlan3_heatmap_bar_bac.pdf --size 4




# fun
rm -rf track*_fun
cat cfg/global_fun.cfg tree2_label_color_fun.txt > track0_fun
cat track*_fun > graphlan_annotate_fun.txt

cut -f 6 -d '.' tree1_backbone_fun.txt > tree1_backbone_fun_id
awk 'BEGIN{OFS=FS="\t"} NR==FNR{a[$1]=$3} NR>FNR {print $1,a[$1]}' filtered_taxonomy_fun.txt tree1_backbone_fun_id > tree1_backbone_fun


cat cfg/ring1.cfg <(awk '$2=="Ascomycota"' tree1_backbone_fun | cut -f 1 | sed 's/$/\tring_shape\t1\tR/') \
    <(awk '$2=="Basidiomycota"' tree1_backbone_fun | cut -f 1 | sed 's/$/\tring_shape\t1\tv/') \
    <(awk '$2=="Mortierellomycota"' tree1_backbone_fun | cut -f 1 | sed 's/$/\tring_shape\t1\t^/') \
    <(awk '$2=="Ascomycota"' tree1_backbone_fun | cut -f 1 | sed 's/$/\tring_color\t1\tm/') \
    <(awk '$2=="Basidiomycota"' tree1_backbone_fun | cut -f 1 | sed 's/$/\tring_color\t1\ty/') \
    <(awk '$2=="Mortierellomycota"' tree1_backbone_fun | cut -f 1 | sed 's/$/\tring_color\t1\tr/') > track1_fun


cat track*_fun > graphlan_annotate_fun.txt

awk 'BEGIN{OFS=FS="\t"} NR==FNR{a[$1]=$3} NR>FNR {print $1,a[$1]}' otutab_mean_fun.txt tree1_backbone_fun_id > tree2_group1_mean_fun
awk 'BEGIN{OFS=FS="\t"} NR==FNR{a[$1]=$4} NR>FNR {print $1,a[$1]}' otutab_mean_fun.txt tree1_backbone_fun_id > tree2_group2_mean_fun
awk 'BEGIN{OFS=FS="\t"} NR==FNR{a[$1]=$5} NR>FNR {print $1,a[$1]}' otutab_mean_fun.txt tree1_backbone_fun_id > tree2_group3_mean_fun
awk 'BEGIN{OFS=FS="\t"} NR==FNR{a[$1]=$6} NR>FNR {print $1,a[$1]}' otutab_mean_fun.txt tree1_backbone_fun_id > tree2_group4_mean_fun
awk 'BEGIN{OFS=FS="\t"} NR==FNR{a[$1]=$7} NR>FNR {print $1,a[$1]}' otutab_mean_fun.txt tree1_backbone_fun_id > tree2_group5_mean_fun
awk 'BEGIN{OFS=FS="\t"} NR==FNR{a[$1]=$8} NR>FNR {print $1,a[$1]}' otutab_mean_fun.txt tree1_backbone_fun_id > tree2_group6_mean_fun


awk 'BEGIN{OFS=FS="\t"} NR==FNR{a[$1]=$0} NR>FNR {print a[$1]}' filtered_taxonomy_fun.txt tree2_group0_mean_fun > filtered_taxonomy_fun_assigned.txt


cat cfg/bar2.cfg <(sed 's/\t/\tring_height\t2\t/g' tree2_group0_mean_fun) > track2_fun
cat cfg/heat3.cfg <(sed 's/\t/\tring_alpha\t3\t/g' tree2_group1_mean_fun) > track3_fun
cat cfg/heat4.cfg <(sed 's/\t/\tring_alpha\t4\t/g' tree2_group2_mean_fun) > track4_fun
cat cfg/heat5.cfg <(sed 's/\t/\tring_alpha\t5\t/g' tree2_group3_mean_fun) > track5_fun
cat cfg/heat6.cfg <(sed 's/\t/\tring_alpha\t6\t/g' tree2_group4_mean_fun) > track6_fun
cat cfg/heat7.cfg <(sed 's/\t/\tring_alpha\t7\t/g' tree2_group5_mean_fun) > track7_fun
cat cfg/heat8.cfg <(sed 's/\t/\tring_alpha\t8\t/g' tree2_group6_mean_fun) > track8_fun

  
cat track*_fun > graphlan_annotate_fun.txt
graphlan_annotate.py --annot graphlan_annotate_fun.txt tree1_backbone_fun.txt graphlan_fun.xml
graphlan.py graphlan_fun.xml graphlan3_heatmap_bar_fun.pdf --size 4
