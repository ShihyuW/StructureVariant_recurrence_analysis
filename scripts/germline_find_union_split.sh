#!/bin/bash
#PBS -q ntu192G
#PBS -l select=1:ncpus=30
#PBS -P MST109178
#PBS -W group_list=MST109178
#PBS -N find_union_SV
#PBS -j oe
#PBS -M s0890003@gmail.com
#PBS -m e

OUTPUT_PATH=/work2/lynn88065/Annotation/20200414_NTUH_Deafness/20200608_Annotation/MDF016/SV_recurrence/test
cd $OUTPUT_PATH

# make the list for union file
#rm variants_files_list_sv.txt variants_files_list_cnv.txt variants_files_list_repeat
while read -r ID;do
echo "$OUTPUT_PATH/NTUH_Deafness_${ID}_sv.sorted.vcf.annotSV.output.tsv" >> variants_files_list_sv.txt
echo "$OUTPUT_PATH/NTUH_Deafness_${ID}_cnv.sorted.vcf.annotSV.output.tsv" >> variants_files_list_cnv.txt
echo "$OUTPUT_PATH/NTUH_Deafness_${ID}_repeats.sorted.vcf.annotSV.output.tsv" >> variants_files_list_repeats.txt
done< $OUTPUT_PATH/family_id.txt

variants_files_sv=$OUTPUT_PATH/variants_files_list_sv.txt
variants_files_cnv=$OUTPUT_PATH/variants_files_list_cnv.txt
variants_files_repeats=$OUTPUT_PATH/variants_files_list_repeats.txt

# Get header from one of the target file
# add variant_type header
head -n 1 $OUTPUT_PATH/NTUH_Deafness_MDD0026_hg38_sv.sorted.vcf.annotSV.output.tsv \
| awk -F "\t" '{print $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $11 "\t" $16 "\t" $17 "\t" $20 "\t" $21 "\t" $32 "\t" $36 "\t" $37 "\t" $77}' \
| awk -F "\t" 'BEGIN {OFS = "\t"} FNR==1{$(NF+1)="Variant_type"}1' > $OUTPUT_PATH/test/union_of_variants_annotsv.txt

# SV / CNV / Repeats
# Write every varints with info to a combine_variants_info.txt
# only export the split variants
while read -r list variant_type; do
#rm combined_variants_info.txt
while read -r annotsv_file; do
grep -P "split" ${annotsv_file} \
| awk -F "\t" '{print $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $11 "\t" $16 "\t" $17 "\t" $20 "\t" $21 "\t" $32 "\t" $36 "\t" $37 "\t" $77}' >> combined_variants_info.txt
done<$list

# add variant_type column and sorting
sort -n -k1 combined_variants_info.txt | uniq |awk -v t="${variant_type}" -F "\t" 'BEGIN {OFS = "\t"} FNR>0{$(NF+1)=t;} 1' > union_of_variants_${variant_type}.txt

# Get the union of all variants with information
less union_of_variants_${variant_type}.txt >> $OUTPUT_PATH/test/union_of_variants_annotsv.txt

done < $OUTPUT_PATH/files_list.txt
