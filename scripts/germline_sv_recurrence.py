# -*- coding: utf-8 -*-
"""
Created on Thu Oct 29 00:47:22 2020

@author: Jacob-Lab
"""
import pandas as pd
import os as os

# set the input/output path 
input_path='D:\\Jacob-Lab\\github\\StructureVariant_recurrence_analysis\\inputs'
output_path='D:\\Jacob-Lab\\github\\StructureVariant_recurrence_analysis\\outputs'

# change the directory
os.chdir('{}'.format(input_path))

# family member AnnotSV file
Sample_list = ["germline_test01","germline_test02"]
union_file_name = 'germline_union_of_variants_annotsv.txt'
union_file = pd.read_csv('{}'.format(union_file_name),sep="\t")

# =============================================================================
# 
# # remove the full rows from raw data
# for i in Sample_list:
#     A = pd.read_csv(f"{i}.tsv", sep = "\t").fillna(0)
#     ary_split = A[A['AnnotSV type'] == 'split']
#     ary_split.to_csv(f"{i}_spilt.csv", sep = "\t" , index = False)
# =============================================================================

# create index
def create_idx(df):
    idx_list=[]
    L=df
    L["SV start"]=L["SV start"].apply(str)
    L["SV end"]=L["SV end"].apply(str)
  
    for i in range(len(L)):
        idx=L.loc[i]["SV chrom"]+"/"+L.loc[i]["SV start"]+"/"+L.loc[i]["SV end"]+"/"+L.loc[i]["Gene name"]
        idx_list.append(idx)
    idx_df=pd.DataFrame(idx_list,columns=["Idx"])
    finaldf=pd.concat([L,idx_df], axis=1)
    return finaldf

# generate info tags for N/T
# FORMAT (there are four format...)
# GT:SO:REPCN:REPCI:ADSP:ADFL:ADIR:LC
# GT:FT:GQ:PL:PR
# GT:FT:GQ:PL:PR:SR
# GT:SM:CN:BC:PE

N_info_tags = ['GT','SO','REPCN','REPCI','ADSP','ADFL','ADIR','LC']
for i in range(8):
    N_info_tags[i]+="_N"
    
# processing columns
for i in Sample_list:
    file = pd.read_csv(f"{i}.tsv", sep = "\t")
    normal_info = file[i]
    
    n_list=[]
    for j in range(len(normal_info)):
        n = normal_info[j].split(":")
        n_list.append(n)   
         
    df1=pd.DataFrame(n_list, columns=N_info_tags)
    final_Annotation_file=pd.concat([file,df1],axis=1)    
# create index on each annotation file
    create_idx(final_Annotation_file).to_excel('{}'.format(output_path) + f"\\Indexed_{i}.xlsx", index=False)

# create index on union file
create_idx(union_file).to_excel('{}'.format(output_path) + "\\Indexed_germline_union_of_variants_annotsv.xlsx", index=False)

# Generate PASSed and 3 filtered array by Idx
U=pd.read_excel('{}'.format(output_path) + "\\Indexed_germline_union_of_variants_annotsv.xlsx", index_col="Idx")[["SV chrom","SV start","SV end","SV length", "SV type", "Gene name" ,"location", "location2", "GD_ID", "GD_AF", "GD_POPMAX_AF","AnnotSV ranking","Variant_type"]]

# Process union file
for i in Sample_list:
    S=pd.read_excel('{}'.format(output_path) + f"\\Indexed_{i}.xlsx", index_col="Idx")
    
    #Filters (Add if needed)
    #filter1=(S['Func.refGene'].isin(["exonic","splicing","exonic;splicing"]))&(S['ExonicFunc.refGene'] != "synonymous SNV")
    #filter2=(S['AF_eas.1'].isin(["0","."]))&(S['TaiwanBiobank-official_Illumina1000-AF'].isin(["0","."]))
    #filter3=(S.AF_T >= 0.05)&(S.AF_T-S.AF_N >= 0.05)
    #GATK_filtered_S=S[(S.Otherinfo10 == "PASS")&(filter1&filter2&filter3)]
    
    ary=U.join(S["GT_N"]).rename(columns={'GT_N':f'{i}'})
    U=ary
  
ary["Counts"]=ary.iloc[:,13:15].count(axis=1)
ary = ary.fillna(".")
# check the count column
#print(ary.iloc[0,13:15])

# add "chr" to SV chrom column
ary["SV chrom"] = 'chr' + ary["SV chrom"].astype(str)

# format SV type column 
# <DEL> to DEL
# <DUP> to DUP
ary["SV type"] = ary["SV type"].str.replace("<DEL>","DEL").str.replace("<DUP>","DUP")

# format AnnotSV ranking column
replace_values = { 1 : "benign" , 2 : "likely-benign" , 3 : "VUS" , 4 : "likely-pathogenic", 5 : "pathogenic" }
ary = ary.replace({"AnnotSV ranking":replace_values})

# export the result to output_path
ary.to_excel('{}'.format(output_path) + "\\germline_sv_VaraintBasedArray.xlsx")
