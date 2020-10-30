#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 11 12:07:24 2020

@author: shihyu
"""


import pandas as pd

CLL_WGS_list=["CLL001","CLL002","CLL003","CLL004","CLL005","CLL006","CLL007","CLL008","CLL009","CLL010","CLL011","CLL012","CLL013","CLL014","CLL015","CLL016","CLL017","CLL018","CLL019","CLL020"]
Targets=pd.read_excel('CLLGeneList_STRING.xlsx',sheet_name='Union', header=None).rename(columns={0:'Genes'})
TargetList=Targets['Genes'].to_list()
A=pd.read_excel('GeneBasedAry_CLLWGS_AFfilter_PASS.xlsx')
targetA=A[A['Gene name'].isin(TargetList)]
targetA.to_excel('recurrenceSNPinSV.xlsx')


# =============================================================================
# IDlist=[]
# Genelist=[]
# 
# for i in CLL_WGS_list:
#     A=pd.read_csv(f"{i}_sv.sorted.vcf.annotSV.output.tsv",sep="\t").fillna(0)
#     ary=A[A['SV type']!='BND']
#     ary_test=A[A['SV length']==0]
#     ary_split=ary[ary['AnnotSV type'] == "split"]
#     ary_full=ary[ary['AnnotSV type'] == "full"]
#     IDlist.extend(ary_full['AnnotSV ID'].tolist())
#     Genelist.extend(ary_split['Gene name'].tolist())
#     
# IDset=pd.DataFrame(sorted(set(IDlist)),columns=['ID']).set_index('ID')
# Geneset=pd.DataFrame(sorted(set(Genelist)),columns=['Gene']).set_index('Gene')
# =============================================================================

# =============================================================================
# Make example for full/split data type
# A=pd.read_csv("CLL020_sv.sorted.vcf.annotSV.output.tsv",sep="\t").set_index('AnnotSV ID')
# ary=A[A['SV type'] != 'BND']
# Test=ary.loc['13_49997960_51050950_DEL_1'][['AnnotSV type', 'Gene name', 'GD_AF', 'GD_POPMAX_AF', 'location', 'location2', 'AnnotSV ranking']]
# =============================================================================


#Make ID based array

#Make union of unique IDs with necessary infos
for i in CLL_WGS_list:
    A=pd.read_csv(f"{i}_sv.sorted.vcf.annotSV.output.tsv",sep="\t").fillna('Nan').set_index('AnnotSV ID', drop=False)
    ary=A[(A['AnnotSV type'] == "full") &(A['SV type']!='BND')]
    if i == 'CLL001':
        Demo=ary[['AnnotSV ID','AnnotSV type', 'Gene name', 'GD_AF', 'GD_POPMAX_AF', 'location', 'location2', 'AnnotSV ranking']]

    else:
        Union_ID=pd.concat([Demo,ary[['AnnotSV ID','AnnotSV type', 'Gene name', 'GD_AF', 'GD_POPMAX_AF', 'location', 'location2', 'AnnotSV ranking']]]).drop_duplicates()
        Demo=Union_ID.copy()

#Create recurrence ID table
for i in CLL_WGS_list:
    S=pd.read_csv(f"{i}_sv.sorted.vcf.annotSV.output.tsv",sep="\t",index_col="AnnotSV ID")
    S_full=S[(S['AnnotSV type'] == "full") &(S['SV type']!='BND')]
    VariantBasedAry=Union_ID.join(S_full['INFO']).rename(columns={'INFO':i})
    Union_ID=VariantBasedAry.copy()

VariantBasedAry['Counts']=VariantBasedAry.iloc[:,8:28].count(axis=1)
VariantBasedAry.to_excel("IDBasedAry_CLLWGS_nofilter_SV.xlsx")





#Make Gene based array

#Make union of unique genes
for i in CLL_WGS_list:
    A=pd.read_csv(f"{i}_sv.sorted.vcf.annotSV.output.tsv",sep="\t").fillna('Nan').set_index('AnnotSV ID', drop=False)
    ary=A[(A['AnnotSV type'] == "split") &(A['SV type']!='BND')]
    if i == 'CLL001':
        Demo=ary[['AnnotSV ID','Gene name']]

    else:
        Union_GN=pd.concat([Demo,ary[['AnnotSV ID','Gene name']]]).drop_duplicates()
        Demo=Union_GN.copy()

#Add column of recurrence IDs to each gene
Geneset=Union_GN.groupby(['Gene name'])['AnnotSV ID'].apply(lambda x: ','.join(x.astype(str))).reset_index()
df_Gene=pd.concat([Geneset,pd.DataFrame(columns=CLL_WGS_list)],axis=1).set_index('Gene name')

for i in CLL_WGS_list:
    S=pd.read_csv(f"{i}_sv.sorted.vcf.annotSV.output.tsv",sep="\t",index_col="Gene name").fillna(0)
    filter1=(S['GD_AF']<0.01)&(S['GD_POPMAX_AF']<0.01)
    filter2=(S['FILTER']=='PASS')
    S_split=S[(S['AnnotSV type'] == "split") & (S['SV type']!='BND')&filter1&filter2] 
    for j in S_split.index:
        df_Gene.at[j,i]=S_split.index.value_counts()[j]

df_Gene['Sample Counts']=df_Gene.iloc[:,1:21].count(axis=1)
df_Gene['Variant Counts']=df_Gene.iloc[:,1:21].sum(axis=1)
df_Gene_final=df_Gene[df_Gene['Variant Counts']>0]
df_Gene_final.to_excel('GeneBasedAry_CLLWGS_AFfilter_PASS.xlsx')









    