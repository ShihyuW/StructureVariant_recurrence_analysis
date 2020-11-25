#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 11 12:07:24 2020

@author: shihyu
"""


import pandas as pd

sample_list=["somatic_test1","somatic_test2","somatic_test3"]

IDlist=[]
for i in sample_list:
    A=pd.read_csv(rf"inputs/{i}.tsv",sep="\t").fillna(0)
    ary=A[A['SV type']!='BND']
    ary_split=ary[ary['AnnotSV type'] == "split"]
    ary_full=ary[ary['AnnotSV type'] == "full"]
    IDlist.extend(ary_full['AnnotSV ID'].tolist())
    
IDset=pd.DataFrame(sorted(set(IDlist)),columns=['ID']).set_index('ID')



#Make ID based array

#Make union of unique IDs with necessary infos
for i in sample_list:
    A=pd.read_csv(rf"inputs/{i}.tsv",sep="\t").fillna('Nan').set_index('AnnotSV ID', drop=False)
    ary=A[(A['AnnotSV type'] == "full") &(A['SV type']!='BND')]
    if i == sample_list[0]:
        tmp=ary[['AnnotSV ID','AnnotSV type', 'Gene name', 'GD_AF', 'GD_POPMAX_AF', 'location', 'location2', 'AnnotSV ranking']]

    else:
        Union_ID=pd.concat([tmp,ary[['AnnotSV ID','AnnotSV type', 'Gene name', 'GD_AF', 'GD_POPMAX_AF', 'location', 'location2', 'AnnotSV ranking']]]).drop_duplicates()
        tmp=Union_ID.copy()

#Create recurrence ID table
for i in sample_list:
    S=pd.read_csv(rf"inputs/{i}.tsv",sep="\t",index_col="AnnotSV ID")
    S_full=S[(S['AnnotSV type'] == "full") &(S['SV type']!='BND')]
    VariantBasedAry=Union_ID.join(S_full['INFO']).rename(columns={'INFO':i})
    Union_ID=VariantBasedAry.copy()

VariantBasedAry['Counts']=VariantBasedAry.loc[:,sample_list].count(axis=1)
VariantBasedAry.to_excel("outputs/somatic_IDBasedAry_SV.xlsx")





#Make Gene based array

#Make union of unique genes
for i in sample_list:
    A=pd.read_csv(rf"inputs/{i}.tsv",sep="\t").fillna('Nan').set_index('AnnotSV ID', drop=False)
    ary=A[(A['AnnotSV type'] == "split") &(A['SV type']!='BND')]
    if i == sample_list[0]:
        tmp=ary[['AnnotSV ID','Gene name']]

    else:
        Union_GN=pd.concat([tmp,ary[['AnnotSV ID','Gene name']]]).drop_duplicates()
        tmp=Union_GN.copy()

#Add column of recurrence IDs to each gene
Geneset=Union_GN.groupby(['Gene name'])['AnnotSV ID'].apply(lambda x: ','.join(x.astype(str))).reset_index()
df_Gene=pd.concat([Geneset,pd.DataFrame(columns=sample_list)],axis=1).set_index('Gene name')

for i in sample_list:
    S=pd.read_csv(rf"inputs/{i}.tsv",sep="\t",index_col="Gene name").fillna(0)
    filter1=(S['GD_AF']<0.01)&(S['GD_POPMAX_AF']<0.01)
    filter2=(S['FILTER']=='PASS')
    S_split=S[(S['AnnotSV type'] == "split") & (S['SV type']!='BND')&filter1&filter2] 
    for j in S_split.index:
        df_Gene.at[j,i]=S_split.index.value_counts()[j]

df_Gene['Sample Counts']=df_Gene.loc[:,sample_list].count(axis=1)
df_Gene['Variant Counts']=df_Gene.loc[:,sample_list].sum(axis=1)
df_Gene_final=df_Gene[df_Gene['Variant Counts']>0]
df_Gene_final.to_excel('outputs/somatic_GeneBasedAry_SV.xlsx')









    