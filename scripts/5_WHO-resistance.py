#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from functools import reduce
plt.rcParams["figure.figsize"] = (20,10)


#Check wheter a subs is Syn
def syn(s):
    import re
    temp = re.compile("([a-zA-Z]+)([0-9]+)([a-zA-Z]+)")
    res = temp.match(s).groups()
    if res[0] == res[2]:
        print(s)

#Check wheter a subs is NoSyn
def nsyn(s):
    import re
    temp = re.compile("([a-zA-Z]+)([0-9]+)([a-zA-Z]+)")
    res = temp.match(s).groups()
    if res[0] != res[2]:
        print(s)
        
def nsyn2(s):
    import re
    temp = re.compile("([a-zA-Z]+)([0-9]+)([a-zA-Z]+)")
    res1 = temp.match(s).group(1)
    res2 = temp.match(s).group(2)
    res3 = temp.match(s).group(3)
    if res1 != res3:
        print(s)
        
# Check whether there is STOP codon       
def nsyn3(s):
    import re
    temp = re.compile("([a-zA-Z_]+)([0-9]+)([a-zA-Z_]+)")
    #check whether s IS NOT white space
    if s and not s.isspace():
        res1 = temp.match(s).group(1)
        res2 = temp.match(s).group(2)
        res3 = temp.match(s).group(3)
        if res1 != res3:
            return(True)
        else:
            return(False)
    else:
        return(False)

# remove content between () in s
def cleanvar(s):
    import re
    return(re.sub(r"\([^()]*\)", "", s))  

# convert tab to space
def tab2space(s):
    import re
    return(s.replace('\t', ' '))

def removespace(s):
    import re
    return(s.replace(' ', ''))

def remove_(s):
    S=s.split('_')
    return(S[0])

def remove1_(s):
    S=s.split('_')
    return(S[1])

def first(s):
    return(s[0])

def last(s):
    return(s[-1])

def one2tree(s):
    from Bio.PDB.Polypeptide import one_to_three
    s=s.rstrip()
    if s.isupper() and not s.isnumeric() :
        try:
            return(one_to_three(s))
        except:
            return(s)
    else:
        return(s)
    
def codon(s):
    import re
    temp = re.findall(r'\d+', s)
    return(temp[0])


myfile = open("head", "r")
head = myfile.read()
head_list=head.split('\t')
#remove \n character
head_list[16]=head_list[16].strip()
head_list.insert(0,"File")
head_list.append('GeneName1')
head_list.append('Gene1')
head_list


df=pd.read_csv("pharma_gene.tab",sep='\t',header=None,names=head_list)


#clean substitution and add a new column SubstClean
df['SubstClean']=df.Subst.apply(cleanvar)
df['SubstClean2'] = df.SubstClean.apply(removespace)
df['ID']=df.File.apply(remove_)
df['NSYN'] = df.SubstClean2.apply(nsyn3)
df.rename({'#Pos': 'Genome position'}, axis=1,inplace=True)
df['genome_index']=df['Genome position']

whoG=pd.read_csv('WHO-UCN-GTB-PCI-2021.7-eng_genome_ind.csv',sep='\t')

whoG.genome_index=whoG.genome_index.str.split(',')
whoG = whoG.explode('genome_index').reset_index(drop=True)
whoG.genome_index=whoG.genome_index.astype(int)
whoG['Allel']=whoG['final_annotation.AlternativeNucleotide'].str.upper()
whoG['Ref']=whoG['final_annotation.ReferenceNucleotide'].str.upper()
whoG['common_name']=whoG['final_annotation.Gene'] + '_' + whoG['final_annotation.TentativeHGVSNucleotidicAnnotation']

df_whoG=pd.merge(df,whoG,on=['genome_index','Allel','Ref'])

df_whoG.sort_values('File',inplace=True)

df_whoG.filter(regex=r'(File|_Conf_Grade)')



# keep columns name by regex as array
drugs_conf=df_whoG.filter(regex=r'(_Conf_Grade)').columns
# create array of drugs
drugs=drugs_conf.str.split('_').str[0]
#drugs=[]
#for i in drugs_conf:
#    t=i.split('_')
#    drugs.append(t[0])

cutoff=10
D = {}
Dsubst = {}
Dsubst1 = {}
Dsubst5 = {}
li=[]
lisubst=[]
lisubst5=[]
lisubst1=[]
for i in drugs_conf:
    t = i.split('_')[0]
    D[t] = df_whoG[((df_whoG[i]=='1) Assoc w R') | (df_whoG[i]=='2) Assoc w R - Interim')) & (df_whoG.Freq>cutoff) & (df_whoG.Qual20>4)][['File',i]]
    Dsubst[t] = df_whoG[((df_whoG[i]=='1) Assoc w R') | (df_whoG[i]=='2) Assoc w R - Interim')) & (df_whoG.Freq>cutoff) & (df_whoG.Qual20>4)][['File','variant','common_name','Freq',i]]
    Dsubst1[t] = df_whoG[((df_whoG[i]=='1) Assoc w R') | (df_whoG[i]=='2) Assoc w R - Interim') | (df_whoG[i]=='3) Uncertain significance')) & (df_whoG.Freq>cutoff) & (df_whoG.Qual20>4)][['File','variant','common_name','Freq',i]]
    Dsubst5[t] = df_whoG[(df_whoG.Freq>cutoff) & (df_whoG.Qual20>4)][['File','variant','common_name','Freq',i]]
    li.append(D[t])
    lisubst.append(Dsubst[t])
    lisubst5.append(Dsubst5[t])
    lisubst1.append(Dsubst1[t])
#    t=df_whoG[df_whoG[i]=='1) Assoc w R'][['File',i]]


T5=pd.concat(lisubst5)

T5.fillna(0,inplace=True)


T5b=T5.groupby(['File','RIF_Conf_Grade','INH_Conf_Grade', 'EMB_Conf_Grade', 'PZA_Conf_Grade', 'LEV_Conf_Grade',       'MXF_Conf_Grade', 'BDQ_Conf_Grade', 'LZD_Conf_Grade', 'CFZ_Conf_Grade',       'DLM_Conf_Grade', 'AMI_Conf_Grade', 'STM_Conf_Grade', 'ETH_Conf_Grade',       'KAN_Conf_Grade', 'CAP_Conf_Grade'])['variant'].apply(', '.join).reset_index()


T5b.insert(1, "ID", T5b.File.str.split('_').str[0], True)
#T5b['ID']=T5b.File.str.split('_').str[0]
T5b['new'] = T5b['ID'].str.extract('(\d+)').astype(int)
#https://stackoverflow.com/questions/66134896/python-pandas-sort-an-alphanumeric-dataframe
T5b = T5b.sort_values(by=['new'], ascending=True).drop('new', axis=1)


T5b.to_excel('qcr26_full.xls',index=False)

myfile = open("head", "r")
head1 = myfile.read()
head_list1=head1.split('\t')
#remove \n character
head_list1[16]=head_list1[16].strip()
head_list1.insert(0,"File")
all=pd.read_csv('all_cf1N.tab',sep='\t',header=None,names=head_list1)


#clean substitution and add a new column SubstClean
all['SubstClean']=all.Subst.apply(cleanvar)
all['SubstClean2'] = all.SubstClean.apply(removespace)
all['ID']=all.File.apply(remove_)
all['NSYN'] = all.SubstClean2.apply(nsyn3)
all.rename({'#Pos': 'Genome position'}, axis=1,inplace=True)
all['genome_index']=all['Genome position']


all_whoG=pd.merge(all,whoG,on=['genome_index','Allel','Ref'])

all_whoG.sort_values('File',inplace=True)


all_whoG.filter(regex=r'(File|_Conf_Grade)')


cutoff=10
A = {}
Asubst = {}
Asubst1 = {}
Asubst5 = {}
Ali=[]
Alisubst=[]
Alisubst5=[]
Alisubst1=[]
for i in drugs_conf:
    t = i.split('_')[0]
    A[t] = all_whoG[((all_whoG[i]=='1) Assoc w R') | (all_whoG[i]=='2) Assoc w R - Interim')) & (all_whoG.Freq>cutoff) & (all_whoG.Qual20>4)][['File',i]]
    Asubst[t] = all_whoG[((all_whoG[i]=='1) Assoc w R') | (all_whoG[i]=='2) Assoc w R - Interim')) & (all_whoG.Freq>cutoff) & (all_whoG.Qual20>4)][['File','variant','common_name','Freq',i]]
    Asubst1[t] = all_whoG[((all_whoG[i]=='1) Assoc w R') | (all_whoG[i]=='2) Assoc w R - Interim') | (all_whoG[i]=='3) Uncertain significance')) & (all_whoG.Freq>cutoff) & (all_whoG.Qual20>4)][['File','variant','common_name','Freq',i]]
    Asubst5[t] = all_whoG[(all_whoG.Freq>cutoff) & (all_whoG.Qual20>4)][['File','variant','common_name','Freq',i]]
    Ali.append(D[t])
    Alisubst.append(Dsubst[t])
    Alisubst5.append(Dsubst5[t])
    Alisubst1.append(Dsubst1[t])
#    t=df_whoG[df_whoG[i]=='1) Assoc w R'][['File',i]]


A5=pd.concat(Alisubst5)

A5.fillna(0,inplace=True)


A5b=A5.groupby(['File','RIF_Conf_Grade','INH_Conf_Grade', 'EMB_Conf_Grade', 'PZA_Conf_Grade', 'LEV_Conf_Grade',       'MXF_Conf_Grade', 'BDQ_Conf_Grade', 'LZD_Conf_Grade', 'CFZ_Conf_Grade',       'DLM_Conf_Grade', 'AMI_Conf_Grade', 'STM_Conf_Grade', 'ETH_Conf_Grade',       'KAN_Conf_Grade', 'CAP_Conf_Grade'])['variant'].apply(', '.join).reset_index()


A5b.insert(1, "ID", A5b.File.str.split('_').str[0], True)
#T5b['ID']=T5b.File.str.split('_').str[0]
A5b['new'] = A5b['ID'].str.extract('(\d+)').astype(int)
#https://stackoverflow.com/questions/66134896/python-pandas-sort-an-alphanumeric-dataframe
A5b = A5b.sort_values(by=['new'], ascending=True).drop('new', axis=1)


A5b.to_excel('qcr26_all_full.xls',index=False)


