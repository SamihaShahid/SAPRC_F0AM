#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 29 11:32:07 2022

@author: samiha
"""
import pandas as pd
import pubchempy as pcp
import numpy as np
import glob
import re
import matplotlib.pyplot as plt

sp=pd.read_excel('/Users/samiha/Desktop/chemical mechanism/Final/specdb.xlsx')
#___
sp['id']=sp['inchi']
sp=sp[sp['id'].notna()]

aa=pd.read_excel('/Users/samiha/Desktop/F0AM simulation/NMOG_dataset/andreae.xlsx')
aa=aa[aa['pollutant category']=='NMOG']

aa=aa.merge(sp[['id','S07','S18B']],on='id',how='left')

# pyrrole
aa.loc[81,'S07']='FURAN'
aa.loc[81,'S18B']='FURNS'

# convert ARO2 to FURAN
ind=list(aa[aa['S07']=='ARO2'].index)

for i in ind:
    aa.loc[i,'S07']='FURAN'


aa=aa[aa['S07'].notna()].reset_index(drop=True)

#S07 furan mass fraction

uu=list(aa['S07'].unique())

adf=pd.DataFrame()
adf['S07']=uu
for i in range(len(adf['S07'])):
    adf.loc[i,'mass']=aa['ef'][aa['S07']==adf['S07'][i]].sum()
    
bb=pd.read_excel('/Users/samiha/Desktop/F0AM simulation/NMOG_dataset/andreae_S07_unassigned_mass_sp.xlsx')

for i in range(len(bb)):
    ind=list(adf[adf['S07']==bb['S07'].iloc[i]].index)
    if len(ind)==1:
        adf.loc[ind[0],'mass']=adf['mass'].iloc[ind[0]]+bb['mass'].iloc[i]
    if len(ind)==0:
        adf.loc[-1,'S07']=bb['S07'].iloc[i]
        adf.loc[-1,'mass']=bb['mass'].iloc[i]
        adf=adf.reset_index(drop=True)

total_mass=adf['mass'].sum()
adf['mf']=adf['mass']/total_mass

#total VOC 91.55 ppb

adf['InitConc_ppb']=adf['mf']*(91.55)

inorganic=pd.read_excel('/Users/samiha/Desktop/F0AM simulation/NMOG_dataset/InitConc_inorganic.xlsx')

inorganic['S07']=inorganic['model_species']
df=inorganic[['S07','InitConc_ppb']].append(adf[['S07','InitConc_ppb']]).reset_index(drop=True)

for i in range(len(df)):
    df.loc[i,'hold me']=0


df.to_excel('/Users/samiha/Desktop/F0AM simulation/NMOG_dataset/aa_new/andreae_S07_furan_InitConc.xlsx',index=False)

# mass fraction S07 base

aa.loc[81,'S07']='ARO2'

uu=list(aa['S07'].unique())

adf=pd.DataFrame()
adf['S07']=uu
for i in range(len(adf['S07'])):
    adf.loc[i,'mass']=aa['ef'][aa['S07']==adf['S07'][i]].sum()
    
bb=pd.read_excel('/Users/samiha/Desktop/F0AM simulation/NMOG_dataset/andreae_S07_unassigned_mass_sp.xlsx')

for i in range(len(bb)):
    ind=list(adf[adf['S07']==bb['S07'].iloc[i]].index)
    if len(ind)==1:
        adf.loc[ind[0],'mass']=adf['mass'].iloc[ind[0]]+bb['mass'].iloc[i]
    if len(ind)==0:
        adf.loc[-1,'S07']=bb['S07'].iloc[i]
        adf.loc[-1,'mass']=bb['mass'].iloc[i]
        adf=adf.reset_index(drop=True)

total_mass=adf['mass'].sum()
adf['mf']=adf['mass']/total_mass

#total VOC 91.55 ppb

adf['InitConc_ppb']=adf['mf']*(91.55)

inorganic=pd.read_excel('/Users/samiha/Desktop/F0AM simulation/NMOG_dataset/InitConc_inorganic.xlsx')

inorganic['S07']=inorganic['model_species']
df=inorganic[['S07','InitConc_ppb']].append(adf[['S07','InitConc_ppb']]).reset_index(drop=True)

for i in range(len(df)):
    df.loc[i,'hold me']=0


df.to_excel('/Users/samiha/Desktop/F0AM simulation/NMOG_dataset/aa_new/andreae_S07_base_InitConc.xlsx',index=False)

# S18B mass fraction

uu=list(aa['S18B'].unique())

adf=pd.DataFrame()
adf['S18B']=uu
for i in range(len(adf['S18B'])):
    adf.loc[i,'mass']=aa['ef'][aa['S18B']==adf['S18B'][i]].sum()

bb=pd.read_excel('/Users/samiha/Desktop/F0AM simulation/NMOG_dataset/andreae_S18_unassigned_mass_sp.xlsx')

for i in range(len(bb)):
    ind=list(adf[adf['S18B']==bb['S18B'].iloc[i]].index)
    if len(ind)==1:
        adf.loc[ind[0],'mass']=adf['mass'].iloc[ind[0]]+bb['mass'].iloc[i]
    if len(ind)==0:
        adf.loc[-1,'S18B']=bb['S18B'].iloc[i]
        adf.loc[-1,'mass']=bb['mass'].iloc[i]
        adf=adf.reset_index(drop=True)


total_mass=adf['mass'].sum()
adf['mf']=adf['mass']/total_mass

adf['InitConc_ppb']=adf['mf']*(91.55)
inorganic['S18B']=inorganic['model_species']
df=inorganic[['S18B','InitConc_ppb']].append(adf[['S18B','InitConc_ppb']]).reset_index(drop=True)
for i in range(len(df)):
    df.loc[i,'hold me']=0

df.to_excel('/Users/samiha/Desktop/F0AM simulation/NMOG_dataset/aa_new/andreae_S18_InitConc.xlsx',index=False)


#NEIVA temperate forest

aa=pd.read_excel('/Users/samiha/Desktop/Supplement/temperate_forest/AVG dataset.xlsx')
aa=aa[aa['pollutant category']=='NMOG']

cm=pd.read_excel('/Users/samiha/Desktop/CM files/all_s07_include_hisomer_unassigned_added_v2.xlsx')

aa=aa.merge(cm[['S07T','id']],how='left')

ind=list(aa[aa['S07T'].isnull()][aa['compound'].str.contains('unknown')].index)
unk_mass=aa['AVG'][aa['S07T'].isnull()][aa['compound'].str.contains('unknown')].sum()

aa=aa.drop(index=ind).reset_index(drop=True)

sp=pd.read_excel('/Users/samiha/Desktop/chemical mechanism/Final/specdb.xlsx')

sp['id']=sp['inchi']

#fixing S07 and S07T
aa=aa.merge(sp[['id','S07']],on='id',how='left')

#check if they exist as S07 species
uu=list(aa['S07T'][aa['S07'].isnull()].unique())

for u in uu:
    if len(sp[sp['S07']==u])==0:
        print(u)
    
#13BDE == OLE2
#MXYL ==  ARO2
#SESQ == TERP

ind=list(aa[aa['S07T']=='13BDE'].index)

for i in ind:
    aa.loc[i,'S07']='OLE2'
    
ind=list(aa[aa['S07T']=='MXYL'].index)

for i in ind:
    aa.loc[i,'S07']='ARO2'
    
ind=list(aa[aa['S07T']=='SESQ'].index)

for i in ind:
    aa.loc[i,'S07']='TERP'

ind=list(aa[aa['S07'].isnull()][aa['S07T'].notna()].index)

for i in ind:
    aa.loc[i,'S07']=aa['S07T'].iloc[i]


unk=aa['AVG'][aa['S07'].isnull()].sum()

aa=aa[aa['S07'].notna()]

#mass fraction

uu=list(aa['S07'].unique())

ndf=pd.DataFrame()
ndf['S07']=uu
for i in range(len(ndf)):
    ndf.loc[i,'mass']=aa['AVG'][aa['S07']==ndf['S07'].iloc[i]].sum()
    
total_mass = ndf['mass'].sum()

ndf['mf']=ndf['mass']/total_mass

#total VOC 91.55 ppb (Muller et al., 2016)

ndf['InitConc']=ndf['mf']*(91.55)

ndf.to_excel('/Users/samiha/Desktop/F0AM simulation/NMOG_dataset/NEIVA_InitConc.xlsx',index=False)

# NEIVA S18 assign

aa=pd.read_excel('/Users/samiha/Desktop/Supplement/temperate_forest/AVG dataset.xlsx')
aa=aa[aa['pollutant category']=='NMOG']

cm=pd.read_excel('/Users/samiha/Desktop/CM files/NMOG_SAPRC_mechanism.xlsx')

aa=aa.merge(cm[['S18B','id']],how='left')

ind=list(aa[aa['S18B'].isnull()][aa['compound'].str.contains('unknown')].index)
unk_mass=aa['AVG'][aa['S18B'].isnull()][aa['compound'].str.contains('unknown')].sum()

aa=aa.drop(index=ind).reset_index(drop=True)

# Manual assignment- Benzofuran, -dimethyl- (isomer)s
aa.loc[474,'S18B']='NAPS'
# Hexenyne
aa.loc[115,'S18B']='OLE2'
# Benzaldehyde, 2+3-methyl-
aa.loc[360,'S18B']='BALD'
# Butene nitrates
aa.loc[341,'S18B']='R1NO3'

aa=aa[aa['S18B'].notna()]
#fix S18B / 
ind=list(aa[aa['S18B'].str.contains('/',na=False)].index)
for i in ind:
    aa.loc[i,'S18B']=aa['S18B'].iloc[i].split('/')[0]
#mass fraction

uu=list(aa['S18B'].unique())

ndf=pd.DataFrame()
ndf['S18B']=uu
for i in range(len(ndf)):
    ndf.loc[i,'mass']=aa['AVG'][aa['S18B']==ndf['S18B'].iloc[i]].sum()
    
total_mass = ndf['mass'].sum()

ndf['mf']=ndf['mass']/total_mass

#total VOC 91.55 ppb

ndf['InitConc']=ndf['mf']*(91.55)

inorganic=pd.read_excel('/Users/samiha/Desktop/F0AM simulation/NMOG_dataset/InitConc_inorganic.xlsx')

inorganic['S18B']=inorganic['model_species']
df=inorganic[['S18B','InitConc_ppb']].append(adf[['S18B','InitConc_ppb']]).reset_index(drop=True)
for i in range(len(df)):
    df.loc[i,'hold me']=0




# Init File analysis

andf=pd.read_excel('/Users/samiha/Desktop/F0AM simulation/NMOG_dataset/aa_new/andreae_S07_base_InitConc.xlsx')
ndf=pd.read_excel('/Users/samiha/Desktop/F0AM simulation/NMOG_dataset/aa_new/neiva_initconc_S07_base.xlsx')

andf=andf[6:]
ndf=ndf[6:]


ndf[~ndf['S07'].isin(andf['S07'])]
andf['InitConc_an']=andf['InitConc_ppb']

ndf=ndf.merge(andf[['S07','InitConc_an']],on='S07',how='left')
ndf=ndf.sort_values(by='InitConc_ppb',ascending=False)


#hlight data
x_axis=np.arange(len(ndf))

plt.figure(figsize=(12,5))
plt.bar(x_axis - 0.2,ndf['InitConc_ppb'],width=0.4,color='blue',label='NEIVA DB')
plt.bar(x_axis + 0.2,ndf['InitConc_an'],width=0.4,color='orange',label='Andreae19 DB')

plt.xticks(x_axis,ndf['S07'],rotation=90,fontsize=10)
plt.ylabel('Emission (ppb)',fontsize=10)
plt.grid(alpha=0.2)
plt.title('VOC emission mapped to SAPRC07 model species')
plt.legend()
#plt.xlim([0,33])
plt.tight_layout()

plt.savefig("VOC_speciation", format='png', dpi=350, bbox_inches='tight',pad_inches=0.1)


