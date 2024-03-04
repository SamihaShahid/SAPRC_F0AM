#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  3 10:37:24 2022

@author: samiha
"""

import pandas as pd
import pubchempy as pcp
import numpy as np
import glob
import re
import matplotlib.pyplot as plt


aa=pd.read_excel('/Users/samiha/Desktop/NEIVA_v_1.1/Data/Recommended_EF.xlsx')
cm=pd.read_excel('/Users/samiha/Desktop/CM files/all_s07_include_hisomer.xlsx')

cc=cm[['S07T','id']]


bb=aa[['compound','AVG temperate_forest','id']][aa['pollutant category']=='NMOG']
bb=bb.rename(columns={'AVG temperate_forest':'ef'})
bb=bb[bb['ef'].notna()]

bb=bb.merge(cc,on='id',how='left')



dd=pd.read_excel('/Users/samiha/Desktop/NEIVA_v_1.1/Data/IntData.xlsx')

nmog=dd[dd['pollutant category']=='NMOG']
nmog=nmog.merge(cc,on='id',how='left')

nmog2=nmog[nmog['S07T'].isnull()][~nmog2['compound'].str.contains('unk')]



#sptool
test=pd.read_excel('/Users/samiha/Desktop/speciation db/shared schema/tbl_mechanism.xlsx')
com=pd.read_excel('/Users/samiha/Desktop/speciation db/sptool_comlist_inchi_out.xlsx')


bb[bb['id'].isin(com['id'])]
bb=bb.merge(com[['SPECIES_ID','id']],on='id',how='left')
bb=bb.rename(columns={'SPECIES_ID':'specie_id'})

            
    
#
sp=list(test['specie_id'][test['mechanism']=='SAPRC07TC_AE7'][test['aqm_poll']=='RCHO'])

mech=test[test['mechanism']=='CB6R3_AE7']

mech['aqm_poll'][mech['specie_id'].isin(sp)].unique()



### match with SPECIATE compounds
test=pd.read_excel('/Users/samiha/Desktop/speciation db/shared schema/tbl_mechanism.xlsx')
com=pd.read_excel('/Users/samiha/Desktop/speciation db/sptool_comlist_inchi_out.xlsx')
aa=pd.read_excel('/Users/samiha/Desktop/NEIVA_v_1.1/Data/Recommended_EF.xlsx')

bb=aa[['compound','AVG temperate_forest','id']][aa['pollutant category']=='NMOG']
bb=bb.rename(columns={'AVG temperate_forest':'ef'})
bb=bb[bb['ef'].notna()]


unkmass=bb['ef'][bb['compound'].str.contains('unk')].sum()
bb=bb[~bb['compound'].str.contains('unk')]

match=bb[bb['id'].isin(com['id'])]

#
match=match.merge(com[['id','SPECIES_ID']],on='id',how='left')
match=match.rename(columns={'SPECIES_ID':'specie_id'})

newmech=mech[mech['specie_id'].isin(match['specie_id'])]
newmech=newmech.reset_index(drop=True)

uu=list(newmech['aqm_poll'].unique())

for i in range(len(match)):
    ss=newmech['moles_per_mole'][newmech['specie_id']==match['specie_id'].iloc[i]].sum()
    match.loc[i,'ef']=match['ef'][i]/ss
    
efsum=[]    
for u in uu:
    spid=list(newmech['specie_id'][newmech['aqm_poll']==u])
    f=list(newmech['moles_per_mole'][newmech['aqm_poll']==u])
    ef=match['ef'][match['specie_id'].isin(spid)]*f
    efsum.append(ef.sum())

df=pd.DataFrame()
df['ms']=uu
df['ef']=efsum


#
unmatch=bb[~bb['id'].isin(com['id'])]

unmatch=unmatch.merge(cm[['id','S07T']],on='id',how='left')

unk=unmatch['ef'][unmatch['S07T'].isnull()].sum()

unmatch=unmatch[unmatch['S07T'].notna()]



df.loc[-1,'ms']='TERP'
df.loc[-1,'ef']=0.3382245497832932

df=df.reset_index(drop=True)

unmatch_mass=unmatch['ef'][unmatch['S07T']!='TERP'].sum()

tef=unmatch_mass+unk

df=df.sort_values(by='ef',ascending=False).reset_index(drop=True)

df['ef']=df['ef'] + (0.54958)


##


