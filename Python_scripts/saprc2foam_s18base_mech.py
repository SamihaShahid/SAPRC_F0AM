#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 12 13:31:27 2022

@author: samiha
"""
import pandas as pd
import pubchempy as pcp
import numpy as np

with open('/Users/samiha/Desktop/desktop_v1/F0AM simulation/SAPRC18/SAPRC18_Rxn.txt') as f:
    rxnfile=f.readlines()

# .RXN line index
line_index = [x for x in range(len(rxnfile)) if ".RXN" in rxnfile[x]]

#
rxn=rxnfile[line_index[0]+1:]

rxn_ind=[]
for i in range(len(rxn)):
    if rxn[i].find(')')!=-1:
        rxn_ind.append(i)


rxn_id=[]
k=[]
rxn_line=[]

for ind in range(len(rxn_ind)):
    i=ind
    if rxn_ind[i]!=rxn_ind[-1]:
        f=ind+1
        aa=rxn[rxn_ind[i]:rxn_ind[f]]
    if rxn_ind[i]==rxn_ind[-1]:
        aa=rxn[rxn_ind[i]:]
    
    if len(aa)==1:
        rxn_line.append(aa[0].split(';')[1].replace('\n','').strip())
        rxn_id.append(aa[0].split(';')[0].split(')')[0])
        k.append(aa[0].split(';')[0].split(')')[1].strip())
    if len(aa)>1:
        if len(aa[1])>29:
            rxn_line.append((' ').join(aa).split(';')[1].replace('  ','').replace('\n',''))
            k.append(aa[0].split(';')[0].strip().split(')')[1].strip())
            rxn_id.append(aa[0].split(';')[0].strip().split(')')[0].strip())
        if len(aa[1])<29:
            k.append(aa[0].split(';')[0].split(')')[1].strip()+','+(',').join(aa[1:]).replace('  ','').replace('\n',''))
            rxn_line.append(aa[0].split(';')[1].strip().replace('\n',''))
            rxn_id.append(aa[0].split(';')[0].split(')')[0])
            
# In order to discard the '.' at the end of the line
rxn_line[-1]='MAPAN_P2 + SumRCO3 = SumRCO3 + #.8 R2CO3 + #.8 HCHO + #.2 PAN2 + #-0.8 XC + #.8 XN + #.8 SumRCO3'

# 
rxn_line[62]='ACETL + OH = #.3 HO2 + #.7 OH + #.3 CO + #.3 HCOOH + #.7 GLY'
k[62]='FALLOFF,5.50e-30 0.000 0.00,8.30e-13 0.000 2.00,0.60 1.00'


            
ii=[]
for i in range(len(rxn_line)):
    if rxn_line[i].find('{')!=-1:
        ii.append(i)
                        
rxn_line[34]='HNO4 + HV = #.8 HO2 + #.8 NO2 + #.2 OH + #0.2 NO3'            
rxn_line[48]='MEOOH + OH = H2O + #.4 HCHO + #0.4 OH + #.6 MEO2 + #0.6 SumRO2'
rxn_line[64]='ETOH + OH = #.95 HO2 + #0.95 MECHO + #.05 ETHEO2 + #0.05 SumRO2'        
rxn_line[68]='MECHO + OH = H2O + #.95 MECO3 + #0.95 SumRCO3 + #.05 HCOMEO2 + #0.05 SumRO2'    
rxn_line[69]='MECHO + HV = HO2 + #.9 CO + #0.9 MEO2 + #0.9 SumRO2 + #.1 MECO3 + #0.1 SumRCO3'        
rxn_line[74]='GLCHO + OH = #.2 HO2 + #.8 HOCCO3 + #0.8 SumRCO3 + #.2 GLY'
rxn_line[75]='GLCHO + NO3 = HNO3 + #.991 HOCCO3 + #0.991 SumRCO3 + #.009 CO + #0.009 HCHO + #0.009 HO2'    
rxn_line[76]='GLCHO + HV = #.93 CO + #.1 MEOH + #.07 OH + #1.66 HO2 + #.83 HCHO + #.07 HCOMEO2 + #0.07 SumRO2'        
rxn_line[79]='GLY + HV = #2 CO + #2 HO2'
rxn_line[81]='GLY + OH = #1.7 CO + #.7 HO2 + #.3 OH + #0.3 CO2'   
rxn_line[82]='GLY + NO3 = HNO3 + #1.7 CO + #.7 HO2 + #.3 OH + #0.3 CO2'
rxn_line[105]='PHOT + HV = #2 HO2 + #2 RO2C + #2 SumRO2 + OTH2 + #1 XC'       
rxn_line[118]='MEO2 + HO2 = #.9 MEOOH + #.1 HCHO + #0.1 H2O + O2'    
rxn_line[121]='MEO2 + SumRCO3 = SumRCO3 + #.9 HCHO + #0.9 HO2 + #.1 HCHO + #0.1 O2'
rxn_line[210]='MECO3 + HO2 = #.37 PAA + #.13 O3 + #0.13 AACID + #.5 OH + #0.5 MEO2 + #0.5 CO2 + #0.5 SumRO2'
rxn_line[212]='MECO3 + SumRO2 = SumRO2 + #.9 MEO2 + #0.9 CO2 + #.1 AACID + #0.9 SumRO2'
rxn_line[215]='PAN + HV = #.7 MECO3 + #0.7 NO2 + #0.7 SumRCO3 + #.3 MEO2 + #0.3 CO2 + #0.3 NO3 + #0.3 SumRO2'
rxn_line[218]='HOCCO3 + HO2 = #.37 PAA + #.13 O3 + #0.13 AACID + #.5 OH + #0.5 HCHO + #0.5 HO2 + #0.5 CO2'
rxn_line[223]='HOPAN + HV = #.6 HOCCO3 + #0.6 NO2 + #0.6 SumRCO3 + #.4 HCHO + #0.4 HO2 + #0.4 CO2 + #0.4 NO3'
rxn_line[226]='ETCO3 + HO2 = #.37 RCOOH + #.13 O3 + #0.13 RCOOH + #.5 OH + #0.5 ETO2 + #0.5 CO2 + #0.5 SumRO2'
rxn_line[231]='PPN + HV = #.6 ETCO3 + #0.6 NO2 + #.4 ETO2 + #0.4 CO2 + #0.4 NO3 + #0.4 SumRO2'
rxn_line[234]='ACO3 + HO2 = #.37 RCOOH + #.13 O3 + #0.13 RCOOH + #.5 OH + #0.5 HO2 + #0.5 CO + #0.5 CO2 + #0.5 HCHO'
rxn_line[240]='MACO3 + HO2 = #.37 RCOOH + #.13 O3 + #0.13 RCOOH + #.5 OH + #0.5 HCHO + #0.5 MECO3 + #0.5 CO2 + #.5 XC + #.5 SumRCO3'
rxn_line[246]='R2CO3 + HO2 = #.37 RCOOH + #.13 O3 + #0.13 RCOOH + #.5 OH + #0.5 C3RO2 + #0.5 CO2 + #-.5 XC + #.5 SumRCO3'
rxn_line[251]='PAN2 + HV = #.6 R2CO3 + #0.6 NO2 + #.4 C3RO2 + #0.4 CO2 + #0.4 NO3 + #-.4 XC + #0.6 SumRCO3'
rxn_line[254]='R2NCO3 + HO2 = #.37 RCOOH + #.13 O3 + #0.13 RCOOH + #.5 OH + #0.5 NO2 + #0.5 HCHO + #0.5 CO2 + #-.5 XC + #.5 XN + #.5 SumRCO3'
rxn_line[259]='PAN2N + HV = #.6 R2NCO3 + #0.6 NO2 + #.4 NO2 + #0.4 HCHO + #0.4 CO2 + #0.4 NO3 + #0.6 SumRCO3'
rxn_line[262]='BZCO3 + HO2 = #.37 RCOOH + #.13 O3 + #0.13 RCOOH + #.5 OH + #0.5 BZO2 + #0.5 CO2 + #2 XC + #5 SumRO2'
rxn_line[267]='PBZN + HV = #.6 BZCO3 + #0.6 NO2 + #.4 CO2 + #0.4 BZO2 + #0.4 NO3 + #0.3 SumRO2'


df=pd.DataFrame()
df['id']=rxn_id
df['k']=k
df['rxn']=rxn_line

# k = A* (T/300)^B * e(-E/R*T)

#Jnames
Jind=list(df['k'][df['k'].str.contains('PF')].index)
for i in Jind:
    if df['k'][i].find('QY')==-1:
        df.loc[i,'k_name']='J'+df['k'][i].split('=')[1].replace('-','_')
    if df['k'][i].find('QY')!=-1:
        df.loc[i,'k_name']='J'+df['k'][i].split('PF=')[1].split(' ')[0].replace('-','_')+'.*'+df['k'][i].split('PF=')[1].split(' ')[1].split('QY=')[1]
        
one_ind=[]
for i in range(len(df)):
    if len(df['k'][i].split(' '))==1:
        one_ind.append(i)

one_ind=set(one_ind)-set(Jind)        
for i in one_ind:
    df.loc[i,'k_name']=df['k'].iloc[i]

samek_ind=set(df[df['k'].str.contains('SAMEK')].index)
K_ind=set(df[df['k'].str.contains('K')].index) - samek_ind - set(Jind)
falloff = set(df[df['k'].str.contains('FALLOFF')].index)

others = set(df[df['k_name'].isnull()].index) - (samek_ind) - (K_ind) - falloff
        
#samek
id_l=[]
for i in samek_ind:
    id_l.append(df['k'][i].split(' ')[1])
        
#others
R=0.0019872
for i in others:
    if len(df['k'][i].split(' '))==2:
        A=df['k'][i].split(' ')[0]
        Ea=-1*float(df['k'][i].split(' ')[1])
        df.loc[i,'k_name'] = A + '.*exp(' + str(round(Ea/R,3)) + './T)'
    if len(df['k'][i].split(' '))==3:
        A=df['k'][i].split(' ')[0]
        Ea=-1*float(df['k'][i].split(' ')[1])
        B=df['k'][i].split(' ')[2]
        df.loc[i,'k_name']= A + '.*(T./300).^'+ B +'.*exp(' + str(round(Ea/R,3)) + './T)'
    
#falloff
for i in falloff:
    df.loc[i,'k_name']='kf_'+df['rxn'][i].split('=')[0].strip().replace('+','_').replace(' ','')
    
for i in K_ind:
    df.loc[i,'k_name']='k_'+df['rxn'][i].split('=')[0].strip().replace('+','_').replace(' ','')

#samek_ind
for i in samek_ind:
    iid=df['k'][i].split(' ')[1]
    ind=list(df[df['id']==iid].index)[0]
    df.loc[i,'k_name']=df['k_name'][ind]

#

df.to_excel('/Users/samiha/Desktop/F0AM simulation/SAPRC18/df_rxn_files/S18_Jia_PHEN_Rxns.xlsx',index=False)

duplines=['CROOH_P3 + HO2 ',
          'CROOH_P3 + NO ',
          'CROOH_P3 + NO3 ',
          'CROOH_P3 + SumRCO3 ',
          'CROOH_P3 + SumRO2 ',
          'CROOH_P3 ',
          'OLEA1_R1 + NO ',
          'OLEA1_R1 ']

# checking XC
for i in range(len(df)):
    df.loc[i,'r']=df['rxn'][i].split('=')[0]
    df.loc[i,'p']=df['rxn'][i].split('=')[1]
    
    
d_ind=list(df['r'][df['rxn'].duplicated()].index)
df=df.drop(index=d_ind)
df=df.reset_index(drop=True)

# XC

for i in range(len(df)):
    df.loc[i,'rxn2']=df['rxn'][i].replace('RO2XC','')

for i in range(len(df)):
    df.loc[i,'p2']=df['rxn2'][i].split('=')[1]

ind=list(df[df['p2'].str.contains('XC')].index)

for i in ind:
    df.loc[i,'fXC']=df['p'].iloc[i].replace('RO2XC','').split('XC')[0].split('+')[-1].replace('#','').strip()

ind=list(df[df['fXC']==''].index)
for i in ind:
    df.loc[i,'fXC']=1

ind=list(df[df['fXC'].notna()].index)
for i in ind:
    df.loc[i,'fXC']=float(df['fXC'][i])

#init model species
m=['CH4','NO','NO2','OH','O3','CO','HCHO','HONO',
   'MEOH','MECHO','PROPE','BENZ','FURNS','ISOP','ACET','ETCHO','HCOOH',
   'AACID','BACL','MVK','MACR','MGLY','OLEA1','KET2','MALAH',
   'PHEN','SVPHE','STYRS','CRES']

m2=['PROPE','BENZ','FURNS','ISOP',
    'ACET','ETCHO','HCOOH',
   'AACID','BACL','MVK','MACR',
   'MGLY','OLEA1','KET2','MALAH',
   'PHEN','SVPHE','STYRS','CRES']

ss=[]
for i in m2:
    ss.append(df['fXC'][df['r'].str.contains(i)].sum())



    