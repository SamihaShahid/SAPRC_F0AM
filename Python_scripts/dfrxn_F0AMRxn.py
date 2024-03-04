#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 16 18:48:28 2022

@author: samiha
"""

import pandas as pd
import pubchempy as pcp
import numpy as np

df=pd.read_excel('/Users/samiha/Desktop/F0AM simulation/SAPRC18/df_rxn_files/S18_r_PHEN_CRES_Rxns.xlsx')


rxn=df.rxn.tolist()
rates=df.k_name.tolist()

reactants = []
for r in rxn[:]:
    reactants.append(r.split('=')[0].strip())

products_old = []
for r in rxn[:]:
    products_old.append(r.split('=')[1].strip())

products = []
for pp in products_old:
    prod = pp.split('+')
    real_p = []
    for p in prod:
        real_p.append(p.split('+')[0].strip().replace(' ','*').replace('#',''))
    products.append(' + '.join(real_p))
    
fvals = []
for pp in products[:]:
    pstr = ''
    aa = pp.split(' ')
    pm = '+'
    for item in aa:
        if '-' in item:
            pm = '-'
        else:
            pm = '+'
        if '*' in item:
            fac = item.split('*')[0]
            pstr += 'f'+item.split('*')[1]+'(i)=f'+item.split('*')[1]+'(i)'+pm+fac+'; '
        if (('+' not in item) and ('-' not in item)) and ('*' not in item):
            pstr += 'f'+item+'(i)=f'+item+'(i)'+pm+'1; '
        pstr = pstr.replace('--','-')
    fvals.append(pstr)

for i in range(len(rxn)):
    l4 = ''
    l5 = ''
    nn = 1
    for r in reactants[i].split('+'):
        rr = r.strip()
        if rr != 'HV':
            l4 += "Gstr{i,"+str(nn)+"} = '"+rr+"'; "
            l5 += 'f'+rr+'(i)=f'+rr+'(i)-1; '
        nn = nn+1
    l1 = 'i=i+1;\n'
    l2 = "Rnames{i} = '"+reactants[i]+" = "+products[i].replace('+ -','- ')+"';\n"
    l3 = 'k(:,i) = '+rates[i]+';\n'
    blocks = l1+l2+l3+l4[:-1]+'\n'+l5+fvals[i][:-1]+'\n\n'
    with open('S18_Jia_PHEN_CRES_Rxns.txt','a') as the_file:
        the_file.write(blocks)

############### Species to add to txt #######################
rr=[]
for r in reactants:
    rsp = r.split('+')
    for rc in rsp:
        rr.append(rc.strip())
# rr = list(set(rr))
pp=[]
for r in products:
    rsp = r.split('+')
    for rc in rsp:
        if '*' in rc:
            pp.append(rc.split('*')[1].strip())
        else:
            pp.append(rc.strip())
# pp = list(set(pp))
def chunks(lst, n):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), n):
        yield lst[i:i + n]        
allsp = list(set(rr+pp))
clist=list(chunks(allsp,10))

spstr=''
for ll in clist:
    val="'; '".join(ll)
    spstr+="'"+val+"'; ...\n"
    
with open('SpeciesToAdd_S18_Jia_PHEN_CRES_Rxns.txt','a') as the_file:
    the_file.write(spstr)
############### Species to add to txt #######################

