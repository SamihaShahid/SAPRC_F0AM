# -*- coding: utf-8 -*-
"""
Created on Sat May 14 17:20:22 2022

@author: Taufiq
"""
import pandas as pd

fzs_dir = "/Users/samiha/Desktop/F0AM simulation/SAPRC18"

# reading as list with each line as an element 
with open(fzs_dir+"/STDZA640.FZS") as file:
    lines = file.readlines()

# data starts at line 15    
data = lines[15:]

## number of photolysis rates round(len(data)/7) = 48
l=0
jnames=[]
A=[]
B=[]
E=[]
F=[]
CZLOW=[]
for i in range(round(len(data)/7)):
    val = data[l]
    ll=list(filter(None, val.split(' ')))
    jnames.append('J'+ll[0].strip().replace('-','_'))
    A.append((ll[1].strip().replace('E','e')))
    B.append((ll[2].strip().replace('E','e')))
    E.append((ll[3].strip().replace('E','e')))
    F.append((ll[4].strip().replace('E','e')))
    CZLOW.append((ll[5].strip().replace('E','e')))
    l = l+7

df = pd.DataFrame()
df['jnames'] = jnames
df['A'] = A
df['B'] = B
df['E'] = E
df['F'] = F
df['CZLOW'] = CZLOW    
# df.to_excel('FZS_out.xlsx',index=False)
for i in range(len(df)):
    l1 = 'i=i+1'
    l2 = 'Jnames{i} = '+"'"+df['jnames'][i]+"'"
    l3 = 'A(i) = '+df['A'][i]
    l4 = 'B(i) = '+df['B'][i]
    l5 = 'E(i) = '+df['E'][i]
    l6 = 'F(i) = '+df['F'][i]
    l7 = 'CZLOW(i) = '+df['CZLOW'][i]
    blocks=l1+';\n'+l2+';\n'+l3+';\n'+l4+';\n'+l5+';\n'+l6+';\n'+l7+';\n\n'
    with open('FZS2matlab_SAPRC18.txt','a') as the_file:
        the_file.write(blocks)
