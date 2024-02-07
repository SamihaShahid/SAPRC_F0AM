
SpeciesToAdd = {...
    'DLIMO';'LIMRXN';};

AddSpecies

i=i+1;
Rnames{i} = 'DLIMO  =  LIMRXN';
k(:,i) = 3.4e-4;
Gstr{i,1} = 'DLIMO';
fDLIMO(i)=fDLIMO(i)-1; fLIMRXN(i)=fLIMRXN(i)+1;
