SpeciesToAdd = {'BZCO3';'NO';'yHFONS';'yBACL';'yRNO3';...
    'BACL';'XC';'RO2C';'MACO3';'MEK';'HFONS';'H3XE25DO';'ROOH';...
    'RNO3';'HO2';'NO3';'x3HXE25DO';'PRD2';'MECO3';...
    'xROOH';'AFG3';'NO3';'yAFG3';'RCO3';'HO2';'MEO2';'zHFONS';...
    'RO2XC';'MGLY';'yMGLY';'xNO3';'xHFONS';};

AddSpecies

i=i+1;
Rnames{i} = 'xNO3 + NO = NO + NO3';
k(:,i) = 2.60e-12.*exp(380./ T);
Gstr{i,1} = 'xNO3'; Gstr{i,2} = 'NO';
fxNO3(i)=fxNO3(i)-1; fNO(i)=fNO(i)-1; fNO(i)=fNO(i)+1; fNO3(i)=fNO3(i)+1;

i=i+1;
Rnames{i} = 'xNO3 + HO2 = HO2';
k(:,i) = 3.80e-13.*exp(900./ T);
Gstr{i,1} = 'xNO3'; Gstr{i,2} = 'HO2';
fxNO3(i)=fxNO3(i)-1; fHO2(i)=fHO2(i)-1; fHO2(i)=fHO2(i)+1;

i=i+1;
Rnames{i} = 'xNO3 + NO3 = NO3 + NO3';
k(:,i) = 2.30e-12;
Gstr{i,1} = 'xNO3'; Gstr{i,2} = 'NO3';
fxNO3(i)=fxNO3(i)-1; fNO3(i)=fNO3(i)-1; fNO3(i)=fNO3(i)+1; fNO3(i)=fNO3(i)+1;

i=i+1;
Rnames{i} = 'xNO3 + MEO2 = NO3 + 0.5*NO3';
k(:,i) = 2.00e-13;
Gstr{i,1} = 'xNO3'; Gstr{i,2} = 'MEO2';
fxNO3(i)=fxNO3(i)-1; fMEO2(i)=fMEO2(i)-1; fNO3(i)=fNO3(i)+1; fNO3(i)=fNO3(i)+0.5;

i=i+1;
Rnames{i} = 'xNO3 + RO2C = RO2C + 0.5*NO3';
k(:,i) = 3.50e-14;
Gstr{i,1} = 'xNO3'; Gstr{i,2} = 'RO2C';
fxNO3(i)=fxNO3(i)-1; fRO2C(i)=fRO2C(i)-1; fRO2C(i)=fRO2C(i)+1; fNO3(i)=fNO3(i)+0.5;

i=i+1;
Rnames{i} = 'xNO3 + RO2XC = RO2XC + 0.5*NO3';
k(:,i) = 3.50e-14;
Gstr{i,1} = 'xNO3'; Gstr{i,2} = 'RO2XC';
fxNO3(i)=fxNO3(i)-1; fRO2XC(i)=fRO2XC(i)-1; fRO2XC(i)=fRO2XC(i)+1; fNO3(i)=fNO3(i)+0.5;

i=i+1;
Rnames{i} = 'xNO3 + MECO3 = MECO3 + NO3';
k(:,i) = 4.40e-13.*exp(1070./ T);
Gstr{i,1} = 'xNO3'; Gstr{i,2} = 'MECO3';
fxNO3(i)=fxNO3(i)-1; fMECO3(i)=fMECO3(i)-1; fMECO3(i)=fMECO3(i)+1; fNO3(i)=fNO3(i)+1;

i=i+1;
Rnames{i} = 'xNO3 + RCO3 = RCO3 + NO3';
k(:,i) = 4.40e-13.*exp(1070./ T);
Gstr{i,1} = 'xNO3'; Gstr{i,2} = 'RCO3';
fxNO3(i)=fxNO3(i)-1; fRCO3(i)=fRCO3(i)-1; fRCO3(i)=fRCO3(i)+1; fNO3(i)=fNO3(i)+1;

i=i+1;
Rnames{i} = 'xNO3 + BZCO3 = BZCO3 + NO3';
k(:,i) = 4.40e-13.*exp(1070./ T);
Gstr{i,1} = 'xNO3'; Gstr{i,2} = 'BZCO3';
fxNO3(i)=fxNO3(i)-1; fBZCO3(i)=fBZCO3(i)-1; fBZCO3(i)=fBZCO3(i)+1; fNO3(i)=fNO3(i)+1;

i=i+1;
Rnames{i} = 'xNO3 + MACO3 = MACO3 + NO3';
k(:,i) = 4.40e-13.*exp(1070./ T);
Gstr{i,1} = 'xNO3'; Gstr{i,2} = 'MACO3';
fxNO3(i)=fxNO3(i)-1; fMACO3(i)=fMACO3(i)-1; fMACO3(i)=fMACO3(i)+1; fNO3(i)=fNO3(i)+1;

i=i+1;
Rnames{i} = 'x3HXE25DO + NO = NO + H3XE25DO';
k(:,i) = 2.60e-12.*exp(380./ T);
Gstr{i,1} = 'x3HXE25DO'; Gstr{i,2} = 'NO';
fx3HXE25DO(i)=fx3HXE25DO(i)-1; fNO(i)=fNO(i)-1; fNO(i)=fNO(i)+1; fH3XE25DO(i)=fH3XE25DO(i)+1;

i=i+1;
Rnames{i} = 'x3HXE25DO + HO2 = HO2 + 6.0*XC';
k(:,i) = 3.80e-13.*exp(900./ T);
Gstr{i,1} = 'x3HXE25DO'; Gstr{i,2} = 'HO2';
fx3HXE25DO(i)=fx3HXE25DO(i)-1; fHO2(i)=fHO2(i)-1; fHO2(i)=fHO2(i)+1; fXC(i)=fXC(i)+6.0;

i=i+1;
Rnames{i} = 'x3HXE25DO + NO3 = NO3 + H3XE25DO';
k(:,i) = 2.30e-12;
Gstr{i,1} = 'x3HXE25DO'; Gstr{i,2} = 'NO3';
fx3HXE25DO(i)=fx3HXE25DO(i)-1; fNO3(i)=fNO3(i)-1; fNO3(i)=fNO3(i)+1; fH3XE25DO(i)=fH3XE25DO(i)+1;

i=i+1;
Rnames{i} = 'x3HXE25DO + MEO2 = H3XE25DO + 0.5*H3XE25DO + 3.0*XC';
k(:,i) = 2.00e-13;
Gstr{i,1} = 'x3HXE25DO'; Gstr{i,2} = 'MEO2';
fx3HXE25DO(i)=fx3HXE25DO(i)-1; fMEO2(i)=fMEO2(i)-1; fH3XE25DO(i)=fH3XE25DO(i)+1; fH3XE25DO(i)=fH3XE25DO(i)+0.5; fXC(i)=fXC(i)+3.0;

i=i+1;
Rnames{i} = 'x3HXE25DO + RO2C = RO2C + 0.5*H3XE25DO + 3.0*XC';
k(:,i) = 3.50e-14;
Gstr{i,1} = 'x3HXE25DO'; Gstr{i,2} = 'RO2C';
fx3HXE25DO(i)=fx3HXE25DO(i)-1; fRO2C(i)=fRO2C(i)-1; fRO2C(i)=fRO2C(i)+1; fH3XE25DO(i)=fH3XE25DO(i)+0.5; fXC(i)=fXC(i)+3.0;

i=i+1;
Rnames{i} = 'x3HXE25DO + RO2XC = RO2XC + 0.5*H3XE25DO + 3.0*XC';
k(:,i) = 3.50e-14;
Gstr{i,1} = 'x3HXE25DO'; Gstr{i,2} = 'RO2XC';
fx3HXE25DO(i)=fx3HXE25DO(i)-1; fRO2XC(i)=fRO2XC(i)-1; fRO2XC(i)=fRO2XC(i)+1; fH3XE25DO(i)=fH3XE25DO(i)+0.5; fXC(i)=fXC(i)+3.0;

i=i+1;
Rnames{i} = 'x3HXE25DO + MECO3 = MECO3 + H3XE25DO';
k(:,i) = 4.40e-13.*exp(1070./ T);
Gstr{i,1} = 'x3HXE25DO'; Gstr{i,2} = 'MECO3';
fx3HXE25DO(i)=fx3HXE25DO(i)-1; fMECO3(i)=fMECO3(i)-1; fMECO3(i)=fMECO3(i)+1; fH3XE25DO(i)=fH3XE25DO(i)+1;

i=i+1;
Rnames{i} = 'x3HXE25DO + RCO3 = RCO3 + H3XE25DO';
k(:,i) = 4.40e-13.*exp(1070./ T);
Gstr{i,1} = 'x3HXE25DO'; Gstr{i,2} = 'RCO3';
fx3HXE25DO(i)=fx3HXE25DO(i)-1; fRCO3(i)=fRCO3(i)-1; fRCO3(i)=fRCO3(i)+1; fH3XE25DO(i)=fH3XE25DO(i)+1;

i=i+1;
Rnames{i} = 'x3HXE25DO + BZCO3 = BZCO3 + H3XE25DO';
k(:,i) = 4.40e-13.*exp(1070./ T);
Gstr{i,1} = 'x3HXE25DO'; Gstr{i,2} = 'BZCO3';
fx3HXE25DO(i)=fx3HXE25DO(i)-1; fBZCO3(i)=fBZCO3(i)-1; fBZCO3(i)=fBZCO3(i)+1; fH3XE25DO(i)=fH3XE25DO(i)+1;

i=i+1;
Rnames{i} = 'x3HXE25DO + MACO3 = MACO3 + H3XE25DO';
k(:,i) = 4.40e-13.*exp(1070./ T);
Gstr{i,1} = 'x3HXE25DO'; Gstr{i,2} = 'MACO3';
fx3HXE25DO(i)=fx3HXE25DO(i)-1; fMACO3(i)=fMACO3(i)-1; fMACO3(i)=fMACO3(i)+1; fH3XE25DO(i)=fH3XE25DO(i)+1;

i=i+1;
Rnames{i} = 'xHFONS + NO = NO + HFONS';
k(:,i) = 2.60e-12.*exp(380./ T);
Gstr{i,1} = 'xHFONS'; Gstr{i,2} = 'NO';
fxHFONS(i)=fxHFONS(i)-1; fNO(i)=fNO(i)-1; fNO(i)=fNO(i)+1; fHFONS(i)=fHFONS(i)+1;

i=i+1;
Rnames{i} = 'xHFONS + HO2 = HO2 + 5.0*XC';
k(:,i) = 3.80e-13.*exp(900./ T);
Gstr{i,1} = 'xHFONS'; Gstr{i,2} = 'HO2';
fxHFONS(i)=fxHFONS(i)-1; fHO2(i)=fHO2(i)-1; fHO2(i)=fHO2(i)+1; fXC(i)=fXC(i)+5.0;

i=i+1;
Rnames{i} = 'xHFONS + NO3 = NO3 + HFONS';
k(:,i) = 2.30e-12;
Gstr{i,1} = 'xHFONS'; Gstr{i,2} = 'NO3';
fxHFONS(i)=fxHFONS(i)-1; fNO3(i)=fNO3(i)-1; fNO3(i)=fNO3(i)+1; fHFONS(i)=fHFONS(i)+1;

i=i+1;
Rnames{i} = 'xHFONS + MEO2 = HFONS + 0.5*HFONS + 2.5*XC';
k(:,i) = 2.00e-13;
Gstr{i,1} = 'xHFONS'; Gstr{i,2} = 'MEO2';
fxHFONS(i)=fxHFONS(i)-1; fMEO2(i)=fMEO2(i)-1; fHFONS(i)=fHFONS(i)+1; fHFONS(i)=fHFONS(i)+0.5; fXC(i)=fXC(i)+2.5;

i=i+1;
Rnames{i} = 'xHFONS + RO2C = RO2C + 0.5*HFONS + 2.5*XC';
k(:,i) = 3.50e-14;
Gstr{i,1} = 'xHFONS'; Gstr{i,2} = 'RO2C';
fxHFONS(i)=fxHFONS(i)-1; fRO2C(i)=fRO2C(i)-1; fRO2C(i)=fRO2C(i)+1; fHFONS(i)=fHFONS(i)+0.5; fXC(i)=fXC(i)+2.5;

i=i+1;
Rnames{i} = 'xHFONS + RO2XC = RO2XC + 0.5*HFONS + 2.5*XC';
k(:,i) = 3.50e-14;
Gstr{i,1} = 'xHFONS'; Gstr{i,2} = 'RO2XC';
fxHFONS(i)=fxHFONS(i)-1; fRO2XC(i)=fRO2XC(i)-1; fRO2XC(i)=fRO2XC(i)+1; fHFONS(i)=fHFONS(i)+0.5; fXC(i)=fXC(i)+2.5;

i=i+1;
Rnames{i} = 'xHFONS + MECO3 = MECO3 + HFONS';
k(:,i) = 4.40e-13.*exp(1070./ T);
Gstr{i,1} = 'xHFONS'; Gstr{i,2} = 'MECO3';
fxHFONS(i)=fxHFONS(i)-1; fMECO3(i)=fMECO3(i)-1; fMECO3(i)=fMECO3(i)+1; fHFONS(i)=fHFONS(i)+1;

i=i+1;
Rnames{i} = 'xHFONS + RCO3 = RCO3 + HFONS';
k(:,i) = 4.40e-13.*exp(1070./ T);
Gstr{i,1} = 'xHFONS'; Gstr{i,2} = 'RCO3';
fxHFONS(i)=fxHFONS(i)-1; fRCO3(i)=fRCO3(i)-1; fRCO3(i)=fRCO3(i)+1; fHFONS(i)=fHFONS(i)+1;

i=i+1;
Rnames{i} = 'xHFONS + BZCO3 = BZCO3 + HFONS';
k(:,i) = 4.40e-13.*exp(1070./ T);
Gstr{i,1} = 'xHFONS'; Gstr{i,2} = 'BZCO3';
fxHFONS(i)=fxHFONS(i)-1; fBZCO3(i)=fBZCO3(i)-1; fBZCO3(i)=fBZCO3(i)+1; fHFONS(i)=fHFONS(i)+1;

i=i+1;
Rnames{i} = 'xHFONS + MACO3 = MACO3 + HFONS';
k(:,i) = 4.40e-13.*exp(1070./ T);
Gstr{i,1} = 'xHFONS'; Gstr{i,2} = 'MACO3';
fxHFONS(i)=fxHFONS(i)-1; fMACO3(i)=fMACO3(i)-1; fMACO3(i)=fMACO3(i)+1; fHFONS(i)=fHFONS(i)+1;

i=i+1;
Rnames{i} = 'xROOH + NO = NO + ROOH';
k(:,i) = 2.60e-12.*exp(380./ T);
Gstr{i,1} = 'xROOH'; Gstr{i,2} = 'NO';
fxROOH(i)=fxROOH(i)-1; fNO(i)=fNO(i)-1; fNO(i)=fNO(i)+1; fROOH(i)=fROOH(i)+1;

i=i+1;
Rnames{i} = 'xROOH + HO2 = HO2 + 3.0*XC';
k(:,i) = 3.80e-13.*exp(900./ T);
Gstr{i,1} = 'xROOH'; Gstr{i,2} = 'HO2';
fxROOH(i)=fxROOH(i)-1; fHO2(i)=fHO2(i)-1; fHO2(i)=fHO2(i)+1; fXC(i)=fXC(i)+3.0;

i=i+1;
Rnames{i} = 'xROOH + NO3 = NO3 + ROOH';
k(:,i) = 2.30e-12;
Gstr{i,1} = 'xROOH'; Gstr{i,2} = 'NO3';
fxROOH(i)=fxROOH(i)-1; fNO3(i)=fNO3(i)-1; fNO3(i)=fNO3(i)+1; fROOH(i)=fROOH(i)+1;

i=i+1;
Rnames{i} = 'xROOH + MEO2 = ROOH + 0.5*ROOH + 1.5*XC';
k(:,i) = 2.00e-13;
Gstr{i,1} = 'xROOH'; Gstr{i,2} = 'MEO2';
fxROOH(i)=fxROOH(i)-1; fMEO2(i)=fMEO2(i)-1; fROOH(i)=fROOH(i)+1; fROOH(i)=fROOH(i)+0.5; fXC(i)=fXC(i)+1.5;

i=i+1;
Rnames{i} = 'xROOH + RO2C = RO2C + 0.5*ROOH + 1.5*XC';
k(:,i) = 3.50e-14;
Gstr{i,1} = 'xROOH'; Gstr{i,2} = 'RO2C';
fxROOH(i)=fxROOH(i)-1; fRO2C(i)=fRO2C(i)-1; fRO2C(i)=fRO2C(i)+1; fROOH(i)=fROOH(i)+0.5; fXC(i)=fXC(i)+1.5;

i=i+1;
Rnames{i} = 'xROOH + RO2XC = RO2XC + 0.5*ROOH + 1.5*XC';
k(:,i) = 3.50e-14;
Gstr{i,1} = 'xROOH'; Gstr{i,2} = 'RO2XC';
fxROOH(i)=fxROOH(i)-1; fRO2XC(i)=fRO2XC(i)-1; fRO2XC(i)=fRO2XC(i)+1; fROOH(i)=fROOH(i)+0.5; fXC(i)=fXC(i)+1.5;

i=i+1;
Rnames{i} = 'xROOH + MECO3 = MECO3 + ROOH';
k(:,i) = 4.40e-13.*exp(1070./ T);
Gstr{i,1} = 'xROOH'; Gstr{i,2} = 'MECO3';
fxROOH(i)=fxROOH(i)-1; fMECO3(i)=fMECO3(i)-1; fMECO3(i)=fMECO3(i)+1; fROOH(i)=fROOH(i)+1;

i=i+1;
Rnames{i} = 'xROOH + RCO3 = RCO3 + ROOH';
k(:,i) = 4.40e-13.*exp(1070./ T);
Gstr{i,1} = 'xROOH'; Gstr{i,2} = 'RCO3';
fxROOH(i)=fxROOH(i)-1; fRCO3(i)=fRCO3(i)-1; fRCO3(i)=fRCO3(i)+1; fROOH(i)=fROOH(i)+1;

i=i+1;
Rnames{i} = 'xROOH + BZCO3 = BZCO3 + ROOH';
k(:,i) = 4.40e-13.*exp(1070./ T);
Gstr{i,1} = 'xROOH'; Gstr{i,2} = 'BZCO3';
fxROOH(i)=fxROOH(i)-1; fBZCO3(i)=fBZCO3(i)-1; fBZCO3(i)=fBZCO3(i)+1; fROOH(i)=fROOH(i)+1;

i=i+1;
Rnames{i} = 'xROOH + MACO3 = MACO3 + ROOH';
k(:,i) = 4.40e-13.*exp(1070./ T);
Gstr{i,1} = 'xROOH'; Gstr{i,2} = 'MACO3';
fxROOH(i)=fxROOH(i)-1; fMACO3(i)=fMACO3(i)-1; fMACO3(i)=fMACO3(i)+1; fROOH(i)=fROOH(i)+1;

i=i+1;
Rnames{i} = 'yHFONS + NO = NO';
k(:,i) = 2.60e-12.*exp(380./ T);
Gstr{i,1} = 'yHFONS'; Gstr{i,2} = 'NO';
fyHFONS(i)=fyHFONS(i)-1; fNO(i)=fNO(i)-1; fNO(i)=fNO(i)+1;

i=i+1;
Rnames{i} = 'yHFONS + HO2 = HO2 + HFONS - 5.0*XC';
k(:,i) = 3.80e-13.*exp(900./ T);
Gstr{i,1} = 'yHFONS'; Gstr{i,2} = 'HO2';
fyHFONS(i)=fyHFONS(i)-1; fHO2(i)=fHO2(i)-1; fHO2(i)=fHO2(i)+1; fHFONS(i)=fHFONS(i)+1; fXC(i)=fXC(i)-5.0;

i=i+1;
Rnames{i} = 'yHFONS + NO3 = NO3';
k(:,i) = 2.30e-12;
Gstr{i,1} = 'yHFONS'; Gstr{i,2} = 'NO3';
fyHFONS(i)=fyHFONS(i)-1; fNO3(i)=fNO3(i)-1; fNO3(i)=fNO3(i)+1;

i=i+1;
Rnames{i} = 'yHFONS + MEO2 = MEO2 + 0.5*PRD2 - 3*XC';
k(:,i) = 2.00e-13;
Gstr{i,1} = 'yHFONS'; Gstr{i,2} = 'MEO2';
fyHFONS(i)=fyHFONS(i)-1; fMEO2(i)=fMEO2(i)-1; fMEO2(i)=fMEO2(i)+1; fPRD2(i)=fPRD2(i)+0.5; fXC(i)=fXC(i)-3;

i=i+1;
Rnames{i} = 'yHFONS + RO2C = RO2C + 0.5*PRD2 - 3*XC';
k(:,i) = 3.50e-14;
Gstr{i,1} = 'yHFONS'; Gstr{i,2} = 'RO2C';
fyHFONS(i)=fyHFONS(i)-1; fRO2C(i)=fRO2C(i)-1; fRO2C(i)=fRO2C(i)+1; fPRD2(i)=fPRD2(i)+0.5; fXC(i)=fXC(i)-3;

i=i+1;
Rnames{i} = 'yHFONS + RO2XC = RO2XC + 0.5*PRD2 - 3*XC';
k(:,i) = 3.50e-14;
Gstr{i,1} = 'yHFONS'; Gstr{i,2} = 'RO2XC';
fyHFONS(i)=fyHFONS(i)-1; fRO2XC(i)=fRO2XC(i)-1; fRO2XC(i)=fRO2XC(i)+1; fPRD2(i)=fPRD2(i)+0.5; fXC(i)=fXC(i)-3;

i=i+1;
Rnames{i} = 'yHFONS + MECO3 = MECO3';
k(:,i) = 4.40e-13.*exp(1070./ T);
Gstr{i,1} = 'yHFONS'; Gstr{i,2} = 'MECO3';
fyHFONS(i)=fyHFONS(i)-1; fMECO3(i)=fMECO3(i)-1; fMECO3(i)=fMECO3(i)+1;

i=i+1;
Rnames{i} = 'yHFONS + RCO3 = RCO3';
k(:,i) = 4.40e-13.*exp(1070./ T);
Gstr{i,1} = 'yHFONS'; Gstr{i,2} = 'RCO3';
fyHFONS(i)=fyHFONS(i)-1; fRCO3(i)=fRCO3(i)-1; fRCO3(i)=fRCO3(i)+1;

i=i+1;
Rnames{i} = 'yHFONS + BZCO3 = BZCO3';
k(:,i) = 4.40e-13.*exp(1070./ T);
Gstr{i,1} = 'yHFONS'; Gstr{i,2} = 'BZCO3';
fyHFONS(i)=fyHFONS(i)-1; fBZCO3(i)=fBZCO3(i)-1; fBZCO3(i)=fBZCO3(i)+1;

i=i+1;
Rnames{i} = 'yHFONS + MACO3 = MACO3';
k(:,i) = 4.40e-13.*exp(1070./ T);
Gstr{i,1} = 'yHFONS'; Gstr{i,2} = 'MACO3';
fyHFONS(i)=fyHFONS(i)-1; fMACO3(i)=fMACO3(i)-1; fMACO3(i)=fMACO3(i)+1;

i=i+1;
Rnames{i} = 'yMGLY + NO = NO';
k(:,i) = 2.60e-12.*exp(380./ T);
Gstr{i,1} = 'yMGLY'; Gstr{i,2} = 'NO';
fyMGLY(i)=fyMGLY(i)-1; fNO(i)=fNO(i)-1; fNO(i)=fNO(i)+1;

i=i+1;
Rnames{i} = 'yMGLY + HO2 = HO2 + MGLY - 3.0*XC';
k(:,i) = 3.80e-13.*exp(900./ T);
Gstr{i,1} = 'yMGLY'; Gstr{i,2} = 'HO2';
fyMGLY(i)=fyMGLY(i)-1; fHO2(i)=fHO2(i)-1; fHO2(i)=fHO2(i)+1; fMGLY(i)=fMGLY(i)+1; fXC(i)=fXC(i)-3.0;

i=i+1;
Rnames{i} = 'yMGLY + NO3 = NO3';
k(:,i) = 2.30e-12;
Gstr{i,1} = 'yMGLY'; Gstr{i,2} = 'NO3';
fyMGLY(i)=fyMGLY(i)-1; fNO3(i)=fNO3(i)-1; fNO3(i)=fNO3(i)+1;

i=i+1;
Rnames{i} = 'yMGLY + MEO2 = MEO2 + 0.5*MEK - 2*XC';
k(:,i) = 2.00e-13;
Gstr{i,1} = 'yMGLY'; Gstr{i,2} = 'MEO2';
fyMGLY(i)=fyMGLY(i)-1; fMEO2(i)=fMEO2(i)-1; fMEO2(i)=fMEO2(i)+1; fMEK(i)=fMEK(i)+0.5; fXC(i)=fXC(i)-2;

i=i+1;
Rnames{i} = 'yMGLY + RO2C = RO2C + 0.5*MEK - 2*XC';
k(:,i) = 3.50e-14;
Gstr{i,1} = 'yMGLY'; Gstr{i,2} = 'RO2C';
fyMGLY(i)=fyMGLY(i)-1; fRO2C(i)=fRO2C(i)-1; fRO2C(i)=fRO2C(i)+1; fMEK(i)=fMEK(i)+0.5; fXC(i)=fXC(i)-2;

i=i+1;
Rnames{i} = 'yMGLY + RO2XC = RO2XC + 0.5*MEK - 2*XC';
k(:,i) = 3.50e-14;
Gstr{i,1} = 'yMGLY'; Gstr{i,2} = 'RO2XC';
fyMGLY(i)=fyMGLY(i)-1; fRO2XC(i)=fRO2XC(i)-1; fRO2XC(i)=fRO2XC(i)+1; fMEK(i)=fMEK(i)+0.5; fXC(i)=fXC(i)-2;

i=i+1;
Rnames{i} = 'yMGLY + MECO3 = MECO3';
k(:,i) = 4.40e-13.*exp(1070./ T);
Gstr{i,1} = 'yMGLY'; Gstr{i,2} = 'MECO3';
fyMGLY(i)=fyMGLY(i)-1; fMECO3(i)=fMECO3(i)-1; fMECO3(i)=fMECO3(i)+1;

i=i+1;
Rnames{i} = 'yMGLY + RCO3 = RCO3';
k(:,i) = 4.40e-13.*exp(1070./ T);
Gstr{i,1} = 'yMGLY'; Gstr{i,2} = 'RCO3';
fyMGLY(i)=fyMGLY(i)-1; fRCO3(i)=fRCO3(i)-1; fRCO3(i)=fRCO3(i)+1;

i=i+1;
Rnames{i} = 'yMGLY + BZCO3 = BZCO3';
k(:,i) = 4.40e-13.*exp(1070./ T);
Gstr{i,1} = 'yMGLY'; Gstr{i,2} = 'BZCO3';
fyMGLY(i)=fyMGLY(i)-1; fBZCO3(i)=fBZCO3(i)-1; fBZCO3(i)=fBZCO3(i)+1;

i=i+1;
Rnames{i} = 'yMGLY + MACO3 = MACO3';
k(:,i) = 4.40e-13.*exp(1070./ T);
Gstr{i,1} = 'yMGLY'; Gstr{i,2} = 'MACO3';
fyMGLY(i)=fyMGLY(i)-1; fMACO3(i)=fMACO3(i)-1; fMACO3(i)=fMACO3(i)+1;

i=i+1;
Rnames{i} = 'yRNO3 + NO = NO';
k(:,i) = 2.60e-12.*exp(380./ T);
Gstr{i,1} = 'yRNO3'; Gstr{i,2} = 'NO';
fyRNO3(i)=fyRNO3(i)-1; fNO(i)=fNO(i)-1; fNO(i)=fNO(i)+1;

i=i+1;
Rnames{i} = 'yRNO3 + HO2 = HO2 + RNO3 - 6.0*XC';
k(:,i) = 3.80e-13.*exp(900./ T);
Gstr{i,1} = 'yRNO3'; Gstr{i,2} = 'HO2';
fyRNO3(i)=fyRNO3(i)-1; fHO2(i)=fHO2(i)-1; fHO2(i)=fHO2(i)+1; fRNO3(i)=fRNO3(i)+1; fXC(i)=fXC(i)-6.0;

i=i+1;
Rnames{i} = 'yRNO3 + NO3 = NO3';
k(:,i) = 2.30e-12;
Gstr{i,1} = 'yRNO3'; Gstr{i,2} = 'NO3';
fyRNO3(i)=fyRNO3(i)-1; fNO3(i)=fNO3(i)-1; fNO3(i)=fNO3(i)+1;

i=i+1;
Rnames{i} = 'yRNO3 + MEO2 = MEO2 + 0.5*PRD2 - 3*XC';
k(:,i) = 2.00e-13;
Gstr{i,1} = 'yRNO3'; Gstr{i,2} = 'MEO2';
fyRNO3(i)=fyRNO3(i)-1; fMEO2(i)=fMEO2(i)-1; fMEO2(i)=fMEO2(i)+1; fPRD2(i)=fPRD2(i)+0.5; fXC(i)=fXC(i)-3;

i=i+1;
Rnames{i} = 'yRNO3 + RO2C = RO2C + 0.5*PRD2 - 3*XC';
k(:,i) = 3.50e-14;
Gstr{i,1} = 'yRNO3'; Gstr{i,2} = 'RO2C';
fyRNO3(i)=fyRNO3(i)-1; fRO2C(i)=fRO2C(i)-1; fRO2C(i)=fRO2C(i)+1; fPRD2(i)=fPRD2(i)+0.5; fXC(i)=fXC(i)-3;

i=i+1;
Rnames{i} = 'yRNO3 + RO2XC = RO2XC + 0.5*PRD2 - 3*XC';
k(:,i) = 3.50e-14;
Gstr{i,1} = 'yRNO3'; Gstr{i,2} = 'RO2XC';
fyRNO3(i)=fyRNO3(i)-1; fRO2XC(i)=fRO2XC(i)-1; fRO2XC(i)=fRO2XC(i)+1; fPRD2(i)=fPRD2(i)+0.5; fXC(i)=fXC(i)-3;

i=i+1;
Rnames{i} = 'yRNO3 + MECO3 = MECO3';
k(:,i) = 4.40e-13.*exp(1070./ T);
Gstr{i,1} = 'yRNO3'; Gstr{i,2} = 'MECO3';
fyRNO3(i)=fyRNO3(i)-1; fMECO3(i)=fMECO3(i)-1; fMECO3(i)=fMECO3(i)+1;

i=i+1;
Rnames{i} = 'yRNO3 + RCO3 = RCO3';
k(:,i) = 4.40e-13.*exp(1070./ T);
Gstr{i,1} = 'yRNO3'; Gstr{i,2} = 'RCO3';
fyRNO3(i)=fyRNO3(i)-1; fRCO3(i)=fRCO3(i)-1; fRCO3(i)=fRCO3(i)+1;

i=i+1;
Rnames{i} = 'yRNO3 + BZCO3 = BZCO3';
k(:,i) = 4.40e-13.*exp(1070./ T);
Gstr{i,1} = 'yRNO3'; Gstr{i,2} = 'BZCO3';
fyRNO3(i)=fyRNO3(i)-1; fBZCO3(i)=fBZCO3(i)-1; fBZCO3(i)=fBZCO3(i)+1;

i=i+1;
Rnames{i} = 'yRNO3 + MACO3 = MACO3';
k(:,i) = 4.40e-13.*exp(1070./ T);
Gstr{i,1} = 'yRNO3'; Gstr{i,2} = 'MACO3';
fyRNO3(i)=fyRNO3(i)-1; fMACO3(i)=fMACO3(i)-1; fMACO3(i)=fMACO3(i)+1;

i=i+1;
Rnames{i} = 'yAFG3 + NO = NO';
k(:,i) = 2.60e-12.*exp(380./ T);
Gstr{i,1} = 'yAFG3'; Gstr{i,2} = 'NO';
fyAFG3(i)=fyAFG3(i)-1; fNO(i)=fNO(i)-1; fNO(i)=fNO(i)+1;

i=i+1;
Rnames{i} = 'yAFG3 + HO2 = HO2 + AFG3 - 7.0*XC';
k(:,i) = 3.80e-13.*exp(900./ T);
Gstr{i,1} = 'yAFG3'; Gstr{i,2} = 'HO2';
fyAFG3(i)=fyAFG3(i)-1; fHO2(i)=fHO2(i)-1; fHO2(i)=fHO2(i)+1; fAFG3(i)=fAFG3(i)+1; fXC(i)=fXC(i)-7.0;

i=i+1;
Rnames{i} = 'yAFG3 + NO3 = NO3';
k(:,i) = 2.30e-12;
Gstr{i,1} = 'yAFG3'; Gstr{i,2} = 'NO3';
fyAFG3(i)=fyAFG3(i)-1; fNO3(i)=fNO3(i)-1; fNO3(i)=fNO3(i)+1;

i=i+1;
Rnames{i} = 'yAFG3 + MEO2 = MEO2 + 0.5*PRD2 - 3*XC';
k(:,i) = 2.00e-13;
Gstr{i,1} = 'yAFG3'; Gstr{i,2} = 'MEO2';
fyAFG3(i)=fyAFG3(i)-1; fMEO2(i)=fMEO2(i)-1; fMEO2(i)=fMEO2(i)+1; fPRD2(i)=fPRD2(i)+0.5; fXC(i)=fXC(i)-3;

i=i+1;
Rnames{i} = 'yAFG3 + RO2C = RO2C + 0.5*PRD2 - 3*XC';
k(:,i) = 3.50e-14;
Gstr{i,1} = 'yAFG3'; Gstr{i,2} = 'RO2C';
fyAFG3(i)=fyAFG3(i)-1; fRO2C(i)=fRO2C(i)-1; fRO2C(i)=fRO2C(i)+1; fPRD2(i)=fPRD2(i)+0.5; fXC(i)=fXC(i)-3;

i=i+1;
Rnames{i} = 'yAFG3 + RO2XC = RO2XC + 0.5*PRD2 - 3*XC';
k(:,i) = 3.50e-14;
Gstr{i,1} = 'yAFG3'; Gstr{i,2} = 'RO2XC';
fyAFG3(i)=fyAFG3(i)-1; fRO2XC(i)=fRO2XC(i)-1; fRO2XC(i)=fRO2XC(i)+1; fPRD2(i)=fPRD2(i)+0.5; fXC(i)=fXC(i)-3;

i=i+1;
Rnames{i} = 'yAFG3 + MECO3 = MECO3';
k(:,i) = 4.40e-13.*exp(1070./ T);
Gstr{i,1} = 'yAFG3'; Gstr{i,2} = 'MECO3';
fyAFG3(i)=fyAFG3(i)-1; fMECO3(i)=fMECO3(i)-1; fMECO3(i)=fMECO3(i)+1;

i=i+1;
Rnames{i} = 'yAFG3 + RCO3 = RCO3';
k(:,i) = 4.40e-13.*exp(1070./ T);
Gstr{i,1} = 'yAFG3'; Gstr{i,2} = 'RCO3';
fyAFG3(i)=fyAFG3(i)-1; fRCO3(i)=fRCO3(i)-1; fRCO3(i)=fRCO3(i)+1;

i=i+1;
Rnames{i} = 'yAFG3 + BZCO3 = BZCO3';
k(:,i) = 4.40e-13.*exp(1070./ T);
Gstr{i,1} = 'yAFG3'; Gstr{i,2} = 'BZCO3';
fyAFG3(i)=fyAFG3(i)-1; fBZCO3(i)=fBZCO3(i)-1; fBZCO3(i)=fBZCO3(i)+1;

i=i+1;
Rnames{i} = 'yAFG3 + MACO3 = MACO3';
k(:,i) = 4.40e-13.*exp(1070./ T);
Gstr{i,1} = 'yAFG3'; Gstr{i,2} = 'MACO3';
fyAFG3(i)=fyAFG3(i)-1; fMACO3(i)=fMACO3(i)-1; fMACO3(i)=fMACO3(i)+1;

i=i+1;
Rnames{i} = 'yBACL + NO = NO';
k(:,i) = 2.60e-12.*exp(380./ T);
Gstr{i,1} = 'yBACL'; Gstr{i,2} = 'NO';
fyBACL(i)=fyBACL(i)-1; fNO(i)=fNO(i)-1; fNO(i)=fNO(i)+1;

i=i+1;
Rnames{i} = 'yBACL + HO2 = HO2 + BACL - 4.0*XC';
k(:,i) = 3.80e-13.*exp(900./ T);
Gstr{i,1} = 'yBACL'; Gstr{i,2} = 'HO2';
fyBACL(i)=fyBACL(i)-1; fHO2(i)=fHO2(i)-1; fHO2(i)=fHO2(i)+1; fBACL(i)=fBACL(i)+1; fXC(i)=fXC(i)-4.0;

i=i+1;
Rnames{i} = 'yBACL + NO3 = NO3';
k(:,i) = 2.30e-12;
Gstr{i,1} = 'yBACL'; Gstr{i,2} = 'NO3';
fyBACL(i)=fyBACL(i)-1; fNO3(i)=fNO3(i)-1; fNO3(i)=fNO3(i)+1;

i=i+1;
Rnames{i} = 'yBACL + MEO2 = MEO2 + 0.5*MEK - 2*XC';
k(:,i) = 2.00e-13;
Gstr{i,1} = 'yBACL'; Gstr{i,2} = 'MEO2';
fyBACL(i)=fyBACL(i)-1; fMEO2(i)=fMEO2(i)-1; fMEO2(i)=fMEO2(i)+1; fMEK(i)=fMEK(i)+0.5; fXC(i)=fXC(i)-2;

i=i+1;
Rnames{i} = 'yBACL + RO2C = RO2C + 0.5*MEK - 2*XC';
k(:,i) = 3.50e-14;
Gstr{i,1} = 'yBACL'; Gstr{i,2} = 'RO2C';
fyBACL(i)=fyBACL(i)-1; fRO2C(i)=fRO2C(i)-1; fRO2C(i)=fRO2C(i)+1; fMEK(i)=fMEK(i)+0.5; fXC(i)=fXC(i)-2;

i=i+1;
Rnames{i} = 'yBACL + RO2XC = RO2XC + 0.5*MEK - 2*XC';
k(:,i) = 3.50e-14;
Gstr{i,1} = 'yBACL'; Gstr{i,2} = 'RO2XC';
fyBACL(i)=fyBACL(i)-1; fRO2XC(i)=fRO2XC(i)-1; fRO2XC(i)=fRO2XC(i)+1; fMEK(i)=fMEK(i)+0.5; fXC(i)=fXC(i)-2;

i=i+1;
Rnames{i} = 'yBACL + MECO3 = MECO3';
k(:,i) = 4.40e-13.*exp(1070./ T);
Gstr{i,1} = 'yBACL'; Gstr{i,2} = 'MECO3';
fyBACL(i)=fyBACL(i)-1; fMECO3(i)=fMECO3(i)-1; fMECO3(i)=fMECO3(i)+1;

i=i+1;
Rnames{i} = 'yBACL + RCO3 = RCO3';
k(:,i) = 4.40e-13.*exp(1070./ T);
Gstr{i,1} = 'yBACL'; Gstr{i,2} = 'RCO3';
fyBACL(i)=fyBACL(i)-1; fRCO3(i)=fRCO3(i)-1; fRCO3(i)=fRCO3(i)+1;

i=i+1;
Rnames{i} = 'yBACL + BZCO3 = BZCO3';
k(:,i) = 4.40e-13.*exp(1070./ T);
Gstr{i,1} = 'yBACL'; Gstr{i,2} = 'BZCO3';
fyBACL(i)=fyBACL(i)-1; fBZCO3(i)=fBZCO3(i)-1; fBZCO3(i)=fBZCO3(i)+1;

i=i+1;
Rnames{i} = 'yBACL + MACO3 = MACO3';
k(:,i) = 4.40e-13.*exp(1070./ T);
Gstr{i,1} = 'yBACL'; Gstr{i,2} = 'MACO3';
fyBACL(i)=fyBACL(i)-1; fMACO3(i)=fMACO3(i)-1; fMACO3(i)=fMACO3(i)+1;

i=i+1;
Rnames{i} = 'zHFONS + NO = NO + HFONS';
k(:,i) = 2.60e-12.*exp(380./ T);
Gstr{i,1} = 'zHFONS'; Gstr{i,2} = 'NO';
fzHFONS(i)=fzHFONS(i)-1; fNO(i)=fNO(i)-1; fNO(i)=fNO(i)+1; fHFONS(i)=fHFONS(i)+1;

i=i+1;
Rnames{i} = 'zHFONS + HO2 = HO2 + 5.0*XC';
k(:,i) = 3.80e-13.*exp(900./ T);
Gstr{i,1} = 'zHFONS'; Gstr{i,2} = 'HO2';
fzHFONS(i)=fzHFONS(i)-1; fHO2(i)=fHO2(i)-1; fHO2(i)=fHO2(i)+1; fXC(i)=fXC(i)+5.0;

i=i+1;
Rnames{i} = 'zHFONS + NO3 = NO3 + PRD2 + HO2';
k(:,i) = 2.30e-12;
Gstr{i,1} = 'zHFONS'; Gstr{i,2} = 'NO3';
fzHFONS(i)=fzHFONS(i)-1; fNO3(i)=fNO3(i)-1; fNO3(i)=fNO3(i)+1; fPRD2(i)=fPRD2(i)+1; fHO2(i)=fHO2(i)+1;

i=i+1;
Rnames{i} = 'zHFONS + MEO2 = MEO2 + 0.5*PRD2 + 0.5*HO2 + 3*XC';
k(:,i) = 2.00e-13;
Gstr{i,1} = 'zHFONS'; Gstr{i,2} = 'MEO2';
fzHFONS(i)=fzHFONS(i)-1; fMEO2(i)=fMEO2(i)-1; fMEO2(i)=fMEO2(i)+1; fPRD2(i)=fPRD2(i)+0.5; fHO2(i)=fHO2(i)+0.5; fXC(i)=fXC(i)+3;

i=i+1;
Rnames{i} = 'zHFONS + RO2C = RO2C + 0.5*PRD2 + 0.5*HO2 + 3*XC';
k(:,i) = 3.50e-14;
Gstr{i,1} = 'zHFONS'; Gstr{i,2} = 'RO2C';
fzHFONS(i)=fzHFONS(i)-1; fRO2C(i)=fRO2C(i)-1; fRO2C(i)=fRO2C(i)+1; fPRD2(i)=fPRD2(i)+0.5; fHO2(i)=fHO2(i)+0.5; fXC(i)=fXC(i)+3;

i=i+1;
Rnames{i} = 'zHFONS + RO2XC = RO2XC + 0.5*PRD2 + 0.5*HO2 + 3*XC';
k(:,i) = 3.50e-14;
Gstr{i,1} = 'zHFONS'; Gstr{i,2} = 'RO2XC';
fzHFONS(i)=fzHFONS(i)-1; fRO2XC(i)=fRO2XC(i)-1; fRO2XC(i)=fRO2XC(i)+1; fPRD2(i)=fPRD2(i)+0.5; fHO2(i)=fHO2(i)+0.5; fXC(i)=fXC(i)+3;

i=i+1;
Rnames{i} = 'zHFONS + MECO3 = MECO3 + PRD2 + HO2';
k(:,i) = 4.40e-13.*exp(1070./ T);
Gstr{i,1} = 'zHFONS'; Gstr{i,2} = 'MECO3';
fzHFONS(i)=fzHFONS(i)-1; fMECO3(i)=fMECO3(i)-1; fMECO3(i)=fMECO3(i)+1; fPRD2(i)=fPRD2(i)+1; fHO2(i)=fHO2(i)+1;

i=i+1;
Rnames{i} = 'zHFONS + RCO3 = RCO3 + PRD2 + HO2';
k(:,i) = 4.40e-13.*exp(1070./ T);
Gstr{i,1} = 'zHFONS'; Gstr{i,2} = 'RCO3';
fzHFONS(i)=fzHFONS(i)-1; fRCO3(i)=fRCO3(i)-1; fRCO3(i)=fRCO3(i)+1; fPRD2(i)=fPRD2(i)+1; fHO2(i)=fHO2(i)+1;

i=i+1;
Rnames{i} = 'zHFONS + BZCO3 = BZCO3 + PRD2 + HO2';
k(:,i) = 4.40e-13.*exp(1070./ T);
Gstr{i,1} = 'zHFONS'; Gstr{i,2} = 'BZCO3';
fzHFONS(i)=fzHFONS(i)-1; fBZCO3(i)=fBZCO3(i)-1; fBZCO3(i)=fBZCO3(i)+1; fPRD2(i)=fPRD2(i)+1; fHO2(i)=fHO2(i)+1;

i=i+1;
Rnames{i} = 'zHFONS + MACO3 = MACO3 + PRD2 + HO2';
k(:,i) = 4.40e-13.*exp(1070./ T);
Gstr{i,1} = 'zHFONS'; Gstr{i,2} = 'MACO3';
fzHFONS(i)=fzHFONS(i)-1; fMACO3(i)=fMACO3(i)-1; fMACO3(i)=fMACO3(i)+1; fPRD2(i)=fPRD2(i)+1; fHO2(i)=fHO2(i)+1;

