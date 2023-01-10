SpeciesToAdd = {...
'RUOOH'; 'RCHO2'; 'H3XE25DO_P3'; 'HFON52_P1'; 'H3F2ONM5_P1'; 'OTH1'; 'H3F2ONM5_P2'; 'OH'; 'M3_FURAN'; ...
'H3F2ONE_R1'; 'O4X2PEAL_R1'; 'H3XE25DO'; 'H3F2ONM5_P3'; 'CROOH'; 'HFON52M3_P1'; 'M3_FURAN_P1'; 'XC'; 'M2BUTDAL_P2'; 'OLEP'; ...
'OLEA1'; 'M25_FUR'; 'HPALD'; 'M2BUTDAL_R2'; 'ROOH'; 'MEO2'; 'RDNO3'; 'RO2XC'; 'M2_FURAN_P2'; 'M25_FUR_P2'; ...
'M25_FUR_P5'; 'XN'; 'HFON52M4'; 'H3XE25DO_P1'; 'SumRCO3'; 'HFON52M5'; 'NO'; 'HFON52M3_P3'; 'M2_FURAN_R1'; 'M25_FUR_R1'; ...
'RCNO3'; 'AFG3'; 'O4X2PEAL'; 'M2BUTDAL_R1'; 'M2_FURAN_P3'; 'M2_FURAN_P4'; 'H3F2ONE_P1'; 'R2CO3'; 'FURAN'; 'M25_FUR_P1'; ...
'HFON52M3_P2'; 'H3F2ONM5_P4'; 'FURAN_P2'; 'NROG'; 'RCHO'; 'CO2'; 'M25_FUR_P3'; 'GLY'; 'M2BUTDAL_P1'; 'HFON52M4_P2'; ...
'CO'; 'M25_FUR_P4'; 'O4X2PEAL_P3'; 'HFON52M3'; 'BUTEDIAL_R1'; 'O4X2PEAL_P4'; 'H3F2ONM4_P3'; 'HO2'; 'M3_FURAN_P2'; 'zRCNO3'; ...
'IEPOX'; 'FURAN_P3'; 'FURAN_P4'; 'H3F2ONE_P2'; 'H3F2ONM4_P2'; 'HFON52'; 'H3F2ONM5'; 'HFON52_P2'; 'M2BUTDAL_R3'; 'zRHNO3'; ...
'HFON52M4_P1'; 'H3F2ONE_P4'; 'HCHO'; 'O3'; 'HFON52M5_P2'; 'HFON52M5_P1'; 'RO2C'; 'SumRO2'; 'FURAN_P1'; 'H3F2ONM4_P1'; ...
'H3XE25DO_P2'; 'KET2'; 'M2_FURAN'; 'M2BUTDAL'; 'O4X2PEAL_P2'; 'H3F2ONE'; 'H3F2ONM5_R1'; 'O4X2PEAL_P1'; 'MGLY'; 'NO3'; ...
'HV'; 'RHNO3'; 'MECO3'; 'BACL'; 'zRDNO3'; 'MALAH'; 'H3F2ONE_P3'; 'M2_FURAN_P5'; 'MACO3'; 'BUTEDIAL'; ...
'M2_FURAN_P1'; 'OTH3'; 'H3F2ONM4';};

AddSpecies

i=i+1;
Rnames{i} = 'FURAN + OH = .76*BUTEDIAL + .76*HO2 + .24*FURAN_P1 + .24*SumRO2';
k(:,i) = 1.30e-11.*(T./300).^0.00.*exp(333.132./T);
Gstr{i,1} = 'FURAN'; Gstr{i,2} = 'OH';
fFURAN(i)=fFURAN(i)-1; fOH(i)=fOH(i)-1; fBUTEDIAL(i)=fBUTEDIAL(i)+.76; fHO2(i)=fHO2(i)+.76; fFURAN_P1(i)=fFURAN_P1(i)+.24; fSumRO2(i)=fSumRO2(i)+.24;

i=i+1;
Rnames{i} = 'FURAN + O3 = RCHO';
k(:,i) = 2.42e-18;
Gstr{i,1} = 'FURAN'; Gstr{i,2} = 'O3';
fFURAN(i)=fFURAN(i)-1; fO3(i)=fO3(i)-1; fRCHO(i)=fRCHO(i)+1;

i=i+1;
Rnames{i} = 'FURAN + NO3 = FURAN_P2 + SumRO2';
k(:,i) = 1.30e-13.*(T./300).^0.00.*exp(699.98./T);
Gstr{i,1} = 'FURAN'; Gstr{i,2} = 'NO3';
fFURAN(i)=fFURAN(i)-1; fNO3(i)=fNO3(i)-1; fFURAN_P2(i)=fFURAN_P2(i)+1; fSumRO2(i)=fSumRO2(i)+1;

i=i+1;
Rnames{i} = 'FURAN_P1 + NO = .974*NO2 + .652*HFON52 + .652*HO2 + .323*FURAN_P3 + .026*RHNO3 - 0.108*XC + .323*SumRO2';
k(:,i) = 9.13e-12;
Gstr{i,1} = 'FURAN_P1'; Gstr{i,2} = 'NO';
fFURAN_P1(i)=fFURAN_P1(i)-1; fNO(i)=fNO(i)-1; fNO2(i)=fNO2(i)+.974; fHFON52(i)=fHFON52(i)+.652; fHO2(i)=fHO2(i)+.652; fFURAN_P3(i)=fFURAN_P3(i)+.323; fRHNO3(i)=fRHNO3(i)+.026; fXC(i)=fXC(i)-0.108; fSumRO2(i)=fSumRO2(i)+.323;

i=i+1;
Rnames{i} = 'FURAN_P1 + NO3 = NO2 + .669*HFON52 + .669*HO2 + .331*FURAN_P3 + .331*SumRO2';
k(:,i) = 2.30e-12;
Gstr{i,1} = 'FURAN_P1'; Gstr{i,2} = 'NO3';
fFURAN_P1(i)=fFURAN_P1(i)-1; fNO3(i)=fNO3(i)-1; fNO2(i)=fNO2(i)+1; fHFON52(i)=fHFON52(i)+.669; fHO2(i)=fHO2(i)+.669; fFURAN_P3(i)=fFURAN_P3(i)+.331; fSumRO2(i)=fSumRO2(i)+.331;

i=i+1;
Rnames{i} = 'FURAN_P1 + HO2 = .9*RUOOH + .1*HFON52 - 1.8*XC';
k(:,i) = 1.61e-11;
Gstr{i,1} = 'FURAN_P1'; Gstr{i,2} = 'HO2';
fFURAN_P1(i)=fFURAN_P1(i)-1; fHO2(i)=fHO2(i)-1; fRUOOH(i)=fRUOOH(i)+.9; fHFON52(i)=fHFON52(i)+.1; fXC(i)=fXC(i)-1.8;

i=i+1;
Rnames{i} = 'FURAN_P1 + SumRO2 = SumRO2 + .585*HFON52 + .335*HO2 + .25*OLEP + .166*FURAN_P3 - 1.004*XC + .166*SumRO2';
k(:,i) = 2.57e-12;
Gstr{i,1} = 'FURAN_P1'; Gstr{i,2} = 'SumRO2';
fFURAN_P1(i)=fFURAN_P1(i)-1; fSumRO2(i)=fSumRO2(i)-1; fSumRO2(i)=fSumRO2(i)+1; fHFON52(i)=fHFON52(i)+.585; fHO2(i)=fHO2(i)+.335; fOLEP(i)=fOLEP(i)+.25; fFURAN_P3(i)=fFURAN_P3(i)+.166; fXC(i)=fXC(i)-1.004; fSumRO2(i)=fSumRO2(i)+.166;

i=i+1;
Rnames{i} = 'FURAN_P1 + SumRCO3 = SumRCO3 + .735*HFON52 + .535*HO2 + .265*FURAN_P3 + .265*SumRO2';
k(:,i) = 1.37e-11;
Gstr{i,1} = 'FURAN_P1'; Gstr{i,2} = 'SumRCO3';
fFURAN_P1(i)=fFURAN_P1(i)-1; fSumRCO3(i)=fSumRCO3(i)-1; fSumRCO3(i)=fSumRCO3(i)+1; fHFON52(i)=fHFON52(i)+.735; fHO2(i)=fHO2(i)+.535; fFURAN_P3(i)=fFURAN_P3(i)+.265; fSumRO2(i)=fSumRO2(i)+.265;

i=i+1;
Rnames{i} = 'FURAN_P2 + NO = .974*NO2 + .653*RCNO3 + .653*HO2 + .322*FURAN_P4 + .026*RDNO3 - 0.108*XC - 0.001*XN + .322*SumRO2';
k(:,i) = 9.13e-12;
Gstr{i,1} = 'FURAN_P2'; Gstr{i,2} = 'NO';
fFURAN_P2(i)=fFURAN_P2(i)-1; fNO(i)=fNO(i)-1; fNO2(i)=fNO2(i)+.974; fRCNO3(i)=fRCNO3(i)+.653; fHO2(i)=fHO2(i)+.653; fFURAN_P4(i)=fFURAN_P4(i)+.322; fRDNO3(i)=fRDNO3(i)+.026; fXC(i)=fXC(i)-0.108; fXN(i)=fXN(i)-0.001; fSumRO2(i)=fSumRO2(i)+.322;

i=i+1;
Rnames{i} = 'FURAN_P2 + NO3 = NO2 + .67*RCNO3 + .67*HO2 + .33*FURAN_P4 + .33*SumRO2';
k(:,i) = 2.30e-12;
Gstr{i,1} = 'FURAN_P2'; Gstr{i,2} = 'NO3';
fFURAN_P2(i)=fFURAN_P2(i)-1; fNO3(i)=fNO3(i)-1; fNO2(i)=fNO2(i)+1; fRCNO3(i)=fRCNO3(i)+.67; fHO2(i)=fHO2(i)+.67; fFURAN_P4(i)=fFURAN_P4(i)+.33; fSumRO2(i)=fSumRO2(i)+.33;

i=i+1;
Rnames{i} = 'FURAN_P2 + HO2 = .9*RHNO3 + .1*RCNO3 - 3.6*XC';
k(:,i) = 1.61e-11;
Gstr{i,1} = 'FURAN_P2'; Gstr{i,2} = 'HO2';
fFURAN_P2(i)=fFURAN_P2(i)-1; fHO2(i)=fHO2(i)-1; fRHNO3(i)=fRHNO3(i)+.9; fRCNO3(i)=fRCNO3(i)+.1; fXC(i)=fXC(i)-3.6;

i=i+1;
Rnames{i} = 'FURAN_P2 + SumRO2 = SumRO2 + .585*RCNO3 + .335*HO2 + .25*RHNO3 + .165*FURAN_P4 - 1*XC + .165*SumRO2';
k(:,i) = 2.57e-12;
Gstr{i,1} = 'FURAN_P2'; Gstr{i,2} = 'SumRO2';
fFURAN_P2(i)=fFURAN_P2(i)-1; fSumRO2(i)=fSumRO2(i)-1; fSumRO2(i)=fSumRO2(i)+1; fRCNO3(i)=fRCNO3(i)+.585; fHO2(i)=fHO2(i)+.335; fRHNO3(i)=fRHNO3(i)+.25; fFURAN_P4(i)=fFURAN_P4(i)+.165; fXC(i)=fXC(i)-1; fSumRO2(i)=fSumRO2(i)+.165;

i=i+1;
Rnames{i} = 'FURAN_P2 + SumRCO3 = SumRCO3 + .736*RCNO3 + .536*HO2 + .264*FURAN_P4 + .264*SumRO2';
k(:,i) = 1.37e-11;
Gstr{i,1} = 'FURAN_P2'; Gstr{i,2} = 'SumRCO3';
fFURAN_P2(i)=fFURAN_P2(i)-1; fSumRCO3(i)=fSumRCO3(i)-1; fSumRCO3(i)=fSumRCO3(i)+1; fRCNO3(i)=fRCNO3(i)+.736; fHO2(i)=fHO2(i)+.536; fFURAN_P4(i)=fFURAN_P4(i)+.264; fSumRO2(i)=fSumRO2(i)+.264;

i=i+1;
Rnames{i} = 'FURAN_P3 + NO = .974*RCHO + .974*HO2 + .974*NO2 + .026*RHNO3 - 0.104*XC';
k(:,i) = 9.13e-12;
Gstr{i,1} = 'FURAN_P3'; Gstr{i,2} = 'NO';
fFURAN_P3(i)=fFURAN_P3(i)-1; fNO(i)=fNO(i)-1; fRCHO(i)=fRCHO(i)+.974; fHO2(i)=fHO2(i)+.974; fNO2(i)=fNO2(i)+.974; fRHNO3(i)=fRHNO3(i)+.026; fXC(i)=fXC(i)-0.104;

i=i+1;
Rnames{i} = 'FURAN_P3 + NO3 = RCHO + HO2 + NO2';
k(:,i) = 2.30e-12;
Gstr{i,1} = 'FURAN_P3'; Gstr{i,2} = 'NO3';
fFURAN_P3(i)=fFURAN_P3(i)-1; fNO3(i)=fNO3(i)-1; fRCHO(i)=fRCHO(i)+1; fHO2(i)=fHO2(i)+1; fNO2(i)=fNO2(i)+1;

i=i+1;
Rnames{i} = 'FURAN_P3 + HO2 = ROOH - 1*XC';
k(:,i) = 1.61e-11;
Gstr{i,1} = 'FURAN_P3'; Gstr{i,2} = 'HO2';
fFURAN_P3(i)=fFURAN_P3(i)-1; fHO2(i)=fHO2(i)-1; fROOH(i)=fROOH(i)+1; fXC(i)=fXC(i)-1;

i=i+1;
Rnames{i} = 'FURAN_P3 + SumRO2 = SumRO2 + .5*RCHO + .5*HO2 + .25*KET2 + .25*IEPOX - 0.75*XC';
k(:,i) = 2.57e-12;
Gstr{i,1} = 'FURAN_P3'; Gstr{i,2} = 'SumRO2';
fFURAN_P3(i)=fFURAN_P3(i)-1; fSumRO2(i)=fSumRO2(i)-1; fSumRO2(i)=fSumRO2(i)+1; fRCHO(i)=fRCHO(i)+.5; fHO2(i)=fHO2(i)+.5; fKET2(i)=fKET2(i)+.25; fIEPOX(i)=fIEPOX(i)+.25; fXC(i)=fXC(i)-0.75;

i=i+1;
Rnames{i} = 'FURAN_P3 + SumRCO3 = SumRCO3 + .8*RCHO + .8*HO2 + .2*KET2 - 0.4*XC';
k(:,i) = 1.37e-11;
Gstr{i,1} = 'FURAN_P3'; Gstr{i,2} = 'SumRCO3';
fFURAN_P3(i)=fFURAN_P3(i)-1; fSumRCO3(i)=fSumRCO3(i)-1; fSumRCO3(i)=fSumRCO3(i)+1; fRCHO(i)=fRCHO(i)+.8; fHO2(i)=fHO2(i)+.8; fKET2(i)=fKET2(i)+.2; fXC(i)=fXC(i)-0.4;

i=i+1;
Rnames{i} = 'FURAN_P4 + NO = 1.949*NO2 + .974*RCHO + .026*RDNO3 - 0.104*XC - 0.001*XN';
k(:,i) = 9.13e-12;
Gstr{i,1} = 'FURAN_P4'; Gstr{i,2} = 'NO';
fFURAN_P4(i)=fFURAN_P4(i)-1; fNO(i)=fNO(i)-1; fNO2(i)=fNO2(i)+1.949; fRCHO(i)=fRCHO(i)+.974; fRDNO3(i)=fRDNO3(i)+.026; fXC(i)=fXC(i)-0.104; fXN(i)=fXN(i)-0.001;

i=i+1;
Rnames{i} = 'FURAN_P4 + NO3 = 2*NO2 + RCHO';
k(:,i) = 2.30e-12;
Gstr{i,1} = 'FURAN_P4'; Gstr{i,2} = 'NO3';
fFURAN_P4(i)=fFURAN_P4(i)-1; fNO3(i)=fNO3(i)-1; fNO2(i)=fNO2(i)+2; fRCHO(i)=fRCHO(i)+1;

i=i+1;
Rnames{i} = 'FURAN_P4 + HO2 = RHNO3 - 4*XC';
k(:,i) = 1.61e-11;
Gstr{i,1} = 'FURAN_P4'; Gstr{i,2} = 'HO2';
fFURAN_P4(i)=fFURAN_P4(i)-1; fHO2(i)=fHO2(i)-1; fRHNO3(i)=fRHNO3(i)+1; fXC(i)=fXC(i)-4;

i=i+1;
Rnames{i} = 'FURAN_P4 + SumRO2 = SumRO2 + .5*RCHO + .5*NO2 + .25*RCNO3 + .25*RHNO3 - 1*XC';
k(:,i) = 2.57e-12;
Gstr{i,1} = 'FURAN_P4'; Gstr{i,2} = 'SumRO2';
fFURAN_P4(i)=fFURAN_P4(i)-1; fSumRO2(i)=fSumRO2(i)-1; fSumRO2(i)=fSumRO2(i)+1; fRCHO(i)=fRCHO(i)+.5; fNO2(i)=fNO2(i)+.5; fRCNO3(i)=fRCNO3(i)+.25; fRHNO3(i)=fRHNO3(i)+.25; fXC(i)=fXC(i)-1;

i=i+1;
Rnames{i} = 'FURAN_P4 + SumRCO3 = SumRCO3 + .8*RCHO + .8*NO2 + .2*RCNO3';
k(:,i) = 1.37e-11;
Gstr{i,1} = 'FURAN_P4'; Gstr{i,2} = 'SumRCO3';
fFURAN_P4(i)=fFURAN_P4(i)-1; fSumRCO3(i)=fSumRCO3(i)-1; fSumRCO3(i)=fSumRCO3(i)+1; fRCHO(i)=fRCHO(i)+.8; fNO2(i)=fNO2(i)+.8; fRCNO3(i)=fRCNO3(i)+.2;

i=i+1;
Rnames{i} = 'M2_FURAN + OH = .56*M2_FURAN_P1 + .44*O4X2PEAL + .44*HO2 + .56*SumRO2';
k(:,i) = 7.31e-11;
Gstr{i,1} = 'M2_FURAN'; Gstr{i,2} = 'OH';
fM2_FURAN(i)=fM2_FURAN(i)-1; fOH(i)=fOH(i)-1; fM2_FURAN_P1(i)=fM2_FURAN_P1(i)+.56; fO4X2PEAL(i)=fO4X2PEAL(i)+.44; fHO2(i)=fHO2(i)+.44; fSumRO2(i)=fSumRO2(i)+.56;

i=i+1;
Rnames{i} = 'M2_FURAN + O3 = .613*OH + .591*RCHO + .204*OLEP + .204*M2_FURAN_R1 - 0.016*XC';
k(:,i) = 2.00e-17;
Gstr{i,1} = 'M2_FURAN'; Gstr{i,2} = 'O3';
fM2_FURAN(i)=fM2_FURAN(i)-1; fO3(i)=fO3(i)-1; fOH(i)=fOH(i)+.613; fRCHO(i)=fRCHO(i)+.591; fOLEP(i)=fOLEP(i)+.204; fM2_FURAN_R1(i)=fM2_FURAN_R1(i)+.204; fXC(i)=fXC(i)-0.016;

i=i+1;
Rnames{i} = 'M2_FURAN + NO3 = M2_FURAN_P2 + SumRO2';
k(:,i) = 2.57e-11;
Gstr{i,1} = 'M2_FURAN'; Gstr{i,2} = 'NO3';
fM2_FURAN(i)=fM2_FURAN(i)-1; fNO3(i)=fNO3(i)-1; fM2_FURAN_P2(i)=fM2_FURAN_P2(i)+1; fSumRO2(i)=fSumRO2(i)+1;

i=i+1;
Rnames{i} = 'M2_FURAN_P1 + NO = .938*NO2 + .462*HFON52M5 + .462*HO2 + .309*HFON52 + .309*MEO2 + .167*M2_FURAN_P3 + .062*RHNO3 - 0.186*XC + .476*SumRO2';
k(:,i) = 9.13e-12;
Gstr{i,1} = 'M2_FURAN_P1'; Gstr{i,2} = 'NO';
fM2_FURAN_P1(i)=fM2_FURAN_P1(i)-1; fNO(i)=fNO(i)-1; fNO2(i)=fNO2(i)+.938; fHFON52M5(i)=fHFON52M5(i)+.462; fHO2(i)=fHO2(i)+.462; fHFON52(i)=fHFON52(i)+.309; fMEO2(i)=fMEO2(i)+.309; fM2_FURAN_P3(i)=fM2_FURAN_P3(i)+.167; fRHNO3(i)=fRHNO3(i)+.062; fXC(i)=fXC(i)-0.186; fSumRO2(i)=fSumRO2(i)+.476;

i=i+1;
Rnames{i} = 'M2_FURAN_P1 + NO3 = NO2 + .493*HFON52M5 + .493*HO2 + .33*HFON52 + .33*MEO2 + .177*M2_FURAN_P3 + .507*SumRO2';
k(:,i) = 2.30e-12;
Gstr{i,1} = 'M2_FURAN_P1'; Gstr{i,2} = 'NO3';
fM2_FURAN_P1(i)=fM2_FURAN_P1(i)-1; fNO3(i)=fNO3(i)-1; fNO2(i)=fNO2(i)+1; fHFON52M5(i)=fHFON52M5(i)+.493; fHO2(i)=fHO2(i)+.493; fHFON52(i)=fHFON52(i)+.33; fMEO2(i)=fMEO2(i)+.33; fM2_FURAN_P3(i)=fM2_FURAN_P3(i)+.177; fSumRO2(i)=fSumRO2(i)+.507;

i=i+1;
Rnames{i} = 'M2_FURAN_P1 + HO2 = .946*RUOOH + .054*HFON52M5 - 0.946*XC';
k(:,i) = 1.78e-11;
Gstr{i,1} = 'M2_FURAN_P1'; Gstr{i,2} = 'HO2';
fM2_FURAN_P1(i)=fM2_FURAN_P1(i)-1; fHO2(i)=fHO2(i)-1; fRUOOH(i)=fRUOOH(i)+.946; fHFON52M5(i)=fHFON52M5(i)+.054; fXC(i)=fXC(i)-0.946;

i=i+1;
Rnames{i} = 'M2_FURAN_P1 + SumRO2 = SumRO2 + .38*HFON52M5 + .366*OLEP + .246*HO2 + .165*HFON52 + .165*MEO2 + .089*M2_FURAN_P3 - 1.098*XC + .254*SumRO2';
k(:,i) = 1.38e-12;
Gstr{i,1} = 'M2_FURAN_P1'; Gstr{i,2} = 'SumRO2';
fM2_FURAN_P1(i)=fM2_FURAN_P1(i)-1; fSumRO2(i)=fSumRO2(i)-1; fSumRO2(i)=fSumRO2(i)+1; fHFON52M5(i)=fHFON52M5(i)+.38; fOLEP(i)=fOLEP(i)+.366; fHO2(i)=fHO2(i)+.246; fHFON52(i)=fHFON52(i)+.165; fMEO2(i)=fMEO2(i)+.165; fM2_FURAN_P3(i)=fM2_FURAN_P3(i)+.089; fXC(i)=fXC(i)-1.098; fSumRO2(i)=fSumRO2(i)+.254;

i=i+1;
Rnames{i} = 'M2_FURAN_P1 + SumRCO3 = SumRCO3 + .501*HFON52M5 + .394*HO2 + .33*HFON52 + .33*MEO2 + .169*M2_FURAN_P3 + .499*SumRO2';
k(:,i) = 1.37e-11;
Gstr{i,1} = 'M2_FURAN_P1'; Gstr{i,2} = 'SumRCO3';
fM2_FURAN_P1(i)=fM2_FURAN_P1(i)-1; fSumRCO3(i)=fSumRCO3(i)-1; fSumRCO3(i)=fSumRCO3(i)+1; fHFON52M5(i)=fHFON52M5(i)+.501; fHO2(i)=fHO2(i)+.394; fHFON52(i)=fHFON52(i)+.33; fMEO2(i)=fMEO2(i)+.33; fM2_FURAN_P3(i)=fM2_FURAN_P3(i)+.169; fSumRO2(i)=fSumRO2(i)+.499;

i=i+1;
Rnames{i} = 'M2_FURAN_R1 = M2_FURAN_P4 + SumRO2';
k(:,i) = 1.64e+11.*(T./300).^0.00.*exp(-7848.229./T);
Gstr{i,1} = 'M2_FURAN_R1';
fM2_FURAN_R1(i)=fM2_FURAN_R1(i)-1; fM2_FURAN_P4(i)=fM2_FURAN_P4(i)+1; fSumRO2(i)=fSumRO2(i)+1;

i=i+1;
Rnames{i} = 'M2_FURAN_R1 + NO = .938*OLEA1 + .938*CO + .938*HO2 + .938*NO2 + .062*RCNO3 - 1.814*XC';
k(:,i) = 9.13e-12;
Gstr{i,1} = 'M2_FURAN_R1'; Gstr{i,2} = 'NO';
fM2_FURAN_R1(i)=fM2_FURAN_R1(i)-1; fNO(i)=fNO(i)-1; fOLEA1(i)=fOLEA1(i)+.938; fCO(i)=fCO(i)+.938; fHO2(i)=fHO2(i)+.938; fNO2(i)=fNO2(i)+.938; fRCNO3(i)=fRCNO3(i)+.062; fXC(i)=fXC(i)-1.814;

i=i+1;
Rnames{i} = 'M2_FURAN_P2 + NO = .938*NO2 + .772*RCNO3 + .463*HO2 + .309*MEO2 + .166*M2_FURAN_P5 + .062*RDNO3 + .277*XC + .475*SumRO2';
k(:,i) = 9.13e-12;
Gstr{i,1} = 'M2_FURAN_P2'; Gstr{i,2} = 'NO';
fM2_FURAN_P2(i)=fM2_FURAN_P2(i)-1; fNO(i)=fNO(i)-1; fNO2(i)=fNO2(i)+.938; fRCNO3(i)=fRCNO3(i)+.772; fHO2(i)=fHO2(i)+.463; fMEO2(i)=fMEO2(i)+.309; fM2_FURAN_P5(i)=fM2_FURAN_P5(i)+.166; fRDNO3(i)=fRDNO3(i)+.062; fXC(i)=fXC(i)+.277; fSumRO2(i)=fSumRO2(i)+.475;

i=i+1;
Rnames{i} = 'M2_FURAN_P2 + NO3 = NO2 + .823*RCNO3 + .493*HO2 + .329*MEO2 + .177*M2_FURAN_P5 + .494*XC + .506*SumRO2';
k(:,i) = 2.30e-12;
Gstr{i,1} = 'M2_FURAN_P2'; Gstr{i,2} = 'NO3';
fM2_FURAN_P2(i)=fM2_FURAN_P2(i)-1; fNO3(i)=fNO3(i)-1; fNO2(i)=fNO2(i)+1; fRCNO3(i)=fRCNO3(i)+.823; fHO2(i)=fHO2(i)+.493; fMEO2(i)=fMEO2(i)+.329; fM2_FURAN_P5(i)=fM2_FURAN_P5(i)+.177; fXC(i)=fXC(i)+.494; fSumRO2(i)=fSumRO2(i)+.506;

i=i+1;
Rnames{i} = 'M2_FURAN_P2 + HO2 = .946*RHNO3 + .054*RCNO3 - 2.784*XC';
k(:,i) = 1.78e-11;
Gstr{i,1} = 'M2_FURAN_P2'; Gstr{i,2} = 'HO2';
fM2_FURAN_P2(i)=fM2_FURAN_P2(i)-1; fHO2(i)=fHO2(i)-1; fRHNO3(i)=fRHNO3(i)+.946; fRCNO3(i)=fRCNO3(i)+.054; fXC(i)=fXC(i)-2.784;

i=i+1;
Rnames{i} = 'M2_FURAN_P2 + SumRO2 = SumRO2 + .545*RCNO3 + .366*RHNO3 + .247*HO2 + .165*MEO2 + .089*M2_FURAN_P5 - 0.718*XC + .254*SumRO2';
k(:,i) = 1.38e-12;
Gstr{i,1} = 'M2_FURAN_P2'; Gstr{i,2} = 'SumRO2';
fM2_FURAN_P2(i)=fM2_FURAN_P2(i)-1; fSumRO2(i)=fSumRO2(i)-1; fSumRO2(i)=fSumRO2(i)+1; fRCNO3(i)=fRCNO3(i)+.545; fRHNO3(i)=fRHNO3(i)+.366; fHO2(i)=fHO2(i)+.247; fMEO2(i)=fMEO2(i)+.165; fM2_FURAN_P5(i)=fM2_FURAN_P5(i)+.089; fXC(i)=fXC(i)-0.718; fSumRO2(i)=fSumRO2(i)+.254;

i=i+1;
Rnames{i} = 'M2_FURAN_P2 + SumRCO3 = SumRCO3 + .831*RCNO3 + .394*HO2 + .329*MEO2 + .169*M2_FURAN_P5 + .502*XC + .498*SumRO2';
k(:,i) = 1.37e-11;
Gstr{i,1} = 'M2_FURAN_P2'; Gstr{i,2} = 'SumRCO3';
fM2_FURAN_P2(i)=fM2_FURAN_P2(i)-1; fSumRCO3(i)=fSumRCO3(i)-1; fSumRCO3(i)=fSumRCO3(i)+1; fRCNO3(i)=fRCNO3(i)+.831; fHO2(i)=fHO2(i)+.394; fMEO2(i)=fMEO2(i)+.329; fM2_FURAN_P5(i)=fM2_FURAN_P5(i)+.169; fXC(i)=fXC(i)+.502; fSumRO2(i)=fSumRO2(i)+.498;

i=i+1;
Rnames{i} = 'M2_FURAN_P3 + NO = .938*RCHO + .938*HO2 + .938*NO2 + .062*RHNO3 + .752*XC';
k(:,i) = 9.13e-12;
Gstr{i,1} = 'M2_FURAN_P3'; Gstr{i,2} = 'NO';
fM2_FURAN_P3(i)=fM2_FURAN_P3(i)-1; fNO(i)=fNO(i)-1; fRCHO(i)=fRCHO(i)+.938; fHO2(i)=fHO2(i)+.938; fNO2(i)=fNO2(i)+.938; fRHNO3(i)=fRHNO3(i)+.062; fXC(i)=fXC(i)+.752;

i=i+1;
Rnames{i} = 'M2_FURAN_P3 + NO3 = RCHO + HO2 + NO2 + XC';
k(:,i) = 2.30e-12;
Gstr{i,1} = 'M2_FURAN_P3'; Gstr{i,2} = 'NO3';
fM2_FURAN_P3(i)=fM2_FURAN_P3(i)-1; fNO3(i)=fNO3(i)-1; fRCHO(i)=fRCHO(i)+1; fHO2(i)=fHO2(i)+1; fNO2(i)=fNO2(i)+1; fXC(i)=fXC(i)+1;

i=i+1;
Rnames{i} = 'M2_FURAN_P3 + HO2 = ROOH';
k(:,i) = 1.78e-11;
Gstr{i,1} = 'M2_FURAN_P3'; Gstr{i,2} = 'HO2';
fM2_FURAN_P3(i)=fM2_FURAN_P3(i)-1; fHO2(i)=fHO2(i)-1; fROOH(i)=fROOH(i)+1;

i=i+1;
Rnames{i} = 'M2_FURAN_P3 + SumRO2 = SumRO2 + .5*RCHO + .5*HO2 + .25*KET2 + .25*IEPOX + .25*XC';
k(:,i) = 2.57e-12;
Gstr{i,1} = 'M2_FURAN_P3'; Gstr{i,2} = 'SumRO2';
fM2_FURAN_P3(i)=fM2_FURAN_P3(i)-1; fSumRO2(i)=fSumRO2(i)-1; fSumRO2(i)=fSumRO2(i)+1; fRCHO(i)=fRCHO(i)+.5; fHO2(i)=fHO2(i)+.5; fKET2(i)=fKET2(i)+.25; fIEPOX(i)=fIEPOX(i)+.25; fXC(i)=fXC(i)+.25;

i=i+1;
Rnames{i} = 'M2_FURAN_P3 + SumRCO3 = SumRCO3 + .8*RCHO + .8*HO2 + .2*KET2 + .6*XC';
k(:,i) = 1.37e-11;
Gstr{i,1} = 'M2_FURAN_P3'; Gstr{i,2} = 'SumRCO3';
fM2_FURAN_P3(i)=fM2_FURAN_P3(i)-1; fSumRCO3(i)=fSumRCO3(i)-1; fSumRCO3(i)=fSumRCO3(i)+1; fRCHO(i)=fRCHO(i)+.8; fHO2(i)=fHO2(i)+.8; fKET2(i)=fKET2(i)+.2; fXC(i)=fXC(i)+.6;

i=i+1;
Rnames{i} = 'M2_FURAN_P4 + NO = .938*OLEA1 + .938*CO2 + .938*OH + .938*NO2 + .062*RCNO3 - 1.814*XC';
k(:,i) = 9.13e-12;
Gstr{i,1} = 'M2_FURAN_P4'; Gstr{i,2} = 'NO';
fM2_FURAN_P4(i)=fM2_FURAN_P4(i)-1; fNO(i)=fNO(i)-1; fOLEA1(i)=fOLEA1(i)+.938; fCO2(i)=fCO2(i)+.938; fOH(i)=fOH(i)+.938; fNO2(i)=fNO2(i)+.938; fRCNO3(i)=fRCNO3(i)+.062; fXC(i)=fXC(i)-1.814;

i=i+1;
Rnames{i} = 'M2_FURAN_P4 + NO3 = OLEA1 + CO2 + OH + NO2 - 2*XC';
k(:,i) = 2.30e-12;
Gstr{i,1} = 'M2_FURAN_P4'; Gstr{i,2} = 'NO3';
fM2_FURAN_P4(i)=fM2_FURAN_P4(i)-1; fNO3(i)=fNO3(i)-1; fOLEA1(i)=fOLEA1(i)+1; fCO2(i)=fCO2(i)+1; fOH(i)=fOH(i)+1; fNO2(i)=fNO2(i)+1; fXC(i)=fXC(i)-2;

i=i+1;
Rnames{i} = 'M2_FURAN_P4 + HO2 = .9*RUOOH + .1*BACL - 0.8*XC';
k(:,i) = 1.78e-11;
Gstr{i,1} = 'M2_FURAN_P4'; Gstr{i,2} = 'HO2';
fM2_FURAN_P4(i)=fM2_FURAN_P4(i)-1; fHO2(i)=fHO2(i)-1; fRUOOH(i)=fRUOOH(i)+.9; fBACL(i)=fBACL(i)+.1; fXC(i)=fXC(i)-0.8;

i=i+1;
Rnames{i} = 'M2_FURAN_P4 + SumRO2 = SumRO2 + .5*OLEA1 + .5*CO2 + .5*OH + .25*BACL + .25*RUOOH - 1*XC';
k(:,i) = 2.57e-12;
Gstr{i,1} = 'M2_FURAN_P4'; Gstr{i,2} = 'SumRO2';
fM2_FURAN_P4(i)=fM2_FURAN_P4(i)-1; fSumRO2(i)=fSumRO2(i)-1; fSumRO2(i)=fSumRO2(i)+1; fOLEA1(i)=fOLEA1(i)+.5; fCO2(i)=fCO2(i)+.5; fOH(i)=fOH(i)+.5; fBACL(i)=fBACL(i)+.25; fRUOOH(i)=fRUOOH(i)+.25; fXC(i)=fXC(i)-1;

i=i+1;
Rnames{i} = 'M2_FURAN_P4 + SumRCO3 = SumRCO3 + .8*OLEA1 + .8*CO2 + .8*OH + .2*BACL - 1.4*XC';
k(:,i) = 1.37e-11;
Gstr{i,1} = 'M2_FURAN_P4'; Gstr{i,2} = 'SumRCO3';
fM2_FURAN_P4(i)=fM2_FURAN_P4(i)-1; fSumRCO3(i)=fSumRCO3(i)-1; fSumRCO3(i)=fSumRCO3(i)+1; fOLEA1(i)=fOLEA1(i)+.8; fCO2(i)=fCO2(i)+.8; fOH(i)=fOH(i)+.8; fBACL(i)=fBACL(i)+.2; fXC(i)=fXC(i)-1.4;

i=i+1;
Rnames{i} = 'M2_FURAN_P5 + NO = 1.82*NO2 + .882*RCHO + .062*RDNO3 + .056*RCNO3 + .056*OH + .752*XC';
k(:,i) = 9.13e-12;
Gstr{i,1} = 'M2_FURAN_P5'; Gstr{i,2} = 'NO';
fM2_FURAN_P5(i)=fM2_FURAN_P5(i)-1; fNO(i)=fNO(i)-1; fNO2(i)=fNO2(i)+1.82; fRCHO(i)=fRCHO(i)+.882; fRDNO3(i)=fRDNO3(i)+.062; fRCNO3(i)=fRCNO3(i)+.056; fOH(i)=fOH(i)+.056; fXC(i)=fXC(i)+.752;

i=i+1;
Rnames{i} = 'M2_FURAN_P5 + NO3 = 1.94*NO2 + .94*RCHO + .06*RCNO3 + .06*OH + XC';
k(:,i) = 2.30e-12;
Gstr{i,1} = 'M2_FURAN_P5'; Gstr{i,2} = 'NO3';
fM2_FURAN_P5(i)=fM2_FURAN_P5(i)-1; fNO3(i)=fNO3(i)-1; fNO2(i)=fNO2(i)+1.94; fRCHO(i)=fRCHO(i)+.94; fRCNO3(i)=fRCNO3(i)+.06; fOH(i)=fOH(i)+.06; fXC(i)=fXC(i)+1;

i=i+1;
Rnames{i} = 'M2_FURAN_P5 + HO2 = RHNO3 - 3*XC';
k(:,i) = 1.78e-11;
Gstr{i,1} = 'M2_FURAN_P5'; Gstr{i,2} = 'HO2';
fM2_FURAN_P5(i)=fM2_FURAN_P5(i)-1; fHO2(i)=fHO2(i)-1; fRHNO3(i)=fRHNO3(i)+1; fXC(i)=fXC(i)-3;

i=i+1;
Rnames{i} = 'M2_FURAN_P5 + SumRO2 = SumRO2 + .47*RCHO + .47*NO2 + .28*RCNO3 + .25*RHNO3 + .03*OH';
k(:,i) = 2.57e-12;
Gstr{i,1} = 'M2_FURAN_P5'; Gstr{i,2} = 'SumRO2';
fM2_FURAN_P5(i)=fM2_FURAN_P5(i)-1; fSumRO2(i)=fSumRO2(i)-1; fSumRO2(i)=fSumRO2(i)+1; fRCHO(i)=fRCHO(i)+.47; fNO2(i)=fNO2(i)+.47; fRCNO3(i)=fRCNO3(i)+.28; fRHNO3(i)=fRHNO3(i)+.25; fOH(i)=fOH(i)+.03;

i=i+1;
Rnames{i} = 'M2_FURAN_P5 + SumRCO3 = SumRCO3 + .752*RCHO + .752*NO2 + .248*RCNO3 + .048*OH + XC';
k(:,i) = 1.37e-11;
Gstr{i,1} = 'M2_FURAN_P5'; Gstr{i,2} = 'SumRCO3';
fM2_FURAN_P5(i)=fM2_FURAN_P5(i)-1; fSumRCO3(i)=fSumRCO3(i)-1; fSumRCO3(i)=fSumRCO3(i)+1; fRCHO(i)=fRCHO(i)+.752; fNO2(i)=fNO2(i)+.752; fRCNO3(i)=fRCNO3(i)+.248; fOH(i)=fOH(i)+.048; fXC(i)=fXC(i)+1;

i=i+1;
Rnames{i} = 'M3_FURAN + OH = .73*M3_FURAN_P1 + .27*M2BUTDAL + .27*HO2 + .73*SumRO2';
k(:,i) = 8.73e-11;
Gstr{i,1} = 'M3_FURAN'; Gstr{i,2} = 'OH';
fM3_FURAN(i)=fM3_FURAN(i)-1; fOH(i)=fOH(i)-1; fM3_FURAN_P1(i)=fM3_FURAN_P1(i)+.73; fM2BUTDAL(i)=fM2BUTDAL(i)+.27; fHO2(i)=fHO2(i)+.27; fSumRO2(i)=fSumRO2(i)+.73;

i=i+1;
Rnames{i} = 'M3_FURAN + O3 = .591*KET2 + .409*RCHO - 0.182*XC';
k(:,i) = 2.05e-17;
Gstr{i,1} = 'M3_FURAN'; Gstr{i,2} = 'O3';
fM3_FURAN(i)=fM3_FURAN(i)-1; fO3(i)=fO3(i)-1; fKET2(i)=fKET2(i)+.591; fRCHO(i)=fRCHO(i)+.409; fXC(i)=fXC(i)-0.182;

i=i+1;
Rnames{i} = 'M3_FURAN + NO3 = M3_FURAN_P2 + SumRO2';
k(:,i) = 1.26e-11;
Gstr{i,1} = 'M3_FURAN'; Gstr{i,2} = 'NO3';
fM3_FURAN(i)=fM3_FURAN(i)-1; fNO3(i)=fNO3(i)-1; fM3_FURAN_P2(i)=fM3_FURAN_P2(i)+1; fSumRO2(i)=fSumRO2(i)+1;

i=i+1;
Rnames{i} = 'M3_FURAN_P1 + NO = .938*NO2 + .937*HO2 + .81*HFON52M4 + .114*HFON52M3 + .062*RHNO3 + .013*RO2C + .013*RCHO + .001*RO2XC + .001*zRHNO3 - 0.176*XC + .014*SumRO2';
k(:,i) = 9.13e-12;
Gstr{i,1} = 'M3_FURAN_P1'; Gstr{i,2} = 'NO';
fM3_FURAN_P1(i)=fM3_FURAN_P1(i)-1; fNO(i)=fNO(i)-1; fNO2(i)=fNO2(i)+.938; fHO2(i)=fHO2(i)+.937; fHFON52M4(i)=fHFON52M4(i)+.81; fHFON52M3(i)=fHFON52M3(i)+.114; fRHNO3(i)=fRHNO3(i)+.062; fRO2C(i)=fRO2C(i)+.013; fRCHO(i)=fRCHO(i)+.013; fRO2XC(i)=fRO2XC(i)+.001; fzRHNO3(i)=fzRHNO3(i)+.001; fXC(i)=fXC(i)-0.176; fSumRO2(i)=fSumRO2(i)+.014;

i=i+1;
Rnames{i} = 'M3_FURAN_P1 + NO3 = NO2 + .999*HO2 + .863*HFON52M4 + .122*HFON52M3 + .014*RO2C + .014*RCHO + .001*RO2XC + .001*zRHNO3 + .011*XC + .015*SumRO2';
k(:,i) = 2.30e-12;
Gstr{i,1} = 'M3_FURAN_P1'; Gstr{i,2} = 'NO3';
fM3_FURAN_P1(i)=fM3_FURAN_P1(i)-1; fNO3(i)=fNO3(i)-1; fNO2(i)=fNO2(i)+1; fHO2(i)=fHO2(i)+.999; fHFON52M4(i)=fHFON52M4(i)+.863; fHFON52M3(i)=fHFON52M3(i)+.122; fRO2C(i)=fRO2C(i)+.014; fRCHO(i)=fRCHO(i)+.014; fRO2XC(i)=fRO2XC(i)+.001; fzRHNO3(i)=fzRHNO3(i)+.001; fXC(i)=fXC(i)+.011; fSumRO2(i)=fSumRO2(i)+.015;

i=i+1;
Rnames{i} = 'M3_FURAN_P1 + HO2 = .9*RUOOH + .086*HFON52M4 + .014*HFON52M3 - 0.9*XC';
k(:,i) = 1.78e-11;
Gstr{i,1} = 'M3_FURAN_P1'; Gstr{i,2} = 'HO2';
fM3_FURAN_P1(i)=fM3_FURAN_P1(i)-1; fHO2(i)=fHO2(i)-1; fRUOOH(i)=fRUOOH(i)+.9; fHFON52M4(i)=fHFON52M4(i)+.086; fHFON52M3(i)=fHFON52M3(i)+.014; fXC(i)=fXC(i)-0.9;

i=i+1;
Rnames{i} = 'M3_FURAN_P1 + SumRO2 = SumRO2 + .647*HFON52M4 + .5*HO2 + .25*OLEP + .095*HFON52M3 + .007*RO2C + .007*RCHO - 0.738*XC + .007*SumRO2';
k(:,i) = 2.57e-12;
Gstr{i,1} = 'M3_FURAN_P1'; Gstr{i,2} = 'SumRO2';
fM3_FURAN_P1(i)=fM3_FURAN_P1(i)-1; fSumRO2(i)=fSumRO2(i)-1; fSumRO2(i)=fSumRO2(i)+1; fHFON52M4(i)=fHFON52M4(i)+.647; fHO2(i)=fHO2(i)+.5; fOLEP(i)=fOLEP(i)+.25; fHFON52M3(i)=fHFON52M3(i)+.095; fRO2C(i)=fRO2C(i)+.007; fRCHO(i)=fRCHO(i)+.007; fXC(i)=fXC(i)-0.738; fSumRO2(i)=fSumRO2(i)+.007;

i=i+1;
Rnames{i} = 'M3_FURAN_P1 + SumRCO3 = SumRCO3 + .863*HFON52M4 + .799*HO2 + .125*HFON52M3 + .011*RO2C + .011*RCHO + .001*RO2XC + .001*zRHNO3 + .008*XC + .012*SumRO2';
k(:,i) = 1.37e-11;
Gstr{i,1} = 'M3_FURAN_P1'; Gstr{i,2} = 'SumRCO3';
fM3_FURAN_P1(i)=fM3_FURAN_P1(i)-1; fSumRCO3(i)=fSumRCO3(i)-1; fSumRCO3(i)=fSumRCO3(i)+1; fHFON52M4(i)=fHFON52M4(i)+.863; fHO2(i)=fHO2(i)+.799; fHFON52M3(i)=fHFON52M3(i)+.125; fRO2C(i)=fRO2C(i)+.011; fRCHO(i)=fRCHO(i)+.011; fRO2XC(i)=fRO2XC(i)+.001; fzRHNO3(i)=fzRHNO3(i)+.001; fXC(i)=fXC(i)+.008; fSumRO2(i)=fSumRO2(i)+.012;

i=i+1;
Rnames{i} = 'M3_FURAN_P2 + NO = .951*NO2 + .924*RCNO3 + .924*HO2 + .278*CO + .062*RDNO3 + .013*RO2C + .013*RCHO + .001*RO2XC + .001*zRDNO3 + .47*XC + .014*SumRO2';
k(:,i) = 9.13e-12;
Gstr{i,1} = 'M3_FURAN_P2'; Gstr{i,2} = 'NO';
fM3_FURAN_P2(i)=fM3_FURAN_P2(i)-1; fNO(i)=fNO(i)-1; fNO2(i)=fNO2(i)+.951; fRCNO3(i)=fRCNO3(i)+.924; fHO2(i)=fHO2(i)+.924; fCO(i)=fCO(i)+.278; fRDNO3(i)=fRDNO3(i)+.062; fRO2C(i)=fRO2C(i)+.013; fRCHO(i)=fRCHO(i)+.013; fRO2XC(i)=fRO2XC(i)+.001; fzRDNO3(i)=fzRDNO3(i)+.001; fXC(i)=fXC(i)+.47; fSumRO2(i)=fSumRO2(i)+.014;

i=i+1;
Rnames{i} = 'M3_FURAN_P2 + NO3 = 1.014*NO2 + .985*RCNO3 + .985*HO2 + .296*CO + .014*RO2C + .014*RCHO + .001*RO2XC + .001*zRDNO3 + .7*XC + .015*SumRO2';
k(:,i) = 2.30e-12;
Gstr{i,1} = 'M3_FURAN_P2'; Gstr{i,2} = 'NO3';
fM3_FURAN_P2(i)=fM3_FURAN_P2(i)-1; fNO3(i)=fNO3(i)-1; fNO2(i)=fNO2(i)+1.014; fRCNO3(i)=fRCNO3(i)+.985; fHO2(i)=fHO2(i)+.985; fCO(i)=fCO(i)+.296; fRO2C(i)=fRO2C(i)+.014; fRCHO(i)=fRCHO(i)+.014; fRO2XC(i)=fRO2XC(i)+.001; fzRDNO3(i)=fzRDNO3(i)+.001; fXC(i)=fXC(i)+.7; fSumRO2(i)=fSumRO2(i)+.015;

i=i+1;
Rnames{i} = 'M3_FURAN_P2 + HO2 = .9*RHNO3 + .1*RCNO3 - 2.6*XC';
k(:,i) = 1.78e-11;
Gstr{i,1} = 'M3_FURAN_P2'; Gstr{i,2} = 'HO2';
fM3_FURAN_P2(i)=fM3_FURAN_P2(i)-1; fHO2(i)=fHO2(i)-1; fRHNO3(i)=fRHNO3(i)+.9; fRCNO3(i)=fRCNO3(i)+.1; fXC(i)=fXC(i)-2.6;

i=i+1;
Rnames{i} = 'M3_FURAN_P2 + SumRO2 = SumRO2 + .742*RCNO3 + .492*HO2 + .25*RHNO3 + .148*CO + .008*NO2 + .007*RO2C + .007*RCHO - 0.144*XC + .007*SumRO2';
k(:,i) = 2.57e-12;
Gstr{i,1} = 'M3_FURAN_P2'; Gstr{i,2} = 'SumRO2';
fM3_FURAN_P2(i)=fM3_FURAN_P2(i)-1; fSumRO2(i)=fSumRO2(i)-1; fSumRO2(i)=fSumRO2(i)+1; fRCNO3(i)=fRCNO3(i)+.742; fHO2(i)=fHO2(i)+.492; fRHNO3(i)=fRHNO3(i)+.25; fCO(i)=fCO(i)+.148; fNO2(i)=fNO2(i)+.008; fRO2C(i)=fRO2C(i)+.007; fRCHO(i)=fRCHO(i)+.007; fXC(i)=fXC(i)-0.144; fSumRO2(i)=fSumRO2(i)+.007;

i=i+1;
Rnames{i} = 'M3_FURAN_P2 + SumRCO3 = SumRCO3 + .988*RCNO3 + .788*HO2 + .237*CO + .011*NO2 + .011*RO2C + .011*RCHO + .001*RO2XC + .001*zRDNO3 + .759*XC + .012*SumRO2';
k(:,i) = 1.37e-11;
Gstr{i,1} = 'M3_FURAN_P2'; Gstr{i,2} = 'SumRCO3';
fM3_FURAN_P2(i)=fM3_FURAN_P2(i)-1; fSumRCO3(i)=fSumRCO3(i)-1; fSumRCO3(i)=fSumRCO3(i)+1; fRCNO3(i)=fRCNO3(i)+.988; fHO2(i)=fHO2(i)+.788; fCO(i)=fCO(i)+.237; fNO2(i)=fNO2(i)+.011; fRO2C(i)=fRO2C(i)+.011; fRCHO(i)=fRCHO(i)+.011; fRO2XC(i)=fRO2XC(i)+.001; fzRDNO3(i)=fzRDNO3(i)+.001; fXC(i)=fXC(i)+.759; fSumRO2(i)=fSumRO2(i)+.012;

i=i+1;
Rnames{i} = 'M25_FUR + OH = .72*M25_FUR_P1 + .28*H3XE25DO + .28*HO2 + .72*SumRO2';
k(:,i) = 1.26e-10;
Gstr{i,1} = 'M25_FUR'; Gstr{i,2} = 'OH';
fM25_FUR(i)=fM25_FUR(i)-1; fOH(i)=fOH(i)-1; fM25_FUR_P1(i)=fM25_FUR_P1(i)+.72; fH3XE25DO(i)=fH3XE25DO(i)+.28; fHO2(i)=fHO2(i)+.28; fSumRO2(i)=fSumRO2(i)+.72;

i=i+1;
Rnames{i} = 'M25_FUR + O3 = 1.5*OH + .5*OLEP + .5*M25_FUR_R1 - 1*XC';
k(:,i) = 4.20e-16;
Gstr{i,1} = 'M25_FUR'; Gstr{i,2} = 'O3';
fM25_FUR(i)=fM25_FUR(i)-1; fO3(i)=fO3(i)-1; fOH(i)=fOH(i)+1.5; fOLEP(i)=fOLEP(i)+.5; fM25_FUR_R1(i)=fM25_FUR_R1(i)+.5; fXC(i)=fXC(i)-1;

i=i+1;
Rnames{i} = 'M25_FUR + NO3 = M25_FUR_P2 + SumRO2';
k(:,i) = 5.80e-11;
Gstr{i,1} = 'M25_FUR'; Gstr{i,2} = 'NO3';
fM25_FUR(i)=fM25_FUR(i)-1; fNO3(i)=fNO3(i)-1; fM25_FUR_P2(i)=fM25_FUR_P2(i)+1; fSumRO2(i)=fSumRO2(i)+1;

i=i+1;
Rnames{i} = 'M25_FUR_P1 + NO = .876*NO2 + .692*HFON52M5 + .692*MEO2 + .184*M25_FUR_P3 + .124*RHNO3 - 0.248*XC + .876*SumRO2';
k(:,i) = 9.13e-12;
Gstr{i,1} = 'M25_FUR_P1'; Gstr{i,2} = 'NO';
fM25_FUR_P1(i)=fM25_FUR_P1(i)-1; fNO(i)=fNO(i)-1; fNO2(i)=fNO2(i)+.876; fHFON52M5(i)=fHFON52M5(i)+.692; fMEO2(i)=fMEO2(i)+.692; fM25_FUR_P3(i)=fM25_FUR_P3(i)+.184; fRHNO3(i)=fRHNO3(i)+.124; fXC(i)=fXC(i)-0.248; fSumRO2(i)=fSumRO2(i)+.876;

i=i+1;
Rnames{i} = 'M25_FUR_P1 + NO3 = NO2 + .79*HFON52M5 + .79*MEO2 + .21*M25_FUR_P3 + SumRO2';
k(:,i) = 2.30e-12;
Gstr{i,1} = 'M25_FUR_P1'; Gstr{i,2} = 'NO3';
fM25_FUR_P1(i)=fM25_FUR_P1(i)-1; fNO3(i)=fNO3(i)-1; fNO2(i)=fNO2(i)+1; fHFON52M5(i)=fHFON52M5(i)+.79; fMEO2(i)=fMEO2(i)+.79; fM25_FUR_P3(i)=fM25_FUR_P3(i)+.21; fSumRO2(i)=fSumRO2(i)+1;

i=i+1;
Rnames{i} = 'M25_FUR_P1 + HO2 = RUOOH';
k(:,i) = 1.95e-11;
Gstr{i,1} = 'M25_FUR_P1'; Gstr{i,2} = 'HO2';
fM25_FUR_P1(i)=fM25_FUR_P1(i)-1; fHO2(i)=fHO2(i)-1; fRUOOH(i)=fRUOOH(i)+1;

i=i+1;
Rnames{i} = 'M25_FUR_P1 + SumRO2 = SumRO2 + .5*OLEP + .395*HFON52M5 + .395*MEO2 + .105*M25_FUR_P3 - 1*XC + .5*SumRO2';
k(:,i) = 3.71e-16;
Gstr{i,1} = 'M25_FUR_P1'; Gstr{i,2} = 'SumRO2';
fM25_FUR_P1(i)=fM25_FUR_P1(i)-1; fSumRO2(i)=fSumRO2(i)-1; fSumRO2(i)=fSumRO2(i)+1; fOLEP(i)=fOLEP(i)+.5; fHFON52M5(i)=fHFON52M5(i)+.395; fMEO2(i)=fMEO2(i)+.395; fM25_FUR_P3(i)=fM25_FUR_P3(i)+.105; fXC(i)=fXC(i)-1; fSumRO2(i)=fSumRO2(i)+.5;

i=i+1;
Rnames{i} = 'M25_FUR_P1 + SumRCO3 = SumRCO3 + .79*HFON52M5 + .79*MEO2 + .21*M25_FUR_P3 + SumRO2';
k(:,i) = 1.37e-11;
Gstr{i,1} = 'M25_FUR_P1'; Gstr{i,2} = 'SumRCO3';
fM25_FUR_P1(i)=fM25_FUR_P1(i)-1; fSumRCO3(i)=fSumRCO3(i)-1; fSumRCO3(i)=fSumRCO3(i)+1; fHFON52M5(i)=fHFON52M5(i)+.79; fMEO2(i)=fMEO2(i)+.79; fM25_FUR_P3(i)=fM25_FUR_P3(i)+.21; fSumRO2(i)=fSumRO2(i)+1;

i=i+1;
Rnames{i} = 'M25_FUR_R1 = M25_FUR_P4 + SumRO2';
k(:,i) = 1.64e+11.*(T./300).^0.00.*exp(-7848.229./T);
Gstr{i,1} = 'M25_FUR_R1';
fM25_FUR_R1(i)=fM25_FUR_R1(i)-1; fM25_FUR_P4(i)=fM25_FUR_P4(i)+1; fSumRO2(i)=fSumRO2(i)+1;

i=i+1;
Rnames{i} = 'M25_FUR_R1 + NO = .876*OLEA1 + .876*CO + .876*HO2 + .876*NO2 + .124*RCNO3 - 0.628*XC';
k(:,i) = 9.13e-12;
Gstr{i,1} = 'M25_FUR_R1'; Gstr{i,2} = 'NO';
fM25_FUR_R1(i)=fM25_FUR_R1(i)-1; fNO(i)=fNO(i)-1; fOLEA1(i)=fOLEA1(i)+.876; fCO(i)=fCO(i)+.876; fHO2(i)=fHO2(i)+.876; fNO2(i)=fNO2(i)+.876; fRCNO3(i)=fRCNO3(i)+.124; fXC(i)=fXC(i)-0.628;

i=i+1;
Rnames{i} = 'M25_FUR_P2 + NO = .876*NO2 + .692*RCNO3 + .692*MEO2 + .184*M25_FUR_P5 + .124*RDNO3 + .444*XC + .876*SumRO2';
k(:,i) = 9.13e-12;
Gstr{i,1} = 'M25_FUR_P2'; Gstr{i,2} = 'NO';
fM25_FUR_P2(i)=fM25_FUR_P2(i)-1; fNO(i)=fNO(i)-1; fNO2(i)=fNO2(i)+.876; fRCNO3(i)=fRCNO3(i)+.692; fMEO2(i)=fMEO2(i)+.692; fM25_FUR_P5(i)=fM25_FUR_P5(i)+.184; fRDNO3(i)=fRDNO3(i)+.124; fXC(i)=fXC(i)+.444; fSumRO2(i)=fSumRO2(i)+.876;

i=i+1;
Rnames{i} = 'M25_FUR_P2 + NO3 = NO2 + .79*RCNO3 + .79*MEO2 + .21*M25_FUR_P5 + .79*XC + SumRO2';
k(:,i) = 2.30e-12;
Gstr{i,1} = 'M25_FUR_P2'; Gstr{i,2} = 'NO3';
fM25_FUR_P2(i)=fM25_FUR_P2(i)-1; fNO3(i)=fNO3(i)-1; fNO2(i)=fNO2(i)+1; fRCNO3(i)=fRCNO3(i)+.79; fMEO2(i)=fMEO2(i)+.79; fM25_FUR_P5(i)=fM25_FUR_P5(i)+.21; fXC(i)=fXC(i)+.79; fSumRO2(i)=fSumRO2(i)+1;

i=i+1;
Rnames{i} = 'M25_FUR_P2 + HO2 = RHNO3 - 2*XC';
k(:,i) = 1.95e-11;
Gstr{i,1} = 'M25_FUR_P2'; Gstr{i,2} = 'HO2';
fM25_FUR_P2(i)=fM25_FUR_P2(i)-1; fHO2(i)=fHO2(i)-1; fRHNO3(i)=fRHNO3(i)+1; fXC(i)=fXC(i)-2;

i=i+1;
Rnames{i} = 'M25_FUR_P2 + SumRO2 = SumRO2 + .5*RHNO3 + .395*RCNO3 + .395*MEO2 + .105*M25_FUR_P5 - 0.605*XC + .5*SumRO2';
k(:,i) = 3.71e-16;
Gstr{i,1} = 'M25_FUR_P2'; Gstr{i,2} = 'SumRO2';
fM25_FUR_P2(i)=fM25_FUR_P2(i)-1; fSumRO2(i)=fSumRO2(i)-1; fSumRO2(i)=fSumRO2(i)+1; fRHNO3(i)=fRHNO3(i)+.5; fRCNO3(i)=fRCNO3(i)+.395; fMEO2(i)=fMEO2(i)+.395; fM25_FUR_P5(i)=fM25_FUR_P5(i)+.105; fXC(i)=fXC(i)-0.605; fSumRO2(i)=fSumRO2(i)+.5;

i=i+1;
Rnames{i} = 'M25_FUR_P2 + SumRCO3 = SumRCO3 + .79*RCNO3 + .79*MEO2 + .21*M25_FUR_P5 + .79*XC + SumRO2';
k(:,i) = 1.37e-11;
Gstr{i,1} = 'M25_FUR_P2'; Gstr{i,2} = 'SumRCO3';
fM25_FUR_P2(i)=fM25_FUR_P2(i)-1; fSumRCO3(i)=fSumRCO3(i)-1; fSumRCO3(i)=fSumRCO3(i)+1; fRCNO3(i)=fRCNO3(i)+.79; fMEO2(i)=fMEO2(i)+.79; fM25_FUR_P5(i)=fM25_FUR_P5(i)+.21; fXC(i)=fXC(i)+.79; fSumRO2(i)=fSumRO2(i)+1;

i=i+1;
Rnames{i} = 'M25_FUR_P3 + NO = .876*RCHO + .876*HO2 + .876*NO2 + .124*RHNO3 + 1.504*XC';
k(:,i) = 9.13e-12;
Gstr{i,1} = 'M25_FUR_P3'; Gstr{i,2} = 'NO';
fM25_FUR_P3(i)=fM25_FUR_P3(i)-1; fNO(i)=fNO(i)-1; fRCHO(i)=fRCHO(i)+.876; fHO2(i)=fHO2(i)+.876; fNO2(i)=fNO2(i)+.876; fRHNO3(i)=fRHNO3(i)+.124; fXC(i)=fXC(i)+1.504;

i=i+1;
Rnames{i} = 'M25_FUR_P3 + NO3 = RCHO + HO2 + NO2 + 2*XC';
k(:,i) = 2.30e-12;
Gstr{i,1} = 'M25_FUR_P3'; Gstr{i,2} = 'NO3';
fM25_FUR_P3(i)=fM25_FUR_P3(i)-1; fNO3(i)=fNO3(i)-1; fRCHO(i)=fRCHO(i)+1; fHO2(i)=fHO2(i)+1; fNO2(i)=fNO2(i)+1; fXC(i)=fXC(i)+2;

i=i+1;
Rnames{i} = 'M25_FUR_P3 + HO2 = ROOH + XC';
k(:,i) = 1.95e-11;
Gstr{i,1} = 'M25_FUR_P3'; Gstr{i,2} = 'HO2';
fM25_FUR_P3(i)=fM25_FUR_P3(i)-1; fHO2(i)=fHO2(i)-1; fROOH(i)=fROOH(i)+1; fXC(i)=fXC(i)+1;

i=i+1;
Rnames{i} = 'M25_FUR_P3 + SumRO2 = SumRO2 + .5*RCHO + .5*HO2 + .25*KET2 + .25*IEPOX + 1.25*XC';
k(:,i) = 2.57e-12;
Gstr{i,1} = 'M25_FUR_P3'; Gstr{i,2} = 'SumRO2';
fM25_FUR_P3(i)=fM25_FUR_P3(i)-1; fSumRO2(i)=fSumRO2(i)-1; fSumRO2(i)=fSumRO2(i)+1; fRCHO(i)=fRCHO(i)+.5; fHO2(i)=fHO2(i)+.5; fKET2(i)=fKET2(i)+.25; fIEPOX(i)=fIEPOX(i)+.25; fXC(i)=fXC(i)+1.25;

i=i+1;
Rnames{i} = 'M25_FUR_P3 + SumRCO3 = SumRCO3 + .8*RCHO + .8*HO2 + .2*KET2 + 1.6*XC';
k(:,i) = 1.37e-11;
Gstr{i,1} = 'M25_FUR_P3'; Gstr{i,2} = 'SumRCO3';
fM25_FUR_P3(i)=fM25_FUR_P3(i)-1; fSumRCO3(i)=fSumRCO3(i)-1; fSumRCO3(i)=fSumRCO3(i)+1; fRCHO(i)=fRCHO(i)+.8; fHO2(i)=fHO2(i)+.8; fKET2(i)=fKET2(i)+.2; fXC(i)=fXC(i)+1.6;

i=i+1;
Rnames{i} = 'M25_FUR_P4 + NO = .876*OLEA1 + .876*CO2 + .876*OH + .876*NO2 + .124*RCNO3 - 0.628*XC';
k(:,i) = 9.13e-12;
Gstr{i,1} = 'M25_FUR_P4'; Gstr{i,2} = 'NO';
fM25_FUR_P4(i)=fM25_FUR_P4(i)-1; fNO(i)=fNO(i)-1; fOLEA1(i)=fOLEA1(i)+.876; fCO2(i)=fCO2(i)+.876; fOH(i)=fOH(i)+.876; fNO2(i)=fNO2(i)+.876; fRCNO3(i)=fRCNO3(i)+.124; fXC(i)=fXC(i)-0.628;

i=i+1;
Rnames{i} = 'M25_FUR_P4 + NO3 = OLEA1 + CO2 + OH + NO2 - 1*XC';
k(:,i) = 2.30e-12;
Gstr{i,1} = 'M25_FUR_P4'; Gstr{i,2} = 'NO3';
fM25_FUR_P4(i)=fM25_FUR_P4(i)-1; fNO3(i)=fNO3(i)-1; fOLEA1(i)=fOLEA1(i)+1; fCO2(i)=fCO2(i)+1; fOH(i)=fOH(i)+1; fNO2(i)=fNO2(i)+1; fXC(i)=fXC(i)-1;

i=i+1;
Rnames{i} = 'M25_FUR_P4 + HO2 = .9*RUOOH + .1*BACL + .2*XC';
k(:,i) = 1.95e-11;
Gstr{i,1} = 'M25_FUR_P4'; Gstr{i,2} = 'HO2';
fM25_FUR_P4(i)=fM25_FUR_P4(i)-1; fHO2(i)=fHO2(i)-1; fRUOOH(i)=fRUOOH(i)+.9; fBACL(i)=fBACL(i)+.1; fXC(i)=fXC(i)+.2;

i=i+1;
Rnames{i} = 'M25_FUR_P4 + SumRO2 = SumRO2 + .5*OLEA1 + .5*CO2 + .5*OH + .25*BACL + .25*RUOOH';
k(:,i) = 2.57e-12;
Gstr{i,1} = 'M25_FUR_P4'; Gstr{i,2} = 'SumRO2';
fM25_FUR_P4(i)=fM25_FUR_P4(i)-1; fSumRO2(i)=fSumRO2(i)-1; fSumRO2(i)=fSumRO2(i)+1; fOLEA1(i)=fOLEA1(i)+.5; fCO2(i)=fCO2(i)+.5; fOH(i)=fOH(i)+.5; fBACL(i)=fBACL(i)+.25; fRUOOH(i)=fRUOOH(i)+.25;

i=i+1;
Rnames{i} = 'M25_FUR_P4 + SumRCO3 = SumRCO3 + .8*OLEA1 + .8*CO2 + .8*OH + .2*BACL - 0.4*XC';
k(:,i) = 1.37e-11;
Gstr{i,1} = 'M25_FUR_P4'; Gstr{i,2} = 'SumRCO3';
fM25_FUR_P4(i)=fM25_FUR_P4(i)-1; fSumRCO3(i)=fSumRCO3(i)-1; fSumRCO3(i)=fSumRCO3(i)+1; fOLEA1(i)=fOLEA1(i)+.8; fCO2(i)=fCO2(i)+.8; fOH(i)=fOH(i)+.8; fBACL(i)=fBACL(i)+.2; fXC(i)=fXC(i)-0.4;

i=i+1;
Rnames{i} = 'M25_FUR_P5 + NO = 1.536*NO2 + .66*RCHO + .216*RCNO3 + .216*OH + .124*RDNO3 + 1.504*XC';
k(:,i) = 9.13e-12;
Gstr{i,1} = 'M25_FUR_P5'; Gstr{i,2} = 'NO';
fM25_FUR_P5(i)=fM25_FUR_P5(i)-1; fNO(i)=fNO(i)-1; fNO2(i)=fNO2(i)+1.536; fRCHO(i)=fRCHO(i)+.66; fRCNO3(i)=fRCNO3(i)+.216; fOH(i)=fOH(i)+.216; fRDNO3(i)=fRDNO3(i)+.124; fXC(i)=fXC(i)+1.504;

i=i+1;
Rnames{i} = 'M25_FUR_P5 + NO3 = 1.753*NO2 + .753*RCHO + .247*RCNO3 + .247*OH + 2*XC';
k(:,i) = 2.30e-12;
Gstr{i,1} = 'M25_FUR_P5'; Gstr{i,2} = 'NO3';
fM25_FUR_P5(i)=fM25_FUR_P5(i)-1; fNO3(i)=fNO3(i)-1; fNO2(i)=fNO2(i)+1.753; fRCHO(i)=fRCHO(i)+.753; fRCNO3(i)=fRCNO3(i)+.247; fOH(i)=fOH(i)+.247; fXC(i)=fXC(i)+2;

i=i+1;
Rnames{i} = 'M25_FUR_P5 + HO2 = RHNO3 - 2*XC';
k(:,i) = 1.95e-11;
Gstr{i,1} = 'M25_FUR_P5'; Gstr{i,2} = 'HO2';
fM25_FUR_P5(i)=fM25_FUR_P5(i)-1; fHO2(i)=fHO2(i)-1; fRHNO3(i)=fRHNO3(i)+1; fXC(i)=fXC(i)-2;

i=i+1;
Rnames{i} = 'M25_FUR_P5 + SumRO2 = SumRO2 + .376*RCHO + .376*NO2 + .374*RCNO3 + .25*RHNO3 + .124*OH + XC';
k(:,i) = 2.57e-12;
Gstr{i,1} = 'M25_FUR_P5'; Gstr{i,2} = 'SumRO2';
fM25_FUR_P5(i)=fM25_FUR_P5(i)-1; fSumRO2(i)=fSumRO2(i)-1; fSumRO2(i)=fSumRO2(i)+1; fRCHO(i)=fRCHO(i)+.376; fNO2(i)=fNO2(i)+.376; fRCNO3(i)=fRCNO3(i)+.374; fRHNO3(i)=fRHNO3(i)+.25; fOH(i)=fOH(i)+.124; fXC(i)=fXC(i)+1;

i=i+1;
Rnames{i} = 'M25_FUR_P5 + SumRCO3 = SumRCO3 + .602*RCHO + .602*NO2 + .398*RCNO3 + .198*OH + 2*XC';
k(:,i) = 1.37e-11;
Gstr{i,1} = 'M25_FUR_P5'; Gstr{i,2} = 'SumRCO3';
fM25_FUR_P5(i)=fM25_FUR_P5(i)-1; fSumRCO3(i)=fSumRCO3(i)-1; fSumRCO3(i)=fSumRCO3(i)+1; fRCHO(i)=fRCHO(i)+.602; fNO2(i)=fNO2(i)+.602; fRCNO3(i)=fRCNO3(i)+.398; fOH(i)=fOH(i)+.198; fXC(i)=fXC(i)+2;

i=i+1;
Rnames{i} = 'BUTEDIAL + OH = .817*OH + .741*MALAH + .183*BUTEDIAL_R1 + .076*ROOH - 0.076*XC';
k(:,i) = 4.22e-11;
Gstr{i,1} = 'BUTEDIAL'; Gstr{i,2} = 'OH';
fBUTEDIAL(i)=fBUTEDIAL(i)-1; fOH(i)=fOH(i)-1; fOH(i)=fOH(i)+.817; fMALAH(i)=fMALAH(i)+.741; fBUTEDIAL_R1(i)=fBUTEDIAL_R1(i)+.183; fROOH(i)=fROOH(i)+.076; fXC(i)=fXC(i)-0.076;

i=i+1;
Rnames{i} = 'BUTEDIAL + O3 = GLY + .668*HO2 + .607*CO2 + .393*RCHO2 + .334*CO + .273*HCHO - 0.393*XC';
k(:,i) = 1.60e-18;
Gstr{i,1} = 'BUTEDIAL'; Gstr{i,2} = 'O3';
fBUTEDIAL(i)=fBUTEDIAL(i)-1; fO3(i)=fO3(i)-1; fGLY(i)=fGLY(i)+1; fHO2(i)=fHO2(i)+.668; fCO2(i)=fCO2(i)+.607; fRCHO2(i)=fRCHO2(i)+.393; fCO(i)=fCO(i)+.334; fHCHO(i)=fHCHO(i)+.273; fXC(i)=fXC(i)-0.393;

i=i+1;
Rnames{i} = 'BUTEDIAL + NO3 = .994*OH + .95*MALAH + .05*RCNO3 + .005*RO2C + .005*CO + .005*HO2 - 0.005*XC + .95*XN + .005*SumRO2';
k(:,i) = 2.09e-15;
Gstr{i,1} = 'BUTEDIAL'; Gstr{i,2} = 'NO3';
fBUTEDIAL(i)=fBUTEDIAL(i)-1; fNO3(i)=fNO3(i)-1; fOH(i)=fOH(i)+.994; fMALAH(i)=fMALAH(i)+.95; fRCNO3(i)=fRCNO3(i)+.05; fRO2C(i)=fRO2C(i)+.005; fCO(i)=fCO(i)+.005; fHO2(i)=fHO2(i)+.005; fXC(i)=fXC(i)-0.005; fXN(i)=fXN(i)+.95; fSumRO2(i)=fSumRO2(i)+.005;

i=i+1;
Rnames{i} = 'BUTEDIAL + HV = .59*MALAH + .59*OH + .59*HO2 + .41*H3F2ONE';
k(:,i) = JAFGS.*2.88e-1;
Gstr{i,1} = 'BUTEDIAL';
fBUTEDIAL(i)=fBUTEDIAL(i)-1; fMALAH(i)=fMALAH(i)+.59; fOH(i)=fOH(i)+.59; fHO2(i)=fHO2(i)+.59; fH3F2ONE(i)=fH3F2ONE(i)+.41;

i=i+1;
Rnames{i} = 'BUTEDIAL_R1 = ROOH + OH - 1*XC';
k(:,i) = 1.64e+11.*(T./300).^0.00.*exp(-7848.229./T);
Gstr{i,1} = 'BUTEDIAL_R1';
fBUTEDIAL_R1(i)=fBUTEDIAL_R1(i)-1; fROOH(i)=fROOH(i)+1; fOH(i)=fOH(i)+1; fXC(i)=fXC(i)-1;

i=i+1;
Rnames{i} = 'BUTEDIAL_R1 + NO = .974*HO2 + .974*NO2 + .884*GLY + .884*MGLY + .09*CROOH + .09*CO + .026*RCNO3 - 1.334*XC';
k(:,i) = 9.13e-12;
Gstr{i,1} = 'BUTEDIAL_R1'; Gstr{i,2} = 'NO';
fBUTEDIAL_R1(i)=fBUTEDIAL_R1(i)-1; fNO(i)=fNO(i)-1; fHO2(i)=fHO2(i)+.974; fNO2(i)=fNO2(i)+.974; fGLY(i)=fGLY(i)+.884; fMGLY(i)=fMGLY(i)+.884; fCROOH(i)=fCROOH(i)+.09; fCO(i)=fCO(i)+.09; fRCNO3(i)=fRCNO3(i)+.026; fXC(i)=fXC(i)-1.334;

i=i+1;
Rnames{i} = 'M2BUTDAL + OH = .767*OH + .67*MALAH + .167*M2BUTDAL_R2 + .097*ROOH + .066*M2BUTDAL_R1 + .67*XC';
k(:,i) = 4.15e-11;
Gstr{i,1} = 'M2BUTDAL'; Gstr{i,2} = 'OH';
fM2BUTDAL(i)=fM2BUTDAL(i)-1; fOH(i)=fOH(i)-1; fOH(i)=fOH(i)+.767; fMALAH(i)=fMALAH(i)+.67; fM2BUTDAL_R2(i)=fM2BUTDAL_R2(i)+.167; fROOH(i)=fROOH(i)+.097; fM2BUTDAL_R1(i)=fM2BUTDAL_R1(i)+.066; fXC(i)=fXC(i)+.67;

i=i+1;
Rnames{i} = 'M2BUTDAL + O3 = .5*GLY + .5*MGLY + .447*RCHO2 + .417*CO + .334*HO2 + .303*CO2 + .25*M2BUTDAL_P1 + .25*OH + .137*HCHO - 0.198*XC + .25*SumRO2';
k(:,i) = 4.80e-18;
Gstr{i,1} = 'M2BUTDAL'; Gstr{i,2} = 'O3';
fM2BUTDAL(i)=fM2BUTDAL(i)-1; fO3(i)=fO3(i)-1; fGLY(i)=fGLY(i)+.5; fMGLY(i)=fMGLY(i)+.5; fRCHO2(i)=fRCHO2(i)+.447; fCO(i)=fCO(i)+.417; fHO2(i)=fHO2(i)+.334; fCO2(i)=fCO2(i)+.303; fM2BUTDAL_P1(i)=fM2BUTDAL_P1(i)+.25; fOH(i)=fOH(i)+.25; fHCHO(i)=fHCHO(i)+.137; fXC(i)=fXC(i)-0.198; fSumRO2(i)=fSumRO2(i)+.25;

i=i+1;
Rnames{i} = 'M2BUTDAL + NO3 = .523*OH + .477*M2BUTDAL_R3 + .325*MALAH + .198*RCNO3 + .523*XC + .325*XN';
k(:,i) = 3.13e-15;
Gstr{i,1} = 'M2BUTDAL'; Gstr{i,2} = 'NO3';
fM2BUTDAL(i)=fM2BUTDAL(i)-1; fNO3(i)=fNO3(i)-1; fOH(i)=fOH(i)+.523; fM2BUTDAL_R3(i)=fM2BUTDAL_R3(i)+.477; fMALAH(i)=fMALAH(i)+.325; fRCNO3(i)=fRCNO3(i)+.198; fXC(i)=fXC(i)+.523; fXN(i)=fXN(i)+.325;

i=i+1;
Rnames{i} = 'M2BUTDAL + HV = .7*M2BUTDAL_P2 + .7*HO2 + .3*H3F2ONM4 + .7*SumRO2';
k(:,i) = JAFGS.*3.42e-1;
Gstr{i,1} = 'M2BUTDAL';
fM2BUTDAL(i)=fM2BUTDAL(i)-1; fM2BUTDAL_P2(i)=fM2BUTDAL_P2(i)+.7; fHO2(i)=fHO2(i)+.7; fH3F2ONM4(i)=fH3F2ONM4(i)+.3; fSumRO2(i)=fSumRO2(i)+.7;

i=i+1;
Rnames{i} = 'M2BUTDAL_R1 = ROOH + OH';
k(:,i) = 1.64e+11.*(T./300).^0.00.*exp(-7848.229./T);
Gstr{i,1} = 'M2BUTDAL_R1';
fM2BUTDAL_R1(i)=fM2BUTDAL_R1(i)-1; fROOH(i)=fROOH(i)+1; fOH(i)=fOH(i)+1;

i=i+1;
Rnames{i} = 'M2BUTDAL_R1 + NO = .938*HO2 + .938*NO2 + .585*BACL + .585*GLY + .353*CROOH + .353*CO + .062*RCNO3 - 1.935*XC';
k(:,i) = 9.13e-12;
Gstr{i,1} = 'M2BUTDAL_R1'; Gstr{i,2} = 'NO';
fM2BUTDAL_R1(i)=fM2BUTDAL_R1(i)-1; fNO(i)=fNO(i)-1; fHO2(i)=fHO2(i)+.938; fNO2(i)=fNO2(i)+.938; fBACL(i)=fBACL(i)+.585; fGLY(i)=fGLY(i)+.585; fCROOH(i)=fCROOH(i)+.353; fCO(i)=fCO(i)+.353; fRCNO3(i)=fRCNO3(i)+.062; fXC(i)=fXC(i)-1.935;

i=i+1;
Rnames{i} = 'M2BUTDAL_R2 = ROOH + OH';
k(:,i) = 1.64e+11.*(T./300).^0.00.*exp(-7848.229./T);
Gstr{i,1} = 'M2BUTDAL_R2';
fM2BUTDAL_R2(i)=fM2BUTDAL_R2(i)-1; fROOH(i)=fROOH(i)+1; fOH(i)=fOH(i)+1;

i=i+1;
Rnames{i} = 'M2BUTDAL_R2 + NO = 1.825*MGLY + .938*HO2 + .938*NO2 + .062*RCNO3 + .025*CROOH + .025*CO - 0.948*XC';
k(:,i) = 9.13e-12;
Gstr{i,1} = 'M2BUTDAL_R2'; Gstr{i,2} = 'NO';
fM2BUTDAL_R2(i)=fM2BUTDAL_R2(i)-1; fNO(i)=fNO(i)-1; fMGLY(i)=fMGLY(i)+1.825; fHO2(i)=fHO2(i)+.938; fNO2(i)=fNO2(i)+.938; fRCNO3(i)=fRCNO3(i)+.062; fCROOH(i)=fCROOH(i)+.025; fCO(i)=fCO(i)+.025; fXC(i)=fXC(i)-0.948;

i=i+1;
Rnames{i} = 'M2BUTDAL_P1 + NO = NO2 + .749*MGLY + .749*HO2 + .251*CO2 + .251*HCHO + .251*OH - 0.749*XC';
k(:,i) = 2.55e-12.*(T./300).^0.00.*exp(379.932./T);
Gstr{i,1} = 'M2BUTDAL_P1'; Gstr{i,2} = 'NO';
fM2BUTDAL_P1(i)=fM2BUTDAL_P1(i)-1; fNO(i)=fNO(i)-1; fNO2(i)=fNO2(i)+1; fMGLY(i)=fMGLY(i)+.749; fHO2(i)=fHO2(i)+.749; fCO2(i)=fCO2(i)+.251; fHCHO(i)=fHCHO(i)+.251; fOH(i)=fOH(i)+.251; fXC(i)=fXC(i)-0.749;

i=i+1;
Rnames{i} = 'M2BUTDAL_P1 + NO3 = NO2 + .749*MGLY + .749*HO2 + .251*CO2 + .251*HCHO + .251*OH - 0.749*XC';
k(:,i) = 2.30e-12;
Gstr{i,1} = 'M2BUTDAL_P1'; Gstr{i,2} = 'NO3';
fM2BUTDAL_P1(i)=fM2BUTDAL_P1(i)-1; fNO3(i)=fNO3(i)-1; fNO2(i)=fNO2(i)+1; fMGLY(i)=fMGLY(i)+.749; fHO2(i)=fHO2(i)+.749; fCO2(i)=fCO2(i)+.251; fHCHO(i)=fHCHO(i)+.251; fOH(i)=fOH(i)+.251; fXC(i)=fXC(i)-0.749;

i=i+1;
Rnames{i} = 'M2BUTDAL_P1 + HO2 = .9*CROOH + .1*MGLY - 5.5*XC';
k(:,i) = 1.27e-11;
Gstr{i,1} = 'M2BUTDAL_P1'; Gstr{i,2} = 'HO2';
fM2BUTDAL_P1(i)=fM2BUTDAL_P1(i)-1; fHO2(i)=fHO2(i)-1; fCROOH(i)=fCROOH(i)+.9; fMGLY(i)=fMGLY(i)+.1; fXC(i)=fXC(i)-5.5;

i=i+1;
Rnames{i} = 'M2BUTDAL_P1 + SumRO2 = SumRO2 + .624*MGLY + .374*HO2 + .25*CROOH + .126*CO2 + .126*HCHO + .126*OH - 2.124*XC';
k(:,i) = 2.03e-13;
Gstr{i,1} = 'M2BUTDAL_P1'; Gstr{i,2} = 'SumRO2';
fM2BUTDAL_P1(i)=fM2BUTDAL_P1(i)-1; fSumRO2(i)=fSumRO2(i)-1; fSumRO2(i)=fSumRO2(i)+1; fMGLY(i)=fMGLY(i)+.624; fHO2(i)=fHO2(i)+.374; fCROOH(i)=fCROOH(i)+.25; fCO2(i)=fCO2(i)+.126; fHCHO(i)=fHCHO(i)+.126; fOH(i)=fOH(i)+.126; fXC(i)=fXC(i)-2.124;

i=i+1;
Rnames{i} = 'M2BUTDAL_P1 + SumRCO3 = SumRCO3 + .799*MGLY + .599*HO2 + .201*CO2 + .201*HCHO + .201*OH - 0.799*XC';
k(:,i) = 1.37e-11;
Gstr{i,1} = 'M2BUTDAL_P1'; Gstr{i,2} = 'SumRCO3';
fM2BUTDAL_P1(i)=fM2BUTDAL_P1(i)-1; fSumRCO3(i)=fSumRCO3(i)-1; fSumRCO3(i)=fSumRCO3(i)+1; fMGLY(i)=fMGLY(i)+.799; fHO2(i)=fHO2(i)+.599; fCO2(i)=fCO2(i)+.201; fHCHO(i)=fHCHO(i)+.201; fOH(i)=fOH(i)+.201; fXC(i)=fXC(i)-0.799;

i=i+1;
Rnames{i} = 'M2BUTDAL_R3 = RCNO3 + OH + XC';
k(:,i) = 1.64e+11.*(T./300).^0.00.*exp(-7848.229./T);
Gstr{i,1} = 'M2BUTDAL_R3';
fM2BUTDAL_R3(i)=fM2BUTDAL_R3(i)-1; fRCNO3(i)=fRCNO3(i)+1; fOH(i)=fOH(i)+1; fXC(i)=fXC(i)+1;

i=i+1;
Rnames{i} = 'M2BUTDAL_R3 + NO = RCNO3 + .938*CO + .938*HO2 + .938*NO2 + .062*XC + .062*XN';
k(:,i) = 9.13e-12;
Gstr{i,1} = 'M2BUTDAL_R3'; Gstr{i,2} = 'NO';
fM2BUTDAL_R3(i)=fM2BUTDAL_R3(i)-1; fNO(i)=fNO(i)-1; fRCNO3(i)=fRCNO3(i)+1; fCO(i)=fCO(i)+.938; fHO2(i)=fHO2(i)+.938; fNO2(i)=fNO2(i)+.938; fXC(i)=fXC(i)+.062; fXN(i)=fXN(i)+.062;

i=i+1;
Rnames{i} = 'M2BUTDAL_P2 + NO = .938*NO2 + .469*MACO3 + .284*MALAH + .284*HO2 + .185*MGLY + .185*MEO2 + .185*CO + .062*RCNO3 + .815*XC + .185*SumRO2 + .469*SumRCO3';
k(:,i) = 9.13e-12;
Gstr{i,1} = 'M2BUTDAL_P2'; Gstr{i,2} = 'NO';
fM2BUTDAL_P2(i)=fM2BUTDAL_P2(i)-1; fNO(i)=fNO(i)-1; fNO2(i)=fNO2(i)+.938; fMACO3(i)=fMACO3(i)+.469; fMALAH(i)=fMALAH(i)+.284; fHO2(i)=fHO2(i)+.284; fMGLY(i)=fMGLY(i)+.185; fMEO2(i)=fMEO2(i)+.185; fCO(i)=fCO(i)+.185; fRCNO3(i)=fRCNO3(i)+.062; fXC(i)=fXC(i)+.815; fSumRO2(i)=fSumRO2(i)+.185; fSumRCO3(i)=fSumRCO3(i)+.469;

i=i+1;
Rnames{i} = 'M2BUTDAL_P2 + NO3 = NO2 + .5*MACO3 + .303*MALAH + .303*HO2 + .197*MGLY + .197*MEO2 + .197*CO + .803*XC + .197*SumRO2 + .5*SumRCO3';
k(:,i) = 2.30e-12;
Gstr{i,1} = 'M2BUTDAL_P2'; Gstr{i,2} = 'NO3';
fM2BUTDAL_P2(i)=fM2BUTDAL_P2(i)-1; fNO3(i)=fNO3(i)-1; fNO2(i)=fNO2(i)+1; fMACO3(i)=fMACO3(i)+.5; fMALAH(i)=fMALAH(i)+.303; fHO2(i)=fHO2(i)+.303; fMGLY(i)=fMGLY(i)+.197; fMEO2(i)=fMEO2(i)+.197; fCO(i)=fCO(i)+.197; fXC(i)=fXC(i)+.803; fSumRO2(i)=fSumRO2(i)+.197; fSumRCO3(i)=fSumRCO3(i)+.5;

i=i+1;
Rnames{i} = 'M2BUTDAL_P2 + HO2 = .45*RUOOH + .45*HPALD + .05*BACL + .05*MALAH - 0.35*XC';
k(:,i) = 1.78e-11;
Gstr{i,1} = 'M2BUTDAL_P2'; Gstr{i,2} = 'HO2';
fM2BUTDAL_P2(i)=fM2BUTDAL_P2(i)-1; fHO2(i)=fHO2(i)-1; fRUOOH(i)=fRUOOH(i)+.45; fHPALD(i)=fHPALD(i)+.45; fBACL(i)=fBACL(i)+.05; fMALAH(i)=fMALAH(i)+.05; fXC(i)=fXC(i)-0.35;

i=i+1;
Rnames{i} = 'M2BUTDAL_P2 + SumRO2 = SumRO2 + .277*MALAH + .25*MACO3 + .152*HO2 + .125*BACL + .125*OLEP + .125*HFON52M4 + .098*MGLY + .098*MEO2 + .098*CO + .277*XC + .098*SumRO2 + .25*SumRCO3';
k(:,i) = 2.57e-12;
Gstr{i,1} = 'M2BUTDAL_P2'; Gstr{i,2} = 'SumRO2';
fM2BUTDAL_P2(i)=fM2BUTDAL_P2(i)-1; fSumRO2(i)=fSumRO2(i)-1; fSumRO2(i)=fSumRO2(i)+1; fMALAH(i)=fMALAH(i)+.277; fMACO3(i)=fMACO3(i)+.25; fHO2(i)=fHO2(i)+.152; fBACL(i)=fBACL(i)+.125; fOLEP(i)=fOLEP(i)+.125; fHFON52M4(i)=fHFON52M4(i)+.125; fMGLY(i)=fMGLY(i)+.098; fMEO2(i)=fMEO2(i)+.098; fCO(i)=fCO(i)+.098; fXC(i)=fXC(i)+.277; fSumRO2(i)=fSumRO2(i)+.098; fSumRCO3(i)=fSumRCO3(i)+.25;

i=i+1;
Rnames{i} = 'M2BUTDAL_P2 + SumRCO3 = SumRCO3 + .4*MACO3 + .342*MALAH + .242*HO2 + .158*MGLY + .158*MEO2 + .158*CO + .1*BACL + .842*XC + .158*SumRO2 + .4*SumRCO3';
k(:,i) = 1.37e-11;
Gstr{i,1} = 'M2BUTDAL_P2'; Gstr{i,2} = 'SumRCO3';
fM2BUTDAL_P2(i)=fM2BUTDAL_P2(i)-1; fSumRCO3(i)=fSumRCO3(i)-1; fSumRCO3(i)=fSumRCO3(i)+1; fMACO3(i)=fMACO3(i)+.4; fMALAH(i)=fMALAH(i)+.342; fHO2(i)=fHO2(i)+.242; fMGLY(i)=fMGLY(i)+.158; fMEO2(i)=fMEO2(i)+.158; fCO(i)=fCO(i)+.158; fBACL(i)=fBACL(i)+.1; fXC(i)=fXC(i)+.842; fSumRO2(i)=fSumRO2(i)+.158; fSumRCO3(i)=fSumRCO3(i)+.4;

i=i+1;
Rnames{i} = 'O4X2PEAL + OH = .381*MACO3 + .368*O4X2PEAL_R1 + .251*O4X2PEAL_P1 + .381*XC + .251*SumRO2 + .381*SumRCO3';
k(:,i) = 5.67e-11;
Gstr{i,1} = 'O4X2PEAL'; Gstr{i,2} = 'OH';
fO4X2PEAL(i)=fO4X2PEAL(i)-1; fOH(i)=fOH(i)-1; fMACO3(i)=fMACO3(i)+.381; fO4X2PEAL_R1(i)=fO4X2PEAL_R1(i)+.368; fO4X2PEAL_P1(i)=fO4X2PEAL_P1(i)+.251; fXC(i)=fXC(i)+.381; fSumRO2(i)=fSumRO2(i)+.251; fSumRCO3(i)=fSumRCO3(i)+.381;

i=i+1;
Rnames{i} = 'O4X2PEAL + O3 = .742*MGLY + .496*HO2 + .451*CO2 + .421*RCHO2 + .377*CO + .258*GLY + .203*HCHO + .129*O4X2PEAL_P2 + .129*OH - 0.294*XC + .129*SumRO2';
k(:,i) = 4.80e-18;
Gstr{i,1} = 'O4X2PEAL'; Gstr{i,2} = 'O3';
fO4X2PEAL(i)=fO4X2PEAL(i)-1; fO3(i)=fO3(i)-1; fMGLY(i)=fMGLY(i)+.742; fHO2(i)=fHO2(i)+.496; fCO2(i)=fCO2(i)+.451; fRCHO2(i)=fRCHO2(i)+.421; fCO(i)=fCO(i)+.377; fGLY(i)=fGLY(i)+.258; fHCHO(i)=fHCHO(i)+.203; fO4X2PEAL_P2(i)=fO4X2PEAL_P2(i)+.129; fOH(i)=fOH(i)+.129; fXC(i)=fXC(i)-0.294; fSumRO2(i)=fSumRO2(i)+.129;

i=i+1;
Rnames{i} = 'O4X2PEAL + NO3 = .924*MACO3 + .065*O4X2PEAL_P3 + .011*RO2C + .011*RCNO3 + .009*CO2 + .009*OH + .002*CO + .002*HO2 + .001*RO2XC + .001*zRCNO3 + .92*XC + .924*XN + .077*SumRO2 + .924*SumRCO3';
k(:,i) = 1.00e-14;
Gstr{i,1} = 'O4X2PEAL'; Gstr{i,2} = 'NO3';
fO4X2PEAL(i)=fO4X2PEAL(i)-1; fNO3(i)=fNO3(i)-1; fMACO3(i)=fMACO3(i)+.924; fO4X2PEAL_P3(i)=fO4X2PEAL_P3(i)+.065; fRO2C(i)=fRO2C(i)+.011; fRCNO3(i)=fRCNO3(i)+.011; fCO2(i)=fCO2(i)+.009; fOH(i)=fOH(i)+.009; fCO(i)=fCO(i)+.002; fHO2(i)=fHO2(i)+.002; fRO2XC(i)=fRO2XC(i)+.001; fzRCNO3(i)=fzRCNO3(i)+.001; fXC(i)=fXC(i)+.92; fXN(i)=fXN(i)+.924; fSumRO2(i)=fSumRO2(i)+.077; fSumRCO3(i)=fSumRCO3(i)+.924;

i=i+1;
Rnames{i} = 'O4X2PEAL + HV = .7*MALAH + .7*OH + .7*MEO2 + .3*H3F2ONM5 + .7*SumRO2';
k(:,i) = JAFGS.*3.42e-1;
Gstr{i,1} = 'O4X2PEAL';
fO4X2PEAL(i)=fO4X2PEAL(i)-1; fMALAH(i)=fMALAH(i)+.7; fOH(i)=fOH(i)+.7; fMEO2(i)=fMEO2(i)+.7; fH3F2ONM5(i)=fH3F2ONM5(i)+.3; fSumRO2(i)=fSumRO2(i)+.7;

i=i+1;
Rnames{i} = 'O4X2PEAL_R1 = O4X2PEAL_P4 + SumRO2';
k(:,i) = 1.64e+11.*(T./300).^0.00.*exp(-7848.229./T);
Gstr{i,1} = 'O4X2PEAL_R1';
fO4X2PEAL_R1(i)=fO4X2PEAL_R1(i)-1; fO4X2PEAL_P4(i)=fO4X2PEAL_P4(i)+1; fSumRO2(i)=fSumRO2(i)+1;

i=i+1;
Rnames{i} = 'O4X2PEAL_R1 + NO = .948*HO2 + .948*NO2 + .86*MGLY + .86*GLY + .088*RCHO + .088*CO + .052*RCNO3 + .052*XC';
k(:,i) = 9.13e-12;
Gstr{i,1} = 'O4X2PEAL_R1'; Gstr{i,2} = 'NO';
fO4X2PEAL_R1(i)=fO4X2PEAL_R1(i)-1; fNO(i)=fNO(i)-1; fHO2(i)=fHO2(i)+.948; fNO2(i)=fNO2(i)+.948; fMGLY(i)=fMGLY(i)+.86; fGLY(i)=fGLY(i)+.86; fRCHO(i)=fRCHO(i)+.088; fCO(i)=fCO(i)+.088; fRCNO3(i)=fRCNO3(i)+.052; fXC(i)=fXC(i)+.052;

i=i+1;
Rnames{i} = 'O4X2PEAL_P1 + NO = .938*NO2 + .75*MECO3 + .75*CROOH + .376*MGLY + .188*HO2 + .062*RCNO3 - 3.876*XC + .75*SumRCO3';
k(:,i) = 9.13e-12;
Gstr{i,1} = 'O4X2PEAL_P1'; Gstr{i,2} = 'NO';
fO4X2PEAL_P1(i)=fO4X2PEAL_P1(i)-1; fNO(i)=fNO(i)-1; fNO2(i)=fNO2(i)+.938; fMECO3(i)=fMECO3(i)+.75; fCROOH(i)=fCROOH(i)+.75; fMGLY(i)=fMGLY(i)+.376; fHO2(i)=fHO2(i)+.188; fRCNO3(i)=fRCNO3(i)+.062; fXC(i)=fXC(i)-3.876; fSumRCO3(i)=fSumRCO3(i)+.75;

i=i+1;
Rnames{i} = 'O4X2PEAL_P1 + NO3 = NO2 + .8*MECO3 + .8*CROOH + .401*MGLY + .2*HO2 - 4.203*XC + .8*SumRCO3';
k(:,i) = 2.30e-12;
Gstr{i,1} = 'O4X2PEAL_P1'; Gstr{i,2} = 'NO3';
fO4X2PEAL_P1(i)=fO4X2PEAL_P1(i)-1; fNO3(i)=fNO3(i)-1; fNO2(i)=fNO2(i)+1; fMECO3(i)=fMECO3(i)+.8; fCROOH(i)=fCROOH(i)+.8; fMGLY(i)=fMGLY(i)+.401; fHO2(i)=fHO2(i)+.2; fXC(i)=fXC(i)-4.203; fSumRCO3(i)=fSumRCO3(i)+.8;

i=i+1;
Rnames{i} = 'O4X2PEAL_P1 + HO2 = .9*CROOH + .1*BACL - 2.6*XC';
k(:,i) = 1.78e-11;
Gstr{i,1} = 'O4X2PEAL_P1'; Gstr{i,2} = 'HO2';
fO4X2PEAL_P1(i)=fO4X2PEAL_P1(i)-1; fHO2(i)=fHO2(i)-1; fCROOH(i)=fCROOH(i)+.9; fBACL(i)=fBACL(i)+.1; fXC(i)=fXC(i)-2.6;

i=i+1;
Rnames{i} = 'O4X2PEAL_P1 + SumRO2 = SumRO2 + .65*CROOH + .4*MECO3 + .25*BACL + .2*MGLY + .1*HO2 - 2.6*XC + .4*SumRCO3';
k(:,i) = 2.57e-12;
Gstr{i,1} = 'O4X2PEAL_P1'; Gstr{i,2} = 'SumRO2';
fO4X2PEAL_P1(i)=fO4X2PEAL_P1(i)-1; fSumRO2(i)=fSumRO2(i)-1; fSumRO2(i)=fSumRO2(i)+1; fCROOH(i)=fCROOH(i)+.65; fMECO3(i)=fMECO3(i)+.4; fBACL(i)=fBACL(i)+.25; fMGLY(i)=fMGLY(i)+.2; fHO2(i)=fHO2(i)+.1; fXC(i)=fXC(i)-2.6; fSumRCO3(i)=fSumRCO3(i)+.4;

i=i+1;
Rnames{i} = 'O4X2PEAL_P1 + SumRCO3 = SumRCO3 + .64*MECO3 + .64*CROOH + .32*MGLY + .2*BACL + .16*HO2 - 3.16*XC + .64*SumRCO3';
k(:,i) = 1.37e-11;
Gstr{i,1} = 'O4X2PEAL_P1'; Gstr{i,2} = 'SumRCO3';
fO4X2PEAL_P1(i)=fO4X2PEAL_P1(i)-1; fSumRCO3(i)=fSumRCO3(i)-1; fSumRCO3(i)=fSumRCO3(i)+1; fMECO3(i)=fMECO3(i)+.64; fCROOH(i)=fCROOH(i)+.64; fMGLY(i)=fMGLY(i)+.32; fBACL(i)=fBACL(i)+.2; fHO2(i)=fHO2(i)+.16; fXC(i)=fXC(i)-3.16; fSumRCO3(i)=fSumRCO3(i)+.64;

i=i+1;
Rnames{i} = 'O4X2PEAL_P2 + NO = NO2 + .749*MGLY + .749*HO2 + .251*CO2 + .251*HCHO + .251*OH - 0.749*XC';
k(:,i) = 2.55e-12.*(T./300).^0.00.*exp(379.932./T);
Gstr{i,1} = 'O4X2PEAL_P2'; Gstr{i,2} = 'NO';
fO4X2PEAL_P2(i)=fO4X2PEAL_P2(i)-1; fNO(i)=fNO(i)-1; fNO2(i)=fNO2(i)+1; fMGLY(i)=fMGLY(i)+.749; fHO2(i)=fHO2(i)+.749; fCO2(i)=fCO2(i)+.251; fHCHO(i)=fHCHO(i)+.251; fOH(i)=fOH(i)+.251; fXC(i)=fXC(i)-0.749;

i=i+1;
Rnames{i} = 'O4X2PEAL_P2 + NO3 = NO2 + .749*MGLY + .749*HO2 + .251*CO2 + .251*HCHO + .251*OH - 0.749*XC';
k(:,i) = 2.30e-12;
Gstr{i,1} = 'O4X2PEAL_P2'; Gstr{i,2} = 'NO3';
fO4X2PEAL_P2(i)=fO4X2PEAL_P2(i)-1; fNO3(i)=fNO3(i)-1; fNO2(i)=fNO2(i)+1; fMGLY(i)=fMGLY(i)+.749; fHO2(i)=fHO2(i)+.749; fCO2(i)=fCO2(i)+.251; fHCHO(i)=fHCHO(i)+.251; fOH(i)=fOH(i)+.251; fXC(i)=fXC(i)-0.749;

i=i+1;
Rnames{i} = 'O4X2PEAL_P2 + HO2 = .9*CROOH + .1*MGLY - 5.5*XC';
k(:,i) = 1.27e-11;
Gstr{i,1} = 'O4X2PEAL_P2'; Gstr{i,2} = 'HO2';
fO4X2PEAL_P2(i)=fO4X2PEAL_P2(i)-1; fHO2(i)=fHO2(i)-1; fCROOH(i)=fCROOH(i)+.9; fMGLY(i)=fMGLY(i)+.1; fXC(i)=fXC(i)-5.5;

i=i+1;
Rnames{i} = 'O4X2PEAL_P2 + SumRO2 = SumRO2 + .624*MGLY + .374*HO2 + .25*CROOH + .126*CO2 + .126*HCHO + .126*OH - 2.124*XC';
k(:,i) = 2.03e-13;
Gstr{i,1} = 'O4X2PEAL_P2'; Gstr{i,2} = 'SumRO2';
fO4X2PEAL_P2(i)=fO4X2PEAL_P2(i)-1; fSumRO2(i)=fSumRO2(i)-1; fSumRO2(i)=fSumRO2(i)+1; fMGLY(i)=fMGLY(i)+.624; fHO2(i)=fHO2(i)+.374; fCROOH(i)=fCROOH(i)+.25; fCO2(i)=fCO2(i)+.126; fHCHO(i)=fHCHO(i)+.126; fOH(i)=fOH(i)+.126; fXC(i)=fXC(i)-2.124;

i=i+1;
Rnames{i} = 'O4X2PEAL_P2 + SumRCO3 = SumRCO3 + .799*MGLY + .599*HO2 + .201*CO2 + .201*HCHO + .201*OH - 0.799*XC';
k(:,i) = 1.37e-11;
Gstr{i,1} = 'O4X2PEAL_P2'; Gstr{i,2} = 'SumRCO3';
fO4X2PEAL_P2(i)=fO4X2PEAL_P2(i)-1; fSumRCO3(i)=fSumRCO3(i)-1; fSumRCO3(i)=fSumRCO3(i)+1; fMGLY(i)=fMGLY(i)+.799; fHO2(i)=fHO2(i)+.599; fCO2(i)=fCO2(i)+.201; fHCHO(i)=fHCHO(i)+.201; fOH(i)=fOH(i)+.201; fXC(i)=fXC(i)-0.799;

i=i+1;
Rnames{i} = 'O4X2PEAL_P3 + NO = RCNO3 + .938*MECO3 + .938*NO2 - 0.876*XC + .062*XN + .938*SumRCO3';
k(:,i) = 9.13e-12;
Gstr{i,1} = 'O4X2PEAL_P3'; Gstr{i,2} = 'NO';
fO4X2PEAL_P3(i)=fO4X2PEAL_P3(i)-1; fNO(i)=fNO(i)-1; fRCNO3(i)=fRCNO3(i)+1; fMECO3(i)=fMECO3(i)+.938; fNO2(i)=fNO2(i)+.938; fXC(i)=fXC(i)-0.876; fXN(i)=fXN(i)+.062; fSumRCO3(i)=fSumRCO3(i)+.938;

i=i+1;
Rnames{i} = 'O4X2PEAL_P3 + NO3 = MECO3 + RCNO3 + NO2 - 1*XC + SumRCO3';
k(:,i) = 2.30e-12;
Gstr{i,1} = 'O4X2PEAL_P3'; Gstr{i,2} = 'NO3';
fO4X2PEAL_P3(i)=fO4X2PEAL_P3(i)-1; fNO3(i)=fNO3(i)-1; fMECO3(i)=fMECO3(i)+1; fRCNO3(i)=fRCNO3(i)+1; fNO2(i)=fNO2(i)+1; fXC(i)=fXC(i)-1; fSumRCO3(i)=fSumRCO3(i)+1;

i=i+1;
Rnames{i} = 'O4X2PEAL_P3 + HO2 = RCNO3 + XC';
k(:,i) = 1.78e-11;
Gstr{i,1} = 'O4X2PEAL_P3'; Gstr{i,2} = 'HO2';
fO4X2PEAL_P3(i)=fO4X2PEAL_P3(i)-1; fHO2(i)=fHO2(i)-1; fRCNO3(i)=fRCNO3(i)+1; fXC(i)=fXC(i)+1;

i=i+1;
Rnames{i} = 'O4X2PEAL_P3 + SumRO2 = SumRO2 + RCNO3 + .5*MECO3 + .5*SumRCO3';
k(:,i) = 2.57e-12;
Gstr{i,1} = 'O4X2PEAL_P3'; Gstr{i,2} = 'SumRO2';
fO4X2PEAL_P3(i)=fO4X2PEAL_P3(i)-1; fSumRO2(i)=fSumRO2(i)-1; fSumRO2(i)=fSumRO2(i)+1; fRCNO3(i)=fRCNO3(i)+1; fMECO3(i)=fMECO3(i)+.5; fSumRCO3(i)=fSumRCO3(i)+.5;

i=i+1;
Rnames{i} = 'O4X2PEAL_P3 + SumRCO3 = SumRCO3 + RCNO3 + .8*MECO3 - 0.6*XC + .8*SumRCO3';
k(:,i) = 1.37e-11;
Gstr{i,1} = 'O4X2PEAL_P3'; Gstr{i,2} = 'SumRCO3';
fO4X2PEAL_P3(i)=fO4X2PEAL_P3(i)-1; fSumRCO3(i)=fSumRCO3(i)-1; fSumRCO3(i)=fSumRCO3(i)+1; fRCNO3(i)=fRCNO3(i)+1; fMECO3(i)=fMECO3(i)+.8; fXC(i)=fXC(i)-0.6; fSumRCO3(i)=fSumRCO3(i)+.8;

i=i+1;
Rnames{i} = 'O4X2PEAL_P4 + NO = .938*NO2 + .899*RCHO + .899*CO2 + .899*OH + .078*MGLY + .062*RCNO3 + .039*HO2 + .023*XC';
k(:,i) = 9.13e-12;
Gstr{i,1} = 'O4X2PEAL_P4'; Gstr{i,2} = 'NO';
fO4X2PEAL_P4(i)=fO4X2PEAL_P4(i)-1; fNO(i)=fNO(i)-1; fNO2(i)=fNO2(i)+.938; fRCHO(i)=fRCHO(i)+.899; fCO2(i)=fCO2(i)+.899; fOH(i)=fOH(i)+.899; fMGLY(i)=fMGLY(i)+.078; fRCNO3(i)=fRCNO3(i)+.062; fHO2(i)=fHO2(i)+.039; fXC(i)=fXC(i)+.023;

i=i+1;
Rnames{i} = 'O4X2PEAL_P4 + NO3 = NO2 + .959*RCHO + .959*CO2 + .959*OH + .083*MGLY + .041*HO2 - 0.044*XC';
k(:,i) = 2.30e-12;
Gstr{i,1} = 'O4X2PEAL_P4'; Gstr{i,2} = 'NO3';
fO4X2PEAL_P4(i)=fO4X2PEAL_P4(i)-1; fNO3(i)=fNO3(i)-1; fNO2(i)=fNO2(i)+1; fRCHO(i)=fRCHO(i)+.959; fCO2(i)=fCO2(i)+.959; fOH(i)=fOH(i)+.959; fMGLY(i)=fMGLY(i)+.083; fHO2(i)=fHO2(i)+.041; fXC(i)=fXC(i)-0.044;

i=i+1;
Rnames{i} = 'O4X2PEAL_P4 + HO2 = .9*CROOH + .1*BACL - 2.6*XC';
k(:,i) = 1.78e-11;
Gstr{i,1} = 'O4X2PEAL_P4'; Gstr{i,2} = 'HO2';
fO4X2PEAL_P4(i)=fO4X2PEAL_P4(i)-1; fHO2(i)=fHO2(i)-1; fCROOH(i)=fCROOH(i)+.9; fBACL(i)=fBACL(i)+.1; fXC(i)=fXC(i)-2.6;

i=i+1;
Rnames{i} = 'O4X2PEAL_P4 + SumRO2 = SumRO2 + .479*RCHO + .479*CO2 + .479*OH + .25*BACL + .25*CROOH + .041*MGLY + .021*HO2 - 0.518*XC';
k(:,i) = 2.57e-12;
Gstr{i,1} = 'O4X2PEAL_P4'; Gstr{i,2} = 'SumRO2';
fO4X2PEAL_P4(i)=fO4X2PEAL_P4(i)-1; fSumRO2(i)=fSumRO2(i)-1; fSumRO2(i)=fSumRO2(i)+1; fRCHO(i)=fRCHO(i)+.479; fCO2(i)=fCO2(i)+.479; fOH(i)=fOH(i)+.479; fBACL(i)=fBACL(i)+.25; fCROOH(i)=fCROOH(i)+.25; fMGLY(i)=fMGLY(i)+.041; fHO2(i)=fHO2(i)+.021; fXC(i)=fXC(i)-0.518;

i=i+1;
Rnames{i} = 'O4X2PEAL_P4 + SumRCO3 = SumRCO3 + .767*RCHO + .767*CO2 + .767*OH + .2*BACL + .066*MGLY + .033*HO2 + .167*XC';
k(:,i) = 1.37e-11;
Gstr{i,1} = 'O4X2PEAL_P4'; Gstr{i,2} = 'SumRCO3';
fO4X2PEAL_P4(i)=fO4X2PEAL_P4(i)-1; fSumRCO3(i)=fSumRCO3(i)-1; fSumRCO3(i)=fSumRCO3(i)+1; fRCHO(i)=fRCHO(i)+.767; fCO2(i)=fCO2(i)+.767; fOH(i)=fOH(i)+.767; fBACL(i)=fBACL(i)+.2; fMGLY(i)=fMGLY(i)+.066; fHO2(i)=fHO2(i)+.033; fXC(i)=fXC(i)+.167;

i=i+1;
Rnames{i} = 'H3XE25DO + OH = H3XE25DO_P1 + SumRO2';
k(:,i) = 4.71e-11;
Gstr{i,1} = 'H3XE25DO'; Gstr{i,2} = 'OH';
fH3XE25DO(i)=fH3XE25DO(i)-1; fOH(i)=fOH(i)-1; fH3XE25DO_P1(i)=fH3XE25DO_P1(i)+1; fSumRO2(i)=fSumRO2(i)+1;

i=i+1;
Rnames{i} = 'H3XE25DO + O3 = MGLY + .5*H3XE25DO_P2 + .5*RCHO2 + .5*OH + .5*CO + .5*SumRO2';
k(:,i) = 1.64e-14;
Gstr{i,1} = 'H3XE25DO'; Gstr{i,2} = 'O3';
fH3XE25DO(i)=fH3XE25DO(i)-1; fO3(i)=fO3(i)-1; fMGLY(i)=fMGLY(i)+1; fH3XE25DO_P2(i)=fH3XE25DO_P2(i)+.5; fRCHO2(i)=fRCHO2(i)+.5; fOH(i)=fOH(i)+.5; fCO(i)=fCO(i)+.5; fSumRO2(i)=fSumRO2(i)+.5;

i=i+1;
Rnames{i} = 'H3XE25DO + NO3 = H3XE25DO_P3 + SumRO2';
k(:,i) = 1.80e-18;
Gstr{i,1} = 'H3XE25DO'; Gstr{i,2} = 'NO3';
fH3XE25DO(i)=fH3XE25DO(i)-1; fNO3(i)=fNO3(i)-1; fH3XE25DO_P3(i)=fH3XE25DO_P3(i)+1; fSumRO2(i)=fSumRO2(i)+1;

i=i+1;
Rnames{i} = 'H3XE25DO_P1 + NO = .93*NO2 + .743*MECO3 + .743*RCHO + .372*MGLY + .186*HO2 + .07*RCNO3 + .146*XC + .743*SumRCO3';
k(:,i) = 9.13e-12;
Gstr{i,1} = 'H3XE25DO_P1'; Gstr{i,2} = 'NO';
fH3XE25DO_P1(i)=fH3XE25DO_P1(i)-1; fNO(i)=fNO(i)-1; fNO2(i)=fNO2(i)+.93; fMECO3(i)=fMECO3(i)+.743; fRCHO(i)=fRCHO(i)+.743; fMGLY(i)=fMGLY(i)+.372; fHO2(i)=fHO2(i)+.186; fRCNO3(i)=fRCNO3(i)+.07; fXC(i)=fXC(i)+.146; fSumRCO3(i)=fSumRCO3(i)+.743;

i=i+1;
Rnames{i} = 'H3XE25DO_P1 + NO3 = NO2 + .8*MECO3 + .8*RCHO + .401*MGLY + .2*HO2 - 0.003*XC + .8*SumRCO3';
k(:,i) = 2.30e-12;
Gstr{i,1} = 'H3XE25DO_P1'; Gstr{i,2} = 'NO3';
fH3XE25DO_P1(i)=fH3XE25DO_P1(i)-1; fNO3(i)=fNO3(i)-1; fNO2(i)=fNO2(i)+1; fMECO3(i)=fMECO3(i)+.8; fRCHO(i)=fRCHO(i)+.8; fMGLY(i)=fMGLY(i)+.401; fHO2(i)=fHO2(i)+.2; fXC(i)=fXC(i)-0.003; fSumRCO3(i)=fSumRCO3(i)+.8;

i=i+1;
Rnames{i} = 'H3XE25DO_P1 + HO2 = .9*CROOH + .1*BACL - 1.6*XC';
k(:,i) = 1.95e-11;
Gstr{i,1} = 'H3XE25DO_P1'; Gstr{i,2} = 'HO2';
fH3XE25DO_P1(i)=fH3XE25DO_P1(i)-1; fHO2(i)=fHO2(i)-1; fCROOH(i)=fCROOH(i)+.9; fBACL(i)=fBACL(i)+.1; fXC(i)=fXC(i)-1.6;

i=i+1;
Rnames{i} = 'H3XE25DO_P1 + SumRO2 = SumRO2 + .4*MECO3 + .4*RCHO + .25*BACL + .25*KET2 + .2*MGLY + .1*HO2 + .5*XC + .4*SumRCO3';
k(:,i) = 2.57e-12;
Gstr{i,1} = 'H3XE25DO_P1'; Gstr{i,2} = 'SumRO2';
fH3XE25DO_P1(i)=fH3XE25DO_P1(i)-1; fSumRO2(i)=fSumRO2(i)-1; fSumRO2(i)=fSumRO2(i)+1; fMECO3(i)=fMECO3(i)+.4; fRCHO(i)=fRCHO(i)+.4; fBACL(i)=fBACL(i)+.25; fKET2(i)=fKET2(i)+.25; fMGLY(i)=fMGLY(i)+.2; fHO2(i)=fHO2(i)+.1; fXC(i)=fXC(i)+.5; fSumRCO3(i)=fSumRCO3(i)+.4;

i=i+1;
Rnames{i} = 'H3XE25DO_P1 + SumRCO3 = SumRCO3 + .64*MECO3 + .64*RCHO + .32*MGLY + .2*BACL + .16*HO2 + .4*XC + .64*SumRCO3';
k(:,i) = 1.37e-11;
Gstr{i,1} = 'H3XE25DO_P1'; Gstr{i,2} = 'SumRCO3';
fH3XE25DO_P1(i)=fH3XE25DO_P1(i)-1; fSumRCO3(i)=fSumRCO3(i)-1; fSumRCO3(i)=fSumRCO3(i)+1; fMECO3(i)=fMECO3(i)+.64; fRCHO(i)=fRCHO(i)+.64; fMGLY(i)=fMGLY(i)+.32; fBACL(i)=fBACL(i)+.2; fHO2(i)=fHO2(i)+.16; fXC(i)=fXC(i)+.4; fSumRCO3(i)=fSumRCO3(i)+.64;

i=i+1;
Rnames{i} = 'H3XE25DO_P2 + NO = NO2 + .749*MGLY + .749*HO2 + .251*CO2 + .251*HCHO + .251*OH - 0.749*XC';
k(:,i) = 2.55e-12.*(T./300).^0.00.*exp(379.932./T);
Gstr{i,1} = 'H3XE25DO_P2'; Gstr{i,2} = 'NO';
fH3XE25DO_P2(i)=fH3XE25DO_P2(i)-1; fNO(i)=fNO(i)-1; fNO2(i)=fNO2(i)+1; fMGLY(i)=fMGLY(i)+.749; fHO2(i)=fHO2(i)+.749; fCO2(i)=fCO2(i)+.251; fHCHO(i)=fHCHO(i)+.251; fOH(i)=fOH(i)+.251; fXC(i)=fXC(i)-0.749;

i=i+1;
Rnames{i} = 'H3XE25DO_P2 + NO3 = NO2 + .749*MGLY + .749*HO2 + .251*CO2 + .251*HCHO + .251*OH - 0.749*XC';
k(:,i) = 2.30e-12;
Gstr{i,1} = 'H3XE25DO_P2'; Gstr{i,2} = 'NO3';
fH3XE25DO_P2(i)=fH3XE25DO_P2(i)-1; fNO3(i)=fNO3(i)-1; fNO2(i)=fNO2(i)+1; fMGLY(i)=fMGLY(i)+.749; fHO2(i)=fHO2(i)+.749; fCO2(i)=fCO2(i)+.251; fHCHO(i)=fHCHO(i)+.251; fOH(i)=fOH(i)+.251; fXC(i)=fXC(i)-0.749;

i=i+1;
Rnames{i} = 'H3XE25DO_P2 + HO2 = .9*CROOH + .1*MGLY - 5.5*XC';
k(:,i) = 1.27e-11;
Gstr{i,1} = 'H3XE25DO_P2'; Gstr{i,2} = 'HO2';
fH3XE25DO_P2(i)=fH3XE25DO_P2(i)-1; fHO2(i)=fHO2(i)-1; fCROOH(i)=fCROOH(i)+.9; fMGLY(i)=fMGLY(i)+.1; fXC(i)=fXC(i)-5.5;

i=i+1;
Rnames{i} = 'H3XE25DO_P2 + SumRO2 = SumRO2 + .624*MGLY + .374*HO2 + .25*CROOH + .126*CO2 + .126*HCHO + .126*OH - 2.124*XC';
k(:,i) = 2.03e-13;
Gstr{i,1} = 'H3XE25DO_P2'; Gstr{i,2} = 'SumRO2';
fH3XE25DO_P2(i)=fH3XE25DO_P2(i)-1; fSumRO2(i)=fSumRO2(i)-1; fSumRO2(i)=fSumRO2(i)+1; fMGLY(i)=fMGLY(i)+.624; fHO2(i)=fHO2(i)+.374; fCROOH(i)=fCROOH(i)+.25; fCO2(i)=fCO2(i)+.126; fHCHO(i)=fHCHO(i)+.126; fOH(i)=fOH(i)+.126; fXC(i)=fXC(i)-2.124;

i=i+1;
Rnames{i} = 'H3XE25DO_P2 + SumRCO3 = SumRCO3 + .799*MGLY + .599*HO2 + .201*CO2 + .201*HCHO + .201*OH - 0.799*XC';
k(:,i) = 1.37e-11;
Gstr{i,1} = 'H3XE25DO_P2'; Gstr{i,2} = 'SumRCO3';
fH3XE25DO_P2(i)=fH3XE25DO_P2(i)-1; fSumRCO3(i)=fSumRCO3(i)-1; fSumRCO3(i)=fSumRCO3(i)+1; fMGLY(i)=fMGLY(i)+.799; fHO2(i)=fHO2(i)+.599; fCO2(i)=fCO2(i)+.201; fHCHO(i)=fHCHO(i)+.201; fOH(i)=fOH(i)+.201; fXC(i)=fXC(i)-0.799;

i=i+1;
Rnames{i} = 'H3XE25DO_P3 + NO = .936*RCNO3 + .847*NO2 + .783*MECO3 + .064*MACO3 + .064*HCHO + .37*XC + .217*XN + .847*SumRCO3';
k(:,i) = 9.13e-12;
Gstr{i,1} = 'H3XE25DO_P3'; Gstr{i,2} = 'NO';
fH3XE25DO_P3(i)=fH3XE25DO_P3(i)-1; fNO(i)=fNO(i)-1; fRCNO3(i)=fRCNO3(i)+.936; fNO2(i)=fNO2(i)+.847; fMECO3(i)=fMECO3(i)+.783; fMACO3(i)=fMACO3(i)+.064; fHCHO(i)=fHCHO(i)+.064; fXC(i)=fXC(i)+.37; fXN(i)=fXN(i)+.217; fSumRCO3(i)=fSumRCO3(i)+.847;

i=i+1;
Rnames{i} = 'H3XE25DO_P3 + NO3 = NO2 + .925*MECO3 + .925*RCNO3 + .075*MACO3 + .075*HCHO + .075*XC + .075*XN + SumRCO3';
k(:,i) = 2.30e-12;
Gstr{i,1} = 'H3XE25DO_P3'; Gstr{i,2} = 'NO3';
fH3XE25DO_P3(i)=fH3XE25DO_P3(i)-1; fNO3(i)=fNO3(i)-1; fNO2(i)=fNO2(i)+1; fMECO3(i)=fMECO3(i)+.925; fRCNO3(i)=fRCNO3(i)+.925; fMACO3(i)=fMACO3(i)+.075; fHCHO(i)=fHCHO(i)+.075; fXC(i)=fXC(i)+.075; fXN(i)=fXN(i)+.075; fSumRCO3(i)=fSumRCO3(i)+1;

i=i+1;
Rnames{i} = 'H3XE25DO_P3 + HO2 = .925*RCNO3 + .068*HPALD + .008*AFG3 + 1.912*XC + .075*XN';
k(:,i) = 1.95e-11;
Gstr{i,1} = 'H3XE25DO_P3'; Gstr{i,2} = 'HO2';
fH3XE25DO_P3(i)=fH3XE25DO_P3(i)-1; fHO2(i)=fHO2(i)-1; fRCNO3(i)=fRCNO3(i)+.925; fHPALD(i)=fHPALD(i)+.068; fAFG3(i)=fAFG3(i)+.008; fXC(i)=fXC(i)+1.912; fXN(i)=fXN(i)+.075;

i=i+1;
Rnames{i} = 'H3XE25DO_P3 + SumRO2 = SumRO2 + .925*RCNO3 + .462*MECO3 + .038*MACO3 + .038*HCHO + .038*AFG3 + .958*XC + .075*XN + .5*SumRCO3';
k(:,i) = 2.39e-12;
Gstr{i,1} = 'H3XE25DO_P3'; Gstr{i,2} = 'SumRO2';
fH3XE25DO_P3(i)=fH3XE25DO_P3(i)-1; fSumRO2(i)=fSumRO2(i)-1; fSumRO2(i)=fSumRO2(i)+1; fRCNO3(i)=fRCNO3(i)+.925; fMECO3(i)=fMECO3(i)+.462; fMACO3(i)=fMACO3(i)+.038; fHCHO(i)=fHCHO(i)+.038; fAFG3(i)=fAFG3(i)+.038; fXC(i)=fXC(i)+.958; fXN(i)=fXN(i)+.075; fSumRCO3(i)=fSumRCO3(i)+.5;

i=i+1;
Rnames{i} = 'H3XE25DO_P3 + SumRCO3 = SumRCO3 + .925*RCNO3 + .74*MECO3 + .06*MACO3 + .06*HCHO + .015*AFG3 + .43*XC + .075*XN + .8*SumRCO3';
k(:,i) = 1.37e-11;
Gstr{i,1} = 'H3XE25DO_P3'; Gstr{i,2} = 'SumRCO3';
fH3XE25DO_P3(i)=fH3XE25DO_P3(i)-1; fSumRCO3(i)=fSumRCO3(i)-1; fSumRCO3(i)=fSumRCO3(i)+1; fRCNO3(i)=fRCNO3(i)+.925; fMECO3(i)=fMECO3(i)+.74; fMACO3(i)=fMACO3(i)+.06; fHCHO(i)=fHCHO(i)+.06; fAFG3(i)=fAFG3(i)+.015; fXC(i)=fXC(i)+.43; fXN(i)=fXN(i)+.075; fSumRCO3(i)=fSumRCO3(i)+.8;

i=i+1;
Rnames{i} = 'HFON52 + OH = .93*HFON52_P1 + .07*MALAH + .07*HO2 + .93*SumRO2';
k(:,i) = 3.83e-11;
Gstr{i,1} = 'HFON52'; Gstr{i,2} = 'OH';
fHFON52(i)=fHFON52(i)-1; fOH(i)=fOH(i)-1; fHFON52_P1(i)=fHFON52_P1(i)+.93; fMALAH(i)=fMALAH(i)+.07; fHO2(i)=fHO2(i)+.07; fSumRO2(i)=fSumRO2(i)+.93;

i=i+1;
Rnames{i} = 'HFON52 + O3 = RCHO2 + XC';
k(:,i) = 5.07e-17;
Gstr{i,1} = 'HFON52'; Gstr{i,2} = 'O3';
fHFON52(i)=fHFON52(i)-1; fO3(i)=fO3(i)-1; fRCHO2(i)=fRCHO2(i)+1; fXC(i)=fXC(i)+1;

i=i+1;
Rnames{i} = 'HFON52 + NO3 = .738*HFON52_P2 + .262*MALAH + .262*HO2 + .262*XN + .738*SumRO2';
k(:,i) = 5.20e-15;
Gstr{i,1} = 'HFON52'; Gstr{i,2} = 'NO3';
fHFON52(i)=fHFON52(i)-1; fNO3(i)=fNO3(i)-1; fHFON52_P2(i)=fHFON52_P2(i)+.738; fMALAH(i)=fMALAH(i)+.262; fHO2(i)=fHO2(i)+.262; fXN(i)=fXN(i)+.262; fSumRO2(i)=fSumRO2(i)+.738;

i=i+1;
Rnames{i} = 'HFON52_P1 + NO = .974*NO2 + .968*HO2 + .69*MGLY + .278*RCHO + .026*RCNO3 + .007*R2CO3 + .693*XC + .007*SumRCO3';
k(:,i) = 9.13e-12;
Gstr{i,1} = 'HFON52_P1'; Gstr{i,2} = 'NO';
fHFON52_P1(i)=fHFON52_P1(i)-1; fNO(i)=fNO(i)-1; fNO2(i)=fNO2(i)+.974; fHO2(i)=fHO2(i)+.968; fMGLY(i)=fMGLY(i)+.69; fRCHO(i)=fRCHO(i)+.278; fRCNO3(i)=fRCNO3(i)+.026; fR2CO3(i)=fR2CO3(i)+.007; fXC(i)=fXC(i)+.693; fSumRCO3(i)=fSumRCO3(i)+.007;

i=i+1;
Rnames{i} = 'HFON52_P1 + NO3 = NO2 + .993*HO2 + .708*MGLY + .285*RCHO + .007*R2CO3 + .715*XC + .007*SumRCO3';
k(:,i) = 2.30e-12;
Gstr{i,1} = 'HFON52_P1'; Gstr{i,2} = 'NO3';
fHFON52_P1(i)=fHFON52_P1(i)-1; fNO3(i)=fNO3(i)-1; fNO2(i)=fNO2(i)+1; fHO2(i)=fHO2(i)+.993; fMGLY(i)=fMGLY(i)+.708; fRCHO(i)=fRCHO(i)+.285; fR2CO3(i)=fR2CO3(i)+.007; fXC(i)=fXC(i)+.715; fSumRCO3(i)=fSumRCO3(i)+.007;

i=i+1;
Rnames{i} = 'HFON52_P1 + HO2 = .935*ROOH + .065*BACL - 0.935*XC';
k(:,i) = 1.61e-11;
Gstr{i,1} = 'HFON52_P1'; Gstr{i,2} = 'HO2';
fHFON52_P1(i)=fHFON52_P1(i)-1; fHO2(i)=fHO2(i)-1; fROOH(i)=fROOH(i)+.935; fBACL(i)=fBACL(i)+.065; fXC(i)=fXC(i)-0.935;

i=i+1;
Rnames{i} = 'HFON52_P1 + SumRO2 = SumRO2 + .497*HO2 + .354*MGLY + .25*OTH3 + .163*BACL + .143*RCHO + .087*KET2 + .003*R2CO3 + .183*XC + .003*SumRCO3';
k(:,i) = 2.57e-12;
Gstr{i,1} = 'HFON52_P1'; Gstr{i,2} = 'SumRO2';
fHFON52_P1(i)=fHFON52_P1(i)-1; fSumRO2(i)=fSumRO2(i)-1; fSumRO2(i)=fSumRO2(i)+1; fHO2(i)=fHO2(i)+.497; fMGLY(i)=fMGLY(i)+.354; fOTH3(i)=fOTH3(i)+.25; fBACL(i)=fBACL(i)+.163; fRCHO(i)=fRCHO(i)+.143; fKET2(i)=fKET2(i)+.087; fR2CO3(i)=fR2CO3(i)+.003; fXC(i)=fXC(i)+.183; fSumRCO3(i)=fSumRCO3(i)+.003;

i=i+1;
Rnames{i} = 'HFON52_P1 + SumRCO3 = SumRCO3 + .795*HO2 + .567*MGLY + .228*RCHO + .131*BACL + .069*KET2 + .005*R2CO3 + .434*XC + .005*SumRCO3';
k(:,i) = 1.37e-11;
Gstr{i,1} = 'HFON52_P1'; Gstr{i,2} = 'SumRCO3';
fHFON52_P1(i)=fHFON52_P1(i)-1; fSumRCO3(i)=fSumRCO3(i)-1; fSumRCO3(i)=fSumRCO3(i)+1; fHO2(i)=fHO2(i)+.795; fMGLY(i)=fMGLY(i)+.567; fRCHO(i)=fRCHO(i)+.228; fBACL(i)=fBACL(i)+.131; fKET2(i)=fKET2(i)+.069; fR2CO3(i)=fR2CO3(i)+.005; fXC(i)=fXC(i)+.434; fSumRCO3(i)=fSumRCO3(i)+.005;

i=i+1;
Rnames{i} = 'HFON52_P2 + NO = 1.462*NO2 + .513*RCNO3 + .487*HO2 + .388*RCHO + .099*MGLY + .099*XC + .025*XN';
k(:,i) = 9.13e-12;
Gstr{i,1} = 'HFON52_P2'; Gstr{i,2} = 'NO';
fHFON52_P2(i)=fHFON52_P2(i)-1; fNO(i)=fNO(i)-1; fNO2(i)=fNO2(i)+1.462; fRCNO3(i)=fRCNO3(i)+.513; fHO2(i)=fHO2(i)+.487; fRCHO(i)=fRCHO(i)+.388; fMGLY(i)=fMGLY(i)+.099; fXC(i)=fXC(i)+.099; fXN(i)=fXN(i)+.025;

i=i+1;
Rnames{i} = 'HFON52_P2 + NO3 = 1.5*NO2 + .5*RCNO3 + .5*HO2 + .399*RCHO + .101*MGLY + .101*XC';
k(:,i) = 2.30e-12;
Gstr{i,1} = 'HFON52_P2'; Gstr{i,2} = 'NO3';
fHFON52_P2(i)=fHFON52_P2(i)-1; fNO3(i)=fNO3(i)-1; fNO2(i)=fNO2(i)+1.5; fRCNO3(i)=fRCNO3(i)+.5; fHO2(i)=fHO2(i)+.5; fRCHO(i)=fRCHO(i)+.399; fMGLY(i)=fMGLY(i)+.101; fXC(i)=fXC(i)+.101;

i=i+1;
Rnames{i} = 'HFON52_P2 + HO2 = RCNO3';
k(:,i) = 1.61e-11;
Gstr{i,1} = 'HFON52_P2'; Gstr{i,2} = 'HO2';
fHFON52_P2(i)=fHFON52_P2(i)-1; fHO2(i)=fHO2(i)-1; fRCNO3(i)=fRCNO3(i)+1;

i=i+1;
Rnames{i} = 'HFON52_P2 + SumRO2 = SumRO2 + .75*RCNO3 + .25*NO2 + .25*HO2 + .199*RCHO + .051*MGLY + .051*XC';
k(:,i) = 2.57e-12;
Gstr{i,1} = 'HFON52_P2'; Gstr{i,2} = 'SumRO2';
fHFON52_P2(i)=fHFON52_P2(i)-1; fSumRO2(i)=fSumRO2(i)-1; fSumRO2(i)=fSumRO2(i)+1; fRCNO3(i)=fRCNO3(i)+.75; fNO2(i)=fNO2(i)+.25; fHO2(i)=fHO2(i)+.25; fRCHO(i)=fRCHO(i)+.199; fMGLY(i)=fMGLY(i)+.051; fXC(i)=fXC(i)+.051;

i=i+1;
Rnames{i} = 'HFON52_P2 + SumRCO3 = SumRCO3 + .6*RCNO3 + .4*NO2 + .4*HO2 + .319*RCHO + .081*MGLY + .081*XC';
k(:,i) = 1.37e-11;
Gstr{i,1} = 'HFON52_P2'; Gstr{i,2} = 'SumRCO3';
fHFON52_P2(i)=fHFON52_P2(i)-1; fSumRCO3(i)=fSumRCO3(i)-1; fSumRCO3(i)=fSumRCO3(i)+1; fRCNO3(i)=fRCNO3(i)+.6; fNO2(i)=fNO2(i)+.4; fHO2(i)=fHO2(i)+.4; fRCHO(i)=fRCHO(i)+.319; fMGLY(i)=fMGLY(i)+.081; fXC(i)=fXC(i)+.081;

i=i+1;
Rnames{i} = 'HFON52M4 + OH = .941*HFON52M4_P1 + .059*MALAH + .059*HO2 + .059*XC + .941*SumRO2';
k(:,i) = 4.61e-11;
Gstr{i,1} = 'HFON52M4'; Gstr{i,2} = 'OH';
fHFON52M4(i)=fHFON52M4(i)-1; fOH(i)=fOH(i)-1; fHFON52M4_P1(i)=fHFON52M4_P1(i)+.941; fMALAH(i)=fMALAH(i)+.059; fHO2(i)=fHO2(i)+.059; fXC(i)=fXC(i)+.059; fSumRO2(i)=fSumRO2(i)+.941;

i=i+1;
Rnames{i} = 'HFON52M4 + O3 = RCHO2 + 2*XC';
k(:,i) = 1.63e-16;
Gstr{i,1} = 'HFON52M4'; Gstr{i,2} = 'O3';
fHFON52M4(i)=fHFON52M4(i)-1; fO3(i)=fO3(i)-1; fRCHO2(i)=fRCHO2(i)+1; fXC(i)=fXC(i)+2;

i=i+1;
Rnames{i} = 'HFON52M4 + NO3 = .984*HFON52M4_P2 + .016*MALAH + .016*HO2 + .016*XC + .016*XN + .984*SumRO2';
k(:,i) = 7.72e-14;
Gstr{i,1} = 'HFON52M4'; Gstr{i,2} = 'NO3';
fHFON52M4(i)=fHFON52M4(i)-1; fNO3(i)=fNO3(i)-1; fHFON52M4_P2(i)=fHFON52M4_P2(i)+.984; fMALAH(i)=fMALAH(i)+.016; fHO2(i)=fHO2(i)+.016; fXC(i)=fXC(i)+.016; fXN(i)=fXN(i)+.016; fSumRO2(i)=fSumRO2(i)+.984;

i=i+1;
Rnames{i} = 'HFON52M4_P1 + NO = .938*HO2 + .938*NO2 + .519*KET2 + .419*MGLY + .062*RCNO3 + .381*XC';
k(:,i) = 9.13e-12;
Gstr{i,1} = 'HFON52M4_P1'; Gstr{i,2} = 'NO';
fHFON52M4_P1(i)=fHFON52M4_P1(i)-1; fNO(i)=fNO(i)-1; fHO2(i)=fHO2(i)+.938; fNO2(i)=fNO2(i)+.938; fKET2(i)=fKET2(i)+.519; fMGLY(i)=fMGLY(i)+.419; fRCNO3(i)=fRCNO3(i)+.062; fXC(i)=fXC(i)+.381;

i=i+1;
Rnames{i} = 'HFON52M4_P1 + NO3 = HO2 + NO2 + .554*KET2 + .446*MGLY + .338*XC';
k(:,i) = 2.30e-12;
Gstr{i,1} = 'HFON52M4_P1'; Gstr{i,2} = 'NO3';
fHFON52M4_P1(i)=fHFON52M4_P1(i)-1; fNO3(i)=fNO3(i)-1; fHO2(i)=fHO2(i)+1; fNO2(i)=fNO2(i)+1; fKET2(i)=fKET2(i)+.554; fMGLY(i)=fMGLY(i)+.446; fXC(i)=fXC(i)+.338;

i=i+1;
Rnames{i} = 'HFON52M4_P1 + HO2 = .958*ROOH + .042*BACL + .042*XC';
k(:,i) = 1.78e-11;
Gstr{i,1} = 'HFON52M4_P1'; Gstr{i,2} = 'HO2';
fHFON52M4_P1(i)=fHFON52M4_P1(i)-1; fHO2(i)=fHO2(i)-1; fROOH(i)=fROOH(i)+.958; fBACL(i)=fBACL(i)+.042; fXC(i)=fXC(i)+.042;

i=i+1;
Rnames{i} = 'HFON52M4_P1 + SumRO2 = SumRO2 + .5*HO2 + .394*OTH3 + .277*KET2 + .223*MGLY + .106*BACL + .669*XC';
k(:,i) = 1.09e-12;
Gstr{i,1} = 'HFON52M4_P1'; Gstr{i,2} = 'SumRO2';
fHFON52M4_P1(i)=fHFON52M4_P1(i)-1; fSumRO2(i)=fSumRO2(i)-1; fSumRO2(i)=fSumRO2(i)+1; fHO2(i)=fHO2(i)+.5; fOTH3(i)=fOTH3(i)+.394; fKET2(i)=fKET2(i)+.277; fMGLY(i)=fMGLY(i)+.223; fBACL(i)=fBACL(i)+.106; fXC(i)=fXC(i)+.669;

i=i+1;
Rnames{i} = 'HFON52M4_P1 + SumRCO3 = SumRCO3 + .915*HO2 + .554*KET2 + .361*MGLY + .085*BACL + .253*XC';
k(:,i) = 1.37e-11;
Gstr{i,1} = 'HFON52M4_P1'; Gstr{i,2} = 'SumRCO3';
fHFON52M4_P1(i)=fHFON52M4_P1(i)-1; fSumRCO3(i)=fSumRCO3(i)-1; fSumRCO3(i)=fSumRCO3(i)+1; fHO2(i)=fHO2(i)+.915; fKET2(i)=fKET2(i)+.554; fMGLY(i)=fMGLY(i)+.361; fBACL(i)=fBACL(i)+.085; fXC(i)=fXC(i)+.253;

i=i+1;
Rnames{i} = 'HFON52M4_P2 + NO = 1.876*NO2 + .909*MGLY + .062*RCNO3 + .029*RCHO + 1.909*XC + .062*XN';
k(:,i) = 9.13e-12;
Gstr{i,1} = 'HFON52M4_P2'; Gstr{i,2} = 'NO';
fHFON52M4_P2(i)=fHFON52M4_P2(i)-1; fNO(i)=fNO(i)-1; fNO2(i)=fNO2(i)+1.876; fMGLY(i)=fMGLY(i)+.909; fRCNO3(i)=fRCNO3(i)+.062; fRCHO(i)=fRCHO(i)+.029; fXC(i)=fXC(i)+1.909; fXN(i)=fXN(i)+.062;

i=i+1;
Rnames{i} = 'HFON52M4_P2 + NO3 = 2*NO2 + .969*MGLY + .031*RCHO + 1.969*XC';
k(:,i) = 2.30e-12;
Gstr{i,1} = 'HFON52M4_P2'; Gstr{i,2} = 'NO3';
fHFON52M4_P2(i)=fHFON52M4_P2(i)-1; fNO3(i)=fNO3(i)-1; fNO2(i)=fNO2(i)+2; fMGLY(i)=fMGLY(i)+.969; fRCHO(i)=fRCHO(i)+.031; fXC(i)=fXC(i)+1.969;

i=i+1;
Rnames{i} = 'HFON52M4_P2 + HO2 = RCNO3 + XC';
k(:,i) = 1.78e-11;
Gstr{i,1} = 'HFON52M4_P2'; Gstr{i,2} = 'HO2';
fHFON52M4_P2(i)=fHFON52M4_P2(i)-1; fHO2(i)=fHO2(i)-1; fRCNO3(i)=fRCNO3(i)+1; fXC(i)=fXC(i)+1;

i=i+1;
Rnames{i} = 'HFON52M4_P2 + SumRO2 = SumRO2 + .5*NO2 + .5*RCNO3 + .485*MGLY + .015*RCHO + 1.485*XC';
k(:,i) = 2.57e-12;
Gstr{i,1} = 'HFON52M4_P2'; Gstr{i,2} = 'SumRO2';
fHFON52M4_P2(i)=fHFON52M4_P2(i)-1; fSumRO2(i)=fSumRO2(i)-1; fSumRO2(i)=fSumRO2(i)+1; fNO2(i)=fNO2(i)+.5; fRCNO3(i)=fRCNO3(i)+.5; fMGLY(i)=fMGLY(i)+.485; fRCHO(i)=fRCHO(i)+.015; fXC(i)=fXC(i)+1.485;

i=i+1;
Rnames{i} = 'HFON52M4_P2 + SumRCO3 = SumRCO3 + .8*NO2 + .775*MGLY + .2*RCNO3 + .025*RCHO + 1.775*XC';
k(:,i) = 1.37e-11;
Gstr{i,1} = 'HFON52M4_P2'; Gstr{i,2} = 'SumRCO3';
fHFON52M4_P2(i)=fHFON52M4_P2(i)-1; fSumRCO3(i)=fSumRCO3(i)-1; fSumRCO3(i)=fSumRCO3(i)+1; fNO2(i)=fNO2(i)+.8; fMGLY(i)=fMGLY(i)+.775; fRCNO3(i)=fRCNO3(i)+.2; fRCHO(i)=fRCHO(i)+.025; fXC(i)=fXC(i)+1.775;

i=i+1;
Rnames{i} = 'HFON52M3 + OH = .955*HFON52M3_P1 + .045*MALAH + .045*HO2 + .045*XC + .955*SumRO2';
k(:,i) = 5.95e-11;
Gstr{i,1} = 'HFON52M3'; Gstr{i,2} = 'OH';
fHFON52M3(i)=fHFON52M3(i)-1; fOH(i)=fOH(i)-1; fHFON52M3_P1(i)=fHFON52M3_P1(i)+.955; fMALAH(i)=fMALAH(i)+.045; fHO2(i)=fHO2(i)+.045; fXC(i)=fXC(i)+.045; fSumRO2(i)=fSumRO2(i)+.955;

i=i+1;
Rnames{i} = 'HFON52M3 + O3 = HFON52M3_P2 + OH + SumRO2';
k(:,i) = 1.63e-16;
Gstr{i,1} = 'HFON52M3'; Gstr{i,2} = 'O3';
fHFON52M3(i)=fHFON52M3(i)-1; fO3(i)=fO3(i)-1; fHFON52M3_P2(i)=fHFON52M3_P2(i)+1; fOH(i)=fOH(i)+1; fSumRO2(i)=fSumRO2(i)+1;

i=i+1;
Rnames{i} = 'HFON52M3 + NO3 = .984*HFON52M3_P3 + .016*MALAH + .016*HO2 + .016*XC + .016*XN + .984*SumRO2';
k(:,i) = 7.72e-14;
Gstr{i,1} = 'HFON52M3'; Gstr{i,2} = 'NO3';
fHFON52M3(i)=fHFON52M3(i)-1; fNO3(i)=fNO3(i)-1; fHFON52M3_P3(i)=fHFON52M3_P3(i)+.984; fMALAH(i)=fMALAH(i)+.016; fHO2(i)=fHO2(i)+.016; fXC(i)=fXC(i)+.016; fXN(i)=fXN(i)+.016; fSumRO2(i)=fSumRO2(i)+.984;

i=i+1;
Rnames{i} = 'HFON52M3_P1 + NO = .938*NO2 + .657*R2CO3 + .281*HO2 + .156*RCHO + .125*BACL + .062*RCNO3 + 1.657*XC + .657*SumRCO3';
k(:,i) = 9.13e-12;
Gstr{i,1} = 'HFON52M3_P1'; Gstr{i,2} = 'NO';
fHFON52M3_P1(i)=fHFON52M3_P1(i)-1; fNO(i)=fNO(i)-1; fNO2(i)=fNO2(i)+.938; fR2CO3(i)=fR2CO3(i)+.657; fHO2(i)=fHO2(i)+.281; fRCHO(i)=fRCHO(i)+.156; fBACL(i)=fBACL(i)+.125; fRCNO3(i)=fRCNO3(i)+.062; fXC(i)=fXC(i)+1.657; fSumRCO3(i)=fSumRCO3(i)+.657;

i=i+1;
Rnames{i} = 'HFON52M3_P1 + NO3 = NO2 + .701*R2CO3 + .299*HO2 + .166*RCHO + .133*BACL + 1.701*XC + .701*SumRCO3';
k(:,i) = 2.30e-12;
Gstr{i,1} = 'HFON52M3_P1'; Gstr{i,2} = 'NO3';
fHFON52M3_P1(i)=fHFON52M3_P1(i)-1; fNO3(i)=fNO3(i)-1; fNO2(i)=fNO2(i)+1; fR2CO3(i)=fR2CO3(i)+.701; fHO2(i)=fHO2(i)+.299; fRCHO(i)=fRCHO(i)+.166; fBACL(i)=fBACL(i)+.133; fXC(i)=fXC(i)+1.701; fSumRCO3(i)=fSumRCO3(i)+.701;

i=i+1;
Rnames{i} = 'HFON52M3_P1 + HO2 = ROOH';
k(:,i) = 1.78e-11;
Gstr{i,1} = 'HFON52M3_P1'; Gstr{i,2} = 'HO2';
fHFON52M3_P1(i)=fHFON52M3_P1(i)-1; fHO2(i)=fHO2(i)-1; fROOH(i)=fROOH(i)+1;

i=i+1;
Rnames{i} = 'HFON52M3_P1 + SumRO2 = SumRO2 + .457*OTH3 + .35*R2CO3 + .15*HO2 + .083*RCHO + .067*BACL + .043*KET2 + 1.264*XC + .35*SumRCO3';
k(:,i) = 4.42e-13;
Gstr{i,1} = 'HFON52M3_P1'; Gstr{i,2} = 'SumRO2';
fHFON52M3_P1(i)=fHFON52M3_P1(i)-1; fSumRO2(i)=fSumRO2(i)-1; fSumRO2(i)=fSumRO2(i)+1; fOTH3(i)=fOTH3(i)+.457; fR2CO3(i)=fR2CO3(i)+.35; fHO2(i)=fHO2(i)+.15; fRCHO(i)=fRCHO(i)+.083; fBACL(i)=fBACL(i)+.067; fKET2(i)=fKET2(i)+.043; fXC(i)=fXC(i)+1.264; fSumRCO3(i)=fSumRCO3(i)+.35;

i=i+1;
Rnames{i} = 'HFON52M3_P1 + SumRCO3 = SumRCO3 + .701*R2CO3 + .265*HO2 + .133*RCHO + .132*BACL + .034*KET2 + 1.633*XC + .701*SumRCO3';
k(:,i) = 1.37e-11;
Gstr{i,1} = 'HFON52M3_P1'; Gstr{i,2} = 'SumRCO3';
fHFON52M3_P1(i)=fHFON52M3_P1(i)-1; fSumRCO3(i)=fSumRCO3(i)-1; fSumRCO3(i)=fSumRCO3(i)+1; fR2CO3(i)=fR2CO3(i)+.701; fHO2(i)=fHO2(i)+.265; fRCHO(i)=fRCHO(i)+.133; fBACL(i)=fBACL(i)+.132; fKET2(i)=fKET2(i)+.034; fXC(i)=fXC(i)+1.633; fSumRCO3(i)=fSumRCO3(i)+.701;

i=i+1;
Rnames{i} = 'HFON52M3_P2 + NO = .938*R2CO3 + .938*HCHO + .938*NO2 + .062*RCNO3 + XC + .938*SumRCO3';
k(:,i) = 9.13e-12;
Gstr{i,1} = 'HFON52M3_P2'; Gstr{i,2} = 'NO';
fHFON52M3_P2(i)=fHFON52M3_P2(i)-1; fNO(i)=fNO(i)-1; fR2CO3(i)=fR2CO3(i)+.938; fHCHO(i)=fHCHO(i)+.938; fNO2(i)=fNO2(i)+.938; fRCNO3(i)=fRCNO3(i)+.062; fXC(i)=fXC(i)+1; fSumRCO3(i)=fSumRCO3(i)+.938;

i=i+1;
Rnames{i} = 'HFON52M3_P2 + NO3 = R2CO3 + HCHO + NO2 + XC + SumRCO3';
k(:,i) = 2.30e-12;
Gstr{i,1} = 'HFON52M3_P2'; Gstr{i,2} = 'NO3';
fHFON52M3_P2(i)=fHFON52M3_P2(i)-1; fNO3(i)=fNO3(i)-1; fR2CO3(i)=fR2CO3(i)+1; fHCHO(i)=fHCHO(i)+1; fNO2(i)=fNO2(i)+1; fXC(i)=fXC(i)+1; fSumRCO3(i)=fSumRCO3(i)+1;

i=i+1;
Rnames{i} = 'HFON52M3_P2 + HO2 = BACL + XC';
k(:,i) = 1.78e-11;
Gstr{i,1} = 'HFON52M3_P2'; Gstr{i,2} = 'HO2';
fHFON52M3_P2(i)=fHFON52M3_P2(i)-1; fHO2(i)=fHO2(i)-1; fBACL(i)=fBACL(i)+1; fXC(i)=fXC(i)+1;

i=i+1;
Rnames{i} = 'HFON52M3_P2 + SumRO2 = SumRO2 + .5*R2CO3 + .5*HCHO + .5*BACL + XC + .5*SumRCO3';
k(:,i) = 2.03e-13;
Gstr{i,1} = 'HFON52M3_P2'; Gstr{i,2} = 'SumRO2';
fHFON52M3_P2(i)=fHFON52M3_P2(i)-1; fSumRO2(i)=fSumRO2(i)-1; fSumRO2(i)=fSumRO2(i)+1; fR2CO3(i)=fR2CO3(i)+.5; fHCHO(i)=fHCHO(i)+.5; fBACL(i)=fBACL(i)+.5; fXC(i)=fXC(i)+1; fSumRCO3(i)=fSumRCO3(i)+.5;

i=i+1;
Rnames{i} = 'HFON52M3_P2 + SumRCO3 = SumRCO3 + .8*R2CO3 + .8*HCHO + .2*BACL + XC + .8*SumRCO3';
k(:,i) = 1.37e-11;
Gstr{i,1} = 'HFON52M3_P2'; Gstr{i,2} = 'SumRCO3';
fHFON52M3_P2(i)=fHFON52M3_P2(i)-1; fSumRCO3(i)=fSumRCO3(i)-1; fSumRCO3(i)=fSumRCO3(i)+1; fR2CO3(i)=fR2CO3(i)+.8; fHCHO(i)=fHCHO(i)+.8; fBACL(i)=fBACL(i)+.2; fXC(i)=fXC(i)+1; fSumRCO3(i)=fSumRCO3(i)+.8;

i=i+1;
Rnames{i} = 'HFON52M3_P3 + NO = RCNO3 + .938*HO2 + .938*NO2 + XC + .062*XN';
k(:,i) = 9.13e-12;
Gstr{i,1} = 'HFON52M3_P3'; Gstr{i,2} = 'NO';
fHFON52M3_P3(i)=fHFON52M3_P3(i)-1; fNO(i)=fNO(i)-1; fRCNO3(i)=fRCNO3(i)+1; fHO2(i)=fHO2(i)+.938; fNO2(i)=fNO2(i)+.938; fXC(i)=fXC(i)+1; fXN(i)=fXN(i)+.062;

i=i+1;
Rnames{i} = 'HFON52M3_P3 + NO3 = RCNO3 + HO2 + NO2 + XC';
k(:,i) = 2.30e-12;
Gstr{i,1} = 'HFON52M3_P3'; Gstr{i,2} = 'NO3';
fHFON52M3_P3(i)=fHFON52M3_P3(i)-1; fNO3(i)=fNO3(i)-1; fRCNO3(i)=fRCNO3(i)+1; fHO2(i)=fHO2(i)+1; fNO2(i)=fNO2(i)+1; fXC(i)=fXC(i)+1;

i=i+1;
Rnames{i} = 'HFON52M3_P3 + HO2 = RCNO3 + XC';
k(:,i) = 1.78e-11;
Gstr{i,1} = 'HFON52M3_P3'; Gstr{i,2} = 'HO2';
fHFON52M3_P3(i)=fHFON52M3_P3(i)-1; fHO2(i)=fHO2(i)-1; fRCNO3(i)=fRCNO3(i)+1; fXC(i)=fXC(i)+1;

i=i+1;
Rnames{i} = 'HFON52M3_P3 + SumRO2 = SumRO2 + RCNO3 + .5*HO2 + XC';
k(:,i) = 2.57e-12;
Gstr{i,1} = 'HFON52M3_P3'; Gstr{i,2} = 'SumRO2';
fHFON52M3_P3(i)=fHFON52M3_P3(i)-1; fSumRO2(i)=fSumRO2(i)-1; fSumRO2(i)=fSumRO2(i)+1; fRCNO3(i)=fRCNO3(i)+1; fHO2(i)=fHO2(i)+.5; fXC(i)=fXC(i)+1;

i=i+1;
Rnames{i} = 'HFON52M3_P3 + SumRCO3 = SumRCO3 + RCNO3 + .8*HO2 + XC';
k(:,i) = 1.37e-11;
Gstr{i,1} = 'HFON52M3_P3'; Gstr{i,2} = 'SumRCO3';
fHFON52M3_P3(i)=fHFON52M3_P3(i)-1; fSumRCO3(i)=fSumRCO3(i)-1; fSumRCO3(i)=fSumRCO3(i)+1; fRCNO3(i)=fRCNO3(i)+1; fHO2(i)=fHO2(i)+.8; fXC(i)=fXC(i)+1;

i=i+1;
Rnames{i} = 'HFON52M5 + OH = HFON52M5_P1 + SumRO2';
k(:,i) = 3.58e-11;
Gstr{i,1} = 'HFON52M5'; Gstr{i,2} = 'OH';
fHFON52M5(i)=fHFON52M5(i)-1; fOH(i)=fOH(i)-1; fHFON52M5_P1(i)=fHFON52M5_P1(i)+1; fSumRO2(i)=fSumRO2(i)+1;

i=i+1;
Rnames{i} = 'HFON52M5 + O3 = RCHO2 + 2*XC';
k(:,i) = 5.07e-17;
Gstr{i,1} = 'HFON52M5'; Gstr{i,2} = 'O3';
fHFON52M5(i)=fHFON52M5(i)-1; fO3(i)=fO3(i)-1; fRCHO2(i)=fRCHO2(i)+1; fXC(i)=fXC(i)+2;

i=i+1;
Rnames{i} = 'HFON52M5 + NO3 = .97*HFON52M5_P2 + .03*MALAH + .03*MEO2 + .03*XN + SumRO2';
k(:,i) = 3.96e-15;
Gstr{i,1} = 'HFON52M5'; Gstr{i,2} = 'NO3';
fHFON52M5(i)=fHFON52M5(i)-1; fNO3(i)=fNO3(i)-1; fHFON52M5_P2(i)=fHFON52M5_P2(i)+.97; fMALAH(i)=fMALAH(i)+.03; fMEO2(i)=fMEO2(i)+.03; fXN(i)=fXN(i)+.03; fSumRO2(i)=fSumRO2(i)+1;

i=i+1;
Rnames{i} = 'HFON52M5_P1 + NO = .938*NO2 + .932*HO2 + .626*MGLY + .306*RCHO + .062*RCNO3 + .006*R2CO3 + 1.632*XC + .006*SumRCO3';
k(:,i) = 9.13e-12;
Gstr{i,1} = 'HFON52M5_P1'; Gstr{i,2} = 'NO';
fHFON52M5_P1(i)=fHFON52M5_P1(i)-1; fNO(i)=fNO(i)-1; fNO2(i)=fNO2(i)+.938; fHO2(i)=fHO2(i)+.932; fMGLY(i)=fMGLY(i)+.626; fRCHO(i)=fRCHO(i)+.306; fRCNO3(i)=fRCNO3(i)+.062; fR2CO3(i)=fR2CO3(i)+.006; fXC(i)=fXC(i)+1.632; fSumRCO3(i)=fSumRCO3(i)+.006;

i=i+1;
Rnames{i} = 'HFON52M5_P1 + NO3 = NO2 + .993*HO2 + .667*MGLY + .326*RCHO + .007*R2CO3 + 1.674*XC + .007*SumRCO3';
k(:,i) = 2.30e-12;
Gstr{i,1} = 'HFON52M5_P1'; Gstr{i,2} = 'NO3';
fHFON52M5_P1(i)=fHFON52M5_P1(i)-1; fNO3(i)=fNO3(i)-1; fNO2(i)=fNO2(i)+1; fHO2(i)=fHO2(i)+.993; fMGLY(i)=fMGLY(i)+.667; fRCHO(i)=fRCHO(i)+.326; fR2CO3(i)=fR2CO3(i)+.007; fXC(i)=fXC(i)+1.674; fSumRCO3(i)=fSumRCO3(i)+.007;

i=i+1;
Rnames{i} = 'HFON52M5_P1 + HO2 = .935*ROOH + .065*BACL + .065*XC';
k(:,i) = 1.78e-11;
Gstr{i,1} = 'HFON52M5_P1'; Gstr{i,2} = 'HO2';
fHFON52M5_P1(i)=fHFON52M5_P1(i)-1; fHO2(i)=fHO2(i)-1; fROOH(i)=fROOH(i)+.935; fBACL(i)=fBACL(i)+.065; fXC(i)=fXC(i)+.065;

i=i+1;
Rnames{i} = 'HFON52M5_P1 + SumRO2 = SumRO2 + .497*HO2 + .333*MGLY + .25*OTH3 + .163*BACL + .163*RCHO + .087*KET2 + .003*R2CO3 + 1.166*XC + .003*SumRCO3';
k(:,i) = 2.57e-12;
Gstr{i,1} = 'HFON52M5_P1'; Gstr{i,2} = 'SumRO2';
fHFON52M5_P1(i)=fHFON52M5_P1(i)-1; fSumRO2(i)=fSumRO2(i)-1; fSumRO2(i)=fSumRO2(i)+1; fHO2(i)=fHO2(i)+.497; fMGLY(i)=fMGLY(i)+.333; fOTH3(i)=fOTH3(i)+.25; fBACL(i)=fBACL(i)+.163; fRCHO(i)=fRCHO(i)+.163; fKET2(i)=fKET2(i)+.087; fR2CO3(i)=fR2CO3(i)+.003; fXC(i)=fXC(i)+1.166; fSumRCO3(i)=fSumRCO3(i)+.003;

i=i+1;
Rnames{i} = 'HFON52M5_P1 + SumRCO3 = SumRCO3 + .795*HO2 + .534*MGLY + .261*RCHO + .131*BACL + .069*KET2 + .005*R2CO3 + 1.401*XC + .005*SumRCO3';
k(:,i) = 1.37e-11;
Gstr{i,1} = 'HFON52M5_P1'; Gstr{i,2} = 'SumRCO3';
fHFON52M5_P1(i)=fHFON52M5_P1(i)-1; fSumRCO3(i)=fSumRCO3(i)-1; fSumRCO3(i)=fSumRCO3(i)+1; fHO2(i)=fHO2(i)+.795; fMGLY(i)=fMGLY(i)+.534; fRCHO(i)=fRCHO(i)+.261; fBACL(i)=fBACL(i)+.131; fKET2(i)=fKET2(i)+.069; fR2CO3(i)=fR2CO3(i)+.005; fXC(i)=fXC(i)+1.401; fSumRCO3(i)=fSumRCO3(i)+.005;

i=i+1;
Rnames{i} = 'HFON52M5_P2 + NO = 1.407*NO2 + .531*RCNO3 + .469*HO2 + .374*RCHO + .095*MGLY + 1.095*XC + .062*XN';
k(:,i) = 9.13e-12;
Gstr{i,1} = 'HFON52M5_P2'; Gstr{i,2} = 'NO';
fHFON52M5_P2(i)=fHFON52M5_P2(i)-1; fNO(i)=fNO(i)-1; fNO2(i)=fNO2(i)+1.407; fRCNO3(i)=fRCNO3(i)+.531; fHO2(i)=fHO2(i)+.469; fRCHO(i)=fRCHO(i)+.374; fMGLY(i)=fMGLY(i)+.095; fXC(i)=fXC(i)+1.095; fXN(i)=fXN(i)+.062;

i=i+1;
Rnames{i} = 'HFON52M5_P2 + NO3 = 1.5*NO2 + .5*RCNO3 + .5*HO2 + .399*RCHO + .101*MGLY + 1.101*XC';
k(:,i) = 2.30e-12;
Gstr{i,1} = 'HFON52M5_P2'; Gstr{i,2} = 'NO3';
fHFON52M5_P2(i)=fHFON52M5_P2(i)-1; fNO3(i)=fNO3(i)-1; fNO2(i)=fNO2(i)+1.5; fRCNO3(i)=fRCNO3(i)+.5; fHO2(i)=fHO2(i)+.5; fRCHO(i)=fRCHO(i)+.399; fMGLY(i)=fMGLY(i)+.101; fXC(i)=fXC(i)+1.101;

i=i+1;
Rnames{i} = 'HFON52M5_P2 + HO2 = RCNO3 + XC';
k(:,i) = 1.78e-11;
Gstr{i,1} = 'HFON52M5_P2'; Gstr{i,2} = 'HO2';
fHFON52M5_P2(i)=fHFON52M5_P2(i)-1; fHO2(i)=fHO2(i)-1; fRCNO3(i)=fRCNO3(i)+1; fXC(i)=fXC(i)+1;

i=i+1;
Rnames{i} = 'HFON52M5_P2 + SumRO2 = SumRO2 + .75*RCNO3 + .25*NO2 + .25*HO2 + .199*RCHO + .051*MGLY + 1.051*XC';
k(:,i) = 2.57e-12;
Gstr{i,1} = 'HFON52M5_P2'; Gstr{i,2} = 'SumRO2';
fHFON52M5_P2(i)=fHFON52M5_P2(i)-1; fSumRO2(i)=fSumRO2(i)-1; fSumRO2(i)=fSumRO2(i)+1; fRCNO3(i)=fRCNO3(i)+.75; fNO2(i)=fNO2(i)+.25; fHO2(i)=fHO2(i)+.25; fRCHO(i)=fRCHO(i)+.199; fMGLY(i)=fMGLY(i)+.051; fXC(i)=fXC(i)+1.051;

i=i+1;
Rnames{i} = 'HFON52M5_P2 + SumRCO3 = SumRCO3 + .6*RCNO3 + .4*NO2 + .4*HO2 + .319*RCHO + .081*MGLY + 1.081*XC';
k(:,i) = 1.37e-11;
Gstr{i,1} = 'HFON52M5_P2'; Gstr{i,2} = 'SumRCO3';
fHFON52M5_P2(i)=fHFON52M5_P2(i)-1; fSumRCO3(i)=fSumRCO3(i)-1; fSumRCO3(i)=fSumRCO3(i)+1; fRCNO3(i)=fRCNO3(i)+.6; fNO2(i)=fNO2(i)+.4; fHO2(i)=fHO2(i)+.4; fRCHO(i)=fRCHO(i)+.319; fMGLY(i)=fMGLY(i)+.081; fXC(i)=fXC(i)+1.081;

i=i+1;
Rnames{i} = 'H3F2ONE + OH = H3F2ONE_P1 + SumRO2';
k(:,i) = 4.45e-11;
Gstr{i,1} = 'H3F2ONE'; Gstr{i,2} = 'OH';
fH3F2ONE(i)=fH3F2ONE(i)-1; fOH(i)=fOH(i)-1; fH3F2ONE_P1(i)=fH3F2ONE_P1(i)+1; fSumRO2(i)=fSumRO2(i)+1;

i=i+1;
Rnames{i} = 'H3F2ONE + O3 = 42.886*NROG + .547*CO2 + .453*H3F2ONE_R1 + .453*OH + .06*H3F2ONE_P2 + .06*HO2 + 1.461*XC + .06*SumRO2';
k(:,i) = 2.52e-17;
Gstr{i,1} = 'H3F2ONE'; Gstr{i,2} = 'O3';
fH3F2ONE(i)=fH3F2ONE(i)-1; fO3(i)=fO3(i)-1; fNROG(i)=fNROG(i)+42.886; fCO2(i)=fCO2(i)+.547; fH3F2ONE_R1(i)=fH3F2ONE_R1(i)+.453; fOH(i)=fOH(i)+.453; fH3F2ONE_P2(i)=fH3F2ONE_P2(i)+.06; fHO2(i)=fHO2(i)+.06; fXC(i)=fXC(i)+1.461; fSumRO2(i)=fSumRO2(i)+.06;

i=i+1;
Rnames{i} = 'H3F2ONE + NO3 = H3F2ONE_P3 + SumRO2';
k(:,i) = 2.82e-13;
Gstr{i,1} = 'H3F2ONE'; Gstr{i,2} = 'NO3';
fH3F2ONE(i)=fH3F2ONE(i)-1; fNO3(i)=fNO3(i)-1; fH3F2ONE_P3(i)=fH3F2ONE_P3(i)+1; fSumRO2(i)=fSumRO2(i)+1;

i=i+1;
Rnames{i} = 'H3F2ONE_P1 + NO = .974*RCHO + .974*HO2 + .974*NO2 + .026*RCNO3';
k(:,i) = 9.13e-12;
Gstr{i,1} = 'H3F2ONE_P1'; Gstr{i,2} = 'NO';
fH3F2ONE_P1(i)=fH3F2ONE_P1(i)-1; fNO(i)=fNO(i)-1; fRCHO(i)=fRCHO(i)+.974; fHO2(i)=fHO2(i)+.974; fNO2(i)=fNO2(i)+.974; fRCNO3(i)=fRCNO3(i)+.026;

i=i+1;
Rnames{i} = 'H3F2ONE_P1 + NO3 = RCHO + HO2 + NO2';
k(:,i) = 2.30e-12;
Gstr{i,1} = 'H3F2ONE_P1'; Gstr{i,2} = 'NO3';
fH3F2ONE_P1(i)=fH3F2ONE_P1(i)-1; fNO3(i)=fNO3(i)-1; fRCHO(i)=fRCHO(i)+1; fHO2(i)=fHO2(i)+1; fNO2(i)=fNO2(i)+1;

i=i+1;
Rnames{i} = 'H3F2ONE_P1 + HO2 = .912*ROOH + .088*OTH1 - 0.824*XC';
k(:,i) = 1.61e-11;
Gstr{i,1} = 'H3F2ONE_P1'; Gstr{i,2} = 'HO2';
fH3F2ONE_P1(i)=fH3F2ONE_P1(i)-1; fHO2(i)=fHO2(i)-1; fROOH(i)=fROOH(i)+.912; fOTH1(i)=fOTH1(i)+.088; fXC(i)=fXC(i)-0.824;

i=i+1;
Rnames{i} = 'H3F2ONE_P1 + SumRO2 = SumRO2 + .5*RCHO + .5*HO2 + .25*OTH3 + .221*OTH1 + .029*KET2 + .163*XC';
k(:,i) = 2.57e-12;
Gstr{i,1} = 'H3F2ONE_P1'; Gstr{i,2} = 'SumRO2';
fH3F2ONE_P1(i)=fH3F2ONE_P1(i)-1; fSumRO2(i)=fSumRO2(i)-1; fSumRO2(i)=fSumRO2(i)+1; fRCHO(i)=fRCHO(i)+.5; fHO2(i)=fHO2(i)+.5; fOTH3(i)=fOTH3(i)+.25; fOTH1(i)=fOTH1(i)+.221; fKET2(i)=fKET2(i)+.029; fXC(i)=fXC(i)+.163;

i=i+1;
Rnames{i} = 'H3F2ONE_P1 + SumRCO3 = SumRCO3 + .8*RCHO + .8*HO2 + .177*OTH1 + .023*KET2 + .131*XC';
k(:,i) = 1.37e-11;
Gstr{i,1} = 'H3F2ONE_P1'; Gstr{i,2} = 'SumRCO3';
fH3F2ONE_P1(i)=fH3F2ONE_P1(i)-1; fSumRCO3(i)=fSumRCO3(i)-1; fSumRCO3(i)=fSumRCO3(i)+1; fRCHO(i)=fRCHO(i)+.8; fHO2(i)=fHO2(i)+.8; fOTH1(i)=fOTH1(i)+.177; fKET2(i)=fKET2(i)+.023; fXC(i)=fXC(i)+.131;

i=i+1;
Rnames{i} = 'H3F2ONE_R1 = H3F2ONE_P4 + SumRO2';
k(:,i) = 1.64e+11.*(T./300).^0.00.*exp(-7848.229./T);
Gstr{i,1} = 'H3F2ONE_R1';
fH3F2ONE_R1(i)=fH3F2ONE_R1(i)-1; fH3F2ONE_P4(i)=fH3F2ONE_P4(i)+1; fSumRO2(i)=fSumRO2(i)+1;

i=i+1;
Rnames{i} = 'H3F2ONE_R1 + NO = .974*HO2 + .974*NO2 + .956*MGLY + .956*CO + .026*RCNO3 + .018*BACL';
k(:,i) = 9.13e-12;
Gstr{i,1} = 'H3F2ONE_R1'; Gstr{i,2} = 'NO';
fH3F2ONE_R1(i)=fH3F2ONE_R1(i)-1; fNO(i)=fNO(i)-1; fHO2(i)=fHO2(i)+.974; fNO2(i)=fNO2(i)+.974; fMGLY(i)=fMGLY(i)+.956; fCO(i)=fCO(i)+.956; fRCNO3(i)=fRCNO3(i)+.026; fBACL(i)=fBACL(i)+.018;

i=i+1;
Rnames{i} = 'H3F2ONE_P2 + NO = .99*NO2 + .742*MGLY + .742*HO2 + .249*R2CO3 + .249*HCHO + .01*RCNO3 - 0.262*XC + .249*SumRCO3';
k(:,i) = 9.13e-12;
Gstr{i,1} = 'H3F2ONE_P2'; Gstr{i,2} = 'NO';
fH3F2ONE_P2(i)=fH3F2ONE_P2(i)-1; fNO(i)=fNO(i)-1; fNO2(i)=fNO2(i)+.99; fMGLY(i)=fMGLY(i)+.742; fHO2(i)=fHO2(i)+.742; fR2CO3(i)=fR2CO3(i)+.249; fHCHO(i)=fHCHO(i)+.249; fRCNO3(i)=fRCNO3(i)+.01; fXC(i)=fXC(i)-0.262; fSumRCO3(i)=fSumRCO3(i)+.249;

i=i+1;
Rnames{i} = 'H3F2ONE_P2 + NO3 = NO2 + .749*MGLY + .749*HO2 + .251*R2CO3 + .251*HCHO - 0.251*XC + .251*SumRCO3';
k(:,i) = 2.30e-12;
Gstr{i,1} = 'H3F2ONE_P2'; Gstr{i,2} = 'NO3';
fH3F2ONE_P2(i)=fH3F2ONE_P2(i)-1; fNO3(i)=fNO3(i)-1; fNO2(i)=fNO2(i)+1; fMGLY(i)=fMGLY(i)+.749; fHO2(i)=fHO2(i)+.749; fR2CO3(i)=fR2CO3(i)+.251; fHCHO(i)=fHCHO(i)+.251; fXC(i)=fXC(i)-0.251; fSumRCO3(i)=fSumRCO3(i)+.251;

i=i+1;
Rnames{i} = 'H3F2ONE_P2 + HO2 = .9*ROOH + .1*MGLY - 1.8*XC';
k(:,i) = 1.44e-11;
Gstr{i,1} = 'H3F2ONE_P2'; Gstr{i,2} = 'HO2';
fH3F2ONE_P2(i)=fH3F2ONE_P2(i)-1; fHO2(i)=fHO2(i)-1; fROOH(i)=fROOH(i)+.9; fMGLY(i)=fMGLY(i)+.1; fXC(i)=fXC(i)-1.8;

i=i+1;
Rnames{i} = 'H3F2ONE_P2 + SumRO2 = SumRO2 + .624*MGLY + .374*HO2 + .25*OTH1 + .126*R2CO3 + .126*HCHO - 0.126*XC + .126*SumRCO3';
k(:,i) = 2.03e-13;
Gstr{i,1} = 'H3F2ONE_P2'; Gstr{i,2} = 'SumRO2';
fH3F2ONE_P2(i)=fH3F2ONE_P2(i)-1; fSumRO2(i)=fSumRO2(i)-1; fSumRO2(i)=fSumRO2(i)+1; fMGLY(i)=fMGLY(i)+.624; fHO2(i)=fHO2(i)+.374; fOTH1(i)=fOTH1(i)+.25; fR2CO3(i)=fR2CO3(i)+.126; fHCHO(i)=fHCHO(i)+.126; fXC(i)=fXC(i)-0.126; fSumRCO3(i)=fSumRCO3(i)+.126;

i=i+1;
Rnames{i} = 'H3F2ONE_P2 + SumRCO3 = SumRCO3 + .799*MGLY + .599*HO2 + .201*R2CO3 + .201*HCHO - 0.201*XC + .201*SumRCO3';
k(:,i) = 1.37e-11;
Gstr{i,1} = 'H3F2ONE_P2'; Gstr{i,2} = 'SumRCO3';
fH3F2ONE_P2(i)=fH3F2ONE_P2(i)-1; fSumRCO3(i)=fSumRCO3(i)-1; fSumRCO3(i)=fSumRCO3(i)+1; fMGLY(i)=fMGLY(i)+.799; fHO2(i)=fHO2(i)+.599; fR2CO3(i)=fR2CO3(i)+.201; fHCHO(i)=fHCHO(i)+.201; fXC(i)=fXC(i)-0.201; fSumRCO3(i)=fSumRCO3(i)+.201;

i=i+1;
Rnames{i} = 'H3F2ONE_P3 + NO = 1.949*NO2 + .974*RCHO + .026*RCNO3 + .025*XN';
k(:,i) = 9.13e-12;
Gstr{i,1} = 'H3F2ONE_P3'; Gstr{i,2} = 'NO';
fH3F2ONE_P3(i)=fH3F2ONE_P3(i)-1; fNO(i)=fNO(i)-1; fNO2(i)=fNO2(i)+1.949; fRCHO(i)=fRCHO(i)+.974; fRCNO3(i)=fRCNO3(i)+.026; fXN(i)=fXN(i)+.025;

i=i+1;
Rnames{i} = 'H3F2ONE_P3 + NO3 = 2*NO2 + RCHO';
k(:,i) = 2.30e-12;
Gstr{i,1} = 'H3F2ONE_P3'; Gstr{i,2} = 'NO3';
fH3F2ONE_P3(i)=fH3F2ONE_P3(i)-1; fNO3(i)=fNO3(i)-1; fNO2(i)=fNO2(i)+2; fRCHO(i)=fRCHO(i)+1;

i=i+1;
Rnames{i} = 'H3F2ONE_P3 + HO2 = RCNO3';
k(:,i) = 1.61e-11;
Gstr{i,1} = 'H3F2ONE_P3'; Gstr{i,2} = 'HO2';
fH3F2ONE_P3(i)=fH3F2ONE_P3(i)-1; fHO2(i)=fHO2(i)-1; fRCNO3(i)=fRCNO3(i)+1;

i=i+1;
Rnames{i} = 'H3F2ONE_P3 + SumRO2 = SumRO2 + .5*RCHO + .5*NO2 + .5*RCNO3';
k(:,i) = 2.57e-12;
Gstr{i,1} = 'H3F2ONE_P3'; Gstr{i,2} = 'SumRO2';
fH3F2ONE_P3(i)=fH3F2ONE_P3(i)-1; fSumRO2(i)=fSumRO2(i)-1; fSumRO2(i)=fSumRO2(i)+1; fRCHO(i)=fRCHO(i)+.5; fNO2(i)=fNO2(i)+.5; fRCNO3(i)=fRCNO3(i)+.5;

i=i+1;
Rnames{i} = 'H3F2ONE_P3 + SumRCO3 = SumRCO3 + .8*RCHO + .8*NO2 + .2*RCNO3';
k(:,i) = 1.37e-11;
Gstr{i,1} = 'H3F2ONE_P3'; Gstr{i,2} = 'SumRCO3';
fH3F2ONE_P3(i)=fH3F2ONE_P3(i)-1; fSumRCO3(i)=fSumRCO3(i)-1; fSumRCO3(i)=fSumRCO3(i)+1; fRCHO(i)=fRCHO(i)+.8; fNO2(i)=fNO2(i)+.8; fRCNO3(i)=fRCNO3(i)+.2;

i=i+1;
Rnames{i} = 'H3F2ONE_P4 + NO = .974*NO2 + .942*BACL + .942*HO2 + .032*MGLY + .026*RCNO3 + .016*R2CO3 + .016*CO2 + .016*OH - 0.032*XC + .016*SumRCO3';
k(:,i) = 9.13e-12;
Gstr{i,1} = 'H3F2ONE_P4'; Gstr{i,2} = 'NO';
fH3F2ONE_P4(i)=fH3F2ONE_P4(i)-1; fNO(i)=fNO(i)-1; fNO2(i)=fNO2(i)+.974; fBACL(i)=fBACL(i)+.942; fHO2(i)=fHO2(i)+.942; fMGLY(i)=fMGLY(i)+.032; fRCNO3(i)=fRCNO3(i)+.026; fR2CO3(i)=fR2CO3(i)+.016; fCO2(i)=fCO2(i)+.016; fOH(i)=fOH(i)+.016; fXC(i)=fXC(i)-0.032; fSumRCO3(i)=fSumRCO3(i)+.016;

i=i+1;
Rnames{i} = 'H3F2ONE_P4 + NO3 = NO2 + .967*BACL + .967*HO2 + .033*MGLY + .017*R2CO3 + .017*CO2 + .017*OH - 0.035*XC + .017*SumRCO3';
k(:,i) = 2.30e-12;
Gstr{i,1} = 'H3F2ONE_P4'; Gstr{i,2} = 'NO3';
fH3F2ONE_P4(i)=fH3F2ONE_P4(i)-1; fNO3(i)=fNO3(i)-1; fNO2(i)=fNO2(i)+1; fBACL(i)=fBACL(i)+.967; fHO2(i)=fHO2(i)+.967; fMGLY(i)=fMGLY(i)+.033; fR2CO3(i)=fR2CO3(i)+.017; fCO2(i)=fCO2(i)+.017; fOH(i)=fOH(i)+.017; fXC(i)=fXC(i)-0.035; fSumRCO3(i)=fSumRCO3(i)+.017;

i=i+1;
Rnames{i} = 'H3F2ONE_P4 + HO2 = .9*ROOH + .1*BACL - 0.9*XC';
k(:,i) = 1.61e-11;
Gstr{i,1} = 'H3F2ONE_P4'; Gstr{i,2} = 'HO2';
fH3F2ONE_P4(i)=fH3F2ONE_P4(i)-1; fHO2(i)=fHO2(i)-1; fROOH(i)=fROOH(i)+.9; fBACL(i)=fBACL(i)+.1; fXC(i)=fXC(i)-0.9;

i=i+1;
Rnames{i} = 'H3F2ONE_P4 + SumRO2 = SumRO2 + .733*BACL + .483*HO2 + .25*ROOH + .017*MGLY + .008*R2CO3 + .008*CO2 + .008*OH - 0.265*XC + .008*SumRCO3';
k(:,i) = 2.57e-12;
Gstr{i,1} = 'H3F2ONE_P4'; Gstr{i,2} = 'SumRO2';
fH3F2ONE_P4(i)=fH3F2ONE_P4(i)-1; fSumRO2(i)=fSumRO2(i)-1; fSumRO2(i)=fSumRO2(i)+1; fBACL(i)=fBACL(i)+.733; fHO2(i)=fHO2(i)+.483; fROOH(i)=fROOH(i)+.25; fMGLY(i)=fMGLY(i)+.017; fR2CO3(i)=fR2CO3(i)+.008; fCO2(i)=fCO2(i)+.008; fOH(i)=fOH(i)+.008; fXC(i)=fXC(i)-0.265; fSumRCO3(i)=fSumRCO3(i)+.008;

i=i+1;
Rnames{i} = 'H3F2ONE_P4 + SumRCO3 = SumRCO3 + .973*BACL + .773*HO2 + .027*MGLY + .013*R2CO3 + .013*CO2 + .013*OH - 0.025*XC + .013*SumRCO3';
k(:,i) = 1.37e-11;
Gstr{i,1} = 'H3F2ONE_P4'; Gstr{i,2} = 'SumRCO3';
fH3F2ONE_P4(i)=fH3F2ONE_P4(i)-1; fSumRCO3(i)=fSumRCO3(i)-1; fSumRCO3(i)=fSumRCO3(i)+1; fBACL(i)=fBACL(i)+.973; fHO2(i)=fHO2(i)+.773; fMGLY(i)=fMGLY(i)+.027; fR2CO3(i)=fR2CO3(i)+.013; fCO2(i)=fCO2(i)+.013; fOH(i)=fOH(i)+.013; fXC(i)=fXC(i)-0.025; fSumRCO3(i)=fSumRCO3(i)+.013;

i=i+1;
Rnames{i} = 'H3F2ONM5 + OH = H3F2ONM5_P1 + SumRO2';
k(:,i) = 7.53e-11;
Gstr{i,1} = 'H3F2ONM5'; Gstr{i,2} = 'OH';
fH3F2ONM5(i)=fH3F2ONM5(i)-1; fOH(i)=fOH(i)-1; fH3F2ONM5_P1(i)=fH3F2ONM5_P1(i)+1; fSumRO2(i)=fSumRO2(i)+1;

i=i+1;
Rnames{i} = 'H3F2ONM5 + O3 = 49.717*NROG + .547*CO2 + .453*H3F2ONM5_R1 + .453*OH + .06*H3F2ONM5_P2 + .06*HO2 + 1.948*XC + .06*SumRO2';
k(:,i) = 8.11e-17;
Gstr{i,1} = 'H3F2ONM5'; Gstr{i,2} = 'O3';
fH3F2ONM5(i)=fH3F2ONM5(i)-1; fO3(i)=fO3(i)-1; fNROG(i)=fNROG(i)+49.717; fCO2(i)=fCO2(i)+.547; fH3F2ONM5_R1(i)=fH3F2ONM5_R1(i)+.453; fOH(i)=fOH(i)+.453; fH3F2ONM5_P2(i)=fH3F2ONM5_P2(i)+.06; fHO2(i)=fHO2(i)+.06; fXC(i)=fXC(i)+1.948; fSumRO2(i)=fSumRO2(i)+.06;

i=i+1;
Rnames{i} = 'H3F2ONM5 + NO3 = H3F2ONM5_P3 + SumRO2';
k(:,i) = 9.85e-12;
Gstr{i,1} = 'H3F2ONM5'; Gstr{i,2} = 'NO3';
fH3F2ONM5(i)=fH3F2ONM5(i)-1; fNO3(i)=fNO3(i)-1; fH3F2ONM5_P3(i)=fH3F2ONM5_P3(i)+1; fSumRO2(i)=fSumRO2(i)+1;

i=i+1;
Rnames{i} = 'H3F2ONM5_P1 + NO = .938*RCHO + .938*HO2 + .938*NO2 + .062*RCNO3 + XC';
k(:,i) = 9.13e-12;
Gstr{i,1} = 'H3F2ONM5_P1'; Gstr{i,2} = 'NO';
fH3F2ONM5_P1(i)=fH3F2ONM5_P1(i)-1; fNO(i)=fNO(i)-1; fRCHO(i)=fRCHO(i)+.938; fHO2(i)=fHO2(i)+.938; fNO2(i)=fNO2(i)+.938; fRCNO3(i)=fRCNO3(i)+.062; fXC(i)=fXC(i)+1;

i=i+1;
Rnames{i} = 'H3F2ONM5_P1 + NO3 = RCHO + HO2 + NO2 + XC';
k(:,i) = 2.30e-12;
Gstr{i,1} = 'H3F2ONM5_P1'; Gstr{i,2} = 'NO3';
fH3F2ONM5_P1(i)=fH3F2ONM5_P1(i)-1; fNO3(i)=fNO3(i)-1; fRCHO(i)=fRCHO(i)+1; fHO2(i)=fHO2(i)+1; fNO2(i)=fNO2(i)+1; fXC(i)=fXC(i)+1;

i=i+1;
Rnames{i} = 'H3F2ONM5_P1 + HO2 = ROOH';
k(:,i) = 1.78e-11;
Gstr{i,1} = 'H3F2ONM5_P1'; Gstr{i,2} = 'HO2';
fH3F2ONM5_P1(i)=fH3F2ONM5_P1(i)-1; fHO2(i)=fHO2(i)-1; fROOH(i)=fROOH(i)+1;

i=i+1;
Rnames{i} = 'H3F2ONM5_P1 + SumRO2 = SumRO2 + .5*RCHO + .5*HO2 + .488*OTH3 + .012*KET2 + .976*XC';
k(:,i) = 1.26e-13;
Gstr{i,1} = 'H3F2ONM5_P1'; Gstr{i,2} = 'SumRO2';
fH3F2ONM5_P1(i)=fH3F2ONM5_P1(i)-1; fSumRO2(i)=fSumRO2(i)-1; fSumRO2(i)=fSumRO2(i)+1; fRCHO(i)=fRCHO(i)+.5; fHO2(i)=fHO2(i)+.5; fOTH3(i)=fOTH3(i)+.488; fKET2(i)=fKET2(i)+.012; fXC(i)=fXC(i)+.976;

i=i+1;
Rnames{i} = 'H3F2ONM5_P1 + SumRCO3 = SumRCO3 + .99*RCHO + .99*HO2 + .01*KET2 + .98*XC';
k(:,i) = 1.37e-11;
Gstr{i,1} = 'H3F2ONM5_P1'; Gstr{i,2} = 'SumRCO3';
fH3F2ONM5_P1(i)=fH3F2ONM5_P1(i)-1; fSumRCO3(i)=fSumRCO3(i)-1; fSumRCO3(i)=fSumRCO3(i)+1; fRCHO(i)=fRCHO(i)+.99; fHO2(i)=fHO2(i)+.99; fKET2(i)=fKET2(i)+.01; fXC(i)=fXC(i)+.98;

i=i+1;
Rnames{i} = 'H3F2ONM5_R1 = H3F2ONM5_P4 + SumRO2';
k(:,i) = 1.64e+11.*(T./300).^0.00.*exp(-7848.229./T);
Gstr{i,1} = 'H3F2ONM5_R1';
fH3F2ONM5_R1(i)=fH3F2ONM5_R1(i)-1; fH3F2ONM5_P4(i)=fH3F2ONM5_P4(i)+1; fSumRO2(i)=fSumRO2(i)+1;

i=i+1;
Rnames{i} = 'H3F2ONM5_R1 + NO = .938*HO2 + .938*NO2 + .921*MGLY + .921*CO + .062*RCNO3 + .018*BACL + .996*XC';
k(:,i) = 9.13e-12;
Gstr{i,1} = 'H3F2ONM5_R1'; Gstr{i,2} = 'NO';
fH3F2ONM5_R1(i)=fH3F2ONM5_R1(i)-1; fNO(i)=fNO(i)-1; fHO2(i)=fHO2(i)+.938; fNO2(i)=fNO2(i)+.938; fMGLY(i)=fMGLY(i)+.921; fCO(i)=fCO(i)+.921; fRCNO3(i)=fRCNO3(i)+.062; fBACL(i)=fBACL(i)+.018; fXC(i)=fXC(i)+.996;

i=i+1;
Rnames{i} = 'H3F2ONM5_P2 + NO = .974*NO2 + .73*MGLY + .73*HO2 + .245*R2CO3 + .245*HCHO + .026*RCNO3 + .726*XC + .245*SumRCO3';
k(:,i) = 9.13e-12;
Gstr{i,1} = 'H3F2ONM5_P2'; Gstr{i,2} = 'NO';
fH3F2ONM5_P2(i)=fH3F2ONM5_P2(i)-1; fNO(i)=fNO(i)-1; fNO2(i)=fNO2(i)+.974; fMGLY(i)=fMGLY(i)+.73; fHO2(i)=fHO2(i)+.73; fR2CO3(i)=fR2CO3(i)+.245; fHCHO(i)=fHCHO(i)+.245; fRCNO3(i)=fRCNO3(i)+.026; fXC(i)=fXC(i)+.726; fSumRCO3(i)=fSumRCO3(i)+.245;

i=i+1;
Rnames{i} = 'H3F2ONM5_P2 + NO3 = NO2 + .749*MGLY + .749*HO2 + .251*R2CO3 + .251*HCHO + .749*XC + .251*SumRCO3';
k(:,i) = 2.30e-12;
Gstr{i,1} = 'H3F2ONM5_P2'; Gstr{i,2} = 'NO3';
fH3F2ONM5_P2(i)=fH3F2ONM5_P2(i)-1; fNO3(i)=fNO3(i)-1; fNO2(i)=fNO2(i)+1; fMGLY(i)=fMGLY(i)+.749; fHO2(i)=fHO2(i)+.749; fR2CO3(i)=fR2CO3(i)+.251; fHCHO(i)=fHCHO(i)+.251; fXC(i)=fXC(i)+.749; fSumRCO3(i)=fSumRCO3(i)+.251;

i=i+1;
Rnames{i} = 'H3F2ONM5_P2 + HO2 = .9*ROOH + .1*MGLY - 0.8*XC';
k(:,i) = 1.61e-11;
Gstr{i,1} = 'H3F2ONM5_P2'; Gstr{i,2} = 'HO2';
fH3F2ONM5_P2(i)=fH3F2ONM5_P2(i)-1; fHO2(i)=fHO2(i)-1; fROOH(i)=fROOH(i)+.9; fMGLY(i)=fMGLY(i)+.1; fXC(i)=fXC(i)-0.8;

i=i+1;
Rnames{i} = 'H3F2ONM5_P2 + SumRO2 = SumRO2 + .624*MGLY + .374*HO2 + .25*OTH1 + .126*R2CO3 + .126*HCHO + .874*XC + .126*SumRCO3';
k(:,i) = 2.03e-13;
Gstr{i,1} = 'H3F2ONM5_P2'; Gstr{i,2} = 'SumRO2';
fH3F2ONM5_P2(i)=fH3F2ONM5_P2(i)-1; fSumRO2(i)=fSumRO2(i)-1; fSumRO2(i)=fSumRO2(i)+1; fMGLY(i)=fMGLY(i)+.624; fHO2(i)=fHO2(i)+.374; fOTH1(i)=fOTH1(i)+.25; fR2CO3(i)=fR2CO3(i)+.126; fHCHO(i)=fHCHO(i)+.126; fXC(i)=fXC(i)+.874; fSumRCO3(i)=fSumRCO3(i)+.126;

i=i+1;
Rnames{i} = 'H3F2ONM5_P2 + SumRCO3 = SumRCO3 + .799*MGLY + .599*HO2 + .201*R2CO3 + .201*HCHO + .799*XC + .201*SumRCO3';
k(:,i) = 1.37e-11;
Gstr{i,1} = 'H3F2ONM5_P2'; Gstr{i,2} = 'SumRCO3';
fH3F2ONM5_P2(i)=fH3F2ONM5_P2(i)-1; fSumRCO3(i)=fSumRCO3(i)-1; fSumRCO3(i)=fSumRCO3(i)+1; fMGLY(i)=fMGLY(i)+.799; fHO2(i)=fHO2(i)+.599; fR2CO3(i)=fR2CO3(i)+.201; fHCHO(i)=fHCHO(i)+.201; fXC(i)=fXC(i)+.799; fSumRCO3(i)=fSumRCO3(i)+.201;

i=i+1;
Rnames{i} = 'H3F2ONM5_P3 + NO = 1.876*NO2 + .938*RCHO + .062*RCNO3 + XC + .062*XN';
k(:,i) = 9.13e-12;
Gstr{i,1} = 'H3F2ONM5_P3'; Gstr{i,2} = 'NO';
fH3F2ONM5_P3(i)=fH3F2ONM5_P3(i)-1; fNO(i)=fNO(i)-1; fNO2(i)=fNO2(i)+1.876; fRCHO(i)=fRCHO(i)+.938; fRCNO3(i)=fRCNO3(i)+.062; fXC(i)=fXC(i)+1; fXN(i)=fXN(i)+.062;

i=i+1;
Rnames{i} = 'H3F2ONM5_P3 + NO3 = 2*NO2 + RCHO + XC';
k(:,i) = 2.30e-12;
Gstr{i,1} = 'H3F2ONM5_P3'; Gstr{i,2} = 'NO3';
fH3F2ONM5_P3(i)=fH3F2ONM5_P3(i)-1; fNO3(i)=fNO3(i)-1; fNO2(i)=fNO2(i)+2; fRCHO(i)=fRCHO(i)+1; fXC(i)=fXC(i)+1;

i=i+1;
Rnames{i} = 'H3F2ONM5_P3 + HO2 = RCNO3 + XC';
k(:,i) = 1.78e-11;
Gstr{i,1} = 'H3F2ONM5_P3'; Gstr{i,2} = 'HO2';
fH3F2ONM5_P3(i)=fH3F2ONM5_P3(i)-1; fHO2(i)=fHO2(i)-1; fRCNO3(i)=fRCNO3(i)+1; fXC(i)=fXC(i)+1;

i=i+1;
Rnames{i} = 'H3F2ONM5_P3 + SumRO2 = SumRO2 + .5*RCHO + .5*NO2 + .5*RCNO3 + XC';
k(:,i) = 2.57e-12;
Gstr{i,1} = 'H3F2ONM5_P3'; Gstr{i,2} = 'SumRO2';
fH3F2ONM5_P3(i)=fH3F2ONM5_P3(i)-1; fSumRO2(i)=fSumRO2(i)-1; fSumRO2(i)=fSumRO2(i)+1; fRCHO(i)=fRCHO(i)+.5; fNO2(i)=fNO2(i)+.5; fRCNO3(i)=fRCNO3(i)+.5; fXC(i)=fXC(i)+1;

i=i+1;
Rnames{i} = 'H3F2ONM5_P3 + SumRCO3 = SumRCO3 + .8*RCHO + .8*NO2 + .2*RCNO3 + XC';
k(:,i) = 1.37e-11;
Gstr{i,1} = 'H3F2ONM5_P3'; Gstr{i,2} = 'SumRCO3';
fH3F2ONM5_P3(i)=fH3F2ONM5_P3(i)-1; fSumRCO3(i)=fSumRCO3(i)-1; fSumRCO3(i)=fSumRCO3(i)+1; fRCHO(i)=fRCHO(i)+.8; fNO2(i)=fNO2(i)+.8; fRCNO3(i)=fRCNO3(i)+.2; fXC(i)=fXC(i)+1;

i=i+1;
Rnames{i} = 'H3F2ONM5_P4 + NO = .938*NO2 + .907*BACL + .907*HO2 + .062*RCNO3 + .031*MGLY + .016*R2CO3 + .016*CO2 + .016*OH + .967*XC + .016*SumRCO3';
k(:,i) = 9.13e-12;
Gstr{i,1} = 'H3F2ONM5_P4'; Gstr{i,2} = 'NO';
fH3F2ONM5_P4(i)=fH3F2ONM5_P4(i)-1; fNO(i)=fNO(i)-1; fNO2(i)=fNO2(i)+.938; fBACL(i)=fBACL(i)+.907; fHO2(i)=fHO2(i)+.907; fRCNO3(i)=fRCNO3(i)+.062; fMGLY(i)=fMGLY(i)+.031; fR2CO3(i)=fR2CO3(i)+.016; fCO2(i)=fCO2(i)+.016; fOH(i)=fOH(i)+.016; fXC(i)=fXC(i)+.967; fSumRCO3(i)=fSumRCO3(i)+.016;

i=i+1;
Rnames{i} = 'H3F2ONM5_P4 + NO3 = NO2 + .967*BACL + .967*HO2 + .033*MGLY + .017*R2CO3 + .017*CO2 + .017*OH + .965*XC + .017*SumRCO3';
k(:,i) = 2.30e-12;
Gstr{i,1} = 'H3F2ONM5_P4'; Gstr{i,2} = 'NO3';
fH3F2ONM5_P4(i)=fH3F2ONM5_P4(i)-1; fNO3(i)=fNO3(i)-1; fNO2(i)=fNO2(i)+1; fBACL(i)=fBACL(i)+.967; fHO2(i)=fHO2(i)+.967; fMGLY(i)=fMGLY(i)+.033; fR2CO3(i)=fR2CO3(i)+.017; fCO2(i)=fCO2(i)+.017; fOH(i)=fOH(i)+.017; fXC(i)=fXC(i)+.965; fSumRCO3(i)=fSumRCO3(i)+.017;

i=i+1;
Rnames{i} = 'H3F2ONM5_P4 + HO2 = .9*ROOH + .1*BACL + .1*XC';
k(:,i) = 1.78e-11;
Gstr{i,1} = 'H3F2ONM5_P4'; Gstr{i,2} = 'HO2';
fH3F2ONM5_P4(i)=fH3F2ONM5_P4(i)-1; fHO2(i)=fHO2(i)-1; fROOH(i)=fROOH(i)+.9; fBACL(i)=fBACL(i)+.1; fXC(i)=fXC(i)+.1;

i=i+1;
Rnames{i} = 'H3F2ONM5_P4 + SumRO2 = SumRO2 + .733*BACL + .483*HO2 + .25*ROOH + .017*MGLY + .008*R2CO3 + .008*CO2 + .008*OH + .735*XC + .008*SumRCO3';
k(:,i) = 2.57e-12;
Gstr{i,1} = 'H3F2ONM5_P4'; Gstr{i,2} = 'SumRO2';
fH3F2ONM5_P4(i)=fH3F2ONM5_P4(i)-1; fSumRO2(i)=fSumRO2(i)-1; fSumRO2(i)=fSumRO2(i)+1; fBACL(i)=fBACL(i)+.733; fHO2(i)=fHO2(i)+.483; fROOH(i)=fROOH(i)+.25; fMGLY(i)=fMGLY(i)+.017; fR2CO3(i)=fR2CO3(i)+.008; fCO2(i)=fCO2(i)+.008; fOH(i)=fOH(i)+.008; fXC(i)=fXC(i)+.735; fSumRCO3(i)=fSumRCO3(i)+.008;

i=i+1;
Rnames{i} = 'H3F2ONM5_P4 + SumRCO3 = SumRCO3 + .973*BACL + .773*HO2 + .027*MGLY + .013*R2CO3 + .013*CO2 + .013*OH + .975*XC + .013*SumRCO3';
k(:,i) = 1.37e-11;
Gstr{i,1} = 'H3F2ONM5_P4'; Gstr{i,2} = 'SumRCO3';
fH3F2ONM5_P4(i)=fH3F2ONM5_P4(i)-1; fSumRCO3(i)=fSumRCO3(i)-1; fSumRCO3(i)=fSumRCO3(i)+1; fBACL(i)=fBACL(i)+.973; fHO2(i)=fHO2(i)+.773; fMGLY(i)=fMGLY(i)+.027; fR2CO3(i)=fR2CO3(i)+.013; fCO2(i)=fCO2(i)+.013; fOH(i)=fOH(i)+.013; fXC(i)=fXC(i)+.975; fSumRCO3(i)=fSumRCO3(i)+.013;

i=i+1;
Rnames{i} = 'H3F2ONM4 + OH = H3F2ONM4_P1 + SumRO2';
k(:,i) = 3.75e-11;
Gstr{i,1} = 'H3F2ONM4'; Gstr{i,2} = 'OH';
fH3F2ONM4(i)=fH3F2ONM4(i)-1; fOH(i)=fOH(i)-1; fH3F2ONM4_P1(i)=fH3F2ONM4_P1(i)+1; fSumRO2(i)=fSumRO2(i)+1;

i=i+1;
Rnames{i} = 'H3F2ONM4 + O3 = H3F2ONM4_P2 + OH + SumRO2';
k(:,i) = 8.11e-17;
Gstr{i,1} = 'H3F2ONM4'; Gstr{i,2} = 'O3';
fH3F2ONM4(i)=fH3F2ONM4(i)-1; fO3(i)=fO3(i)-1; fH3F2ONM4_P2(i)=fH3F2ONM4_P2(i)+1; fOH(i)=fOH(i)+1; fSumRO2(i)=fSumRO2(i)+1;

i=i+1;
Rnames{i} = 'H3F2ONM4 + NO3 = H3F2ONM4_P3 + SumRO2';
k(:,i) = 1.30e-12;
Gstr{i,1} = 'H3F2ONM4'; Gstr{i,2} = 'NO3';
fH3F2ONM4(i)=fH3F2ONM4(i)-1; fNO3(i)=fNO3(i)-1; fH3F2ONM4_P3(i)=fH3F2ONM4_P3(i)+1; fSumRO2(i)=fSumRO2(i)+1;

i=i+1;
Rnames{i} = 'H3F2ONM4_P1 + NO = .938*KET2 + .938*HO2 + .938*NO2 + .062*RCNO3 - 0.876*XC';
k(:,i) = 9.13e-12;
Gstr{i,1} = 'H3F2ONM4_P1'; Gstr{i,2} = 'NO';
fH3F2ONM4_P1(i)=fH3F2ONM4_P1(i)-1; fNO(i)=fNO(i)-1; fKET2(i)=fKET2(i)+.938; fHO2(i)=fHO2(i)+.938; fNO2(i)=fNO2(i)+.938; fRCNO3(i)=fRCNO3(i)+.062; fXC(i)=fXC(i)-0.876;

i=i+1;
Rnames{i} = 'H3F2ONM4_P1 + NO3 = KET2 + HO2 + NO2 - 1*XC';
k(:,i) = 2.30e-12;
Gstr{i,1} = 'H3F2ONM4_P1'; Gstr{i,2} = 'NO3';
fH3F2ONM4_P1(i)=fH3F2ONM4_P1(i)-1; fNO3(i)=fNO3(i)-1; fKET2(i)=fKET2(i)+1; fHO2(i)=fHO2(i)+1; fNO2(i)=fNO2(i)+1; fXC(i)=fXC(i)-1;

i=i+1;
Rnames{i} = 'H3F2ONM4_P1 + HO2 = .925*ROOH + .075*OTH1 + .15*XC';
k(:,i) = 1.78e-11;
Gstr{i,1} = 'H3F2ONM4_P1'; Gstr{i,2} = 'HO2';
fH3F2ONM4_P1(i)=fH3F2ONM4_P1(i)-1; fHO2(i)=fHO2(i)-1; fROOH(i)=fROOH(i)+.925; fOTH1(i)=fOTH1(i)+.075; fXC(i)=fXC(i)+.15;

i=i+1;
Rnames{i} = 'H3F2ONM4_P1 + SumRO2 = SumRO2 + .5*KET2 + .5*HO2 + .313*OTH3 + .187*OTH1 + .187*XC';
k(:,i) = 1.92e-12;
Gstr{i,1} = 'H3F2ONM4_P1'; Gstr{i,2} = 'SumRO2';
fH3F2ONM4_P1(i)=fH3F2ONM4_P1(i)-1; fSumRO2(i)=fSumRO2(i)-1; fSumRO2(i)=fSumRO2(i)+1; fKET2(i)=fKET2(i)+.5; fHO2(i)=fHO2(i)+.5; fOTH3(i)=fOTH3(i)+.313; fOTH1(i)=fOTH1(i)+.187; fXC(i)=fXC(i)+.187;

i=i+1;
Rnames{i} = 'H3F2ONM4_P1 + SumRCO3 = SumRCO3 + .85*KET2 + .85*HO2 + .15*OTH1 - 0.55*XC';
k(:,i) = 1.37e-11;
Gstr{i,1} = 'H3F2ONM4_P1'; Gstr{i,2} = 'SumRCO3';
fH3F2ONM4_P1(i)=fH3F2ONM4_P1(i)-1; fSumRCO3(i)=fSumRCO3(i)-1; fSumRCO3(i)=fSumRCO3(i)+1; fKET2(i)=fKET2(i)+.85; fHO2(i)=fHO2(i)+.85; fOTH1(i)=fOTH1(i)+.15; fXC(i)=fXC(i)-0.55;

i=i+1;
Rnames{i} = 'H3F2ONM4_P2 + NO = .938*NO2 + .469*R2CO3 + .469*HCHO + .469*MECO3 + .469*MGLY + .062*RCNO3 + .531*XC + .938*SumRCO3';
k(:,i) = 9.13e-12;
Gstr{i,1} = 'H3F2ONM4_P2'; Gstr{i,2} = 'NO';
fH3F2ONM4_P2(i)=fH3F2ONM4_P2(i)-1; fNO(i)=fNO(i)-1; fNO2(i)=fNO2(i)+.938; fR2CO3(i)=fR2CO3(i)+.469; fHCHO(i)=fHCHO(i)+.469; fMECO3(i)=fMECO3(i)+.469; fMGLY(i)=fMGLY(i)+.469; fRCNO3(i)=fRCNO3(i)+.062; fXC(i)=fXC(i)+.531; fSumRCO3(i)=fSumRCO3(i)+.938;

i=i+1;
Rnames{i} = 'H3F2ONM4_P2 + NO3 = NO2 + .5*R2CO3 + .5*HCHO + .5*MECO3 + .5*MGLY + .5*XC + SumRCO3';
k(:,i) = 2.30e-12;
Gstr{i,1} = 'H3F2ONM4_P2'; Gstr{i,2} = 'NO3';
fH3F2ONM4_P2(i)=fH3F2ONM4_P2(i)-1; fNO3(i)=fNO3(i)-1; fNO2(i)=fNO2(i)+1; fR2CO3(i)=fR2CO3(i)+.5; fHCHO(i)=fHCHO(i)+.5; fMECO3(i)=fMECO3(i)+.5; fMGLY(i)=fMGLY(i)+.5; fXC(i)=fXC(i)+.5; fSumRCO3(i)=fSumRCO3(i)+1;

i=i+1;
Rnames{i} = 'H3F2ONM4_P2 + HO2 = .9*ROOH + .05*MGLY + .05*BACL + .15*XC';
k(:,i) = 1.78e-11;
Gstr{i,1} = 'H3F2ONM4_P2'; Gstr{i,2} = 'HO2';
fH3F2ONM4_P2(i)=fH3F2ONM4_P2(i)-1; fHO2(i)=fHO2(i)-1; fROOH(i)=fROOH(i)+.9; fMGLY(i)=fMGLY(i)+.05; fBACL(i)=fBACL(i)+.05; fXC(i)=fXC(i)+.15;

i=i+1;
Rnames{i} = 'H3F2ONM4_P2 + SumRO2 = SumRO2 + .375*MGLY + .25*R2CO3 + .25*HCHO + .25*KET2 + .25*MECO3 + .125*BACL + .375*XC + .5*SumRCO3';
k(:,i) = 1.39e-12;
Gstr{i,1} = 'H3F2ONM4_P2'; Gstr{i,2} = 'SumRO2';
fH3F2ONM4_P2(i)=fH3F2ONM4_P2(i)-1; fSumRO2(i)=fSumRO2(i)-1; fSumRO2(i)=fSumRO2(i)+1; fMGLY(i)=fMGLY(i)+.375; fR2CO3(i)=fR2CO3(i)+.25; fHCHO(i)=fHCHO(i)+.25; fKET2(i)=fKET2(i)+.25; fMECO3(i)=fMECO3(i)+.25; fBACL(i)=fBACL(i)+.125; fXC(i)=fXC(i)+.375; fSumRCO3(i)=fSumRCO3(i)+.5;

i=i+1;
Rnames{i} = 'H3F2ONM4_P2 + SumRCO3 = SumRCO3 + .5*MGLY + .4*R2CO3 + .4*HCHO + .4*MECO3 + .1*BACL + .7*XC + .8*SumRCO3';
k(:,i) = 1.37e-11;
Gstr{i,1} = 'H3F2ONM4_P2'; Gstr{i,2} = 'SumRCO3';
fH3F2ONM4_P2(i)=fH3F2ONM4_P2(i)-1; fSumRCO3(i)=fSumRCO3(i)-1; fSumRCO3(i)=fSumRCO3(i)+1; fMGLY(i)=fMGLY(i)+.5; fR2CO3(i)=fR2CO3(i)+.4; fHCHO(i)=fHCHO(i)+.4; fMECO3(i)=fMECO3(i)+.4; fBACL(i)=fBACL(i)+.1; fXC(i)=fXC(i)+.7; fSumRCO3(i)=fSumRCO3(i)+.8;

i=i+1;
Rnames{i} = 'H3F2ONM4_P3 + NO = 1.876*NO2 + .938*KET2 + .062*RCNO3 - 0.876*XC + .062*XN';
k(:,i) = 9.13e-12;
Gstr{i,1} = 'H3F2ONM4_P3'; Gstr{i,2} = 'NO';
fH3F2ONM4_P3(i)=fH3F2ONM4_P3(i)-1; fNO(i)=fNO(i)-1; fNO2(i)=fNO2(i)+1.876; fKET2(i)=fKET2(i)+.938; fRCNO3(i)=fRCNO3(i)+.062; fXC(i)=fXC(i)-0.876; fXN(i)=fXN(i)+.062;

i=i+1;
Rnames{i} = 'H3F2ONM4_P3 + NO3 = 2*NO2 + KET2 - 1*XC';
k(:,i) = 2.30e-12;
Gstr{i,1} = 'H3F2ONM4_P3'; Gstr{i,2} = 'NO3';
fH3F2ONM4_P3(i)=fH3F2ONM4_P3(i)-1; fNO3(i)=fNO3(i)-1; fNO2(i)=fNO2(i)+2; fKET2(i)=fKET2(i)+1; fXC(i)=fXC(i)-1;

i=i+1;
Rnames{i} = 'H3F2ONM4_P3 + HO2 = RCNO3 + XC';
k(:,i) = 1.78e-11;
Gstr{i,1} = 'H3F2ONM4_P3'; Gstr{i,2} = 'HO2';
fH3F2ONM4_P3(i)=fH3F2ONM4_P3(i)-1; fHO2(i)=fHO2(i)-1; fRCNO3(i)=fRCNO3(i)+1; fXC(i)=fXC(i)+1;

i=i+1;
Rnames{i} = 'H3F2ONM4_P3 + SumRO2 = SumRO2 + .5*KET2 + .5*NO2 + .5*RCNO3';
k(:,i) = 2.57e-12;
Gstr{i,1} = 'H3F2ONM4_P3'; Gstr{i,2} = 'SumRO2';
fH3F2ONM4_P3(i)=fH3F2ONM4_P3(i)-1; fSumRO2(i)=fSumRO2(i)-1; fSumRO2(i)=fSumRO2(i)+1; fKET2(i)=fKET2(i)+.5; fNO2(i)=fNO2(i)+.5; fRCNO3(i)=fRCNO3(i)+.5;

i=i+1;
Rnames{i} = 'H3F2ONM4_P3 + SumRCO3 = SumRCO3 + .8*KET2 + .8*NO2 + .2*RCNO3 - 0.6*XC';
k(:,i) = 1.37e-11;
Gstr{i,1} = 'H3F2ONM4_P3'; Gstr{i,2} = 'SumRCO3';
fH3F2ONM4_P3(i)=fH3F2ONM4_P3(i)-1; fSumRCO3(i)=fSumRCO3(i)-1; fSumRCO3(i)=fSumRCO3(i)+1; fKET2(i)=fKET2(i)+.8; fNO2(i)=fNO2(i)+.8; fRCNO3(i)=fRCNO3(i)+.2; fXC(i)=fXC(i)-0.6;

% End


