%% PHEN CRES CATEHOL RXNS

SpeciesToAdd = {...
'PHEN'; 'PHEN_A1'; 'PHEN_P1'; 'SVPHE'; 'HCHDO'; 'PHEN_P2'; 'HCHDO_P1'; 'SVPHE_P1'; 'CRES'; 'CRES_P1';'AFG4A';}

AddSpecies

i=i+1;
Rnames{i} = 'PHEN + OH = .95*PHEN_A1 + .05*BZO + .00*HO2';
k(:,i) = 4.70e-13.*(T./300).^0.00.*exp(1219.807./T);
Gstr{i,1} = 'PHEN'; Gstr{i,2} = 'OH';
fPHEN(i)=fPHEN(i)-1; fOH(i)=fOH(i)-1; fPHEN_A1(i)=fPHEN_A1(i)+0.95; fBZO(i)=fBZO(i)+.05; fHO2(i)=fHO2(i)+0.00;

i=i+1;
Rnames{i} = 'PHEN + NO3 = BZO + HNO3';
k(:,i) = 1.10e-11;
Gstr{i,1} = 'PHEN'; Gstr{i,2} = 'NO3';
fPHEN(i)=fPHEN(i)-1; fNO3(i)=fNO3(i)-1; fBZO(i)=fBZO(i)+1; fHNO3(i)=fHNO3(i)+1;

i=i+1;
Rnames{i} = 'PHEN_A1 + O2 = .25*PHEN_P1 + .75*SVPHE + .75*HO2';
k(:,i) = 5.10e-13;
Gstr{i,1} = 'PHEN_A1'; Gstr{i,2} = 'O2';
fPHEN_A1(i)=fPHEN_A1(i)-1; fO2(i)=fO2(i)-1; fPHEN_P1(i)=fPHEN_P1(i)+.25; fSVPHE(i)=fSVPHE(i)+.75; fHO2(i)=fHO2(i)+.75;

i=i+1;
Rnames{i} = 'PHEN_A1 + NO2 = SVPHE + HONO';
k(:,i) = 3.00e-11;
Gstr{i,1} = 'PHEN_A1'; Gstr{i,2} = 'NO2';
fPHEN_A1(i)=fPHEN_A1(i)-1; fNO2(i)=fNO2(i)-1; fSVPHE(i)=fSVPHE(i)+1; fHONO(i)=fHONO(i)+1;

i=i+1;
Rnames{i} = 'PHEN_P1 = .996*HCHDO + .996*HO2 + .004*PHEN_P2 - 1.992*XC + .004*SumRO2';
k(:,i) = 1.10e+6;
Gstr{i,1} = 'PHEN_P1';
fPHEN_P1(i)=fPHEN_P1(i)-1; fHCHDO(i)=fHCHDO(i)+.996; fHO2(i)=fHO2(i)+.996; fPHEN_P2(i)=fPHEN_P2(i)+.004; fXC(i)=fXC(i)-1.992; fSumRO2(i)=fSumRO2(i)+.004;

i=i+1;
Rnames{i} = 'PHEN_P1 + NO = OLEA1 + HO2 + NO2';
k(:,i) = 9.05e-12;
Gstr{i,1} = 'PHEN_P1'; Gstr{i,2} = 'NO';
fPHEN_P1(i)=fPHEN_P1(i)-1; fNO(i)=fNO(i)-1; fOLEA1(i)=fOLEA1(i)+1; fHO2(i)=fHO2(i)+1; fNO2(i)=fNO2(i)+1;

i=i+1;
Rnames{i} = 'PHEN_P1 + HO2 = .5*OLEA1 + .5*HO2 + .5*OH + .5*RUOOH';
k(:,i) = 6.82e-12;
Gstr{i,1} = 'PHEN_P1'; Gstr{i,2} = 'HO2';
fPHEN_P1(i)=fPHEN_P1(i)-1; fHO2(i)=fHO2(i)-1; fOLEA1(i)=fOLEA1(i)+.5; fHO2(i)=fHO2(i)+.5; fOH(i)=fOH(i)+.5; fRUOOH(i)=fRUOOH(i)+.5;

i=i+1;
Rnames{i} = 'PHEN_P2 + NO = .876*HO2 + .876*NO2 + .438*AFG2A + .438*GLY + .438*BUDAL + .438*MGLY + .124*RPNO3 - 1.124*XC';
k(:,i) = 9.13e-12;
Gstr{i,1} = 'PHEN_P2'; Gstr{i,2} = 'NO';
fPHEN_P2(i)=fPHEN_P2(i)-1; fNO(i)=fNO(i)-1; fHO2(i)=fHO2(i)+.876; fNO2(i)=fNO2(i)+.876; fAFG2A(i)=fAFG2A(i)+.438; fGLY(i)=fGLY(i)+.438; fBUDAL(i)=fBUDAL(i)+.438; fMGLY(i)=fMGLY(i)+.438; fRPNO3(i)=fRPNO3(i)+.124; fXC(i)=fXC(i)-1.124;

i=i+1;
Rnames{i} = 'PHEN_P2 + NO3 = HO2 + NO2 + .5*AFG2A + .5*GLY + .5*BUDAL + .5*MGLY - 1*XC';
k(:,i) = 2.30e-12;
Gstr{i,1} = 'PHEN_P2'; Gstr{i,2} = 'NO3';
fPHEN_P2(i)=fPHEN_P2(i)-1; fNO3(i)=fNO3(i)-1; fHO2(i)=fHO2(i)+1; fNO2(i)=fNO2(i)+1; fAFG2A(i)=fAFG2A(i)+.5; fGLY(i)=fGLY(i)+.5; fBUDAL(i)=fBUDAL(i)+.5; fMGLY(i)=fMGLY(i)+.5; fXC(i)=fXC(i)-1;

i=i+1;
Rnames{i} = 'PHEN_P2 + HO2 = RAOOH - 2*XC';
k(:,i) = 1.95e-11;
Gstr{i,1} = 'PHEN_P2'; Gstr{i,2} = 'HO2';
fPHEN_P2(i)=fPHEN_P2(i)-1; fHO2(i)=fHO2(i)-1; fRAOOH(i)=fRAOOH(i)+1; fXC(i)=fXC(i)-2;

i=i+1;
Rnames{i} = 'PHEN_P2 + SumRO2 = SumRO2 + .5*HO2 + .5*OLEP + .25*AFG2A + .25*GLY + .25*BUDAL + .25*MGLY - 2.5*XC';
k(:,i) = 2.57e-12;
Gstr{i,1} = 'PHEN_P2'; Gstr{i,2} = 'SumRO2';
fPHEN_P2(i)=fPHEN_P2(i)-1; fSumRO2(i)=fSumRO2(i)-1; fSumRO2(i)=fSumRO2(i)+1; fHO2(i)=fHO2(i)+.5; fOLEP(i)=fOLEP(i)+.5; fAFG2A(i)=fAFG2A(i)+.25; fGLY(i)=fGLY(i)+.25; fBUDAL(i)=fBUDAL(i)+.25; fMGLY(i)=fMGLY(i)+.25; fXC(i)=fXC(i)-2.5;

i=i+1;
Rnames{i} = 'PHEN_P2 + SumRCO3 = SumRCO3 + .8*HO2 + .4*AFG2A + .4*GLY + .4*BUDAL + .4*MGLY + .2*OLEP - 1.6*XC';
k(:,i) = 1.37e-11;
Gstr{i,1} = 'PHEN_P2'; Gstr{i,2} = 'SumRCO3';
fPHEN_P2(i)=fPHEN_P2(i)-1; fSumRCO3(i)=fSumRCO3(i)-1; fSumRCO3(i)=fSumRCO3(i)+1; fHO2(i)=fHO2(i)+.8; fAFG2A(i)=fAFG2A(i)+.4; fGLY(i)=fGLY(i)+.4; fBUDAL(i)=fBUDAL(i)+.4; fMGLY(i)=fMGLY(i)+.4; fOLEP(i)=fOLEP(i)+.2; fXC(i)=fXC(i)-1.6;

i=i+1;
Rnames{i} = 'HCHDO + OH = .989*HCHDO_P1 + .011*BACL + .011*HO2 + .022*XC + .989*SumRO2';
k(:,i) = 1.73e-10;
Gstr{i,1} = 'HCHDO'; Gstr{i,2} = 'OH';
fHCHDO(i)=fHCHDO(i)-1; fOH(i)=fOH(i)-1; fHCHDO_P1(i)=fHCHDO_P1(i)+.989; fBACL(i)=fBACL(i)+.011; fHO2(i)=fHO2(i)+.011; fXC(i)=fXC(i)+.022; fSumRO2(i)=fSumRO2(i)+.989;

i=i+1;
Rnames{i} = 'HCHDO + O3 = .75*AFG2A + .643*HO2 + .613*OH + .25*RCHO + .137*CO2 + 1.113*XC';
k(:,i) = 3.77e-17;
Gstr{i,1} = 'HCHDO'; Gstr{i,2} = 'O3';
fHCHDO(i)=fHCHDO(i)-1; fO3(i)=fO3(i)-1; fAFG2A(i)=fAFG2A(i)+.75; fHO2(i)=fHO2(i)+.643; fOH(i)=fOH(i)+.613; fRCHO(i)=fRCHO(i)+.25; fCO2(i)=fCO2(i)+.137; fXC(i)=fXC(i)+1.113;

i=i+1;
Rnames{i} = 'HCHDO + HV = OLEP + CO - 5*XC';
k(:,i) = JMVK_16;
Gstr{i,1} = 'HCHDO';
fHCHDO(i)=fHCHDO(i)-1; fOLEP(i)=fOLEP(i)+1; fCO(i)=fCO(i)+1; fXC(i)=fXC(i)-5;

i=i+1;
Rnames{i} = 'HCHDO_P1 + NO = .931*NO2 + .747*HO2 + .542*MGLY + .204*AFG2A + .185*MACO3 + .069*RCNO3 + 2.269*XC + .185*SumRCO3';
k(:,i) = 9.05e-12;
Gstr{i,1} = 'HCHDO_P1'; Gstr{i,2} = 'NO';
fHCHDO_P1(i)=fHCHDO_P1(i)-1; fNO(i)=fNO(i)-1; fNO2(i)=fNO2(i)+.931; fHO2(i)=fHO2(i)+.747; fMGLY(i)=fMGLY(i)+.542; fAFG2A(i)=fAFG2A(i)+.204; fMACO3(i)=fMACO3(i)+.185; fRCNO3(i)=fRCNO3(i)+.069; fXC(i)=fXC(i)+2.269; fSumRCO3(i)=fSumRCO3(i)+.185;

i=i+1;
Rnames{i} = 'HCHDO_P1 + NO3 = NO2 + .802*HO2 + .582*MGLY + .219*AFG2A + .198*MACO3 + 2.367*XC + .198*SumRCO3';
k(:,i) = 2.30e-12;
Gstr{i,1} = 'HCHDO_P1'; Gstr{i,2} = 'NO3';
fHCHDO_P1(i)=fHCHDO_P1(i)-1; fNO3(i)=fNO3(i)-1; fNO2(i)=fNO2(i)+1; fHO2(i)=fHO2(i)+.802; fMGLY(i)=fMGLY(i)+.582; fAFG2A(i)=fAFG2A(i)+.219; fMACO3(i)=fMACO3(i)+.198; fXC(i)=fXC(i)+2.367; fSumRCO3(i)=fSumRCO3(i)+.198;

i=i+1;
Rnames{i} = 'HCHDO_P1 + HO2 = .759*RUOOH + .219*OTHN + .021*BACL - 1.266*XC';
k(:,i) = 1.95e-11;
Gstr{i,1} = 'HCHDO_P1'; Gstr{i,2} = 'HO2';
fHCHDO_P1(i)=fHCHDO_P1(i)-1; fHO2(i)=fHO2(i)-1; fRUOOH(i)=fRUOOH(i)+.759; fOTHN(i)=fOTHN(i)+.219; fBACL(i)=fBACL(i)+.021; fXC(i)=fXC(i)-1.266;

i=i+1;
Rnames{i} = 'HCHDO_P1 + SumRO2 = SumRO2 + .401*HO2 + .291*MGLY + .202*LVKS + .195*OLEP + .11*AFG2A + .099*MACO3 + .053*BACL + .05*AFG3 + .103*XC + .099*SumRCO3';
k(:,i) = 2.57e-12;
Gstr{i,1} = 'HCHDO_P1'; Gstr{i,2} = 'SumRO2';
fHCHDO_P1(i)=fHCHDO_P1(i)-1; fSumRO2(i)=fSumRO2(i)-1; fSumRO2(i)=fSumRO2(i)+1; fHO2(i)=fHO2(i)+.401; fMGLY(i)=fMGLY(i)+.291; fLVKS(i)=fLVKS(i)+.202; fOLEP(i)=fOLEP(i)+.195; fAFG2A(i)=fAFG2A(i)+.11; fMACO3(i)=fMACO3(i)+.099; fBACL(i)=fBACL(i)+.053; fAFG3(i)=fAFG3(i)+.05; fXC(i)=fXC(i)+.103; fSumRCO3(i)=fSumRCO3(i)+.099;

i=i+1;
Rnames{i} = 'HCHDO_P1 + SumRCO3 = SumRCO3 + .641*HO2 + .466*MGLY + .175*AFG2A + .159*MACO3 + .118*LVKS + .043*BACL + .04*AFG3 + 1.735*XC + .159*SumRCO3';
k(:,i) = 1.37e-11;
Gstr{i,1} = 'HCHDO_P1'; Gstr{i,2} = 'SumRCO3';
fHCHDO_P1(i)=fHCHDO_P1(i)-1; fSumRCO3(i)=fSumRCO3(i)-1; fSumRCO3(i)=fSumRCO3(i)+1; fHO2(i)=fHO2(i)+.641; fMGLY(i)=fMGLY(i)+.466; fAFG2A(i)=fAFG2A(i)+.175; fMACO3(i)=fMACO3(i)+.159; fLVKS(i)=fLVKS(i)+.118; fBACL(i)=fBACL(i)+.043; fAFG3(i)=fAFG3(i)+.04; fXC(i)=fXC(i)+1.735; fSumRCO3(i)=fSumRCO3(i)+.159;

i=i+1;
Rnames{i} = 'SVPHE + OH = .774*HO2 + .386*OLEP + .00*HO2 + .197*SVPHE_P1 + .132*OLEA1 + .123*XYNL + .087*OLEA2 + .045*LVKS + .03*BZO - 2.228*XC + .197*SumRO2';
k(:,i) = 1.04e-10;
Gstr{i,1} = 'SVPHE'; Gstr{i,2} = 'OH';
fSVPHE(i)=fSVPHE(i)-1; fOH(i)=fOH(i)-1; fHO2(i)=fHO2(i)+.774; fOLEP(i)=fOLEP(i)+.386; fHO2(i)=fHO2(i)+.00; fSVPHE_P1(i)=fSVPHE_P1(i)+.197; fOLEA1(i)=fOLEA1(i)+.132; fXYNL(i)=fXYNL(i)+.123; fOLEA2(i)=fOLEA2(i)+.087; fLVKS(i)=fLVKS(i)+.045; fBZO(i)=fBZO(i)+.03; fXC(i)=fXC(i)-2.228; fSumRO2(i)=fSumRO2(i)+.197;

i=i+1;
Rnames{i} = 'SVPHE + NO3 = BZO + XN';
k(:,i) = 0.98e-10;
Gstr{i,1} = 'SVPHE'; Gstr{i,2} = 'NO3';
fSVPHE(i)=fSVPHE(i)-1; fNO3(i)=fNO3(i)-1; fBZO(i)=fBZO(i)+1; fXN(i)=fXN(i)+1;

i=i+1;
Rnames{i} = 'SVPHE_P1 + NO = .876*HO2 + .876*NO2 + .438*AFG4A + .438*MGLY + .438*BUDAL + .438*BACL + .124*RPNO3 - 1.562*XC';
k(:,i) = 9.13e-12;
Gstr{i,1} = 'SVPHE_P1'; Gstr{i,2} = 'NO';
fSVPHE_P1(i)=fSVPHE_P1(i)-1; fNO(i)=fNO(i)-1; fHO2(i)=fHO2(i)+.876; fNO2(i)=fNO2(i)+.876; fAFG4A(i)=fAFG4A(i)+.438; fMGLY(i)=fMGLY(i)+.438; fBUDAL(i)=fBUDAL(i)+.438; fBACL(i)=fBACL(i)+.438; fRPNO3(i)=fRPNO3(i)+.124; fXC(i)=fXC(i)-1.562;

i=i+1;
Rnames{i} = 'SVPHE_P1 + NO3 = HO2 + NO2 + .5*AFG4A + .5*MGLY + .5*BUDAL + .5*BACL - 1.5*XC';
k(:,i) = 2.30e-12;
Gstr{i,1} = 'SVPHE_P1'; Gstr{i,2} = 'NO3';
fSVPHE_P1(i)=fSVPHE_P1(i)-1; fNO3(i)=fNO3(i)-1; fHO2(i)=fHO2(i)+1; fNO2(i)=fNO2(i)+1; fAFG4A(i)=fAFG4A(i)+.5; fMGLY(i)=fMGLY(i)+.5; fBUDAL(i)=fBUDAL(i)+.5; fBACL(i)=fBACL(i)+.5; fXC(i)=fXC(i)-1.5;

i=i+1;
Rnames{i} = 'SVPHE_P1 + HO2 = OTHN - 6*XC';
k(:,i) = 1.95e-11;
Gstr{i,1} = 'SVPHE_P1'; Gstr{i,2} = 'HO2';
fSVPHE_P1(i)=fSVPHE_P1(i)-1; fHO2(i)=fHO2(i)-1; fOTHN(i)=fOTHN(i)+1; fXC(i)=fXC(i)-6;

i=i+1;
Rnames{i} = 'SVPHE_P1 + SumRO2 = SumRO2 + .5*HO2 + .5*OLEP + .25*AFG4A + .25*MGLY + .25*BUDAL + .25*BACL - 2.75*XC';
k(:,i) = 2.57e-12;
Gstr{i,1} = 'SVPHE_P1'; Gstr{i,2} = 'SumRO2';
fSVPHE_P1(i)=fSVPHE_P1(i)-1; fSumRO2(i)=fSumRO2(i)-1; fSumRO2(i)=fSumRO2(i)+1; fHO2(i)=fHO2(i)+.5; fOLEP(i)=fOLEP(i)+.5; fAFG4A(i)=fAFG4A(i)+.25; fMGLY(i)=fMGLY(i)+.25; fBUDAL(i)=fBUDAL(i)+.25; fBACL(i)=fBACL(i)+.25; fXC(i)=fXC(i)-2.75;

i=i+1;
Rnames{i} = 'SVPHE_P1 + SumRCO3 = SumRCO3 + .8*HO2 + .4*AFG4A + .4*MGLY + .4*BUDAL + .4*BACL + .2*OLEP - 2*XC';
k(:,i) = 1.37e-11;
Gstr{i,1} = 'SVPHE_P1'; Gstr{i,2} = 'SumRCO3';
fSVPHE_P1(i)=fSVPHE_P1(i)-1; fSumRCO3(i)=fSumRCO3(i)-1; fSumRCO3(i)=fSumRCO3(i)+1; fHO2(i)=fHO2(i)+.8; fAFG4A(i)=fAFG4A(i)+.4; fMGLY(i)=fMGLY(i)+.4; fBUDAL(i)=fBUDAL(i)+.4; fBACL(i)=fBACL(i)+.4; fOLEP(i)=fOLEP(i)+.2; fXC(i)=fXC(i)-2;

i=i+1;
Rnames{i} = 'CRES + OH = .631*BZO + .215*CRES_P1 + .154*HO2 + .053*OLEA1 + .036*LVKS + .026*OLEP + .021*XYNL + .018*OLEA2 + .495*XC + .215*SumRO2';
k(:,i) = 1.60e-12.*(T./300).^0.00.*exp(998.893./T);
Gstr{i,1} = 'CRES'; Gstr{i,2} = 'OH';
fCRES(i)=fCRES(i)-1; fOH(i)=fOH(i)-1; fBZO(i)=fBZO(i)+.631; fCRES_P1(i)=fCRES_P1(i)+.215; fHO2(i)=fHO2(i)+.154; fOLEA1(i)=fOLEA1(i)+.053; fLVKS(i)=fLVKS(i)+.036; fOLEP(i)=fOLEP(i)+.026; fXYNL(i)=fXYNL(i)+.021; fOLEA2(i)=fOLEA2(i)+.018; fXC(i)=fXC(i)+.495; fSumRO2(i)=fSumRO2(i)+.215;

i=i+1;
Rnames{i} = 'CRES + NO3 = BZO + XC + XN';
k(:,i) = 1.40e-11;
Gstr{i,1} = 'CRES'; Gstr{i,2} = 'NO3';
fCRES(i)=fCRES(i)-1; fNO3(i)=fNO3(i)-1; fBZO(i)=fBZO(i)+1; fXC(i)=fXC(i)+1; fXN(i)=fXN(i)+1;

i=i+1;
Rnames{i} = 'CRES_P1 + NO = .811*HO2 + .811*NO2 + .385*MGLY + .385*AFG2A + .291*BACL + .291*BUDAL + .186*RPNO3 + .095*AFG2B + .095*GLY + .04*BALD + .004*RHNO3 - 0.968*XC - 0.001*XN';
k(:,i) = 9.13e-12;
Gstr{i,1} = 'CRES_P1'; Gstr{i,2} = 'NO';
fCRES_P1(i)=fCRES_P1(i)-1; fNO(i)=fNO(i)-1; fHO2(i)=fHO2(i)+.811; fNO2(i)=fNO2(i)+.811; fMGLY(i)=fMGLY(i)+.385; fAFG2A(i)=fAFG2A(i)+.385; fBACL(i)=fBACL(i)+.291; fBUDAL(i)=fBUDAL(i)+.291; fRPNO3(i)=fRPNO3(i)+.186; fAFG2B(i)=fAFG2B(i)+.095; fGLY(i)=fGLY(i)+.095; fBALD(i)=fBALD(i)+.04; fRHNO3(i)=fRHNO3(i)+.004; fXC(i)=fXC(i)-0.968; fXN(i)=fXN(i)-0.001;

i=i+1;
Rnames{i} = 'CRES_P1 + NO3 = HO2 + NO2 + .478*MGLY + .478*AFG2A + .361*BACL + .361*BUDAL + .117*AFG2B + .117*GLY + .044*BALD - 0.956*XC';
k(:,i) = 2.30e-12;
Gstr{i,1} = 'CRES_P1'; Gstr{i,2} = 'NO3';
fCRES_P1(i)=fCRES_P1(i)-1; fNO3(i)=fNO3(i)-1; fHO2(i)=fHO2(i)+1; fNO2(i)=fNO2(i)+1; fMGLY(i)=fMGLY(i)+.478; fAFG2A(i)=fAFG2A(i)+.478; fBACL(i)=fBACL(i)+.361; fBUDAL(i)=fBUDAL(i)+.361; fAFG2B(i)=fAFG2B(i)+.117; fGLY(i)=fGLY(i)+.117; fBALD(i)=fBALD(i)+.044; fXC(i)=fXC(i)-0.956;

i=i+1;
Rnames{i} = 'CRES_P1 + HO2 = .839*RAOOH + .117*OTHN + .044*ROOH - 1.336*XC';
k(:,i) = 2.11e-11;
Gstr{i,1} = 'CRES_P1'; Gstr{i,2} = 'HO2';
fCRES_P1(i)=fCRES_P1(i)-1; fHO2(i)=fHO2(i)-1; fRAOOH(i)=fRAOOH(i)+.839; fOTHN(i)=fOTHN(i)+.117; fROOH(i)=fROOH(i)+.044; fXC(i)=fXC(i)-1.336;

i=i+1;
Rnames{i} = 'CRES_P1 + SumRO2 = SumRO2 + .5*HO2 + .478*OLEP + .239*MGLY + .239*AFG2A + .18*BACL + .18*BUDAL + .059*AFG2B + .059*GLY + .033*BALD + .011*XYNL - 1.923*XC';
k(:,i) = 2.16e-12;
Gstr{i,1} = 'CRES_P1'; Gstr{i,2} = 'SumRO2';
fCRES_P1(i)=fCRES_P1(i)-1; fSumRO2(i)=fSumRO2(i)-1; fSumRO2(i)=fSumRO2(i)+1; fHO2(i)=fHO2(i)+.5; fOLEP(i)=fOLEP(i)+.478; fMGLY(i)=fMGLY(i)+.239; fAFG2A(i)=fAFG2A(i)+.239; fBACL(i)=fBACL(i)+.18; fBUDAL(i)=fBUDAL(i)+.18; fAFG2B(i)=fAFG2B(i)+.059; fGLY(i)=fGLY(i)+.059; fBALD(i)=fBALD(i)+.033; fXYNL(i)=fXYNL(i)+.011; fXC(i)=fXC(i)-1.923;

i=i+1;
Rnames{i} = 'CRES_P1 + SumRCO3 = SumRCO3 + .823*HO2 + .406*MGLY + .406*AFG2A + .289*BACL + .289*BUDAL + .168*OLEP + .094*AFG2B + .094*GLY + .044*BALD - 1.3*XC';
k(:,i) = 1.37e-11;
Gstr{i,1} = 'CRES_P1'; Gstr{i,2} = 'SumRCO3';
fCRES_P1(i)=fCRES_P1(i)-1; fSumRCO3(i)=fSumRCO3(i)-1; fSumRCO3(i)=fSumRCO3(i)+1; fHO2(i)=fHO2(i)+.823; fMGLY(i)=fMGLY(i)+.406; fAFG2A(i)=fAFG2A(i)+.406; fBACL(i)=fBACL(i)+.289; fBUDAL(i)=fBUDAL(i)+.289; fOLEP(i)=fOLEP(i)+.168; fAFG2B(i)=fAFG2B(i)+.094; fGLY(i)=fGLY(i)+.094; fBALD(i)=fBALD(i)+.044; fXC(i)=fXC(i)-1.3;
