
SpeciesToAdd = {'RNO3';'xOH';'yBACL';'xMACR';'xCO';...
    'xRCO3';'zHFONS';'MEO2';'RO2C';'BUTEDIAL';'xGLY';'O3';...
    'xBACL';'CO2';'xRCHO';'zRNO3';'xHCHO';'HNO3';'yROOH';...
    'OH';'PROD2';'xMACO3';'yAFG3';'OLE2';'HCHO';'HO2';'RCOOH';...
    'xAFG3';'M2BUTDAL';'yRNO3';'M2FURAN';'xHFONS';'GLY';'xRNO3';'M25FUR';...
    'RCHO';'ROOH';'MACO3';'AFG3';'NO3';'M3FURAN';'CO';'RO2XC';'xMEO2';...
    'xROOH';'xNO2';'yMGLY';'xHO2';'xMGLY';'ALK3';'MGLY';'xMECO3';...
    'yHFONS';'FURAN';'HFONS';'O4X2PEAL';'H3XE25DO';'MALANHY';};

AddSpecies

i=i+1;
Rnames{i} = 'FURAN + OH = 0.309*RO2C + 0.008*RO2XC + 0.760*HO2 + 0.232*xHO2 + 0.760*BUTEDIAL + 0.156*xHFONS + 0.075*xRCHO + 0.024*yHFONS + 0.008*zRNO3 + 0.293*yROOH';
k(:,i) = 3.974e-11;
Gstr{i,1} = 'FURAN'; Gstr{i,2} = 'OH';
fFURAN(i)=fFURAN(i)-1; fOH(i)=fOH(i)-1; fRO2C(i)=fRO2C(i)+0.309; fRO2XC(i)=fRO2XC(i)+0.008; fHO2(i)=fHO2(i)+0.760; fxHO2(i)=fxHO2(i)+0.232; fBUTEDIAL(i)=fBUTEDIAL(i)+0.760; fxHFONS(i)=fxHFONS(i)+0.156; fxRCHO(i)=fxRCHO(i)+0.075; fyHFONS(i)=fyHFONS(i)+0.024; fzRNO3(i)=fzRNO3(i)+0.008; fyROOH(i)=fyROOH(i)+0.293;

i=i+1;
Rnames{i} = 'FURAN + NO3 = 1.288*RO2C + 0.034*RO2XC + 0.653*xHO2 + 0.313*xNO2 + 0.653*xHFONS + 0.313*xRCHO + 0.100*yHFONS + 0.026*zRNO3 + 0.900*yRNO3';
k(:,i) = 1.360e-12;
Gstr{i,1} = 'FURAN'; Gstr{i,2} = 'NO3';
fFURAN(i)=fFURAN(i)-1; fNO3(i)=fNO3(i)-1; fRO2C(i)=fRO2C(i)+1.288; fRO2XC(i)=fRO2XC(i)+0.034; fxHO2(i)=fxHO2(i)+0.653; fxNO2(i)=fxNO2(i)+0.313; fxHFONS(i)=fxHFONS(i)+0.653; fxRCHO(i)=fxRCHO(i)+0.313; fyHFONS(i)=fyHFONS(i)+0.100; fzRNO3(i)=fzRNO3(i)+0.026; fyRNO3(i)=fyRNO3(i)+0.900;

i=i+1;
Rnames{i} = 'FURAN + O3 = RCHO';
k(:,i) = 2.410e-18;
Gstr{i,1} = 'FURAN'; Gstr{i,2} = 'O3';
fFURAN(i)=fFURAN(i)-1; fO3(i)=fO3(i)-1; fRCHO(i)=fRCHO(i)+1;

i=i+1;
Rnames{i} = 'M2FURAN + OH = 0.612*RO2C + 0.040*RO2XC + 0.440*HO2 + 0.346*xHO2 + 0.440*O4X2PEAL + 0.432*xHFONS + 0.087*xRCHO + 0.173*xMEO2 + 0.030*yHFONS + 0.040*zRNO3 + 0.623*yROOH';
k(:,i) = 7.310e-11;
Gstr{i,1} = 'M2FURAN'; Gstr{i,2} = 'OH';
fM2FURAN(i)=fM2FURAN(i)-1; fOH(i)=fOH(i)-1; fRO2C(i)=fRO2C(i)+0.612; fRO2XC(i)=fRO2XC(i)+0.040; fHO2(i)=fHO2(i)+0.440; fxHO2(i)=fxHO2(i)+0.346; fO4X2PEAL(i)=fO4X2PEAL(i)+0.440; fxHFONS(i)=fxHFONS(i)+0.432; fxRCHO(i)=fxRCHO(i)+0.087; fxMEO2(i)=fxMEO2(i)+0.173; fyHFONS(i)=fyHFONS(i)+0.030; fzRNO3(i)=fzRNO3(i)+0.040; fyROOH(i)=fyROOH(i)+0.623;

i=i+1;
Rnames{i} = 'M2FURAN + NO3 = 1.094*RO2C + 0.072*RO2XC + 0.463*xHO2 + 0.147*xNO2 + 0.009*xOH + 0.772*xHFONS + 0.147*xRCHO + 0.308*xMEO2 + 0.054*yHFONS + 0.062*zRNO3 + 0.946*yRNO3';
k(:,i) = 2.570e-11;
Gstr{i,1} = 'M2FURAN'; Gstr{i,2} = 'NO3';
fM2FURAN(i)=fM2FURAN(i)-1; fNO3(i)=fNO3(i)-1; fRO2C(i)=fRO2C(i)+1.094; fRO2XC(i)=fRO2XC(i)+0.072; fxHO2(i)=fxHO2(i)+0.463; fxNO2(i)=fxNO2(i)+0.147; fxOH(i)=fxOH(i)+0.009; fxHFONS(i)=fxHFONS(i)+0.772; fxRCHO(i)=fxRCHO(i)+0.147; fxMEO2(i)=fxMEO2(i)+0.308; fyHFONS(i)=fyHFONS(i)+0.054; fzRNO3(i)=fzRNO3(i)+0.062; fyRNO3(i)=fyRNO3(i)+0.946;

i=i+1;
Rnames{i} = 'M2FURAN + O3 = 0.192*RO2C + 0.013*RO2XC + 0.613*OH + 0.052*xHO2 + 0.140*xOH + 0.591*RCHO + 0.204*OLE2 + 0.052*xCO + 0.192*xMACR + 0.013*zRNO3 + 0.015*yBACL + 0.134*yROOH + 0.140*CO2';
k(:,i) = 1.995e-17;
Gstr{i,1} = 'M2FURAN'; Gstr{i,2} = 'O3';
fM2FURAN(i)=fM2FURAN(i)-1; fO3(i)=fO3(i)-1; fRO2C(i)=fRO2C(i)+0.192; fRO2XC(i)=fRO2XC(i)+0.013; fOH(i)=fOH(i)+0.613; fxHO2(i)=fxHO2(i)+0.052; fxOH(i)=fxOH(i)+0.140; fRCHO(i)=fRCHO(i)+0.591; fOLE2(i)=fOLE2(i)+0.204; fxCO(i)=fxCO(i)+0.052; fxMACR(i)=fxMACR(i)+0.192; fzRNO3(i)=fzRNO3(i)+0.013; fyBACL(i)=fyBACL(i)+0.015; fyROOH(i)=fyROOH(i)+0.134; fCO2(i)=fCO2(i)+0.140;

i=i+1;
Rnames{i} = 'M3FURAN + OH = 0.695*RO2C + 0.046*RO2XC + 0.270*HO2 + 0.684*xHO2 + 0.270*M2BUTDAL + 0.675*xHFONS + 0.010*xRCHO + 0.073*yHFONS + 0.046*zRNO3 + 0.667*yROOH';
k(:,i) = 8.729e-11;
Gstr{i,1} = 'M3FURAN'; Gstr{i,2} = 'OH';
fM3FURAN(i)=fM3FURAN(i)-1; fOH(i)=fOH(i)-1; fRO2C(i)=fRO2C(i)+0.695; fRO2XC(i)=fRO2XC(i)+0.046; fHO2(i)=fHO2(i)+0.270; fxHO2(i)=fxHO2(i)+0.684; fM2BUTDAL(i)=fM2BUTDAL(i)+0.270; fxHFONS(i)=fxHFONS(i)+0.675; fxRCHO(i)=fxRCHO(i)+0.010; fyHFONS(i)=fyHFONS(i)+0.073; fzRNO3(i)=fzRNO3(i)+0.046; fyROOH(i)=fyROOH(i)+0.667;

i=i+1;
Rnames{i} = 'M3FURAN + NO3 = 0.951*RO2C + 0.063*RO2XC + 0.924*xHO2 + 0.013*xNO2 + 0.810*xHFONS + 0.013*xRCHO + 0.100*yHFONS + 0.062*zRNO3 + 0.900*yRNO3';
k(:,i) = 1.263e-11;
Gstr{i,1} = 'M3FURAN'; Gstr{i,2} = 'NO3';
fM3FURAN(i)=fM3FURAN(i)-1; fNO3(i)=fNO3(i)-1; fRO2C(i)=fRO2C(i)+0.951; fRO2XC(i)=fRO2XC(i)+0.063; fxHO2(i)=fxHO2(i)+0.924; fxNO2(i)=fxNO2(i)+0.013; fxHFONS(i)=fxHFONS(i)+0.810; fxRCHO(i)=fxRCHO(i)+0.013; fyHFONS(i)=fyHFONS(i)+0.100; fzRNO3(i)=fzRNO3(i)+0.062; fyRNO3(i)=fyRNO3(i)+0.900;

i=i+1;
Rnames{i} = 'M3FURAN + O3 = 0.409*RCHO + 0.591*PROD2';
k(:,i) = 2.046e-17;
Gstr{i,1} = 'M3FURAN'; Gstr{i,2} = 'O3';
fM3FURAN(i)=fM3FURAN(i)-1; fO3(i)=fO3(i)-1; fRCHO(i)=fRCHO(i)+0.409; fPROD2(i)=fPROD2(i)+0.591;

i=i+1;
Rnames{i} = 'M25FUR + OH = 0.747*RO2C + 0.106*RO2XC + 0.280*HO2 + 0.116*xHO2 + 0.280*H3XE25DO + 0.498*xHFONS + 0.116*xRCHO + 0.498*xMEO2 + 0.106*zRNO3 + 0.852*yROOH';
k(:,i) = 1.260e-10;
Gstr{i,1} = 'M25FUR'; Gstr{i,2} = 'OH';
fM25FUR(i)=fM25FUR(i)-1; fOH(i)=fOH(i)-1; fRO2C(i)=fRO2C(i)+0.747; fRO2XC(i)=fRO2XC(i)+0.106; fHO2(i)=fHO2(i)+0.280; fxHO2(i)=fxHO2(i)+0.116; fH3XE25DO(i)=fH3XE25DO(i)+0.280; fxHFONS(i)=fxHFONS(i)+0.498; fxRCHO(i)=fxRCHO(i)+0.116; fxMEO2(i)=fxMEO2(i)+0.498; fzRNO3(i)=fzRNO3(i)+0.106; fyROOH(i)=fyROOH(i)+0.852;

i=i+1;
Rnames{i} = 'M25FUR + NO3 = 1.037*RO2C + 0.147*RO2XC + 0.121*xNO2 + 0.040*xOH + 0.692*xHFONS + 0.121*xRCHO + 0.692*xMEO2 + 0.124*zRNO3 + yRNO3';
k(:,i) = 5.800e-11;
Gstr{i,1} = 'M25FUR'; Gstr{i,2} = 'NO3';
fM25FUR(i)=fM25FUR(i)-1; fNO3(i)=fNO3(i)-1; fRO2C(i)=fRO2C(i)+1.037; fRO2XC(i)=fRO2XC(i)+0.147; fxNO2(i)=fxNO2(i)+0.121; fxOH(i)=fxOH(i)+0.040; fxHFONS(i)=fxHFONS(i)+0.692; fxRCHO(i)=fxRCHO(i)+0.121; fxMEO2(i)=fxMEO2(i)+0.692; fzRNO3(i)=fzRNO3(i)+0.124; fyRNO3(i)=fyRNO3(i)+1;

i=i+1;
Rnames{i} = 'M25FUR + O3 = 0.438*RO2C + 0.062*RO2XC + 1.500*OH + 0.119*xHO2 + 0.319*xOH + 0.500*OLE2 + 0.119*xCO + 0.438*xMACR + 0.062*zRNO3 + 0.037*yBACL + 0.328*yROOH + 0.319*CO2';
k(:,i) = 4.190e-16;
Gstr{i,1} = 'M25FUR'; Gstr{i,2} = 'O3';
fM25FUR(i)=fM25FUR(i)-1; fO3(i)=fO3(i)-1; fRO2C(i)=fRO2C(i)+0.438; fRO2XC(i)=fRO2XC(i)+0.062; fOH(i)=fOH(i)+1.500; fxHO2(i)=fxHO2(i)+0.119; fxOH(i)=fxOH(i)+0.319; fOLE2(i)=fOLE2(i)+0.500; fxCO(i)=fxCO(i)+0.119; fxMACR(i)=fxMACR(i)+0.438; fzRNO3(i)=fzRNO3(i)+0.062; fyBACL(i)=fyBACL(i)+0.037; fyROOH(i)=fyROOH(i)+0.328; fCO2(i)=fCO2(i)+0.319;

i=i+1;
Rnames{i} = 'BUTEDIAL + HV = 0.590*OH + 0.590*HO2 + 0.410*HFONS + 0.590*MALANHY';
k(:,i) = 2.88e-1.*JAFG1;
Gstr{i,1} = 'BUTEDIAL';
fBUTEDIAL(i)=fBUTEDIAL(i)-1; fOH(i)=fOH(i)+0.590; fHO2(i)=fHO2(i)+0.590; fHFONS(i)=fHFONS(i)+0.410; fMALANHY(i)=fMALANHY(i)+0.590;

i=i+1;
Rnames{i} = 'BUTEDIAL + OH = 0.045*RO2C + 0.001*RO2XC + 0.954*OH + 0.045*xHO2 + 0.758*MALANHY + 0.195*ROOH + 0.004*xCO + 0.041*xGLY + 0.041*xMGLY + 0.004*xROOH + 0.001*zRNO3';
k(:,i) = 4.220e-11;
Gstr{i,1} = 'BUTEDIAL'; Gstr{i,2} = 'OH';
fBUTEDIAL(i)=fBUTEDIAL(i)-1; fOH(i)=fOH(i)-1; fRO2C(i)=fRO2C(i)+0.045; fRO2XC(i)=fRO2XC(i)+0.001; fOH(i)=fOH(i)+0.954; fxHO2(i)=fxHO2(i)+0.045; fALK3(i)=fALK3(i)+0.758; fROOH(i)=fROOH(i)+0.195; fxCO(i)=fxCO(i)+0.004; fxGLY(i)=fxGLY(i)+0.041; fxMGLY(i)=fxMGLY(i)+0.041; fxROOH(i)=fxROOH(i)+0.004; fzRNO3(i)=fzRNO3(i)+0.001;

i=i+1;
Rnames{i} = 'BUTEDIAL + NO3 = 0.008*RO2C + 0.991*OH + 0.008*xHO2 + 0.955*HNO3 + 0.955*ALK3 + 0.037*RNO3 + 0.008*xCO + 0.008*xRNO3';
k(:,i) = 2.074e-15;
Gstr{i,1} = 'BUTEDIAL'; Gstr{i,2} = 'NO3';
fBUTEDIAL(i)=fBUTEDIAL(i)-1; fNO3(i)=fNO3(i)-1; fRO2C(i)=fRO2C(i)+0.008; fOH(i)=fOH(i)+0.991; fxHO2(i)=fxHO2(i)+0.008; fHNO3(i)=fHNO3(i)+0.955; fALK3(i)=fALK3(i)+0.955; fRNO3(i)=fRNO3(i)+0.037; fxCO(i)=fxCO(i)+0.008; fxRNO3(i)=fxRNO3(i)+0.008;

i=i+1;
Rnames{i} = 'BUTEDIAL + O3 = 0.667*HO2 + 0.273*HCHO + 0.334*CO + GLY + 0.394*RCOOH + 0.606*CO2';
k(:,i) = 1.600e-18;
Gstr{i,1} = 'BUTEDIAL'; Gstr{i,2} = 'O3';
fBUTEDIAL(i)=fBUTEDIAL(i)-1; fO3(i)=fO3(i)-1; fHO2(i)=fHO2(i)+0.667; fHCHO(i)=fHCHO(i)+0.273; fCO(i)=fCO(i)+0.334; fGLY(i)=fGLY(i)+1; fRCOOH(i)=fRCOOH(i)+0.394; fCO2(i)=fCO2(i)+0.606;

i=i+1;
Rnames{i} = 'M2BUTDAL + HV = 0.656*RO2C + 0.043*RO2XC + 0.699*HO2 + 0.328*xHO2 + 0.328*xAFG3 + 0.328*xMACO3 + 0.035*yAFG3 + 0.043*zHFONS + 0.035*yBACL + 0.629*yROOH';
k(:,i) = 3.42e-1.*JAFG1;
Gstr{i,1} = 'M2BUTDAL';
fM2BUTDAL(i)=fM2BUTDAL(i)-1; fRO2C(i)=fRO2C(i)+0.656; fRO2XC(i)=fRO2XC(i)+0.043; fHO2(i)=fHO2(i)+0.699; fxHO2(i)=fxHO2(i)+0.328; fxAFG3(i)=fxAFG3(i)+0.328; fxMACO3(i)=fxMACO3(i)+0.328; fyAFG3(i)=fyAFG3(i)+0.035; fzHFONS(i)=fzHFONS(i)+0.043; fyBACL(i)=fyBACL(i)+0.035; fyROOH(i)=fyROOH(i)+0.629;

i=i+1;
Rnames{i} = 'M2BUTDAL + OH = 0.055*RO2C + 0.004*RO2XC + 0.941*OH + 0.055*xHO2 + 0.250*ROOH + 0.691*AFG3 + 0.007*xCO + 0.010*xGLY + 0.078*xMGLY + 0.010*xBACL + 0.007*xROOH + 0.004*zRNO3';
k(:,i) = 4.138e-11;
Gstr{i,1} = 'M2BUTDAL'; Gstr{i,2} = 'OH';
fM2BUTDAL(i)=fM2BUTDAL(i)-1; fOH(i)=fOH(i)-1; fRO2C(i)=fRO2C(i)+0.055; fRO2XC(i)=fRO2XC(i)+0.004; fOH(i)=fOH(i)+0.941; fxHO2(i)=fxHO2(i)+0.055; fROOH(i)=fROOH(i)+0.250; fAFG3(i)=fAFG3(i)+0.691; fxCO(i)=fxCO(i)+0.007; fxGLY(i)=fxGLY(i)+0.010; fxMGLY(i)=fxMGLY(i)+0.078; fxBACL(i)=fxBACL(i)+0.010; fxROOH(i)=fxROOH(i)+0.007; fzRNO3(i)=fzRNO3(i)+0.004;

i=i+1;
Rnames{i} = 'M2BUTDAL + NO3 = 0.106*RO2C + 0.007*RO2XC + 0.887*OH + 0.106*xHO2 + 0.408*HNO3 + 0.478*RNO3 + 0.408*AFG3 + 0.106*xCO + 0.106*xRNO3 + 0.007*zRNO3';
k(:,i) = 3.126e-15;
Gstr{i,1} = 'M2BUTDAL'; Gstr{i,2} = 'NO3';
fM2BUTDAL(i)=fM2BUTDAL(i)-1; fNO3(i)=fNO3(i)-1; fRO2C(i)=fRO2C(i)+0.106; fRO2XC(i)=fRO2XC(i)+0.007; fOH(i)=fOH(i)+0.887; fxHO2(i)=fxHO2(i)+0.106; fHNO3(i)=fHNO3(i)+0.408; fRNO3(i)=fRNO3(i)+0.478; fAFG3(i)=fAFG3(i)+0.408; fxCO(i)=fxCO(i)+0.106; fxRNO3(i)=fxRNO3(i)+0.106; fzRNO3(i)=fzRNO3(i)+0.007;

i=i+1;
Rnames{i} = 'M2BUTDAL + O3 = 0.250*RO2C + 0.250*OH + 0.334*HO2 + 0.250*xHO2 + 0.136*HCHO + 0.417*CO + 0.500*GLY + 0.500*MGLY + 0.447*RCOOH + 0.250*xMGLY + 0.025*yMGLY + 0.225*yROOH + 0.303*CO2';
k(:,i) = 4.800e-18;
Gstr{i,1} = 'M2BUTDAL'; Gstr{i,2} = 'O3';
fM2BUTDAL(i)=fM2BUTDAL(i)-1; fO3(i)=fO3(i)-1; fRO2C(i)=fRO2C(i)+0.250; fOH(i)=fOH(i)+0.250; fHO2(i)=fHO2(i)+0.334; fxHO2(i)=fxHO2(i)+0.250; fHCHO(i)=fHCHO(i)+0.136; fCO(i)=fCO(i)+0.417; fGLY(i)=fGLY(i)+0.500; fMGLY(i)=fMGLY(i)+0.500; fRCOOH(i)=fRCOOH(i)+0.447; fxMGLY(i)=fxMGLY(i)+0.250; fyMGLY(i)=fyMGLY(i)+0.025; fyROOH(i)=fyROOH(i)+0.225; fCO2(i)=fCO2(i)+0.303;

i=i+1;
Rnames{i} = 'O4X2PEAL + HV = 0.699*OH + 0.699*ALK3 + 0.699*MEO2';
k(:,i) = 3.42e-1.*JAFG1;
Gstr{i,1} = 'O4X2PEAL';
fO4X2PEAL(i)=fO4X2PEAL(i)-1; fOH(i)=fOH(i)+0.699; fALK3(i)=fALK3(i)+0.699; fMEO2(i)=fMEO2(i)+0.699;

i=i+1;
Rnames{i} = 'O4X2PEAL + OH = 0.585*RO2C + 0.038*RO2XC + 0.144*xHO2 + 0.253*xOH + 0.377*MACO3 + 0.009*xCO + 0.090*xGLY + 0.190*xMGLY + 0.262*xRCHO + 0.178*xROOH + 0.178*xMECO3 + 0.038*zRNO3 + 0.052*yBACL + 0.467*yROOH + 0.253*CO2';
k(:,i) = 5.650e-11;
Gstr{i,1} = 'O4X2PEAL'; Gstr{i,2} = 'OH';
fO4X2PEAL(i)=fO4X2PEAL(i)-1; fOH(i)=fOH(i)-1; fRO2C(i)=fRO2C(i)+0.585; fRO2XC(i)=fRO2XC(i)+0.038; fxHO2(i)=fxHO2(i)+0.144; fxOH(i)=fxOH(i)+0.253; fMACO3(i)=fMACO3(i)+0.377; fxCO(i)=fxCO(i)+0.009; fxGLY(i)=fxGLY(i)+0.090; fxMGLY(i)=fxMGLY(i)+0.190; fxRCHO(i)=fxRCHO(i)+0.262; fxROOH(i)=fxROOH(i)+0.178; fxMECO3(i)=fxMECO3(i)+0.178; fzRNO3(i)=fzRNO3(i)+0.038; fyBACL(i)=fyBACL(i)+0.052; fyROOH(i)=fyROOH(i)+0.467; fCO2(i)=fCO2(i)+0.253;

i=i+1;
Rnames{i} = 'O4X2PEAL + NO3 = 0.064*RO2C + 0.005*RO2XC + 0.009*xOH + 0.931*HNO3 + 0.931*MACO3 + 0.064*xRNO3 + 0.055*xMECO3 + 0.004*zRNO3 + 0.068*yRNO3 + 0.009*CO2';
k(:,i) = 9.985e-15;
Gstr{i,1} = 'O4X2PEAL'; Gstr{i,2} = 'NO3';
fO4X2PEAL(i)=fO4X2PEAL(i)-1; fNO3(i)=fNO3(i)-1; fRO2C(i)=fRO2C(i)+0.064; fRO2XC(i)=fRO2XC(i)+0.005; fxOH(i)=fxOH(i)+0.009; fHNO3(i)=fHNO3(i)+0.931; fMACO3(i)=fMACO3(i)+0.931; fxRNO3(i)=fxRNO3(i)+0.064; fxMECO3(i)=fxMECO3(i)+0.055; fzRNO3(i)=fzRNO3(i)+0.004; fyRNO3(i)=fyRNO3(i)+0.068; fCO2(i)=fCO2(i)+0.009;

i=i+1;
Rnames{i} = 'O4X2PEAL + O3 = 0.129*RO2C + 0.129*OH + 0.495*HO2 + 0.129*xHO2 + 0.202*HCHO + 0.377*CO + 0.258*GLY + 0.742*MGLY + 0.421*RCOOH + 0.129*xMGLY + 0.013*yMGLY + 0.116*yROOH + 0.450*CO2';
k(:,i) = 4.800e-18;
Gstr{i,1} = 'O4X2PEAL'; Gstr{i,2} = 'O3';
fO4X2PEAL(i)=fO4X2PEAL(i)-1; fO3(i)=fO3(i)-1; fRO2C(i)=fRO2C(i)+0.129; fOH(i)=fOH(i)+0.129; fHO2(i)=fHO2(i)+0.495; fxHO2(i)=fxHO2(i)+0.129; fHCHO(i)=fHCHO(i)+0.202; fCO(i)=fCO(i)+0.377; fGLY(i)=fGLY(i)+0.258; fMGLY(i)=fMGLY(i)+0.742; fRCOOH(i)=fRCOOH(i)+0.421; fxMGLY(i)=fxMGLY(i)+0.129; fyMGLY(i)=fyMGLY(i)+0.013; fyROOH(i)=fyROOH(i)+0.116; fCO2(i)=fCO2(i)+0.450;

i=i+1;
Rnames{i} = 'H3XE25DO + OH = 0.930*RO2C + 0.070*RO2XC + 0.186*xHO2 + 0.372*xMGLY + 0.743*xRCHO + 0.743*xMECO3 + 0.070*zRNO3 + 0.100*yBACL + 0.900*yROOH';
k(:,i) = 4.690e-11;
Gstr{i,1} = 'H3XE25DO'; Gstr{i,2} = 'OH';
fH3XE25DO(i)=fH3XE25DO(i)-1; fOH(i)=fOH(i)-1; fRO2C(i)=fRO2C(i)+0.930; fRO2XC(i)=fRO2XC(i)+0.070; fxHO2(i)=fxHO2(i)+0.186; fxMGLY(i)=fxMGLY(i)+0.372; fxRCHO(i)=fxRCHO(i)+0.743; fxMECO3(i)=fxMECO3(i)+0.743; fzRNO3(i)=fzRNO3(i)+0.070; fyBACL(i)=fyBACL(i)+0.100; fyROOH(i)=fyROOH(i)+0.900;

i=i+1;
Rnames{i} = 'H3XE25DO + NO3 = 0.846*RO2C + 0.154*RO2XC + 0.076*HNO3 + 0.065*xHCHO + 0.782*xRNO3 + 0.782*xMECO3 + 0.065*xMACO3 + 0.008*yAFG3 + 0.154*zRNO3 + 0.069*yROOH + 0.924*yRNO3';
k(:,i) = 1.797e-18;
Gstr{i,1} = 'H3XE25DO'; Gstr{i,2} = 'NO3';
fH3XE25DO(i)=fH3XE25DO(i)-1; fNO3(i)=fNO3(i)-1; fRO2C(i)=fRO2C(i)+0.846; fRO2XC(i)=fRO2XC(i)+0.154; fHNO3(i)=fHNO3(i)+0.076; fxHCHO(i)=fxHCHO(i)+0.065; fxRNO3(i)=fxRNO3(i)+0.782; fxMECO3(i)=fxMECO3(i)+0.782; fxMACO3(i)=fxMACO3(i)+0.065; fyAFG3(i)=fyAFG3(i)+0.008; fzRNO3(i)=fzRNO3(i)+0.154; fyROOH(i)=fyROOH(i)+0.069; fyRNO3(i)=fyRNO3(i)+0.924;

i=i+1;
Rnames{i} = 'H3XE25DO + O3 = 0.500*RO2C + 0.500*OH + 0.500*xHO2 + 0.500*CO + MGLY + 0.500*RCOOH + 0.500*xMGLY + 0.050*yMGLY + 0.450*yROOH';
k(:,i) = 1.640e-14;
Gstr{i,1} = 'H3XE25DO'; Gstr{i,2} = 'O3';
fH3XE25DO(i)=fH3XE25DO(i)-1; fO3(i)=fO3(i)-1; fRO2C(i)=fRO2C(i)+0.500; fOH(i)=fOH(i)+0.500; fxHO2(i)=fxHO2(i)+0.500; fCO(i)=fCO(i)+0.500; fMGLY(i)=fMGLY(i)+1; fRCOOH(i)=fRCOOH(i)+0.500; fxMGLY(i)=fxMGLY(i)+0.500; fyMGLY(i)=fyMGLY(i)+0.050; fyROOH(i)=fyROOH(i)+0.450;

i=i+1;
Rnames{i} = 'HFONS + OH = 0.909*RO2C + 0.024*RO2XC + 0.067*HO2 + 0.903*xHO2 + 0.067*ALK3 + 0.645*xMGLY + 0.258*xRCHO + 0.006*xRCO3 + 0.024*zRNO3 + 0.061*yBACL + 0.872*yROOH';
k(:,i) = 3.828e-11;
Gstr{i,1} = 'HFONS'; Gstr{i,2} = 'OH';
fHFONS(i)=fHFONS(i)-1; fOH(i)=fOH(i)-1; fRO2C(i)=fRO2C(i)+0.909; fRO2XC(i)=fRO2XC(i)+0.024; fHO2(i)=fHO2(i)+0.067; fxHO2(i)=fxHO2(i)+0.903; fALK3(i)=fALK3(i)+0.067; fxMGLY(i)=fxMGLY(i)+0.645; fxRCHO(i)=fxRCHO(i)+0.258; fxRCO3(i)=fxRCO3(i)+0.006; fzRNO3(i)=fzRNO3(i)+0.024; fyBACL(i)=fyBACL(i)+0.061; fyROOH(i)=fyROOH(i)+0.872;

i=i+1;
Rnames{i} = 'HFONS + NO3 = 0.842*RO2C + 0.022*RO2XC + 0.136*HO2 + 0.421*xHO2 + 0.421*xNO2 + 0.136*HNO3 + 0.136*ALK3 + 0.085*xMGLY + 0.336*xRCHO + 0.421*xRNO3 + 0.022*zRNO3 + 0.864*yRNO3';
k(:,i) = 4.075e-15;
Gstr{i,1} = 'HFONS'; Gstr{i,2} = 'NO3';
fHFONS(i)=fHFONS(i)-1; fNO3(i)=fNO3(i)-1; fRO2C(i)=fRO2C(i)+0.842; fRO2XC(i)=fRO2XC(i)+0.022; fHO2(i)=fHO2(i)+0.136; fxHO2(i)=fxHO2(i)+0.421; fxNO2(i)=fxNO2(i)+0.421; fHNO3(i)=fHNO3(i)+0.136; fALK3(i)=fALK3(i)+0.136; fxMGLY(i)=fxMGLY(i)+0.085; fxRCHO(i)=fxRCHO(i)+0.336; fxRNO3(i)=fxRNO3(i)+0.421; fzRNO3(i)=fzRNO3(i)+0.022; fyRNO3(i)=fyRNO3(i)+0.864;

i=i+1;
Rnames{i} = 'HFONS + O3 = RCOOH';
k(:,i) = 5.070e-17;
Gstr{i,1} = 'HFONS'; Gstr{i,2} = 'O3'; 
fHFONS(i)=fHFONS(i)-1; fNO3(i)=fO3(i)-1; fRCOOH(i)=fRCOOH(i)+1;

