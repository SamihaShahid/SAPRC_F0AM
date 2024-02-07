function J = FZS_2_Jvalues(Met)

% INPUTS
struct2var(Met);

%%%%%INITIALIAZATION%%%%%
nj = 48; %number of rate constants
A = nan(1,nj);
B = nan(1,nj);
E = nan(1,nj);
F = nan(1,nj);
CZLOW = nan(1,nj);
Jnames = cell(nj,1);
i=0;

i=i+1;
Jnames{i} = 'JACET_06';
A(i) = 5.139e-04;
B(i) = 2.073e+00;
E(i) = 5.748e-02;
F(i) = 3.941e+00;
CZLOW(i) = 3.256e-01;

i=i+1;
Jnames{i} = 'JACRO_09';
A(i) = 4.151e-04;
B(i) = 8.526e-01;
E(i) = 0.000e+00;
F(i) = 2.665e+00;
CZLOW(i) = 3.256e-01;

i=i+1;
Jnames{i} = 'JAFG1';
A(i) = 7.667e-01;
B(i) = 6.813e-01;
E(i) = 0.000e+00;
F(i) = 2.241e+00;
CZLOW(i) = 3.090e-01;

i=i+1;
Jnames{i} = 'JBACL_07';
A(i) = 4.190e-02;
B(i) = 4.383e-01;
E(i) = 0.000e+00;
F(i) = 2.151e+00;
CZLOW(i) = 2.588e-01;

i=i+1;
Jnames{i} = 'JBALD_06';
A(i) = 1.043e-01;
B(i) = 7.105e-01;
E(i) = 0.000e+00;
F(i) = 2.724e+00;
CZLOW(i) = 2.588e-01;

i=i+1;
Jnames{i} = 'JBRNO2';
A(i) = 6.433e-01;
B(i) = 4.476e-01;
E(i) = 0.000e+00;
F(i) = 1.999e+00;
CZLOW(i) = 2.079e-01;

i=i+1;
Jnames{i} = 'JBRONO2';
A(i) = 1.903e-01;
B(i) = 5.343e-01;
E(i) = 0.000e+00;
F(i) = 1.915e+00;
CZLOW(i) = 2.756e-01;

i=i+1;
Jnames{i} = 'JC2CHO';
A(i) = 5.673e-03;
B(i) = 1.399e+00;
E(i) = 2.661e-02;
F(i) = 3.512e+00;
CZLOW(i) = 2.250e-01;

i=i+1;
Jnames{i} = 'JCCHO_R';
A(i) = 2.047e-03;
B(i) = 1.595e+00;
E(i) = 2.917e-02;
F(i) = 3.771e+00;
CZLOW(i) = 2.588e-01;

i=i+1;
Jnames{i} = 'JCH3I';
A(i) = 1.953e-03;
B(i) = 1.331e+00;
E(i) = 3.288e-02;
F(i) = 3.346e+00;
CZLOW(i) = 2.588e-01;

i=i+1;
Jnames{i} = 'JCL2';
A(i) = 3.283e-01;
B(i) = 6.319e-01;
E(i) = 0.000e+00;
F(i) = 2.381e+00;
CZLOW(i) = 2.588e-01;

i=i+1;
Jnames{i} = 'JCLACET';
A(i) = 2.086e-02;
B(i) = 1.124e+00;
E(i) = 1.203e-02;
F(i) = 3.014e+00;
CZLOW(i) = 3.090e-01;

i=i+1;
Jnames{i} = 'JCLCCHO';
A(i) = 1.745e-02;
B(i) = 1.012e+00;
E(i) = 5.044e-03;
F(i) = 2.900e+00;
CZLOW(i) = 3.256e-01;

i=i+1;
Jnames{i} = 'JCLNO_06';
A(i) = 3.933e-01;
B(i) = 4.610e-01;
E(i) = 0.000e+00;
F(i) = 1.640e+00;
CZLOW(i) = 2.250e-01;

i=i+1;
Jnames{i} = 'JCLNO2';
A(i) = 6.271e-02;
B(i) = 8.338e-01;
E(i) = 0.000e+00;
F(i) = 2.626e+00;
CZLOW(i) = 3.256e-01;

i=i+1;
Jnames{i} = 'JCLONO';
A(i) = 6.789e-01;
B(i) = 7.798e-01;
E(i) = 0.000e+00;
F(i) = 2.573e+00;
CZLOW(i) = 3.090e-01;

i=i+1;
Jnames{i} = 'JCLONO2_1';
A(i) = 1.947e-03;
B(i) = 1.112e+00;
E(i) = 1.817e-02;
F(i) = 2.918e+00;
CZLOW(i) = 2.756e-01;

i=i+1;
Jnames{i} = 'JCLONO2_2';
A(i) = 6.201e-03;
B(i) = 6.701e-01;
E(i) = 0.000e+00;
F(i) = 2.202e+00;
CZLOW(i) = 3.090e-01;

i=i+1;
Jnames{i} = 'JCLPICERI';
A(i) = 1.145e-02;
B(i) = 9.073e-01;
E(i) = 9.706e-03;
F(i) = 2.617e+00;
CZLOW(i) = 2.924e-01;

i=i+1;
Jnames{i} = 'JCONO';
A(i) = 2.975e-01;
B(i) = 6.227e-01;
E(i) = 0.000e+00;
F(i) = 3.009e+00;
CZLOW(i) = 2.079e-01;

i=i+1;
Jnames{i} = 'JCOOH';
A(i) = 9.865e-04;
B(i) = 9.177e-01;
E(i) = 5.770e-03;
F(i) = 2.691e+00;
CZLOW(i) = 3.090e-01;

i=i+1;
Jnames{i} = 'JCS2';
A(i) = 2.947e-02;
B(i) = 9.532e-01;
E(i) = 0.000e+00;
F(i) = 2.901e+00;
CZLOW(i) = 3.256e-01;

i=i+1;
Jnames{i} = 'JGLY_07M';
A(i) = 7.679e-03;
B(i) = 8.815e-01;
E(i) = 9.343e-03;
F(i) = 2.440e+00;
CZLOW(i) = 3.256e-01;

i=i+1;
Jnames{i} = 'JGLY_07R';
A(i) = 1.529e-02;
B(i) = 5.173e-01;
E(i) = 0.000e+00;
F(i) = 1.998e+00;
CZLOW(i) = 2.588e-01;

i=i+1;
Jnames{i} = 'JH2O2';
A(i) = 1.494e-03;
B(i) = 9.747e-01;
E(i) = 7.387e-03;
F(i) = 2.792e+00;
CZLOW(i) = 3.090e-01;

i=i+1;
Jnames{i} = 'JHCHOM_06';
A(i) = 7.064e-03;
B(i) = 8.130e-01;
E(i) = 0.000e+00;
F(i) = 2.804e+00;
CZLOW(i) = 2.924e-01;

i=i+1;
Jnames{i} = 'JHCHOR_06';
A(i) = 7.699e-03;
B(i) = 1.024e+00;
E(i) = 2.869e-03;
F(i) = 2.976e+00;
CZLOW(i) = 3.256e-01;

i=i+1;
Jnames{i} = 'JHI';
A(i) = 3.025e-03;
B(i) = 1.228e+00;
E(i) = 2.196e-02;
F(i) = 3.116e+00;
CZLOW(i) = 2.588e-01;

i=i+1;
Jnames{i} = 'JHNO3';
A(i) = 2.418e-04;
B(i) = 1.501e+00;
E(i) = 3.842e-02;
F(i) = 3.462e+00;
CZLOW(i) = 2.756e-01;

i=i+1;
Jnames{i} = 'JHNO4_06';
A(i) = 1.784e-03;
B(i) = 1.193e+00;
E(i) = 2.499e-02;
F(i) = 3.343e+00;
CZLOW(i) = 2.250e-01;

i=i+1;
Jnames{i} = 'JHOCCHO';
A(i) = 2.453e-03;
B(i) = 1.450e+00;
E(i) = 2.930e-02;
F(i) = 3.493e+00;
CZLOW(i) = 2.588e-01;

i=i+1;
Jnames{i} = 'JHOCL_06';
A(i) = 4.239e-02;
B(i) = 7.159e-01;
E(i) = 0.000e+00;
F(i) = 2.397e+00;
CZLOW(i) = 2.924e-01;

i=i+1;
Jnames{i} = 'JHOI';
A(i) = 1.159e+01;
B(i) = 4.536e-01;
E(i) = 0.000e+00;
F(i) = 2.114e+00;
CZLOW(i) = 2.588e-01;

i=i+1;
Jnames{i} = 'JHONO_06';
A(i) = 2.070e-01;
B(i) = 5.861e-01;
E(i) = 0.000e+00;
F(i) = 2.929e+00;
CZLOW(i) = 2.419e-01;

i=i+1;
Jnames{i} = 'JIC3ONO2';
A(i) = 9.246e-04;
B(i) = 1.370e+00;
E(i) = 2.754e-02;
F(i) = 3.470e+00;
CZLOW(i) = 2.419e-01;

i=i+1;
Jnames{i} = 'JMACR_06';
A(i) = 4.779e-04;
B(i) = 8.860e-01;
E(i) = 0.000e+00;
F(i) = 2.698e+00;
CZLOW(i) = 3.256e-01;

i=i+1;
Jnames{i} = 'JMEK_06';
A(i) = 4.044e-03;
B(i) = 1.433e+00;
E(i) = 2.695e-02;
F(i) = 3.543e+00;
CZLOW(i) = 2.250e-01;

i=i+1;
Jnames{i} = 'JMGLY_06';
A(i) = 2.510e-02;
B(i) = 4.662e-01;
E(i) = 0.000e+00;
F(i) = 2.012e+00;
CZLOW(i) = 2.250e-01;

i=i+1;
Jnames{i} = 'JMITC';
A(i) = 1.954e-03;
B(i) = 1.313e+00;
E(i) = 2.305e-02;
F(i) = 3.254e+00;
CZLOW(i) = 2.756e-01;

i=i+1;
Jnames{i} = 'JMVK_06';
A(i) = 1.834e-04;
B(i) = 8.936e-01;
E(i) = 0.000e+00;
F(i) = 2.658e+00;
CZLOW(i) = 2.924e-01;

i=i+1;
Jnames{i} = 'JNO2_06';
A(i) = 1.237e+00;
B(i) = 5.269e-01;
E(i) = 0.000e+00;
F(i) = 2.679e+00;
CZLOW(i) = 2.419e-01;

i=i+1;
Jnames{i} = 'JNO2ex';
A(i) = 3.023e+00;
B(i) = 2.912e-01;
E(i) = 0.000e+00;
F(i) = 1.212e+00;
CZLOW(i) = 2.079e-01;

i=i+1;
Jnames{i} = 'JNO3NO_06';
A(i) = 2.349e+00;
B(i) = 1.870e-01;
E(i) = 0.000e+00;
F(i) = 6.208e-01;
CZLOW(i) = 2.079e-01;

i=i+1;
Jnames{i} = 'JNO3NO2_6';
A(i) = 1.956e+01;
B(i) = 2.195e-01;
E(i) = 0.000e+00;
F(i) = 8.125e-01;
CZLOW(i) = 2.079e-01;

i=i+1;
Jnames{i} = 'JO3O1D_06';
A(i) = 2.602e-02;
B(i) = 2.143e+00;
E(i) = 5.916e-02;
F(i) = 4.011e+00;
CZLOW(i) = 2.924e-01;

i=i+1;
Jnames{i} = 'JO3O3P_06';
A(i) = 4.622e-02;
B(i) = 2.247e-01;
E(i) = 0.000e+00;
F(i) = 5.762e-01;
CZLOW(i) = 2.079e-01;

i=i+1;
Jnames{i} = 'JPAA';
A(i) = 1.874e-04;
B(i) = 1.138e+00;
E(i) = 1.714e-02;
F(i) = 2.988e+00;
CZLOW(i) = 2.924e-01;

i=i+1;
Jnames{i} = 'JPAN';
A(i) = 2.079e-04;
B(i) = 1.225e+00;
E(i) = 2.662e-02;
F(i) = 3.252e+00;
CZLOW(i) = 2.419e-01;

%%%%%CALCULATE J-VALUES%%%%%
ns = length(SZA);
SZAbig = repmat(SZA,1,nj);
Abig = repmat(A,ns,1);
Bbig = repmat(B,ns,1);
Ebig = repmat(E,ns,1);
Fbig = repmat(F,ns,1);
CZLOWbig = repmat(CZLOW,ns,1);

% Better version
% j = Abig.*exp(-Bbig.*sqrt((1+Ebig)./(cosd(SZAbig).^2+Ebig)))./60.; %
%
j = nan(ns,nj);
for ik=1:length(Jnames)
    for jk=1:length(SZA)
        if (cosd(SZA(jk))>=CZLOW(jk))
            j(jk,ik) = A(ik).*exp(-B(ik).*sqrt((1+E(ik))./(cosd(SZA(jk)).^2+E(ik))))./60.;
        elseif ((cosd(SZA(jk))<CZLOW(jk)) && (cosd(SZA(jk))>0.))
            j(jk,ik) = A(jk).*exp(-B(jk).*sqrt((1+E(jk))./(CZLOW(jk).^2+E(jk)))).*(cosd(SZA(jk))./CZLOW(jk)).^F(jk)./60.;
        else
            disp('Y')
            j(jk,ik) = 0.0;
        end
    end
end
warning on
j(SZA>=90,:) = 0;

% accumulate
J = struct;
for i=1:length(Jnames)
    J.(Jnames{i}) = j(:,i);
end

end



