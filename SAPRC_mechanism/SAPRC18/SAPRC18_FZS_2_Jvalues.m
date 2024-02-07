function J = SAPRC18_FZS_2_Jvalues(Met)

% INPUTS
struct2var(Met);

%%%%%INITIALIAZATION%%%%%
nj = 33; %number of rate constants
A = nan(1,nj);
B = nan(1,nj);
E = nan(1,nj);
F = nan(1,nj);
CZLOW = nan(1,nj);
Jnames = cell(nj,1);
i=0;

i=i+1;
Jnames{i} = 'JNO2_06';
A(i) = 1.237e+00;
B(i) = 5.269e-01;
E(i) = 0.000e+00;
F(i) = 2.679e+00;
CZLOW(i) = 2.419e-01;

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
F(i) = 8.124e-01;
CZLOW(i) = 2.079e-01;

i=i+1;
Jnames{i} = 'JO3O1D_06';
A(i) = 2.602e-02;
B(i) = 2.143e+00;
E(i) = 5.915e-02;
F(i) = 3.879e+00;
CZLOW(i) = 3.256e-01;

i=i+1;
Jnames{i} = 'JO3O3P_06';
A(i) = 4.622e-02;
B(i) = 2.247e-01;
E(i) = 0.000e+00;
F(i) = 5.762e-01;
CZLOW(i) = 2.079e-01;

i=i+1;
Jnames{i} = 'JHONO_06';
A(i) = 2.070e-01;
B(i) = 5.861e-01;
E(i) = 0.000e+00;
F(i) = 2.929e+00;
CZLOW(i) = 2.419e-01;

i=i+1;
Jnames{i} = 'JHNO3';
A(i) = 2.418e-04;
B(i) = 1.501e+00;
E(i) = 3.841e-02;
F(i) = 3.462e+00;
CZLOW(i) = 2.756e-01;

i=i+1;
Jnames{i} = 'JHNO4_06';
A(i) = 1.784e-03;
B(i) = 1.193e+00;
E(i) = 2.499e-02;
F(i) = 3.344e+00;
CZLOW(i) = 2.250e-01;

i=i+1;
Jnames{i} = 'JH2O2';
A(i) = 1.494e-03;
B(i) = 9.747e-01;
E(i) = 7.391e-03;
F(i) = 2.792e+00;
CZLOW(i) = 3.090e-01;

i=i+1;
Jnames{i} = 'JCOOH';
A(i) = 9.865e-04;
B(i) = 9.177e-01;
E(i) = 5.770e-03;
F(i) = 2.691e+00;
CZLOW(i) = 3.090e-01;

i=i+1;
Jnames{i} = 'JHCHOR_13';
A(i) = 7.470e-03;
B(i) = 1.033e+00;
E(i) = 3.105e-03;
F(i) = 2.992e+00;
CZLOW(i) = 3.256e-01;

i=i+1;
Jnames{i} = 'JHCHOM_13';
A(i) = 7.820e-03;
B(i) = 8.027e-01;
E(i) = 0.000e+00;
F(i) = 2.797e+00;
CZLOW(i) = 2.924e-01;

i=i+1;
Jnames{i} = 'JCCHOR_13';
A(i) = 2.313e-03;
B(i) = 1.606e+00;
E(i) = 2.990e-02;
F(i) = 3.752e+00;
CZLOW(i) = 2.756e-01;

i=i+1;
Jnames{i} = 'JC2CHO';
A(i) = 5.673e-03;
B(i) = 1.399e+00;
E(i) = 2.661e-02;
F(i) = 3.512e+00;
CZLOW(i) = 2.250e-01;

i=i+1;
Jnames{i} = 'JGLALD_14';
A(i) = 4.086e-03;
B(i) = 1.488e+00;
E(i) = 3.021e-02;
F(i) = 3.522e+00;
CZLOW(i) = 2.756e-01;

i=i+1;
Jnames{i} = 'JPAA';
A(i) = 1.874e-04;
B(i) = 1.138e+00;
E(i) = 1.714e-02;
F(i) = 2.988e+00;
CZLOW(i) = 2.924e-01;

i=i+1;
Jnames{i} = 'JGLY_I13R';
A(i) = 1.654e-02;
B(i) = 5.774e-01;
E(i) = 0.000e+00;
F(i) = 2.215e+00;
CZLOW(i) = 2.588e-01;

i=i+1;
Jnames{i} = 'JGLY_I13M';
A(i) = 4.625e-03;
B(i) = 1.027e+00;
E(i) = 1.003e-02;
F(i) = 2.812e+00;
CZLOW(i) = 3.256e-01;

i=i+1;
Jnames{i} = 'JACET_06';
A(i) = 5.139e-04;
B(i) = 2.074e+00;
E(i) = 5.751e-02;
F(i) = 4.011e+00;
CZLOW(i) = 3.090e-01;

i=i+1;
Jnames{i} = 'JMGLY_13';
A(i) = 2.686e-02;
B(i) = 4.623e-01;
E(i) = 0.000e+00;
F(i) = 2.028e+00;
CZLOW(i) = 2.250e-01;

i=i+1;
Jnames{i} = 'JBACL_11';
A(i) = 4.143e-02;
B(i) = 4.297e-01;
E(i) = 0.000e+00;
F(i) = 2.160e+00;
CZLOW(i) = 2.588e-01;

i=i+1;
Jnames{i} = 'JBALD_11';
A(i) = 6.955e-02;
B(i) = 7.304e-01;
E(i) = 0.000e+00;
F(i) = 2.715e+00;
CZLOW(i) = 2.756e-01;

i=i+1;
Jnames{i} = 'JPAN_11';
A(i) = 2.079e-04;
B(i) = 1.225e+00;
E(i) = 2.662e-02;
F(i) = 3.252e+00;
CZLOW(i) = 2.419e-01;

i=i+1;
Jnames{i} = 'JPPN_11';
A(i) = 3.764e-04;
B(i) = 1.191e+00;
E(i) = 2.039e-02;
F(i) = 3.066e+00;
CZLOW(i) = 2.756e-01;

i=i+1;
Jnames{i} = 'JACROL_16';
A(i) = 5.605e-04;
B(i) = 8.756e-01;
E(i) = 0.000e+00;
F(i) = 2.719e+00;
CZLOW(i) = 3.256e-01;

i=i+1;
Jnames{i} = 'JMEK_06';
A(i) = 4.044e-03;
B(i) = 1.433e+00;
E(i) = 2.695e-02;
F(i) = 3.543e+00;
CZLOW(i) = 2.250e-01;

i=i+1;
Jnames{i} = 'JMACR_06';
A(i) = 4.779e-04;
B(i) = 8.860e-01;
E(i) = 0.000e+00;
F(i) = 2.698e+00;
CZLOW(i) = 3.256e-01;

i=i+1;
Jnames{i} = 'JMVK_16';
A(i) = 7.785e-04;
B(i) = 8.907e-01;
E(i) = 1.140e-03;
F(i) = 2.696e+00;
CZLOW(i) = 3.256e-01;

i=i+1;
Jnames{i} = 'JAFGS';
A(i) = 7.667e-01;
B(i) = 6.813e-01;
E(i) = 0.000e+00;
F(i) = 2.241e+00;
CZLOW(i) = 3.090e-01;

i=i+1;
Jnames{i} = 'JIC3ONO2';
A(i) = 9.246e-04;
B(i) = 1.370e+00;
E(i) = 2.754e-02;
F(i) = 3.470e+00;
CZLOW(i) = 2.419e-01;

i=i+1;
Jnames{i} = 'JCRBNIT';
A(i) = 9.777e-03;
B(i) = 1.173e+00;
E(i) = 1.780e-02;
F(i) = 3.031e+00;
CZLOW(i) = 2.924e-01;

i=i+1;
Jnames{i} = 'JDIONO2';
A(i) = 2.117e-03;
B(i) = 1.101e+00;
E(i) = 1.443e-02;
F(i) = 2.951e+00;
CZLOW(i) = 3.090e-01;

i=i+1;
Jnames{i} = 'JHPALDS';
A(i) = 7.956e-02;
B(i) = 6.929e-01;
E(i) = 0.000e+00;
F(i) = 2.691e+00;
CZLOW(i) = 2.588e-01;

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

%End
