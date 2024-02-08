% ExampleSetup_LagrangianPlume.m
% This example shows a model setup for simulation of an evolving fire plume.
% Model setup is the same as that described in Mu(:)ller et al., ACP (2016).
% Thanks to Markus Mu(:)ller for providing the model setup and data.
% Read comments in each section for a guided tour.
%
% 20160701 GMW

clear

%% DATA
load LagrangianPlumeData_matt.mat
% this structure "DAQ" contains data from an agricultural fire observed from the NASA P3 
% during the DISCOVER-AQ flight on 29 Sept 2013.
% Data are averaged into 1km bins (distance from the source).
% DAQ.TIME is a lagrangian time calculated from observed wind speed and downwind sample distance.

%% METEOROLOGY
%{
Pressure, temperature, and either RH or H2O are required Met inputs.
Dilution is calculated using the change in CO over each model step.
All calculated J-values will be scaled to J4.
%}

%kdil calculation using CO decay and inversion of dilution equation
% dX/dt = -kdil*(X - Xb)
dCOdt        = diff(DAQ.CO)./diff(DAQ.TIME); %loss rate
dCOdt(end+1) = dCOdt(end);
COmid        = (DAQ.CO + [DAQ.CO(2:end);DAQ.CO(end)])/2; %CO in middle of step
kdil         = -dCOdt./(COmid-95);

% similar calculation for tgauss
%dX/dt = (-1/(tgauss + 2t))*(X - Xb)
tgauss = 1./kdil - 2.*DAQ.TIME;
% doesn't work great (gives negative numbers, probably b/c not in center of plume)
% instead just fit it to kdil for best guess
tgauss = fminsearch(@(x) sum((kdil - 1./(x + 2*DAQ.TIME)).^2),300);

% USER: try running with either tgauss or kdil and see how results compare!

Met = {...
    'P'             DAQ.P;... %Pressure, mbar
    'T'             DAQ.T;...    %Temperature, K
    'RH'            DAQ.RH;...    %Relative Humidity, percent    
    'kdil'          kdil;...    %dilution, /s
%     'tgauss'        316; %alternative Gaussian dilution initial timescale
    'SZA'           DAQ.SZA;... %solar zenith angle
    'J4'            DAQ.JNO2;... %NO2 photolysis frequency
%     'jcorr'         'J4';... %correction factor. Muller did not use this.
    };

clear dCOdt COmid kdil tgauss

%% CHEMICAL CONCENTRATIONS
%{
The first observational point is used as initial inputs. 
All concentrations will be calculated "free running," meaning no constraints along the transect.
%}

data=readtable('/Users/samiha/Desktop/NMOG_paper/Data/F0AM_sim_initconc/neiva_S07T.xlsx');
InitConc=table2cell(data);

%% CHEMISTRY
%{
The ChemFiles input is a cell array of strings specifying functions and scripts for the chemical mechanism.
THE FIRST CELL is always a function for generic K-values.
THE SECOND CELL is always a function for J-values (photolysis frequencies).
All other inputs are scripts for mechanisms and sub-mechanisms.
Here we give an example using MCMv3.3.1.  Note that this mechanism was extracted from the MCM website for
the specific set of initial species included above.
"FURFURAL_FURAN" is a very simple set of reactions for initial oxidation of these species,
which are not included in MCM. For extra fun, try toggling this on and off to compare results.
%}

ChemFiles = {...
    'SAPRC07B_K(Met)';
    'FZS_2_Jvalues(Met)'; %Jmethod flag of 0 specifies "MCM" J-value method.
    'SAPRC07B_AllRxns';
    };

%% DILUTION CONCENTRATIONS
% Background concentrations are taken from observations just outside the plume.

BkgdConc = {...
%   names               values
    'DEFAULT'           0         ;... %0 for all zeros, 1 to use InitConc
    'CH4'               1897      ;...
    'NO'                0.031     ;...
    'NO2'               0.112     ;...
    'OH'                0         ;... %OH = 2.75 e6 cm-3
    'O3'                35.4      ;...
    'CO'                95        ;...
    'HCHO'              1.345     ;... %Formaldehyde
    'HONO'              0         ;... %HONO
    'MEOH'             1.937     ;... %Methanol
    'CCHO'             0.389     ;... %acetaldehyde
    'OLE1'             0.084     ;... %propene
    'BENZ'             0.008     ;... %benzene
    'ARO2'             0         ;... %furan; only decay!
    'ISOP'             0.219     ;... %isoprene
    'ACET'             1.491     ;... %acetone
    'RCHO'             0         ;... %propanal
    'FACD'             0         ;... %formic acid
    'AACD'             0.073     ;... %acetic acid
    'BACL'             0         ;... %2,3 butanedione
    'MVK'              0.183.*0.5;... %MVK
    'MACR'             0.183.*0.5;... %MACR
    'MGLY'             0         ;... %methylglyoxal
    'ARO1'             0         ;... %2-furfural; only decay!
    'MEK'              0         ;... %hydroxyacetone
    };

%% OPTIONS
%{
"Verbose" can be set from 0-3; this just affects the level of detail printed to the command
  window regarding model progress.
"EndPointsOnly" is set to 0 because we want output to include all concentrations along each model step.
"LinkSteps" is set to 1 because each step is connected.
"IntTime" is the integration time for each step. The average for constraints is 250s.
"SavePath" is just a filename, which will be saved by default in the \Runs\ folder.
%}

ModelOptions.Verbose = 1;         %flag for verbose output
ModelOptions.EndPointsOnly = 0;   %flag for concentration and rate outputs
ModelOptions.LinkSteps = 1;       %flag for using end-points of one run to initialize next run
ModelOptions.IntTime = 250;       %integration time for each step, seconds
ModelOptions.SavePath = 'LGPlumeResults';

%% MODEL RUN
% Now we call the model.
% Output will be saved in the "SavePath" above and will also be written to the structure S.
% Let's also throw away the inputs (don't worry, they are saved in the output structure).

S = F0AM_ModelCore(Met,InitConc,ChemFiles,BkgdConc,ModelOptions);
clear Met InitConc ChemFiles BkgdConc ModelOptions

%% FIGURES AND ANALYSIS

% calculate normalized excess mixing ratios: delta_X_dil = delta_X * delta_CO_source/delta_CO
% this is standard procedure for biomass burning plumes.
% We will also put these into a new structure compatible with F0AM plotting routines.
delta_CO = S.Conc.CO - S.BkgdConc.CO(1);
fCO = delta_CO(1)./delta_CO;
Sd.Met = S.Met; Sd.Cnames = S.Cnames; Sd.Time = S.Time;
for i=1:length(S.Cnames)
    name = S.Cnames{i};
    if isfield(S.BkgdConc,name), b = S.BkgdConc.(name)(1);
    else b = 0;
    end
    Sd.Conc.(name) = (S.Conc.(name) - b).*fCO;
end

% Next, lets look at ozone production.
figure
plot(Sd.Time,Sd.Conc.O3,'b-',...
    DAQ.TIME,(DAQ.O3-35.4).*(DAQ.CO(1)-95)./(DAQ.CO-95),'ro')
xlabel('Reaction Time (s)')
ylabel('\Delta_{dil}O_3 (ppb)')
legend('Model','Obs')
%% Import the data 
data = struct2table(Sd.Conc);
writetable(data, '/Users/samiha/Desktop/NMOG_paper/Data/F0AM_sim_result/neiva_S07T/conc.xlsx');

time_data = num2cell(S.Time);
time_data=cell2table(time_data,"VariableNames",["t"]);
writetable(time_data, '/Users/samiha/Desktop/NMOG_paper/Data/F0AM_sim_result/neiva_S07T/conc.xlsx/t.xlsx');

data = cell2table(S.Chem.Rnames, "VariableNames",["Rxn"]);
writetable(data, '/Users/samiha/Desktop/NMOG_paper/Data/F0AM_sim_result/neiva_S07T/conc.xlsx/rnames.xlsx');

data = num2cell(S.Chem.Rates);
data = cell2table(data);
writetable(data, '/Users/samiha/Desktop/NMOG_paper/Data/F0AM_sim_result/neiva_S07T/conc.xlsx/rates.xlsx');

data = num2cell(S.Chem.f);
data = cell2table(data);
writetable(data, '/Users/samiha/Desktop/NMOG_paper/Data/F0AM_sim_result/neiva_S07T/conc.xlsx/f.xlsx');

data = num2cell(S.Chem.k);
data = cell2table(data);
writetable(data, '/Users/samiha/Desktop/NMOG_paper/Data/F0AM_sim_result/neiva_S07T/conc.xlsx/k.xlsx');

data = num2cell(S.Chem.iG);
data = cell2table(data);
writetable(data, '/Users/samiha/Desktop/NMOG_paper/Data/F0AM_sim_result/neiva_S07T/conc.xlsx/iG.xlsx');


