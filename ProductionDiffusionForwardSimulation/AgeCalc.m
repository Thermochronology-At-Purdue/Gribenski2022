%Caculate 3He and 10Be exposure age directly from concentrations using 
%the non-time-dependent Stone scaling scheme

%Rq. need to improve error calculation for 3He and 10be ages (l. 54 and 58: for now I put
%roughly 10 % and add instrumental measureemnt error)

function out = AgeCalcFunc(data)

% Independently parameters to provide:
%parameter 10Be
Pref10Be=4.01;%Borchers et al., 2016 (Table 7, St model)
L=0.0000005;%decay constant

%Parameter 3He
Pref3He=data.Pref;
% Pref3He=118.2;%Borchers et al., 2016 (Table 7, St model)

SampleE=0.0000001;% wonder if i can do better

%% Calculation prod. rate at site
% Deal with production rate constants
% (only take in account geographic location and sample thickness)
    
pressure = NCEPatm_2(data.lat,data.long,data.elv); %looks up surface pressure and 1000 mb
    %from NCEP reanalysis. Need this for site specific production rate calc
thickcorr = thickness(data.thick,data.Lsp,data.rho); %calculate the correction needed for
    %sample thickness.
    %calculate geographic scaling factor and multiply by the reference
    %production rate, shielding correction factor, and thickness correction
    %factor to determine the local cosmogenic nulcide production rate
P010Be= Pref10Be.*data.shield.*thickcorr.*stone2000(data.lat,pressure,data.fsp);
P03He= Pref3He.*data.shield.*thickcorr.*stone2000(data.lat,pressure,data.fsp);

%% Age calculation
%non time dependent age (-> prod rate constant)

A10Be= L + data.rho * SampleE./data.Lsp;
Age10Be=(-1/A10Be)*log(1-(data.BeConc * A10Be / P010Be));
Age10BeUnc =(data.BeConcUnc./data.BeConc)*Age10Be+0.1*Age10Be; %Still to fixe that, for teh moment only instrumental unc taken in account

A3He= data.rho * SampleE./data.Lsp;
Age3He=(-1/A3He)*log(1-(data.HeConc * A3He / P03He));
Age3HeUnc =(data.HeConcUnc./data.HeConc)*Age3He+0.1*Age3He; %Still to fixe that, for teh moment only instrumental unc taken in account

RetCalc=Age3He/Age10Be;
RetCalcUnc=sqrt(((Age10BeUnc/Age10Be)^2)+((Age3HeUnc/Age3He)^2))*RetCalc;

%% export results

out.Age3He=Age3He;
out.Age3HeUnc=Age3HeUnc;
out.Age10Be=Age10Be;
out.Age10BeUnc=Age10BeUnc;
out.RetCalc=RetCalc;
out.RetCalcUnc=RetCalcUnc;



% %%
% % % Load sample specific info
% % dataSpe=readtable('SampleSpe_MBTP9.xlsx');%USER INPUT (ID, lat, long, elev, Pref, shield, thick, rho, fsp, Lsp, edot,ztoday)
% % dataSpe=table2struct(dataSpe); 
% % 
% % %Load measured 3He and 10Be data
% % dataObs=readtable('SampleObs_MBTP9AgeCalc.xlsx');%USER INPUT (RetObs, RetObsUnc, HeObs, HeObsUnc, BeObs, BeObsUnc)
% % dataObs=table2struct(dataObs);
% % 
% % data = [fieldnames(dataSpe)' fieldnames(dataObs)'; struct2cell(dataSpe)' struct2cell(dataObs)'];
% % data=struct(data{:});
% 
% %% Measured concentration and prod rate info per nuclides
% Conc10Be=data.BeConc;
% unc10Be=data.BeConcUnc;
% Pref10Be=4.7;
% L=0.0000005;%decay constant
% 
% Conc3He=data.HeConc;
% Unc3He=data.HeConcUnc;
% Pref3He=data.Pref;
% 
% SampleE=0.0000001;
% 
% %% Calculation prod. rate at site
% % Deal with production rate constants
% % (take in account geographic location and sample thickness)
%     
% pressure = NCEPatm_2(data.lat,data.long,data.elv); %looks up surface pressure and 1000 mb
%     %from NCEP reanalysis. Need this for site specific production rate calc
% thickcorr = thickness(data.thick,data.Lsp,data.rho); %calculate the correction needed for
%     %sample thickness.
%     %calculate geographic scaling factor and multiply by the reference
%     %production rate, shielding correction factor, and thickness correction
%     %factor to determine the local cosmogenic nulcide production rate
% ResultsAgeCalc.P010Be= Pref10Be.*data.shield.*thickcorr.*stone2000(data.lat,pressure,data.fsp);
% ResultsAgeCalc.P03He= Pref3He.*data.shield.*thickcorr.*stone2000(data.lat,pressure,data.fsp);
% 
% %% Age calculation
% %non time dependent age (-> prdo rate constant)
% 
% A10Be= L + data.rho * SampleE./data.Lsp;
% ResultsAgeCalc.Age10Be=(-1/A10Be)*log(1-(Conc10Be * A10Be / ResultsAgeCalc.P010Be))
% %%Age10Be_uncint = sqrt(dtdN.^2 * unc10Be.^2); %Still to fixe that
% 
% A3He= data.rho * SampleE./data.Lsp;
% ResultsAgeCalc.Age3He=(-1/A3He)*log(1-(Conc3He * A3He / ResultsAgeCalc.P03He));
% 
% ResultsAgeCalc.RetCalc=ResultsAgeCalc.Age3He/ResultsAgeCalc.Age10Be;
% 
% 
