%%% This is wrapper for the production-diffusion code, applied to multiple diffusion domains. 
%%%Everything the user defines is input into this wrapper.
%%%Last modified by Marissa M. Tremblay, 2017.03.07.

%% NATACHA CHANGES: run forward sim for fixed LGM scenario
%1) Load specific sample info directly from excel file
%2) Load observation(measured) data sample info directly from excel file
%!NEED to enter manually dataObs.r!%
%3) Load Diffusion kinetics data (dataKD) directly from excel file
%-> + l.50: added function for normalisation to dataObs.r (function
%FormatNormLnD0aa)
%->l.53-56: merged all info (sample Spe, samples Obs, sample KD)
%->l.58-70: possibility to calculate internally ratio 3He and 10Be age (function
%AgeCalc)
%4) Load EDT scenarios directly from text file (from visit SUERC)and run
%simulation
%-> l.78-89, for multiple EDT scenarios, loop for automatic storage of data 
%5) Plot results simulation and observed concentrations


%%
clear all; close all;

%% 1) Provide sample specific information
%%%Production parameters

% Load sample specific info
dataSpe=readtable('./Data/SampleSpe_MBTP1.xlsx');%USER INPUT (ID, lat, long, elev, Pref, shield, thick, rho, fsp, Lsp, edot,ztoday)
dataSpe=table2struct(dataSpe); 


%% 2) Provide observational data
%%%3He retention, 3He age, 10Be age

%Load measured 3He and 10Be data
dataObs=readtable('./Data/SampleObs_MBTP1.xlsx');%USER INPUT (RetObs, RetObsUnc, HeObs, HeObsUnc, BeObs, BeObsUnc, HeConc and unc, BeConc and unc)
dataObs=table2struct(dataObs);

%Grain size on which measurement conducted
dataObs.r = 0.045;%USER INPUT; cosmogenic grain size, in cm (radius);

% %age and retention calculated outsidely or internally (1)
% AgeCalcFunc=0;

%% 3) Load Diffusion kinetics data 

% Load raw data diffusion kinetics
dataKD=readtable('./Data/KD_BestFit_MBTP9.xlsx');%USER INPUT (radius, Domin, Ea, LnDoaa, fraction)
%Normalise LnD0aa and format
dataKD = FormatNormLnD0aa(dataKD,dataObs.r);


%% merge all data ( 1),2),3)) in same struct array "data"

data = [fieldnames(dataSpe)' fieldnames(dataObs)' fieldnames(dataKD)'; struct2cell(dataSpe)' struct2cell(dataObs)' struct2cell(dataKD)'];
data=struct(data{:});

% %in case we want 10Be and 3He ages and ratio calculated internally
% %(directly from measured concentrations) using non dependent time model
% 
% if AgeCalcFunc==1;  
%     out=AgeCalc(data);
%     
%     data.HeObs=out.Age3He;
%     data.HeObsUnc=out.Age3HeUnc;
%     data.BeObs=out.Age10Be;
%     data.BeObsUnc=out.Age10BeUnc;
%     data.RetObs=out.RetCalc;
%     data.RetObsUnc=out.RetCalcUnc;
% end

%% 4) Load IsoEDT scenarios
D=importdata('FixedLGMmin15_MBTP1.txt'); %USER INPUT: iso modEDT,isopaloeEDT, chage EDT10ka sharp and since 24 ka

data.T0 = 0; %if you want to impose a surface temperature amplitude


c=2;
    for j=2:1:size(D,2);
    
        data.fname(:,1)=D(:,1);
        data.fname(:,2)=D(:,j);
   
     results = MDD_retention_EDTScenario(data);
     TotHestep(:,1)=data.fname(:,1);%time
     TotHestep(:,c)=results.TC;%IsoEDT
     TotHestep(:,c+1)=results.actTotHestep;%concentration variation
     c=c+2;
end


%% 5) Plot results and obs Conc

%1) plot IsoEDT He Conc
%make consistent plotting colors between plots
colors=colormap(parula(size(TotHestep,2)-1));
%colors=flipud(colors);

for i=2:2:size(TotHestep,2);
    time=flip(TotHestep(:,1));
     
    figure(1);
    plot(time,TotHestep(:,i),'-','LineWidth',1.5,'Color',colors(i,:)); hold on;
    xlabel('Time (yr)'); ylabel('Temperature (ºC)');
    set(gca,'xlim',[0 max(TotHestep(:,1))],'XDir','reverse');
      %set(gca,'xlim',[0 max(TotHestep(:,1))],'XDir','reverse');
    axis square
    
    figure(2);
    plot(time,TotHestep(:,i+1),'-','LineWidth',1.5,'Color',colors(i,:)); hold on;
%     plot(input.maxt-output.tv,output.totProducedConcstep,'--k','LineWidth',1.5); hold on;
    xlabel('Time (yr)'); ylabel('3He concentration');
    set(gca,'xlim',[0 max(TotHestep(:,1))],'XDir','reverse');%'ylim',[0 1e9]);
    %set(gca,'xlim',[0 max(TotHestep(:,1))],'XDir','reverse');%'ylim',[0 1e9]);
    axis square
end
   

%plot observed HeConc

expBe=(time(1)-round(data.BeObs,-3));

figure(2);
errorbar(expBe,data.HeConc,data.BeObsUnc,'horizontal','-','LineWidth',1.5,'Color','k');hold on;
errorbar(expBe,data.HeConc,data.HeConcUnc,'vertical','-','LineWidth',1.5,'Color','k');hold on;
plot(expBe,data.HeConc,'o','LineWidth',1.5,'MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',10);