%%
%This is wrapper for the production-diffusion code, applied to multiple diffusion domains. 
%Everything the user defines is input into this wrapper.

%Written by Marissa Tremblay.
%Contact: tremblam@purdue.edu
%Last modified: 2021.12.28

%Copyright 2021, Marissa Tremblay
%All rights reserved
%Developed in part with funding from the National Science Foundation and
%the American Association for the Advancement of Science.

%This program is free software: you can redistribute it and/or modify
%it under the terms of the GNU General Public License as published by
%the Free Software Foundation, either version 3 of the License, or
%(at your option) any later version.

%This program is distributed in the hope that it will be useful,
%but WITHOUT ANY WARRANTY; without even the implied warranty of
%MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%GNU General Public License for more details.

%You should have received a copy of the GNU General Public License
%along with this program.  If not, see <https://www.gnu.org/licenses/>.


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

%% 3) Load Diffusion kinetics data 

% Load raw data diffusion kinetics
dataKD=readtable('./Data/KD_BestFit_MBTP9.xlsx');%USER INPUT (radius, Domin, Ea, LnDoaa, fraction)
%Normalise LnD0aa and format
dataKD = FormatNormLnD0aa(dataKD,dataObs.r);


%% merge all data from parts 1?3 in same structure called "data"

data = [fieldnames(dataSpe)' fieldnames(dataObs)' fieldnames(dataKD)'; struct2cell(dataSpe)' struct2cell(dataObs)' struct2cell(dataKD)'];
data=struct(data{:});

%% 4) Load EDT scenarios
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


%% 5) Plot results versus observed concentration

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