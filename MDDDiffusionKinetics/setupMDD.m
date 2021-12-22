
%This code is used to fit an MDD model to stepwise degassing diffusion
%experiment data. It is currently set up for only one isotope. The number
%of domains is allowed to vary. The activation energy is assumed to be the
%same across domains while the pre-exponential factor (D0/a^2) and the
%fraction of gas in each domain varies. Needs the companion functions
%D0calc_MonteCarloErros.m and TremblayMDD.m.

%Written by Marissa Tremblay. 
%Contact: tremblam@purdue.edu
%Last modified 2021.12.17

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

%Load data from diffusion experiment and calculate diffusivities

%indicate the experiment of interest here. changed on 2018.05.02 to open a
%file selection dialog box, such that code is not modified every time a new
%dataset is opened
[filename,filepath] = uigetfile('*.txt');
expdata = importdata(fullfile(filepath,filename));
namesplit = strsplit(filename,'.'); sample = namesplit{1};

%specify what diffusion geometry is appropriate to use for this experiment.
%with the geometry flag. 0 = spherical, 1 = plane sheet.
geomflag = 0;

%open a log file in which to write the best fit parameters and
%corresponding residual from each model run
if exist('../MDD_results','dir')~=7
    mkdir ../MDD_results;
end

logfilename = strcat('../MDD_results/',sample,'_MDDresults.txt');
logfile = fopen(logfilename,'a+');
keeptrack = load(logfilename); runnum = size(keeptrack,1) + 1;

%specify plotting colors for MDD model, if desired
facecolor = [1,102/255,102/255];
linecolor =  [1,0,0];

%calculate diffusivities from the diffusion experiment results, save for
%later

diffresults = D0calc_MonteCarloErrors(expdata,geomflag);

%% Define the model inputs

%experiment parameters
input.TC = expdata(:,2);         %temperature from heating step in degrees Celsius
input.thr = expdata(:,3);        %duration of heating step in hours

%modification 2016.11.01. Adjust D/a^2 values determined from diffusion
%experiment to grain size analyzed in cosmogenic measurements
cosmosize = 1; %cm %nat comment, here it doesn't chnage anything teh actual size of our grain used for diff (as long s we keep teh same size for cosmo, but it matters in froward simulations...anyway these two lines can be left at 1 to anlaysi diff kinetics measurements
diffsize = 1; %cm
input.lnDaa = log(exp(diffresults.lnDaa).*(diffsize.^2)./(cosmosize).^2);
%input.lnDaa = diffresults.lnDaa; %D/a^2 values determined frome experiment
input.Fi = diffresults.Fi;

%MDD paramaters - to be changed in each simulation
numdomains = 3; %number of domains 

%range of activation energies to be explored, in kJ/mol. for these models,
%we hold Ea constant across domains
Ea =85.9;%linspace(97.5,98.5,11);David sugesst we keep that fixed (derived from D0_oneIsotop_MC code)

%range of ln(D0/a^2) values to be explored, in ln(s^-1). assign a different
%range for each domain (corresponding to different intercepts). intercepts
%are set to zero for domains not included in the simulation (i.e., for
%domains 5 and 6 if running a four-domain model). the gas fractions for
%the domains not included will also be set to zero below.
lnD0aa1 = linspace(9,11,10);%linspace(10.2,10.7,6);
if numdomains > 1
    lnD0aa2 = linspace(14,16,10);%linspace(12.2,12.7,6);
else
    lnD0aa2 = 0;
end
if numdomains > 2
    lnD0aa3 = linspace(6,8,10);%linspace(15.5,16,6);
else
    lnD0aa3 = 0;
end
if numdomains > 3
    lnD0aa4 = linspace(20,25,3);
else
    lnD0aa4 = 0;
end;
if numdomains > 4
    lnD0aa5 = linspace(0,10,6);
else
    lnD0aa5 = 0;
end
if numdomains > 5
    lnD0aa6 = linspace(0,10,6);
else
    lnD0aa6 = 0;
end

df = 0.01; %increments between gas fractions to be tested
minf = 0.02; %0.05 %minimum gas fraction in any domain %David said this can be tuned to gain time
space = int16(1+(1-minf)/df);

%assign a range of gas fractions to be explored for each domain. if a given
%domain is not included in a simulation (e.g., domains 5 and 6 in a
%4-domain model), the gas fraction will simply be set to zero. 
if numdomains < 2
    fdom1 = ones(1,space);
else
fdom1 = linspace(minf,1,space);
end
if numdomains > 1
   fdom2 = fdom1-df;
else
    fdom2 = 0;
end
if numdomains > 2;
    fdom3 = fdom1-2*df;
    fdom3 = fdom3(2:end);
else
    fdom3 = 0;
end
if numdomains > 3;
    fdom4 = fdom1-3*df;
    fdom4 = fdom4(3:end);
else
    fdom4 = 0;
end
if numdomains > 4;
    fdom5 = fdom1-4*df;
    fdom5 = fdom5(4:end);
else
    fdom5 = 0;
end
if numdomains > 5;
    fdom6 = fdom1-5*df;
    fdom6 = fdom6(5:end);
else
    fdom6 = 0;
end

%this massive for if loop ensures that each combination of gas fractions
%assigned to the various domains is explored. 
index = 1;
if numdomains > 1
    for i = 1:length(fdom1)
        for j = 1:length(fdom2)
            for k = 1:length(fdom3)
                for l = 1:length(fdom4)
                    for m = 1:length(fdom5)
                        for n = 1:length(fdom6)
                            if  fdom1(i) + fdom2(j) + fdom3(k) + ...
                                    fdom4(l) +fdom5(m) + fdom6(n) == 1;
                                if fdom1(i) >=fdom2(j);
                                    if fdom2(j) >= fdom3(k);
                                        if fdom3(k) >= fdom4(l);
                                            if fdom4(l) >= fdom5(m);
                                                if fdom5(m) >= fdom6(n);
                                                    fracs(index,1) = fdom1(i);
                                                    fracs(index,2) = fdom2(j);
                                                    fracs(index,3) = fdom3(k);
                                                    fracs(index,4) = fdom4(l);
                                                    fracs(index,5) = fdom5(m);
                                                    fracs(index,6) = fdom6(n);
                                                    index = index + 1;
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
else
    for i = 1:length(fdom1)
        fracs(index,1) = fdom1(i);
        fracs(index,2) = 0;
        fracs(index,3) = 0;
        fracs(index,4) = 0;
        fracs(index,5) = 0;
        fracs(index,6) = 0;
        index = index + 1;
    end
end

%% Massive MDD model. Be patient!

%allocate matrix of diffusion parameters applied in each iteration
allparams = zeros(length(Ea)*((length(lnD0aa1))^numdomains)...
    *size(fracs,1),13);

%tell the user how many iterations this simulation containts
disp('Number of iterations this simulation:'); disp(size(allparams,1));

%initialize the residuals vector
residuals = zeros(size(allparams,1),1);

%loop through all iterations with the MDD model
index2 = 1;
for a = 1:length(fracs)
     for b = 1:length(Ea)
       for c = 1:length(lnD0aa1)
           for d = 1:length(lnD0aa2)
               for f = 1:length(lnD0aa3)
                   for g = 1:length(lnD0aa4)
                       for h = 1:length(lnD0aa5)
                           for i = 1:length(lnD0aa6)
                               input.fracs = fracs(a,:);
                               input.Ea = [Ea(b) Ea(b) Ea(b) Ea(b)...
                                   Ea(b) Ea(b)];
                               input.lnD0aa = [lnD0aa1(c) lnD0aa2(d) ...
                                   lnD0aa3(f) lnD0aa4(g) lnD0aa5(h) lnD0aa6(i)];
                               
                               MDD = TremblayMDD(input,geomflag);
                               
                               allparams(index2,1) = input.Ea(1);
                               allparams(index2,2:7) = input.lnD0aa;
                               allparams(index2,8:13) = input.fracs;
                               residuals(index2,1) = MDD.residual;
                               index2 = index2 + 1;
%%turn on this if statement if you want to print the current iteration                               
%                                if rem(index2,10) == 0
%                                   disp('iteration #:'); disp(index2);
%                                end
                           end
                       end
                   end
               end
           end
       end
     end
end
%%
input2.TC = input.TC;        %temperature from heating step in degrees Celsius
input2.thr = input.thr;      %duration of heating step in hours
input2.lnDaa = input.lnDaa;  %D/a^2 values determined frome experiment
input2.Fi = input.Fi;

%get uncertainties in diffusivity for each step of the diffusion experiment 
%from the Monte Carlo simulation
%modified on 2016.11.01 to be adjusted for the cosmogenic grain size
%analyzed
lnDaaneg = real(diffresults.lnDaaneg).*log((diffsize.^2)./(cosmosize.^2)); 
lnDaapos = real(diffresults.lnDaapos).*log((diffsize.^2)./(cosmosize.^2));


%now find the set of parameters that minimized the residual in the MDD
%search over all parameters
[val, run] = min(residuals);

%assign these parameters that minimize the residual to be the new
%parameters
newparams = allparams(run,:);

%write these parameters and the corresponding residual to the log file
writetofile = [runnum newparams val];
fprintf(logfile,'%g %3g %3g %3g %3g %3g %3g %3g %3g %3g %3g %3g %3g %3g %4g \n',writetofile);

 
%redo the MDD model for the residual-minimizing parameters for plotting
%purposes
input2.Ea = [newparams(1) newparams(1) newparams(1) newparams(1) ...
    newparams(1) newparams(1)];
input2.lnD0aa = newparams(2:7);
input2.fracs = newparams(8:13);
MDD2 = TremblayMDD(input2,geomflag);

%calculate variables x and y to plot lines corresponding to each domain
R = 83.14; %gas constant 
m = 1*newparams(1)/R;
xplot = linspace(5,32,10);
y1 = -m.*xplot + newparams(2); y2 = -m.*xplot + newparams(3);
y3 = -m.*xplot + newparams(4); y4 = -m.*xplot + newparams(5);
y5 = -m.*xplot + newparams(6); y6 = -m.*xplot + newparams(7);

%begin plotting. first do the Arrhenius figure
figure;
%plot lines corresponding to each domain
plot(xplot,y1,'Color',linecolor,'LineWidth',2); hold on;
if numdomains > 1
    plot(xplot,y2,'Color',linecolor,'LineWidth',2); hold on;
end
if numdomains > 2
    plot(xplot,y3,'Color',linecolor,'LineWidth',2); hold on;
end
if numdomains > 3
    plot(xplot,y4,'Color',linecolor,'LineWidth',2); hold on;
end
if numdomains > 4
    plot(xplot,y5,'Color',linecolor,'LineWidth',2); hold on;
end
if numdomains > 5
    plot(xplot,y6,'Color',linecolor,'LineWidth',2); hold on;
end
%plot uncertainties in diffusivities from the diffusion experiment results
for a = 1:length(diffresults.Tplot);
    xx = [diffresults.Tplot(a) diffresults.Tplot(a)];
    yy = [lnDaaneg(a) lnDaapos(a)];
    plot(xx,yy,'k'); hold on;
end;
%plot diffusivities from diffusion experiment results
plot(diffresults.Tplot,input.lnDaa,'ko','MarkerSize',10,...
     'MarkerFaceColor',[0.88 0.88 0.88]); hold on;
%plot MDD model diffusivities if the number of domains is > 1
%if numdomains > 1
    plot(diffresults.Tplot,MDD2.lnDaa_MDD,'ko','MarkerSize',5,...
        'MarkerFaceColor',facecolor,'Color',linecolor); hold on;
%end
%plot temperatures on top axis in linear space
Tempplot = [100 200 300 400 500 600 800 1000];
Tempplot = (10^4)./(Tempplot + 273.15);
Dplot = -5.*ones(size(Tempplot));
plot(Tempplot,Dplot,'x'); hold on;
set(gca,'xlim',[5 32],'ylim',[-25 -5]);
xlabel('10000/T (K^-^1)');
ylabel('ln(D/a^2)');
plotlabel = strcat(sample,', run #:  ',num2str(runnum));
title(plotlabel);
axis square;

%next plot the residuals with respect to temperature and cumulative gas
%release fraction. use the largest domain (domain 1) for the fit
residplot = input.lnDaa - (-m.*diffresults.Tplot + newparams(2));
modelresidplot = MDD2.lnDaa_MDD - (-m.*diffresults.Tplot + newparams(2));

%first, release fraction
figure;
for a = 1:length(lnDaapos)
    yypos(a) = lnDaapos(a) - (-m.*diffresults.Tplot(a) + newparams(2));
    yyneg(a) = lnDaaneg(a) - (-m.*diffresults.Tplot(a) + newparams(2));
    xx2 = [diffresults.Fi(a) diffresults.Fi(a)];
    yy2 = [yypos(a) yyneg(a)];
    plot(xx2,yy2,'k'); hold on;
end
plot(diffresults.Fi,residplot,'ko','MarkerSize',10,'MarkerFaceColor',...
    [0.88 0.88 0.88]); hold on;
plot(MDD2.f_MDD,modelresidplot,'ko','MarkerSize',5,'MarkerFaceColor',...
    facecolor,'Color',linecolor); hold on;
set(gca,'xlim',[-0.05 1.05],'ylim',[-14 6]);
xlabel('cumulative release fraction'); ylabel('residual');
title(plotlabel);
axis square;

%now temperature
figure;
for a = 1:length(lnDaapos)
    xx3 = [input.TC(a) input.TC(a)];
    yy3 = [yypos(a) yyneg(a)];
    plot(xx3,yy3,'k'); hold on;
end
plot(input.TC,residplot,'ko','MarkerSize',10,'MarkerFaceColor',...
    [0.88 0.88 0.88]); hold on;
plot(input.TC,modelresidplot,'ko','MarkerSize',5,'MarkerFaceColor',...
    facecolor,'Color',linecolor); hold on;
set(gca,'xlim',[0 500],'ylim',[-14 6]);
xlabel('temperature (ºC)'); ylabel('residual');
title(plotlabel);
axis square;

%plot temperature history of experiment
figure;
step = linspace(1,length(input.TC),length(input.TC));
plot(step,input.TC,'ko','MarkerSize',10,'MarkerFaceColor',[0.88 0.88 0.88]); hold on;
set(gca,'xlim',[0 max(step)+1],'ylim',[min(input.TC)-25 max(input.TC)+25]);
xlabel('step'); ylabel('temperature (ºC)');
axis square;

%Add NAt to assess:
fprintf('%g %3g %3g %3g %3g %3g %3g %3g %3g %3g %3g %3g %3g %3g %4g \n',writetofile);