%Setup to do production and diffusion with MDD model.
%Takes user defined inputs from the wrapper code and packages them to be
%passed to the production-diffusion code.

%Written by Marissa Tremblay.
%Contact: tremblam@purdue.edu
%Last modified: 2021.12.17

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

function results = MDD_retention_EDTScenario(data)

input.lat = data.lat; %sample latitude in decimal degrees
input.long = data.long; %sample longitude in decimal degrees; west = negative
input.elv = data.elv; %sample elevation in meters
input.Pref = data.Pref; %sea level high latitude production rate in atoms/g/yr
input.shield = data.shield; %correction factor for shielding
input.thick = data.thick; %sample thickness in cm
input.rho = data.rho; %mineral density in g/cm3
input.fsp = data.fsp; %fraction of spallogenic cosmogenic nuclide production
input.Lsp = data.Lsp; %effective attenuation length for spallation in g/cm^2
input.ee = data.edot.*(1e-6).*100.*data.rho; %erosion rate in g/cm^2/yr
input.ztoday = data.ztoday.*data.rho; %if sample is not at surface today, depth in g/cm^2
input.r = data.r; %average radius of fragments analyzed, in cm 

RetObs = data.RetObs; %observed noble gas retention
RetObsUnc = data.RetObsUnc; %1sd in observed noble gas average retention
HeObs = data.HeObs; %apparent 3He exposure age in years
HeObsUnc = data.HeObsUnc; %1sd in apparent 3He exposure age
BeObs = data.BeObs; %observed 10Be exposure age in years
BeObsUnc = data.BeObsUnc; %1sd in 10Be exposure age

input.a = data.a; % radius of proton-irradiated grain analyzed (cm) 
input.ndom = data.ndom; %number of domains from the best-fit MDD model
input.Ea = ones(1,input.ndom).*data.Ea; %vector of activation energies for each domain.
input.lnD0aa = data.lnD0aa; %vector of ln(D0/a^2) for each domain. technically unitless (normalized to 1/s)
input.fracs = data.fracs; %fraction of gas apportioned to each diffusion domain

input.n = 512; %number of mesh points over which production-diffusion calculation takes place

Temphist = data.fname; %load(data.fname); %load in the user-defined time-temperature histories
times = Temphist(:,1); %the first column in this file should correspond to the time steps in years
temps = Temphist(:,2:end);
input.maxt = max(times); %total duration of the time-temperature history, in years
input.Tt = times';%0:input.dt:input.maxt; 
input.dt = input.maxt./(length(input.Tt)-1); %size of the time step, in years


%run the production-diffusion function for each temperature history the user wants to test.
for i = 1:size(temps,2) 
  %temperatures corresponding to each time step for a particular thermal history in  ?C
    input.TC = temps(:,i)'; %interp1(times,temps(:,i),input.Tt);
    input.TC = fliplr(input.TC);
    input.TZ = fliplr(input.Tt).*input.ee + input.ztoday; %calculate sample depth at each time step in g/cm^2.
    
    output = ProdDiff_EDTScenario(input);
    
    results.actRtstep= output.actRtstep;
    results.romRtstep=output.romRtstep;
    results.TC=input.TC;
    
    %results.actHeatomsgtotal(i) = output.actHeatomsgtotal;
    results.totProduced(i) = output.totProducedConcstep(end);
    results.scale = output.P0./input.Pref;
    results.totProduced = output.totProducedConcstep;
    results.actTotHestep = output.actTotHestep;
    results.romTotHestep = output.romTotHestep;
    %results.depth = output.depth';
    %%%careful each results step correspond to the last EDT-t history
    %%%scenario, while the final R and HE values are saved for all the EDT-t
    %%%scenarios given
    results.P0=output.P0; %3He production rate at site
   
end
