function out = ProdDiff_EDTScenario(data)

% Spherical-geometry production-diffusion model

%This has spatially homogeneous production that can vary with time. 
%Code was originally written by Greg Balco & David Shuster for modeling 
%apatite 4He/3He datasets. Modified here by Marissa Tremblay to model 
%production and diffusion of cosmogenic noble gases.

%Contact: tremblam@purdue.edu
%Last modified: 2021.12.22

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

%% Unpack input data

lat = data.lat; %sample latitude in decimal degrees
long = data.long; %sample longitude in decimal degrees; west = negative
elv = data.elv; %sample elevation in meters
Pref = data.Pref; %sea level high latitude production rate in atoms/g/yr
shield = data.shield; %correction factor for shielding
thick = data.thick; %sample thickness in cm
rho = data.rho; %mineral density in g/cm3
fsp = data.fsp; %fraction of spallogenic cosmogenic nuclide production
Lsp = data.Lsp; %effective attenuation length for spallation in g/cm^2
ee = data.ee; %erosion rate in g/cm^2/yr
r = data.r; %average radius of fragments analyzed, in cm

a = data.a; % radius of proton-irradiated grain analyzed (cm) 
ndom = data.ndom; %number of domains from the best-fit MDD model
allEa = data.Ea; %vector of activation energies for each domain.
alllnD0aa = data.lnD0aa; %vector of ln(D0/a^2) for each domain. technically unitless (normalized to 1/s)
fracs = data.fracs; %fraction of gas apportioned to each diffusion domain

n = data.n; %number of mesh points over which production-diffusion calculation takes place

Tt = data.Tt; %time steps in years
maxt = data.maxt; %total duration of the time-temperature history, in years
dt = data.dt; %size of the time step, in years
TZ = data.TZ; % sample depth at each time step in g/cm^2
TC = data.TC; %temperatures in ?C


%% Done unpacking data. Set up the production-diffusion calculation. Loop through each diffusion domain.

for dom = 1:1:ndom;
    thisdom = dom;
    Ea = allEa(thisdom);
    lnD0aa = alllnD0aa(thisdom);
    
    % Make the calculation mesh
    dx = (r/n); %node dimension, according to Cassata units and values are irrelevant
    x = ((dx/2):dx:(r-dx/2))'; %radial distance of the center of each node
    xlo = x-(dx./2);  %radial distance of the inside of each node
    xhi = x+(dx./2);  %radial distance of the outside of each node
    shellvol = (4/3).*pi.*(xhi.^3 - xlo.^3); %volume of each node
    
    % Weight of shells
    shellwt = shellvol.*rho; % g
    totwt = (4/3).*pi.*(r.^3).*rho;
    
    T = TC+273; % convert to Kelvin
    
    % Deal with diffusion constants
    R = 8.3144; % j/K/mol
    % convert Ea to J/mol
    Ea = Ea.*1e3;
    D0 = (60*60*24*365.25).*((exp(lnD0aa)).*(r.^2)); %frequency factor (cm^2/yr)
    
    % Deal with production rate constants
    pressure = NCEPatm_2(lat,long,elv); %looks up surface pressure and 1000 mb
    %from NCEP reanalysis. Need this for
    %site specific production rate calc
    thickcorr = thickness(thick,Lsp,rho); %calculate the correction needed for
    %sample thickness.
    %calculate geographic scaling factor and multiply by the reference
    %production rate, shielding correction factor, and thickness correction
    %factor to determine the local cosmogenic nulcide production rate
    P0 = Pref.*shield.*thickcorr.*stone2000(lat,pressure,fsp);
    out.P0 = P0; %site specific surface production rate
    
    %% Done setting up. Initalize variables. 
    
    oldC = zeros(size(x)); %initialize concentration vector
    tv = 0:dt:maxt; %initialize time vector;
    Ps = zeros(size(tv)); %initialize production vector
    totHe = zeros(size(tv)); %initialize vector for total He produced
    
 %% Done initializing variables. Start the solver loop.
    %production in time step 1
    for a = 1:length(tv); % Step 1 is all zeros, start at step 2
        %changed from a = 2:length(tv); need production in the first time
        %step
        % Change of variables
        oldU = oldC.*x;
        % Obtain K
        thisT = interp1(Tt,T,tv(a)); %temperature vector, interpolated between
        D = D0*exp(-Ea/R*(1./(thisT))); %diffusivity, in cm^2/yr
        K = D;
        
        % Solver setup
        beta = 2.*(dx.^2)./(K.*dt);
        A1 = diag(ones(1,n)).*(-2-beta) + diag(ones(1,n-1),-1) + diag(ones(1,n-1),1);
        A2 = diag(ones(1,n)).*(2-beta) - diag(ones(1,n-1),-1) - diag(ones(1,n-1),1);
        
        % No-flux LH boundary
        A1(1,1) = A1(1,1)-1;
        A2(1,1) = A2(1,1)+1;
        
        % He source
        % This is due to cosmic ray production
        % Obtain depth
        thisz = interp1(Tt,TZ,tv(a));
        depth(a) = thisz;
        %thisz = (maxt - dt.*(a-1)).*ee;
        %not sure which is the correct/preferred way to calculate thisz
        % Compute production at depth
        Pofx = P0.*exp(-thisz./Lsp); % Atoms/g/yr
        % store for tot-up calculation later
        Ps(a) = Pofx;
        % Add source to totHe -- remember totHe is total He production, not
        % how much He is in grain now. totHe should be exact always - ??
        newHe = Pofx.*dt.*shellwt; % atoms
        
        if a == 1
            totHe(a) = sum(newHe); %atoms
        else
            totHe(a) = totHe(a-1) + sum(newHe); % atoms
        end
        
        % build RHS
        b = A2*oldU - Pofx.*x.*beta.*dt;
        
        % solve
        if isnan(rcond(A1));
            disp(['K = ' num2str(K)]);
        end;
        newU = A1\b;
        
        % un-change of variables
        newC = newU./x;
        
        % Update
        oldC = newC;      
    %end;
    
% Sum up total He left in domain. 
% Use trapezoidal integration
actHe = sum(oldC.*shellwt); % atoms

% This code taken from rombint.m
% Figure number of iterations
decdigs = 1+round(log2(data.n-1 ));
rom=zeros(2,decdigs);
romall=oldC.*4.*pi.*(x.^2).*rho; 
romall(end+1) = 0;
h = r;
rom(1,1)=h*(romall(1)+romall(end))/2;
for i=2:decdigs
  st=2^(decdigs-i+1);
  rom(2,1)=(rom(1,1)+h*sum(romall((st/2)+1:st:2^(decdigs-1))))/2;
  for k=1:i-1
     rom(2,k+1)=((4^k)*rom(2,k)-rom(1,k))/((4^k)-1);
  end
  rom(1,1:i)=rom(2,1:i);
  h=h/2;
end
romHe=rom(1,decdigs); % This should also yield atoms
%Number of atoms 
out.actHe(thisdom) = actHe; out.romHe(thisdom) = romHe; out.totwt(thisdom) = totwt;
%Number of atoms per gram
out.actHeatomsg(thisdom) = actHe./totwt;
out.romHeatomsg(thisdom) = romHe./totwt;

% Compute total produced
out.totProduced(thisdom) = sum(Ps.*dt.*totwt); % atoms
out.totHe(thisdom) = totHe(end);

out.checkTotal(thisdom) = mean(Ps).*maxt;

out.actR(thisdom) = out.actHe(thisdom)./out.totProduced(thisdom);
out.romR(thisdom) = out.romHe(thisdom)./out.totProduced(thisdom);
out.checkR(thisdom) = out.actHe(thisdom)./out.checkTotal(thisdom)./totwt;

%a is the time step, so here we have R and He for each time step for the
%given domain
actRsave(a,dom) = out.actR(thisdom);
romRsave(a,dom) = out.romR(thisdom);
%Add Nat
actTotHesave(a,dom) = out.actHeatomsg(thisdom);
romTotHesave(a,dom) = out.romHeatomsg(thisdom);
totProducedConc(a,thisdom) = sum(Ps.*dt);%.*totwt); % atoms

    end
    
    actRsave(:,dom) = actRsave(:,dom).*fracs(dom);
    romRsave(:,dom) = romRsave(:,dom).*fracs(dom);
    %Add Nat
    actTotHesave(:,dom) = actTotHesave(:,dom).*fracs(dom);
    romTotHesave(:,dom) = romTotHesave(:,dom).*fracs(dom);
    totProducedConc(:,dom) = totProducedConc(:,dom).*fracs(dom);
    
end;

%% sum up retention over all domains
out.actRtotal = sum(out.actR.*fracs);
out.romRtotal = sum(out.romR.*fracs);
out.actHeatomsgtotal = sum(out.actHeatomsg.*fracs);

out.actRtstep = sum(actRsave,2);
out.romRtstep = sum(romRsave,2);
out.actTotHestep = sum(actTotHesave,2); 
out.romTotHestep = sum(romTotHesave,2);
out.totProducedConcstep = sum(totProducedConc,2)
out.tv = tv;
out.depth = depth;

