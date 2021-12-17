function out = D0calc_MonteCarloErrors(expdata,geomflag)

%Code written by Marissa Tremblay, modified from earlier code written by 
%Greg Balco. Calculates diffusivities from step-heating degasing 
%experiments using equations of Fechtig and Kalbitzer (1966). Propagates 
%uncertainty using a Monte-Carlo approach. Future versions will have an 
%iterative scheme that determines how many data points to include in a 
%linear fit (currently this is chosen by the user).

%% Take the data from the diffusion experiment and uses it to calculate ln(D/a2) values using Fechtig and Kalbitzer (1966).

TC = expdata(:,2);      %temperature in degrees Celsius
thr = expdata(:,3);     %duration in hours
M = expdata(:,4);       %gas amount measured in each step
delM = expdata(:,5);    %measurement uncertainty 

%rationalize everything to column vectors, if necessary
if size(TC,1) < size(TC,2);
    TC = TC';
end;
if size(thr,1) < size(thr,2);
    thr = thr';
end;
if size(M,1) < size(M,2);
    M = M';
end;
if size(delM,1) < size(delM,2);
    delM = delM';
end;

%unit conversions
TK = 273.15+TC;                %convert temperature to Kelvin
tsec = thr.*60.*60;            %convert time to seconds
cumtsec = cumsum(tsec);        %cumulative time
Tplot = 1e4./TK;               %Arrhenius plot temperature

nstep = length(M);            %number of steps in experiment

%compute cumulative gas release fraction
Si = cumsum(M);
S = max(Si);
Fi = Si./S;

%initialize diffusivity vectors fore each Fechtig and Kalbitzer equation
DR2_a = zeros(nstep,1); DR2_b = DR2_a; DR2_c = DR2_a;

diffti = cumtsec(2:end) - cumtsec(1:(end-1)); %time difference between steps
diffFi = Fi(2:end) - Fi(1:(end-1)); %gas fraction difference between steps

%if diffusion geometry flag is for spherical geometry (i.e., geomflag = 0),
%use equations 5a through c from Fechtig and Kalbitzer.
if geomflag < 1;
    %Fechtig and Kalbitzer Equation 5a, for cumulative gas fractions up to 10%
    %special case when i = 1; need to insert 0 for previous amount released
    DR2_a(1) = ( (Fi(1)).^2 - 0.^2 ).*pi./(36.*(cumtsec(1)));
    %Equation 5a for all other steps
    DR2_a(2:end) = ( (Fi(2:end)).^2 - (Fi(1:(end-1))).^2 ).*pi./(36.*(diffti));
    
    % Fechtig and Kalbitzer Equation 5b, for cumulative gas fractions between 10 and 90%
    DR2_b(2:end) = (1./((pi.^2).*diffti)).*( -(pi.*pi./3).*diffFi - (2.*pi).*...
        ( sqrt(1-(pi./3).*Fi(2:end)) - sqrt(1 - (pi./3).*Fi(1:(end-1))) ));
    
    % Fechtig and Kalbitzer Equation 5c, for cumulative gas fractions greater than 90%
    DR2_c(2:end) = (1./(pi.*pi.*diffti)).*(log((1-Fi(1:(end-1)))./(1-Fi(2:end))));
    
    % Decide what equation to use based on the cumulative gas fraction of each step
    use_a = (Fi <= 0.1 & Fi > 0.0001);
    use_b = (Fi > 0.1 & Fi <= 0.9);
    use_c = (Fi > 0.9 & Fi <= 1);
    
    %make vector with correct diffusivites for each step
    DR2 = use_a.*DR2_a + use_b.*DR2_b + use_c.*DR2_c;
else
    %special case when i = 1; need to insert zero for previous amount
    %released
    DR2_a(1) = ((((Fi(1).^2) - 0.^2)).*pi)/(4.*cumtsec(1));
    %equation for f < 0.6
    DR2_a(2:end) = ((((Fi(2:end).^2)-((Fi(1:(end-1)))).^2)).*pi)./(4.*diffti);
    
    %equation for f > 0.6
    DR2_b(2:end) = (4./((pi().^2).*diffti)).*log((1-Fi(1:(end-1)))./(1-Fi(2:end)));
    
    %decide what to use
    use_a = (Fi > 0 & Fi < 0.6);
    use_b = (Fi >= 0.6 & Fi <= 1);
    
    DR2 = use_a.*DR2_a + use_b.*DR2_b;
end
    
%% Compute uncertainties in diffusivity using a Monte Carlo simulation
%Generates simulated step degassing datasets, such that each step of the 
%experiment has a Gaussian distribution centered at M and with 1s.d. of 
%delM across the simulated datasets.Then recomputes diffusivities for each 
%simulated dataset and uses the range of diffusivities for each step across
%all simulated datasets to estimate uncertainty. 

n_sim = 30000;              %number of simulations in the Monte Carlo scheme
MCsim = zeros(nstep,n_sim); %initialize matrix for simulated measurements

for m = 1:nstep
    for n = 1:n_sim
        %generate the simulated measurements 
        MCsim(m,n) = randn.*delM(m) + M(m); 
    end
end

%compute cumulative gas release fraction for each simulation
MCSi = cumsum(MCsim,1);
MCS = max(MCSi);
MCFi = zeros(nstep,n_sim);
for m = 1:nstep
    for n = 1:n_sim
        MCFi(m,n) = MCSi(m,n)./MCS(n);
    end
end

%initialize diffusivity vectors fore each Fechtig and Kalbitzer equation
MCDR2_a = zeros(nstep,n_sim); MCDR2_b = MCDR2_a; MCDR2_c = MCDR2_a;

%gas fraction difference between steps in each simulated step degassing experiment
MCdiffFi = zeros(size(DR2_a));
for m = 2:nstep
    for n = 1:n_sim
        MCdiffFi(m,n) = MCFi(m,n) - MCFi(m-1,n);
    end
end

if geomflag < 1 %spherical geometry
    %Fechtig and Kalbitzer Equation 5a, for cumulative gas fractions up to 10%
    %special case when i = 1; need to insert 0 for previous amount released
    for n = 1:n_sim
        MCDR2_a(1,n) = ( (MCFi(1,n)).^2 - 0.^2 ).*pi./(36.*(cumtsec(1)));
    end
    %all other steps
    for m = 2:nstep
        for n = 1:n_sim
            MCDR2_a(m,n) = ( (MCFi(m,n)).^2 - (MCFi(m-1,n)).^2 ).*pi./(36.*(diffti(m-1)));
        end
    end
    
    % Fechtig and Kalbitzer Equation 5b, for cumulative gas fractions between 10 and 90%
    for m = 2:nstep
        for n = 1:n_sim
            MCDR2_b(m,n) = (1./((pi.^2).*diffti(m-1))).*( -(pi.*pi./3).*...
                MCdiffFi(m,n) - (2.*pi).*( sqrt(1-(pi./3).*MCFi(m,n)) - ...
                sqrt(1 - (pi./3).*MCFi(m-1,n)) ));
        end
    end
    
    % Fechtig and Kalbitzer Equation 5c, for cumulative gas fractions greater than 90%
    for m = 2:nstep
        for n = 1:n_sim
            MCDR2_c(m,n) = (1./(pi.*pi.*diffti(m-1))).*(log((1-MCFi(m-1,n))./(1-MCFi(m,n))));
        end
    end
    
    % Decide what equation to use based on the cumulative gas fraction of each step
    MCuse_a = (MCFi <= 0.1 & MCFi > 0.0001);
    MCuse_b = (MCFi > 0.1 & MCFi <= 0.9);
    MCuse_c = (MCFi > 0.9 & MCFi <= 1);
    
    %make vector with correct diffusivites for each step
    MCDR2 = MCuse_a.*MCDR2_a + MCuse_b.*MCDR2_b + MCuse_c.*MCDR2_c;
else % diffusion geometry is plane sheet    
    %equation for f < 0.6
    % Special case when i = 1 -- need to insert 0 for previous amount released
    for n = 1:n_sim
        MCDR2_a(1,n) = ((MCFi(1,n)).^2 - 0.^2 ).*pi./(4.*(cumtsec(1)));
    end
    %equation for f < 0.6 for all other steps
    for m = 2:nstep
        for n = 1:n_sim
            MCDR2_a(m,n) = ((MCFi(m,n)).^2 - (MCFi(m-1,n)).^2 )...
                .*pi./(4.*(diffti(m-1)));
        end
    end
    %equation for f > 0.6
    for m = 2:nstep
        for n = 1:n_sim
            MCDR2_b(m,n) = (4./((pi.^2).*diffti(m-1)))...
                .*log((1-MCFi(m-1,n))./(1-MCFi(m,n)));
        end
    end
    % Decide what to use
    MCuse_a = (MCFi < 0.6 & MCFi > 0);
    MCuse_b = (MCFi >= 0.6 & MCFi <= 1);
    MCDR2 = MCuse_a.*MCDR2_a + MCuse_b.*MCDR2_b;
end

%estimate the uncertainty in diffusivity for each step
MCDR2_uncert = zeros(nstep,1);                              %initialize uncertainty vector
%estimate uncertainty in diffisuvity for one step as as 1s.d. of the diffusivities across all simulations
for m = 1:nstep
    MCDR2_uncert(m,1) = std(MCDR2(m,:));           
end

% This code dumps calculated diffusivities to command window
disp('Step 1/T ln(D0/a2) -1SD +1SD');
for a = 1:nstep;
    str2 = sprintf('%0.0f\t %0.2f\t %0.2f\t %0.2f\t %0.0f',[a Tplot(a) ...
        log(DR2(a)) log(DR2(a)-MCDR2_uncert(a)) log(DR2(a)+MCDR2_uncert(a))]);
    disp(str2);
end;

out.a = a;
out.Tplot = Tplot;
out.Fi = Fi;
out.lnDaa = log(DR2);
out.lnDaaneg = log(DR2 - MCDR2_uncert);
out.lnDaapos = log(DR2 + MCDR2_uncert);
