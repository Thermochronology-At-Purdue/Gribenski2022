function MDD = TremblayMDD(input,geomflag)
%Code written by Marissa Tremblay. Last modified 2016.11.01.

%This function calculates the fraction of gas released from each domain
%in an MDD model during the heating schedule used in the diffusion 
%experiment. Then the fractions released from each domain are combined in
%proportion to one another as specified by the MDD model, and the
%diffusivity of each step is calculated. A residual is calculated as the
%sum of absolute differences between the observed and modeled release
%fractions over all steps.

R = 0.008314; %gas constant

% load inputs

fracs = input.fracs;
Ea = input.Ea;
lnD0aa = input.lnD0aa;
TC = input.TC;
thr = input.thr;
lnDaa = input.lnDaa;
Fi = input.Fi;

% get variables in correct units
tsec = thr.*60.*60;
cumtsec = cumsum(tsec);
TK = 273.15 + TC;
tempPlot = 1e4./TK;

%calculate a D/a^2 for each T in heating schedule using kinetics for each
%domain
for i = 1:length(Ea)
    for j = 1:length(TK)
        Daa(j,i) = exp(lnD0aa(i)).*exp(-Ea(i)./(R.*TK(j)));
    end
end

% Pre-allocate
f = zeros(size(Daa)); Dtaa = f;

% Calculate the Dt/a^2 for each temperature step. This needs to be a
% cumulative value.

for i = 1:size(Daa,2)
    Dtaa(1,i) = Daa(1,i).*tsec(1);
end
for i = 1:size(Daa,2)
    for j = 2:size(Daa,1)
        Dtaa(j,i) = Dtaa(j-1,i) + Daa(j,i).*(cumtsec(j)-cumtsec(j-1));
    end
end

%% Use equations from Fechtig and Kalbitzer to calculate release fractions.
if geomflag < 1 %spherical geometry
    for m = 1:size(Dtaa,2)
        for n = 1:size(Dtaa,1)
            f(n,m) = (6./(pi^(3/2))).*sqrt((pi.^2).*Dtaa(n,m));
            if f(n,m) >= 0.1
                f(n,m) = (6./(pi^(3/2))).*sqrt((pi.^2).*Dtaa(n,m))-(3./(pi.^2))...
                    .*((pi.^2)*Dtaa(n,m));
                if f(n,m) >= 0.9
                    f(n,m) = 1 - (6./(pi.^2)).*exp(-(pi.^2)*Dtaa(n,m));
                end
            end
        end
    end
else
    for m = 1:size(Dtaa,2)
        for n = 1:size(Dtaa,1)
            f(n,m) = (2./pi).*(Dtaa(n,m)).^0.5;
            if f(n,m) > 0.6
                f(n,m) = 1 - (8./(pi.^2)).*exp(-1.*(pi.^2).*Dtaa(n,m)./4);
            end
        end
    end
end
for m = 1:size(Dtaa,2)
    for n = 2:size(Dtaa,1)
        if f(n,m) < f(n-1,m) || f(n-1,m) == 1
            f(n,m) = 1;
        end
    end
end

MDD.f = f;

f_MDD = zeros(size(f));
for m = 1:length(fracs)
    for n = 1:length(f_MDD)
        f_MDD(n,m) = f(n,m).*fracs(m);
    end
end

sumf_MDD = zeros(size(f_MDD,1),1);
for n = 1:size(f_MDD,1)
    sumf_MDD(n) = sum(f_MDD(n,:));
end

MDD.f_MDD = sumf_MDD;

Daa_MDD = zeros(size(sumf_MDD));

for m = 1:size(sumf_MDD,1)
    if geomflag < 1
        if m == 1
            Daa_MDD(m) = (sumf_MDD(m).^2.*pi)./(36.*cumtsec(m));
        else
            Daa_MDD(m) = ((sumf_MDD(m)).^2-(sumf_MDD(m-1)).^2).*pi./...
                (36.*(cumtsec(m)-cumtsec(m-1)));
            if sumf_MDD(m) >= 0.1
                Daa_MDD(m) = (1./((pi.^2).*(cumtsec(m)-cumtsec(m-1)))).*...
                    (-(pi.*pi./3).*(sumf_MDD(m)-sumf_MDD(m-1))-(2.*pi).*...
                    (sqrt(1-(pi./3).*sumf_MDD(m))-sqrt(1-(pi./3).*...
                    sumf_MDD(m-1))));
                if sumf_MDD(m) >= 0.9
                    Daa_MDD(m) = (1./(pi.*pi.*(cumtsec(m)-cumtsec(m-1)))).*...
                        (log((1-sumf_MDD(m-1))./(1-sumf_MDD(m))));
                end
            end
        end
    else
        if m == 1
            %special case when i = 1; need to insert zero for previous amount
            %released
            Daa_MDD(m) = ((((sumf_MDD(m).^2) - 0.^2)).*pi)/(4.*cumtsec(m));
        else
            %equation for f < 0.6
            Daa_MDD(m) = ((sumf_MDD(m).^2-sumf_MDD(m-1).^2).*pi)...
                ./(4.*(cumtsec(m)-cumtsec(m-1)));
        end
        if sumf_MDD(m) > 0.6
            %equation for f > 0.6
            Daa_MDD(m) = (4./((pi.^2).*(cumtsec(m)-cumtsec(m-1))))...
                .*log((1-sumf_MDD(m-1))./(1-sumf_MDD(m)));
        end
    end
end

MDD.Daa_MDD = Daa_MDD;
lnDaa_MDD = log(Daa_MDD);
MDD.lnDaa_MDD = lnDaa_MDD;
%modified on 2016.11.01 to calculate misfit in terms of diffusivities,
%rather than release fractions
%residuals = abs(log(Daa_MDD) - lnDaa);
residuals = abs(sumf_MDD - Fi);

%remove infinities and NaNs that arise from earlier Arrehnius calcs
cut1 = ~isnan(residuals); cut2 = ~isinf(residuals);
residuals = residuals.*cut1.*cut2;
MDD.residuals = residuals(~isnan(residuals));
MDD.residual = sum(MDD.residuals);
