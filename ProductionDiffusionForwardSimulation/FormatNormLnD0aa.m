
function results = FormatNormLnD0aa(dataKD,r);

%Format structure array
NKD.a=dataKD.radius(1);
NKD.ndom=dataKD.Domain(1);
NKD.Ea=dataKD.Ea(1);

%Normalisation to data.a
NKD.lnD0aa= log(exp(dataKD.LnDoaa).*(dataKD.radius.^2)./(r.^2))';
NKD.fracs=dataKD.fraction';

results=NKD;

end

