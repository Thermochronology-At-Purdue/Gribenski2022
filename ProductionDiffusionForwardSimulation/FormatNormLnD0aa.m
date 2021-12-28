
function results = FormatNormLnD0aa(dataKD,r);

%%
%This function transforms values of ln(D0/a^2) from diffusion experiment
%grain sizes to the size appropriate for cosmogenic noble gas measurements.

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


%Format structure array
NKD.a=dataKD.radius(1);
NKD.ndom=dataKD.Domain(1);
NKD.Ea=dataKD.Ea(1);

%Normalisation to data.a
NKD.lnD0aa= log(exp(dataKD.LnDoaa).*(dataKD.radius.^2)./(r.^2))';
NKD.fracs=dataKD.fraction';

results=NKD;

end

