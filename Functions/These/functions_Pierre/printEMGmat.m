function [] = printEMGmat(folder,time,EMGsData,EMGsLabels, tag)
% Print normalized EMG linear envelopes in a .sto or mot file

% The file is part of matlab MOtion data elaboration TOolbox for
% NeuroMusculoSkeletal applications (MOtoNMS). 
% Copyright (C) 2012-2014 Alice Mantoan, Monica Reggiani
%
% MOtoNMS is free software: you can redistribute it and/or modify it under 
% the terms of the GNU General Public License as published by the Free 
% Software Foundation, either version 3 of the License, or (at your option)
% any later version.
%
% Matlab MOtion data elaboration TOolbox for NeuroMusculoSkeletal applications
% is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
% without even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
% PARTICULAR PURPOSE.  See the GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License along 
% with MOtoNMS.  If not, see <http://www.gnu.org/licenses/>.
%
% Alice Mantoan, Monica Reggiani
% <ali.mantoan@gmail.com>, <monica.reggiani@gmail.com>

%%
Envelope_EMG=struct();
for ii=1:length(EMGsLabels)
    Envelope_EMG(ii).name=EMGsLabels{ii};
    Envelope_EMG(ii).value=EMGsData(:,ii);
    Envelope_EMG(ii).time=time;
end
save(fullfile(folder,'Envelope_EMG'),'Envelope_EMG')
end