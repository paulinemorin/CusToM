function [COP_Xp] = COP_Xp_row(filename)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
file=filename;
load('AnalysisParameters.mat')
filename=file;
AnalysisParameters.ExternalForces.Options{2,1}='RFoot';
AnalysisParameters.ExternalForces.Options{1,1}='LFoot';

% Get the solid names on which the forces are applied
Solids = AnalysisParameters.ExternalForces.Options;

% Loading the Biomechanicalmodel file
load('BiomechanicalModel.mat')
Human_model = BiomechanicalModel.OsteoArticularModel;

% Solid list extracted from the OsteoarticularModel
Solid_list = {BiomechanicalModel.OsteoArticularModel.name}';

% Get the numbers of solids on which the forces are applied
[~,num_s]=intersect(Solid_list,Solids);

% Loading external forces computation
load(fullfile(filename,'ExternalForcesComputationResults.mat'))

% Load experimental forces and forces from prediction algorithm.
GRF_Xp=ExternalForcesComputationResults.ExternalForcesExperiments;

q = GRF_Xp;
nbframe=size(q,2);

%% Determination cOP Xp 

for ii=1:numel(num_s)
    COP_Xp.(Solids{ii})=zeros(nbframe,3);
    for i=1:nbframe
        COP_Xp.(Solids{ii})(i,:)=GRF_Xp(i).Visual(1:3,ii);
    end
end

end

