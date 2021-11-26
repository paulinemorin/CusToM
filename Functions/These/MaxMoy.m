function [Max, Moy] = MaxMoy(sujet,filename,ID)

Max=zeros(1,6);
Moy=zeros(1,6);



AnalysisParameters.ExternalForces.Options(1,1)={'LFoot'};
AnalysisParameters.ExternalForces.Options(2,1)={'RFoot'};

% Get the solid names on which the forces are applied
Solids = AnalysisParameters.ExternalForces.Options;

% Loading the Biomechanicalmodel file
load('BiomechanicalModel.mat')

%load('Event_fence.mat')

% Solid list extracted from the OsteoarticularModel
Solid_list = {BiomechanicalModel.OsteoArticularModel.name}';

% Get the numbers of solids on which the forces are applied
[~,num_s]=intersect(Solid_list,Solids);

% Loading external forces computation
load(fullfile(filename,'ExternalForcesComputationResults.mat'))

% Load experimental forces and forces from prediction algorithm.
%ID = [1, length(ExternalForcesComputationResults.ExternalForcesExperiments)];
GRF_Xp=ExternalForcesComputationResults.ExternalForcesExperiments(ID(1):ID(2));

% Number of frames
Nb_frames = numel(GRF_Xp);

% Get the forces applied on the solids
for ii=1:numel(num_s)
    cur_s=num_s(ii); %LFoot and RFoot
    for jj_f=1:Nb_frames         
        % Experimental results
        F_Xp.(Solids{ii})(jj_f,:) = GRF_Xp(jj_f).fext(cur_s).fext(:,1)';
        M_Xp.(Solids{ii})(jj_f,:) = GRF_Xp(jj_f).fext(cur_s).fext(:,2)';
    end
end

%% Max

F =[ F_Xp.(Solids{1});  F_Xp.(Solids{2})];
M =[ M_Xp.(Solids{1});  M_Xp.(Solids{2})];

GRF = [F, M];

Max = max(GRF);

Moy = mean(GRF);


end
