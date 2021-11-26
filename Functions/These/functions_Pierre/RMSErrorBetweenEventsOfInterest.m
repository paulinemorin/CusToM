function [rms_error]=RMSErrorBetweenEventsOfInterest(BiomechanicalModel,AnalysisParameters,filename,FrameOfInterest)
 %Post Processing "SideStep" example
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This example aimed at comparing the consistency of the computed External
% Forces with the Prediction algorithms on cut-off movement of a single 
% subject.
%________________________________________________________
%
% Licence
% Toolbox distributed under GPL 3.0 Licence
%________________________________________________________
%
% Authors : Antoine Muller, Charles Pontonnier, Pierre Puchaud and
% Georges Dumont
%_______________________________________________________


% This function give the RMS between two strike foot for run tests

% Get the solid names on which the forces are applied
Solids = AnalysisParameters.ExternalForces.Options;

% Solid list extracted from the OsteoarticularModel
Solid_list = {BiomechanicalModel.OsteoArticularModel.name}';

% Get the numbers of solids on which the forces are applied
[~,num_s]=intersect(Solid_list,Solids);

% Loading external forces computation
load(fullfile(filename,'ExternalForcesComputationResults.mat'))

% Load experimental forces and forces from prediction algorithm.
GRF_Xp=ExternalForcesComputationResults.ExternalForcesExperiments;
GRF_Pred=ExternalForcesComputationResults.ExternalForcesPrediction;

% Number of frames
Nb_frames = numel(GRF_Xp);


% Get the forces applied on the solids
for ii=1:numel(num_s)
    cur_s=num_s(ii); %LFoot and RFoot
    for jj_f=1:Nb_frames
        % Prediction algorithms
        F_Pred.(Solids{ii})(jj_f,:) = GRF_Pred(jj_f).fext(cur_s).fext(:,1)';
        M_Pred.(Solids{ii})(jj_f,:) = GRF_Pred(jj_f).fext(cur_s).fext(:,2)';
        % Experimental results
        F_Xp.(Solids{ii})(jj_f,:) = GRF_Xp(jj_f).fext(cur_s).fext(:,1)';
        M_Xp.(Solids{ii})(jj_f,:) = GRF_Xp(jj_f).fext(cur_s).fext(:,2)';
        
        % error
        F_error.(Solids{ii})(jj_f,:) = F_Pred.(Solids{ii})(jj_f,:)-F_Xp.(Solids{ii})(jj_f,:);
        M_error.(Solids{ii})(jj_f,:) = M_Pred.(Solids{ii})(jj_f,:)-M_Xp.(Solids{ii})(jj_f,:);
        
    end
end


%% Events of interest

if contains(filename,'Course')
    start = FrameOfInterest.Event(1).frame;
    stop = FrameOfInterest.Event(3).frame;
    if start>stop
        temp = start;
        start = stop;
        stop=temp;
    end
end


%% Plot the result
rms_error=zeros(2,6);

for ii=1:3
    % LFoot
    rms_error(1,ii) = rms(F_error.(Solids{1})(start:stop,ii));
    
    % RFoot
    rms_error(2,ii) = rms(F_error.(Solids{2})(start:stop,ii));
end

for ii=4:6
    % LFoot
    rms_error(1,ii) = rms(M_error.(Solids{1})(start:stop,ii-3));

    % RFoot
    rms_error(2,ii) = rms(M_error.(Solids{2})(start:stop,ii-3));
end

ExternalForcesComputationResults.RMS.ExternalForcesPrediction_plus=rms_error_plus;
ExternalForcesComputationResults.RMS.ExternalForcesPrediction=rms_error;
save([filename '/ExternalForcesComputationResults'],'ExternalForcesComputationResults');

end