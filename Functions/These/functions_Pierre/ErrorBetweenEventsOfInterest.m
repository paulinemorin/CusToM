function []=ErrorBetweenEventsOfInterest(BiomechanicalModel,AnalysisParameters,filename,FrameOfInterest)
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
rms_error=0;
Directions={'X','Y','Z'};
for ii=1:3
    subplot(2,3,ii)
    % LFoot
    plot(F_error.(Solids{1})(start:stop,ii),'b-','LineWidth',2)
    hold on
    % RFoot
    plot(F_error.(Solids{2})(start:stop,ii),'r-','LineWidth',2)

%    xlim([1 length(F_Xp.(Solids{2})(start:stop,ii))])
    xlabel('Frames')
    ylabel('Error force (N)')
    title(['Error force on the ' Directions{ii} '-direction'])
end

for ii=4:6
    subplot(2,3,ii)
    % LFoot
    plot(M_error.(Solids{1})(start:stop,ii-3),'b-','LineWidth',2)
    hold on
    % RFoot
    plot(M_error.(Solids{2})(start:stop,ii-3),'r-','LineWidth',2)
    
%    xlim([1 length(M_Xp.(Solids{2})(start:stop,ii-3))])
    xlabel('Frames')
    ylabel('Error moment (Nm)')
    title(['Moment error on the ' Directions{ii-3} '-direction'])
    
end

sgtitle('Error GRF&M prediction')

legend('Predicted on Left Foot','Predicted on Right Foot')

% On the graph, we can see that the predicted results are comparable to the
% experimental forces measured. However there is artefact when the foot is
% not on the plateforms. More work is needed on filters and algorithms to
% avoid it.

end
