% Post Processing "SideStep" example
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
function [RMSE_Pred, RMSE_COP] = RMSE_PostProcess(filename)


% Loading the Analysis file
%load('AnalysisParameters.mat')
AnalysisParameters.ExternalForces.Options(1,1)={'LFoot'};
AnalysisParameters.ExternalForces.Options(2,1)={'RFoot'};

% Get the solid names on which the forces are applied
Solids = AnalysisParameters.ExternalForces.Options;

% Loading the Biomechanicalmodel file
load('BiomechanicalModel.mat')

% Solid list extracted from the OsteoarticularModel
Solid_list = {BiomechanicalModel.OsteoArticularModel.name}';

% Get the numbers of solids on which the forces are applied
[~,num_s]=intersect(Solid_list,Solids);

% Name of the trial treated
% filename='ChangmntDirectionVinyl0001';
% filename='ChangmntDirectionVinyl0001';
%filename='Step0001';
%filename = 'Trial001';

% Loading external forces computation
load(fullfile(filename,'ExternalForcesComputationResults.mat'))
load(fullfile(filename,'ExperimentalData.mat'))

[start,stop] = croper(ExternalForcesComputationResults,AnalysisParameters, filename);

% Load experimental forces and forces from prediction algorithm.
GRF_Xp=ExternalForcesComputationResults.ExternalForcesExperiments(start:stop);
GRF_Pred=ExternalForcesComputationResults.ExternalForcesPrediction(start:stop);
%GRF_COP=ExternalForcesComputationResults.ExternalForcesPredictionCOP(start:stop);
GRF_COP=ExternalForcesComputationResults.ExternalForcesPrediction(start:stop);


% Number of frames
Nb_frames = numel(GRF_Xp);

% Get the forces applied on the solids
for ii=1:numel(num_s)
    cur_s=num_s(ii); %LFoot and RFoot
    for jj_f=1:Nb_frames
        % Prediction algorithms
        F_Pred.(Solids{ii})(jj_f,:) = GRF_Pred(jj_f).fext(cur_s).fext(:,1)';
        M_Pred.(Solids{ii})(jj_f,:) = GRF_Pred(jj_f).fext(cur_s).fext(:,2);
                
        % Prediction algorithms
        F_COP.(Solids{ii})(jj_f,:) = GRF_COP(jj_f).fext(cur_s).fext(:,1)';
        M_COP.(Solids{ii})(jj_f,:) = GRF_COP(jj_f).fext(cur_s).fext(:,2);

        % Experimental results
        F_Xp.(Solids{ii})(jj_f,:) = GRF_Xp(jj_f).fext(cur_s).fext(:,1)';
        M_Xp.(Solids{ii})(jj_f,:) = GRF_Xp(jj_f).fext(cur_s).fext(:,2)';
    end
end

%% RMSE

ErrorLFoot=[F_Pred.LFoot(:,1:3),M_Pred.LFoot(:,1:3)]-[F_Xp.LFoot(:,1:3),M_Xp.LFoot(:,1:3)];
RMSLFoot=rms(ErrorLFoot);

ErrorRFoot=[F_Pred.RFoot(:,1:3),M_Pred.RFoot(:,1:3)]-[F_Xp.RFoot(:,1:3),M_Xp.RFoot(:,1:3)];
RMSRFoot=rms(ErrorRFoot);

ErrorLFootCOP=[F_COP.LFoot(:,1:3),M_COP.LFoot(:,1:3)]-[F_Xp.LFoot(:,1:3),M_Xp.LFoot(:,1:3)];
RMSLFootCOP=rms(ErrorLFootCOP);

ErrorRFootCOP=[F_COP.RFoot(:,1:3),M_COP.RFoot(:,1:3)]-[F_Xp.RFoot(:,1:3),M_Xp.RFoot(:,1:3)];
RMSRFootCOP=rms(ErrorRFootCOP);

%% Plot the result
figure
set(gcf,'color','w')
suptitle(['External Forces Prediction: ' filename ,newline, newline])
Directions={'X','Y','Z'};
for ii=1:3
    subplot(2,3,ii)
    % LFoot
    plot(F_Pred.(Solids{1})(:,ii),'b-','LineWidth',2)
    hold on
    plot(F_COP.(Solids{1})(:,ii),'b-','LineWidth',1)
    plot(F_Xp.(Solids{1})(:,ii),'b--','LineWidth',2)
    % RFoot
    plot(F_Pred.(Solids{2})(:,ii),'r-','LineWidth',2)
    plot(F_COP.(Solids{2})(:,ii),'r-','LineWidth',1)
    plot(F_Xp.(Solids{2})(:,ii),'r--','LineWidth',2)
    
    %xlim([1 203])
    xlabel('Frames')
    ylabel('Force (N)')
    title({['Force applied on Feets on the ' Directions{ii} '-direction'] ,...
        ['RF RMSE : ' num2str(RMSRFoot(ii)) ' LF RMSE : ' num2str(RMSLFoot(ii)) 'RF RMSE COP : ' num2str(RMSRFootCOP(ii)) ' LF RMSE COP : ' num2str(RMSLFootCOP(ii)) ]})
end

for ii=4:6
    subplot(2,3,ii)
    % LFoot
    plot(M_Pred.(Solids{1})(:,ii-3),'b-','LineWidth',2)
    hold on
    plot(M_COP.(Solids{1})(:,ii-3),'b-','LineWidth',1)
    plot(M_Xp.(Solids{1})(:,ii-3),'b--','LineWidth',2)
    % RFoot
    plot(M_Pred.(Solids{2})(:,ii-3),'r-','LineWidth',2)
    plot(M_COP.(Solids{2})(:,ii-3),'r-','LineWidth',1)
    plot(M_Xp.(Solids{2})(:,ii-3),'r--','LineWidth',2)
    
    %xlim([1 203])
    xlabel('Frames')
    ylabel('Moment (Nm)')
    title({['Moment applied on Feets on the ' Directions{ii-3} '-direction'],...
        ['RF RMSE : ' num2str(RMSRFoot(ii)) ' LF RMSE : ' num2str(RMSLFoot(ii)) 'RF RMSE COP : ' num2str(RMSRFootCOP(ii)) ' LF RMSE COP : ' num2str(RMSLFootCOP(ii))]})
end

legend('Predicted on Left Foot','Experimental on Left Foot',...
    'Predicted on Right Foot','Experimental on Right Foot')

% On the graph, we can see that the predicted results are comparable to the
% experimental forces measured. However there is artefact when the foot is
% not on the plateforms. More work is needed on filters and algorithms to
% avoid it.

RMSE_Pred = rms([ErrorRFoot;ErrorLFoot]);
RMSE_COP = rms([ErrorRFootCOP;ErrorLFootCOP]);

end
