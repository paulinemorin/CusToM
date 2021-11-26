function [RMS_Threshold, RMS_InsoleDetection] = PostProcessISB(sujet,filename,ID)

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

% Name of the trial treated
%filename='ChgtDirection01';
%filename='SQA005';
%filename = 'Trial002';

% Loading external forces computation
load(fullfile(filename,'ExternalForcesComputationResults.mat'))

% Load experimental forces and forces from prediction algorithm.
GRF_Xp=ExternalForcesComputationResults.ExternalForcesExperiments(ID(1):ID(2));
GRF_Pred_Threshold=ExternalForcesComputationResults.ExternalForcesPrediction(ID(1):ID(2));
GRF_Pred_InsoleDetection=ExternalForcesComputationResults.ExternalForcesPredictionInsoleDetection(ID(1):ID(2));
% GRF_Pred_75_25=ExternalForcesComputationResults.ExternalForcesPrediction75_25;
% GRF_Pred_50_50=ExternalForcesComputationResults.ExternalForcesPrediction50_50;
% GRF_Pred_0_1=ExternalForcesComputationResults.ExternalForcesPrediction0_1;


% Number of frames
Nb_frames = numel(GRF_Xp);

% Get the forces applied on the solids
for ii=1:numel(num_s)
    cur_s=num_s(ii); %LFoot and RFoot
    for jj_f=1:Nb_frames
        % Prediction algorithms
        F_Pred_Threshold.(Solids{ii})(jj_f,:) = GRF_Pred_Threshold(jj_f).fext(cur_s).fext(:,1)';
        M_Pred_Threshold.(Solids{ii})(jj_f,:) = GRF_Pred_Threshold(jj_f).fext(cur_s).fext(:,2);
        
        % Prediction algorithms
        F_Pred_InsoleDetection.(Solids{ii})(jj_f,:) = GRF_Pred_InsoleDetection(jj_f).fext(cur_s).fext(:,1)';
        M_Pred_InsoleDetection.(Solids{ii})(jj_f,:) = GRF_Pred_InsoleDetection(jj_f).fext(cur_s).fext(:,2);
        
%         % Prediction algorithms
%         F_Pred_50_50.(Solids{ii})(jj_f,:) = GRF_Pred_50_50(jj_f).fext(cur_s).fext(:,1)';
%         M_Pred_50_50.(Solids{ii})(jj_f,:) = GRF_Pred_50_50(jj_f).fext(cur_s).fext(:,2);
%         
%         % Prediction algorithms
%         F_Pred_25_75.(Solids{ii})(jj_f,:) = GRF_Pred_25_75(jj_f).fext(cur_s).fext(:,1)';
%         M_Pred_25_75.(Solids{ii})(jj_f,:) = GRF_Pred_25_75(jj_f).fext(cur_s).fext(:,2);
%         
%         % Prediction algorithms
%         F_Pred_0_1.(Solids{ii})(jj_f,:) = GRF_Pred_0_1(jj_f).fext(cur_s).fext(:,1)';
%         M_Pred_0_1.(Solids{ii})(jj_f,:) = GRF_Pred_0_1(jj_f).fext(cur_s).fext(:,2);
%         
        % Experimental results
        F_Xp.(Solids{ii})(jj_f,:) = GRF_Xp(jj_f).fext(cur_s).fext(:,1)';
        M_Xp.(Solids{ii})(jj_f,:) = GRF_Xp(jj_f).fext(cur_s).fext(:,2)';
    end
end

%% RMSE 
ErrorLFoot_Threshold=[F_Pred_Threshold.LFoot(:,1:3),M_Pred_Threshold.LFoot(:,1:3)]-[F_Xp.LFoot(:,1:3),M_Xp.LFoot(:,1:3)];
RMSLFoot_Threshold=rms(ErrorLFoot_Threshold);
ErrorRFoot_Threshold=[F_Pred_Threshold.RFoot(:,1:3),M_Pred_Threshold.RFoot(:,1:3)]-[F_Xp.RFoot(:,1:3),M_Xp.RFoot(:,1:3)];
RMSRFoot_Threshold=rms(ErrorRFoot_Threshold);

RMS_Threshold=mean([RMSLFoot_Threshold;RMSRFoot_Threshold],1);

ErrorLFoot_InsoleDetection=[F_Pred_InsoleDetection.LFoot(:,1:3),M_Pred_InsoleDetection.LFoot(:,1:3)]-[F_Xp.LFoot(:,1:3),M_Xp.LFoot(:,1:3)];
RMSLFoot_InsoleDetection=rms(ErrorLFoot_InsoleDetection);
ErrorRFoot_InsoleDetection=[F_Pred_InsoleDetection.RFoot(:,1:3),M_Pred_InsoleDetection.RFoot(:,1:3)]-[F_Xp.RFoot(:,1:3),M_Xp.RFoot(:,1:3)];
RMSRFoot_InsoleDetection=rms(ErrorRFoot_InsoleDetection);

RMS_InsoleDetection=mean([RMSLFoot_InsoleDetection;RMSRFoot_InsoleDetection],1);


% 
% ErrorLFoot_50_50=[F_Pred_50_50.LFoot(:,1:3),M_Pred_50_50.LFoot(:,1:3)]-[F_Xp.LFoot(:,1:3),M_Xp.LFoot(:,1:3)];
% RMSLFoot_50_50=rms(ErrorLFoot_50_50);
% ErrorRFoot_50_50=[F_Pred_50_50.RFoot(:,1:3),M_Pred_50_50.RFoot(:,1:3)]-[F_Xp.RFoot(:,1:3),M_Xp.RFoot(:,1:3)];
% RMSRFoot_50_50=rms(ErrorRFoot_50_50);
% 
% RMS_50_50=mean([RMSLFoot_50_50;RMSRFoot_50_50],1);
% 
% ErrorLFoot_25_75=[F_Pred_25_75.LFoot(:,1:3),M_Pred_25_75.LFoot(:,1:3)]-[F_Xp.LFoot(:,1:3),M_Xp.LFoot(:,1:3)];
% RMSLFoot_25_75=rms(ErrorLFoot_25_75);
% ErrorRFoot_25_75=[F_Pred_25_75.RFoot(:,1:3),M_Pred_25_75.RFoot(:,1:3)]-[F_Xp.RFoot(:,1:3),M_Xp.RFoot(:,1:3)];
% RMSRFoot_25_75=rms(ErrorRFoot_25_75);
% 
% RMS_25_75=mean([RMSLFoot_25_75;RMSRFoot_25_75],1);
% 
% ErrorLFoot_0_1=[F_Pred_0_1.LFoot(:,1:3),M_Pred_0_1.LFoot(:,1:3)]-[F_Xp.LFoot(:,1:3),M_Xp.LFoot(:,1:3)];
% RMSLFoot_0_1=rms(ErrorLFoot_0_1);
% ErrorRFoot_0_1=[F_Pred_0_1.RFoot(:,1:3),M_Pred_0_1.RFoot(:,1:3)]-[F_Xp.RFoot(:,1:3),M_Xp.RFoot(:,1:3)];
% RMSRFoot_0_1=rms(ErrorRFoot_0_1);
% 
% RMS_0_1=mean([RMSLFoot_0_1;RMSRFoot_0_1],1);
% 
% RMS=[RMS_1_0;
%     RMS_75_25;
%     RMS_50_50;
%     RMS_25_75;
%     RMS_0_1];




%% RMSE dynamic
% Event_fence(sujet,2)=length(F_Pred_1_0.LFoot);
% 
% ErrorLFoot_1_0_dynamic=[F_Pred_1_0.LFoot(Event_fence(sujet,1):Event_fence(sujet,2),1:3),M_Pred_0_1.LFoot(Event_fence(sujet,1):Event_fence(sujet,2),1:3)]-[F_Xp.LFoot(Event_fence(sujet,1):Event_fence(sujet,2),1:3),M_Xp.LFoot(Event_fence(sujet,1):Event_fence(sujet,2),1:3)];
% RMSLFoot_1_0_dynamic=rms(ErrorLFoot_1_0_dynamic);
% ErrorRFoot_1_0_dynamic=[F_Pred_1_0.RFoot(Event_fence(sujet,1):Event_fence(sujet,2),1:3),M_Pred_1_0.RFoot(Event_fence(sujet,1):Event_fence(sujet,2),1:3)]-[F_Xp.RFoot(Event_fence(sujet,1):Event_fence(sujet,2),1:3),M_Xp.RFoot(Event_fence(sujet,1):Event_fence(sujet,2),1:3)];
% RMSRFoot_1_0_dynamic=rms(ErrorRFoot_1_0_dynamic);
% 
% RMS_1_0_dynamic=mean([RMSLFoot_1_0_dynamic;RMSRFoot_1_0_dynamic],1);
% 
% ErrorLFoot_75_25_dynamic=[F_Pred_75_25.LFoot(Event_fence(sujet,1):Event_fence(sujet,2),1:3),M_Pred_75_25.LFoot(Event_fence(sujet,1):Event_fence(sujet,2),1:3)]-[F_Xp.LFoot(Event_fence(sujet,1):Event_fence(sujet,2),1:3),M_Xp.LFoot(Event_fence(sujet,1):Event_fence(sujet,2),1:3)];
% RMSLFoot_75_25_dynamic=rms(ErrorLFoot_75_25_dynamic);
% ErrorRFoot_75_25_dynamic=[F_Pred_75_25.RFoot(Event_fence(sujet,1):Event_fence(sujet,2),1:3),M_Pred_75_25.RFoot(Event_fence(sujet,1):Event_fence(sujet,2),1:3)]-[F_Xp.RFoot(Event_fence(sujet,1):Event_fence(sujet,2),1:3),M_Xp.RFoot(Event_fence(sujet,1):Event_fence(sujet,2),1:3)];
% RMSRFoot_75_25_dynamic=rms(ErrorRFoot_75_25_dynamic);
% 
% RMS_75_25_dynamic=mean([RMSLFoot_75_25_dynamic;RMSRFoot_75_25_dynamic],1);
% 
% ErrorLFoot_50_50_dynamic=[F_Pred_50_50.LFoot(Event_fence(sujet,1):Event_fence(sujet,2),1:3),M_Pred_50_50.LFoot(Event_fence(sujet,1):Event_fence(sujet,2),1:3)]-[F_Xp.LFoot(Event_fence(sujet,1):Event_fence(sujet,2),1:3),M_Xp.LFoot(Event_fence(sujet,1):Event_fence(sujet,2),1:3)];
% RMSLFoot_50_50_dynamic=rms(ErrorLFoot_50_50_dynamic);
% ErrorRFoot_50_50_dynamic=[F_Pred_50_50.RFoot(Event_fence(sujet,1):Event_fence(sujet,2),1:3),M_Pred_50_50.RFoot(Event_fence(sujet,1):Event_fence(sujet,2),1:3)]-[F_Xp.RFoot(Event_fence(sujet,1):Event_fence(sujet,2),1:3),M_Xp.RFoot(Event_fence(sujet,1):Event_fence(sujet,2),1:3)];
% RMSRFoot_50_50_dynamic=rms(ErrorRFoot_50_50_dynamic);
% 
% RMS_50_50_dynamic=mean([RMSLFoot_50_50_dynamic;RMSRFoot_50_50_dynamic],1);
% 
% ErrorLFoot_25_75_dynamic=[F_Pred_25_75.LFoot(Event_fence(sujet,1):Event_fence(sujet,2),1:3),M_Pred_25_75.LFoot(Event_fence(sujet,1):Event_fence(sujet,2),1:3)]-[F_Xp.LFoot(Event_fence(sujet,1):Event_fence(sujet,2),1:3),M_Xp.LFoot(Event_fence(sujet,1):Event_fence(sujet,2),1:3)];
% RMSLFoot_25_75_dynamic=rms(ErrorLFoot_25_75_dynamic);
% ErrorRFoot_25_75_dynamic=[F_Pred_25_75.RFoot(Event_fence(sujet,1):Event_fence(sujet,2),1:3),M_Pred_25_75.RFoot(Event_fence(sujet,1):Event_fence(sujet,2),1:3)]-[F_Xp.RFoot(Event_fence(sujet,1):Event_fence(sujet,2),1:3),M_Xp.RFoot(Event_fence(sujet,1):Event_fence(sujet,2),1:3)];
% RMSRFoot_25_75_dynamic=rms(ErrorRFoot_25_75_dynamic);
% 
% RMS_25_75_dynamic=mean([RMSLFoot_25_75_dynamic;RMSRFoot_25_75_dynamic],1);
% 
% ErrorLFoot_0_1_dynamic=[F_Pred_0_1.LFoot(Event_fence(sujet,1):Event_fence(sujet,2),1:3),M_Pred_0_1.LFoot(Event_fence(sujet,1):Event_fence(sujet,2),1:3)]-[F_Xp.LFoot(Event_fence(sujet,1):Event_fence(sujet,2),1:3),M_Xp.LFoot(Event_fence(sujet,1):Event_fence(sujet,2),1:3)];
% RMSLFoot_0_1_dynamic=rms(ErrorLFoot_0_1_dynamic);
% ErrorRFoot_0_1_dynamic=[F_Pred_0_1.RFoot(Event_fence(sujet,1):Event_fence(sujet,2),1:3),M_Pred_0_1.RFoot(Event_fence(sujet,1):Event_fence(sujet,2),1:3)]-[F_Xp.RFoot(Event_fence(sujet,1):Event_fence(sujet,2),1:3),M_Xp.RFoot(Event_fence(sujet,1):Event_fence(sujet,2),1:3)];
% RMSRFoot_0_1_dynamic=rms(ErrorRFoot_0_1_dynamic);
% 
% RMS_0_1_dynamic=mean([RMSLFoot_0_1_dynamic;RMSRFoot_0_1_dynamic],1);
% 
% RMS_dynamic=[RMS_1_0_dynamic;
%     RMS_75_25_dynamic;
%     RMS_50_50_dynamic;
%     RMS_25_75_dynamic;
%     RMS_0_1_dynamic];
% 
% %% RMSE static
% 
% RMSLFoot_1_0_static=rms(ErrorLFoot_1_0(1:Event_fence(sujet,1),:));
% RMSRFoot_1_0_static=rms(ErrorRFoot_1_0(1:Event_fence(sujet,1),:));
% 
% RMS_1_0_static=mean([RMSLFoot_1_0_static;RMSRFoot_1_0_static],1);
% 
% RMSLFoot_75_25_static=rms(ErrorLFoot_75_25(1:Event_fence(sujet,1),:));
% RMSRFoot_75_25_static=rms(ErrorRFoot_75_25(1:Event_fence(sujet,1),:));
% 
% RMS_75_25_static=mean([RMSLFoot_75_25_static;RMSRFoot_75_25_static],1);
% 
% RMSLFoot_50_50_static=rms(ErrorLFoot_50_50(1:Event_fence(sujet,1),:));
% RMSRFoot_50_50_static=rms(ErrorRFoot_50_50(1:Event_fence(sujet,1),:));
% 
% RMS_50_50_static=mean([RMSLFoot_50_50_static;RMSRFoot_50_50_static],1);
% 
% RMSLFoot_25_75_static=rms(ErrorLFoot_25_75(1:Event_fence(sujet,1),:));
% RMSRFoot_25_75_static=rms(ErrorRFoot_25_75(1:Event_fence(sujet,1),:));
% 
% RMS_25_75_static=mean([RMSLFoot_25_75_static;RMSRFoot_25_75_static],1);
% 
% RMSLFoot_0_1_static=rms(ErrorLFoot_0_1(1:Event_fence(sujet,1),:));
% RMSRFoot_0_1_static=rms(ErrorRFoot_0_1(1:Event_fence(sujet,1),:));
% 
% RMS_0_1_static=mean([RMSLFoot_0_1_static;RMSRFoot_0_1_static],1);
% 
% RMS_static=[RMS_1_0_static;
%     RMS_75_25_static;
%     RMS_50_50_static;
%     RMS_25_75_static;
%     RMS_0_1_static];
% 
% ErrorLFoot=[F_Pred.LFoot(:,1:3),M_Pred.LFoot(:,1:3)]-[F_Xp.LFoot(:,1:3),M_Xp.LFoot(:,1:3)];
% RMSLFoot=rms(ErrorLFoot);
% 
% ErrorRFoot=[F_Pred.RFoot(:,1:3),M_Pred.RFoot(:,1:3)]-[F_Xp.RFoot(:,1:3),M_Xp.RFoot(:,1:3)];
% RMSRFoot=rms(ErrorRFoot);

% %% Plot the result
% 
% bleu1=1/255*[151;255;255];
% bleu2=1/255*[141;238;238];
% bleu3=1/255*[121;205;205];
% bleu4=1/255*[82;139;139];
% bleu5=1/255*[47;79;79];
% 
% rouge1=1/255*[255;127;80];
% rouge2=1/255*[255;114;86];
% rouge3=1/255*[238;106;80];
% rouge4=1/255*[205;91;69];
% rouge5=1/255*[139;62;47];
% 
% figure
% set(gcf,'color','w')
% suptitle(['External Forces Prediction: ' filename ,newline, newline])
% %suptitle(['External Forces Prediction: ' ,newline, newline])
% Directions={'X','Y','Z'};
% for ii=1:3
%     subplot(1,6,ii)
%     % LFoot
%     plot(F_Pred_Threshold.(Solids{1})(:,ii),'LineWidth',2)
%     hold on
%     plot(F_Pred_InsoleDetection.(Solids{1})(:,ii),'LineWidth',2)
%     
%     plot(F_Xp.(Solids{1})(:,ii),'b--','LineWidth',2)
%     % RFoot
%     plot(F_Pred_Threshold.(Solids{2})(:,ii),'LineWidth',2)
%     plot(F_Pred_InsoleDetection.(Solids{2})(:,ii),'LineWidth',2)
%    
%     plot(F_Xp.(Solids{2})(:,ii),'r--','LineWidth',2)
%     
%     %xlim([1 203])
%     xlabel('Frames')
%     ylabel('Force (N)')
% 
%  %    title({['Force applied on Feets on the ' Directions{ii} '-direction'] ,...
%  %    ['RF RMSE : ' num2str(RMSRFoot(ii)) ' LF RMSE : ' num2str(RMSLFoot(ii)) ]})
% 
%     
%  %   title({['Force on the ' Directions{ii} '-direction']})
% end
% 
% for ii=4:6
%     subplot(1,6,ii)
%     % LFoot
%     plot(M_Pred_Threshold.(Solids{1})(:,ii-3),'LineWidth',2)
%     hold on
%     plot(M_Pred_InsoleDetection.(Solids{1})(:,ii-3),'LineWidth',2)
%    
%     plot(M_Xp.(Solids{1})(:,ii-3),'b--','LineWidth',2)
%     % RFoot
%     plot(M_Pred_Threshold.(Solids{2})(:,ii-3),'LineWidth',2)
%     plot(M_Pred_InsoleDetection.(Solids{2})(:,ii-3),'LineWidth',2)
%   
%     plot(M_Xp.(Solids{2})(:,ii-3),'r--','LineWidth',2)
%     
%     %xlim([1 203])
%     xlabel('Frames')
%     ylabel('Moment (Nm)')
%  %    title({['Moment applied on Feets on the ' Directions{ii-3} '-direction'],...
%   %       ['RF RMSE : ' num2str(RMSRFoot(ii)) ' LF RMSE : ' num2str(RMSLFoot(ii))]})
% 
%  %   title({['Moment on the ' Directions{ii-3} '-direction']})
% 
% end
% 
% legend('Left Foot Threshold','Left Foot Insole Detection',...
%     'Plateform Left Foot',...
%     'Right Foot Threshold','Right Foot Insole Detection',...
%     'Plateform Right Foot')

% On the graph, we can see that the predicted results are comparable to the
% experimental forces measured. However there is artefact when the foot is
% not on the plateforms. More work is needed on filters and algorithms to
% avoid it.
end

