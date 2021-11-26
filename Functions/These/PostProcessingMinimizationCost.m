% PostProcessing minimization cost

% Loading the Analysis file
load('AnalysisParameters.mat')
AnalysisParameters.ExternalForces.Options{1,1}='RFoot';
AnalysisParameters.ExternalForces.Options{2,1}='LFoot';

% Get the solid names on which the forces are applied
Solids = AnalysisParameters.ExternalForces.Options;
%Solids = ['RFoot','LFoot'];

% Loading the Biomechanicalmodel file
load('BiomechanicalModel.mat')

% Solid list extracted from the OsteoarticularModel
Solid_list = {BiomechanicalModel.OsteoArticularModel.name}';

% Get the numbers of solids on which the forces are applied
[~,num_s]=intersect(Solid_list,Solids);

% Name of the trial treated
filename='ChangmtDirection0001';

load(fullfile(filename,'ExperimentalData.mat'))
time = ExperimentalData.Time;
dt=time(2);

% Loading inverse dynamics results
load(fullfile(filename,'InverseKinematicsResults.mat')) 
JointCoo = InverseKinematicsResults.JointCoordinates;
JointVelo = derivee2(dt,JointCoo);

% Loading external forces computation
load(fullfile(filename,'ExternalForcesComputationResults.mat'))

% % Loading inverse dynamics results
% load(fullfile(filename,'InverseDynamicsResults_1.mat')) 
% torques_Xp = InverseDynamicsResults.JointTorques;
% 
% load(fullfile(filename,'InverseDynamicsResults_2.mat')) 
% torques_Pred = InverseDynamicsResults.JointTorques;

% Load experimental forces and forces from prediction algorithm.
GRF_Xp=ExternalForcesComputationResults.ExternalForcesExperiments;
GRF_Pred=ExternalForcesComputationResults.ExternalForcesPrediction;

% 
% load(fullfile(filename,'ExternalForcesComputationResults.mat'))
% GRF_Pred_int=ExternalForcesComputationResults.ExternalForcesPrediction;

% Number of frames
Nb_frames = numel(GRF_Xp);

% Get the forces applied on the solids
for ii=1:numel(num_s)
    cur_s=num_s(ii); %LFoot and RFoot
    for jj_f=1:Nb_frames
        % Prediction algorithms
        F_Pred.(Solids{ii})(jj_f,:) = GRF_Pred(jj_f).fext(cur_s).fext(:,1)';
        M_Pred.(Solids{ii})(jj_f,:) = GRF_Pred(jj_f).fext(cur_s).fext(:,2);
        
%         F_Pred_int.(Solids{ii})(jj_f,:) = GRF_Pred_int(jj_f).fext(cur_s).fext(:,1)';
%         M_Pred_int.(Solids{ii})(jj_f,:) = GRF_Pred_int(jj_f).fext(cur_s).fext(:,2);
%         
        % Experimental results
        F_Xp.(Solids{ii})(jj_f,:) = GRF_Xp(jj_f).fext(cur_s).fext(:,1)';
        M_Xp.(Solids{ii})(jj_f,:) = GRF_Xp(jj_f).fext(cur_s).fext(:,2)';
    end
end

%% Calcul des fonctions de coût
% 
% % external forces term 
% ext_Pred = zeros(1, Nb_frames);
% ext_Xp = zeros(1, Nb_frames);
% 
% for ii=1:numel(num_s)
%     cur_s=num_s(ii); %LFoot and RFoot
%     for jj_f=1:Nb_frames
%         % Prediction algorithms
%         ext_Pred (jj_f) =  ext_Pred (jj_f) + sum(F_Pred.(Solids{ii})(jj_f,:).^2);
% 
%         % Experimental results
%         ext_Xp  (jj_f) = ext_Xp  (jj_f) + sum(F_Xp.(Solids{ii})(jj_f,:).^2);
%     end
% end
% 
% % internal forces term
% int_Pred = zeros(1, Nb_frames);
% int_Xp = zeros(1, Nb_frames);
% 
% for ii=1:numel(num_s)
%     for jj_f=1:Nb_frames
%         % Prediction algorithms
%         int_Pred (jj_f) = sum(torques_Pred(:,jj_f).^2);
% 
%         % Experimental results
%         int_Xp (jj_f) = sum(torques_Xp(:,jj_f).^2);
%     end
% end
% 
% % internal forces power term
% int_power_Pred = zeros(1, Nb_frames);
% int_power_Xp = zeros(1, Nb_frames);
% 
% for ii=1:numel(num_s)
%     for jj_f=1:Nb_frames
%         % Prediction algorithms
%         int_power_Pred (jj_f) = abs(torques_Pred(:,jj_f)).'*abs(JointVelo(:,jj_f));
% 
%         % Experimental results
%         int_power_Xp (jj_f) = abs(torques_Xp(:,jj_f)).'*abs(JointVelo(:,jj_f));
%     end
% end
% 
% % internal forces jerk term
% int_jerk_Pred = zeros(1, Nb_frames);
% int_jerk_Xp = zeros(1, Nb_frames);
% 
% 
% Acc_Pred = derivee2(dt,torques_Pred);
% Acc_Xp = derivee2(dt,torques_Xp);
% 
% for ii=1:numel(num_s)
%     for jj_f=1:Nb_frames
%         % Prediction algorithms
%         int_jerk_Pred (jj_f) = sum(Acc_Pred(:,jj_f).^2);
% 
%         % Experimental results
%         int_jerk_Xp (jj_f) = sum(Acc_Xp(:,jj_f).^2);
%     end
% end

%% Plot results

% figure
% % External forces term
% subplot(1,2,1)
% plot(ext_Pred,'b-','LineWidth',2)
% hold on
% plot(ext_Xp,'b--','LineWidth',2)
%     
% xlim([1 Nb_frames])
% xlabel('Frames')
% ylabel('Coût associé aux efforts extérieurs (N^2)')
% 
% % Internal forces term
% subplot(1,2,2)
% plot(int_Pred,'b-','LineWidth',2)
% hold on
% plot(int_Xp,'b--','LineWidth',2)
%     
% xlim([1 Nb_frames])
% xlabel('Frames')
% ylabel('Coût associé aux efforts intérieurs ((N.m)^2)')
%    

%figure
% External forces term
%plot(ext_Pred,'b-','LineWidth',2)
%hold on
%plot(ext_Xp,'b--','LineWidth',2)
    
%xlim([1 Nb_frames])
%xlabel('Frames')

% Internal forces term
%plot(int_Pred,'r-','LineWidth',2)
%plot(int_Xp,'r--','LineWidth',2)
    

% Internal forces power term
%plot(int_power_Pred,'g-','LineWidth',2)
%plot(int_power_Xp,'g--','LineWidth',2)

% Internal forces power term
%plot(int_jerk_Pred,'c-','LineWidth',2)
%plot(int_jerk_Xp,'c--','LineWidth',2)

% 
% legend('Coût associé aux efforts extérieurs prédits (N^2)',...
%     'Coût associé aux efforts extérieurs expérimentaux (N^2)',...
%     'Coût associé aux efforts intérieurs prédits ((N.m)^2)',...
%     'Coût associé aux efforts intérieurs expérimentaux ((N.m)^2)',...
%     'Coût associé aux puissances des efforts intérieurs prédits (N.m/rad.s-1)',...
%     'Coût associé aux puissances des efforts intérieurs expérimentaux (N.m/rad.s-1)',...
%     'Coût associé aux jerk des efforts intérieurs prédits (N.m/rad.s-1)',...
%     'Coût associé aux jerk des efforts intérieurs expérimentaux (N.m/rad.s-1)')
% 

% Erreur 
%% Plot the result
figure
Directions={'X','Y','Z'};
for ii=1:3
    subplot(2,3,ii)
    % LFoot
    plot(F_Pred.(Solids{1})(:,ii),'b-','LineWidth',2)
    hold on
    plot(F_Pred_int.(Solids{1})(:,ii),'g-','LineWidth',2)
    plot(F_Xp.(Solids{1})(:,ii),'b--','LineWidth',2)
    % RFoot
    plot(F_Pred.(Solids{2})(:,ii),'r-','LineWidth',2)
    plot(F_Pred_int.(Solids{2})(:,ii),'g-','LineWidth',2)
    plot(F_Xp.(Solids{2})(:,ii),'r--','LineWidth',2)
    
    %xlim([1 203])
    xlabel('Frames')
    ylabel('Force (N)')
    title(['Force applied on Feets on the ' Directions{ii} '-direction'])
end

for ii=4:6
    subplot(2,3,ii)
    % LFoot
    plot(M_Pred.(Solids{1})(:,ii-3),'b-','LineWidth',2)
    hold on
    plot(M_Pred_int.(Solids{1})(:,ii-3),'g-','LineWidth',2)
    plot(M_Xp.(Solids{1})(:,ii-3),'b--','LineWidth',2)
    % RFoot
    plot(M_Pred.(Solids{2})(:,ii-3),'r-','LineWidth',2)
    plot(M_Pred_int.(Solids{2})(:,ii-3),'g-','LineWidth',2)
    plot(M_Xp.(Solids{2})(:,ii-3),'r--','LineWidth',2)
    
    %xlim([1 203])
    xlabel('Frames')
    ylabel('Moment (Nm)')
    title(['Moment applied on Feets on the ' Directions{ii-3} '-direction'])
end

legend('Predicted on Left Foot','Experimental on Left Foot',...
    'Predicted on Right Foot','Experimental on Right Foot')

% On the graph, we can see that the predicted results are comparable to the
% experimental forces measured. However there is artefact when the foot is
% not on the plateforms. More work is needed on filters and algorithms to
% avoid it.








figure
Directions={'X','Y','Z'};
for ii=1:3
    subplot(2,3,ii)
    % LFoot
    plot(F_Pred.(Solids{1})(:,ii)-F_Xp.(Solids{1})(:,ii),'b-','LineWidth',2)
    hold on
    plot(F_Pred_int.(Solids{1})(:,ii)-F_Xp.(Solids{1})(:,ii),'b--','LineWidth',2)
    % RFoot
    plot(F_Pred.(Solids{2})(:,ii)-F_Xp.(Solids{2})(:,ii),'r-','LineWidth',2)
    plot(F_Pred_int.(Solids{2})(:,ii)-F_Xp.(Solids{2})(:,ii),'r--','LineWidth',2)
    
    %xlim([1 203])
    xlabel('Frames')
    ylabel('Force (N)')
    title(['Force applied on Feets on the ' Directions{ii} '-direction'])
end

for ii=4:6
    subplot(2,3,ii)
    % LFoot
    plot(M_Pred.(Solids{1})(:,ii-3)-M_Xp.(Solids{1})(:,ii-3),'b-','LineWidth',2)
    hold on
    plot(M_Pred_int.(Solids{1})(:,ii-3)-M_Xp.(Solids{1})(:,ii-3),'b--','LineWidth',2)
    % RFoot
    plot(M_Pred.(Solids{2})(:,ii-3)-M_Xp.(Solids{2})(:,ii-3),'r-','LineWidth',2)
    plot(M_Pred_int.(Solids{2})(:,ii-3)-M_Xp.(Solids{2})(:,ii-3),'r--','LineWidth',2)
    
    %xlim([1 203])
    xlabel('Frames')
    ylabel('Moment (Nm)')
    title(['Moment applied on Feets on the ' Directions{ii-3} '-direction'])
end

legend('Error on Left Foot (external force minimization)','Error on Left Foot (hybride minimization)',...
    'Error on Right Foot (external force minimization)','Error on Right Foot (hybride minimization)')

% On the graph, we can see that the predicted results are comparable to the
% experimental forces measured. However there is artefact when the foot is
% not on the plateforms. More work is needed on filters and algorithms to
% avoid it.






figure
Directions={'X','Y','Z'};
for ii=1:3
    subplot(1,3,ii)
    % Both Foot
    plot(F_Pred.(Solids{1})(:,ii)+F_Pred.(Solids{2})(:,ii),'b-','LineWidth',2)
    hold on
    plot(F_Pred_int.(Solids{1})(:,ii)+F_Pred_int.(Solids{2})(:,ii),'r-','LineWidth',2)
    plot(F_Xp.(Solids{1})(:,ii)+F_Xp.(Solids{2})(:,ii),'g-','LineWidth',2)
    
    %xlim([1 203])
    xlabel('Frames')
    ylabel('Force (N)')
    title(['Force applied on Feets on the ' Directions{ii} '-direction'])
end

legend('Error on Left Foot (external force minimization)','Error on Left Foot (hybride minimization)',...
    'Error on Right Foot (external force minimization)','Error on Right Foot (hybride minimization)')

% On the graph, we can see that the predicted results are comparable to the
% experimental forces measured. However there is artefact when the foot is
% not on the plateforms. More work is needed on filters and algorithms to
% avoid it.








figure
set(gcf,'color','w')
Directions={'X','Y','Z'};
for ii=1:3
    subplot(1,3,ii)
    % Both Foot
    plot(abs(F_Pred.(Solids{1})(:,ii)+F_Pred.(Solids{2})(:,ii)-F_Xp.(Solids{1})(:,ii)-F_Xp.(Solids{2})(:,ii)),'b-','LineWidth',2)
    hold on
    plot(abs(F_Pred_int.(Solids{1})(:,ii)+F_Pred_int.(Solids{2})(:,ii)-F_Xp.(Solids{1})(:,ii)-F_Xp.(Solids{2})(:,ii)),'r-','LineWidth',2)
    
    %xlim([1 203])
    xlabel('Frames')
    ylabel('Force (N)')
    title(['Force applied on Feets on the ' Directions{ii} '-direction'])
end

legend('Error on Global Equilibrium (external force minimization)','Error on Global Equilibrium (hybride minimization)')

% On the graph, we can see that the predicted results are comparable to the
% experimental forces measured. However there is artefact when the foot is
% not on the plateforms. More work is needed on filters and algorithms to
% avoid it.





figure
set(gcf,'color','w')
Directions={'X','Y','Z'};
for ii=1:3
    subplot(1,3,ii)
    % Both Foot
    plot(abs(F_Pred.(Solids{1})(:,ii)-F_Xp.(Solids{1})(:,ii))+abs(F_Pred.(Solids{2})(:,ii)-F_Xp.(Solids{2})(:,ii)),'b-','LineWidth',2)
    hold on
    plot(abs(F_Pred_int.(Solids{1})(:,ii)-F_Xp.(Solids{1})(:,ii))+abs(F_Pred_int.(Solids{2})(:,ii)-F_Xp.(Solids{2})(:,ii)),'r-','LineWidth',2)
    
    %xlim([1 203])
    xlabel('Frames')
    ylabel('Force (N)')
    title(['Force applied on Feets on the ' Directions{ii} '-direction'])
end

legend('Error totale (external force minimization)','Error totale (hybride minimization)')

% On the graph, we can see that the predicted results are comparable to the
% experimental forces measured. However there is artefact when the foot is
% not on the plateforms. More work is needed on filters and algorithms to
% avoid it.


figure
Directions={'X','Y','Z'};
set(gcf,'color','w')
for ii=1:3
    subplot(1,3,ii)
    % Both Foot
    plot(abs(F_Pred.(Solids{1})(:,ii)-F_Xp.(Solids{1})(:,ii))+abs(F_Pred.(Solids{2})(:,ii)-F_Xp.(Solids{2})(:,ii))...
        -abs(abs(F_Pred.(Solids{1})(:,ii)+F_Pred.(Solids{2})(:,ii))-abs(F_Xp.(Solids{1})(:,ii)+F_Xp.(Solids{2})(:,ii))),'b-','LineWidth',2)
    hold on
    plot(abs(F_Pred_int.(Solids{1})(:,ii)-F_Xp.(Solids{1})(:,ii))+abs(F_Pred_int.(Solids{2})(:,ii)-F_Xp.(Solids{2})(:,ii))...
        -abs(abs(F_Pred_int.(Solids{1})(:,ii)+F_Pred_int.(Solids{2})(:,ii))-abs(F_Xp.(Solids{1})(:,ii)+F_Xp.(Solids{2})(:,ii))),'r-','LineWidth',2)
    
    %xlim([1 203])
    xlabel('Frames')
    ylabel('Force (N)')
    title(['Force applied on Feets on the ' Directions{ii} '-direction'])
end

legend('Error repartition (external force minimization)','Error repartition (hybride minimization)')

% On the graph, we can see that the predicted results are comparable to the
% experimental forces measured. However there is artefact when the foot is
% not on the plateforms. More work is needed on filters and algorithms to
% avoid it.






