%function erreur_COP = Comparer_COP(ExternalForcesComputationResults)

%% Charger les données

% Loading the Analysis file
load('AnalysisParameters.mat')

% Get the solid names on which the forces are applied
%Solids = AnalysisParameters.ExternalForces.Options;
Solids = {'LFoot';'RFoot'};

% Loading the Biomechanicalmodel file
load('BiomechanicalModel.mat')

% Solid list extracted from the OsteoarticularModel
Solid_list = {BiomechanicalModel.OsteoArticularModel.name}';

% Get the numbers of solids on which the forces are applied
[~,num_s]=intersect(Solid_list,Solids);


% Name of the trial treated
filename='ChgtDirection04';

% Loading external forces computation
load(fullfile(filename,'ExternalForcesComputationResults.mat'))

% Load experimental forces and forces from prediction algorithm.
GRF_Xp=ExternalForcesComputationResults.ExternalForcesExperiments;
GRF_Pred=ExternalForcesComputationResults.ExternalForcesPrediction;
%GRF_COP=ExternalForcesComputationResults.ExternalForcesPrediction_plus_COP;

% Number of frames
Nb_frames = numel(GRF_Xp);

% Get the forces applied on the solids
for ii=1:numel(num_s)
    cur_s=num_s(ii); %LFoot and RFoot
    for jj_f=1:Nb_frames
        % Prediction algorithms
        F_Pred.(Solids{ii})(jj_f,:) = GRF_Pred(jj_f).fext(cur_s).fext(:,1)';
        M_Pred.(Solids{ii})(jj_f,:) = GRF_Pred(jj_f).fext(cur_s).fext(:,2);
        % Experimental results
        F_Xp.(Solids{ii})(jj_f,:) = GRF_Xp(jj_f).fext(cur_s).fext(:,1)';
        M_Xp.(Solids{ii})(jj_f,:) = GRF_Xp(jj_f).fext(cur_s).fext(:,2)';
%         % Prediction COP
%         F_COP.(Solids{ii})(jj_f,:) = GRF_COP(jj_f).fext(cur_s).fext(:,1)';
%         M_COP.(Solids{ii})(jj_f,:) = GRF_COP(jj_f).fext(cur_s).fext(:,2);
    end
end

%% pour un contact, pied droit  

COP_exp_1 = zeros( Nb_frames, 3);
COP_pred_1 = zeros( Nb_frames, 3);
erreur_COP_1 = zeros(Nb_frames , 1);

for i=1:Nb_frames
     if norm(ExternalForcesComputationResults.ExternalForcesExperiments(i).Visual(4:6,1)) > 50
        COP_exp_1(i,:)=ExternalForcesComputationResults.ExternalForcesExperiments(i).Visual(1:3,1).';
        COP_pred_1(i,:)=ExternalForcesComputationResults.ExternalForcesPrediction(i).Visual(1:3,1).';
        %COP_impose(i,:)=ExternalForcesComputationResults.ExternalForcesPrediction_plus_COP(i).Visual(1:3,1).';
        %ExternalForcesComputationResults.ExternalForcesExperiments(i).Visual(6,1)
%     if ExternalForcesComputationResults.ExternalForcesPrediction(i).Contact_Pt_Nb==0
%         erreur_COP(i)=0;
%         erreur_COP_impo(i)=0;
%     else
         erreur_COP_1(i)=norm(COP_exp_1(i,:)-COP_pred_1(i,:));
%         erreur_COP_impo(i)=norm(COP_exp(i,:)-COP_impose(i,:));
%    end
     end
end

plot(erreur_COP_1(:),'b+')
hold on
%plot(erreur_COP_impo(:),'r+')
xlabel('Frames')
ylabel('COP location error(m mais ça me parait énorme)')
title('COP location error for the left foot')
legend('erreur COP prédit', 'erreur COP poutant imposé')


%% pour un contact, pied gauche  

COP_exp_2 = zeros( Nb_frames, 3);
COP_pred_2 = zeros( Nb_frames, 3);
erreur_COP_2 = zeros(Nb_frames , 1);

for i=1:Nb_frames
     if norm(ExternalForcesComputationResults.ExternalForcesExperiments(i).Visual(4:6,2)) > 50
        COP_exp_2(i,:)=ExternalForcesComputationResults.ExternalForcesExperiments(i).Visual(1:3,2).';
        COP_pred_2(i,:)=ExternalForcesComputationResults.ExternalForcesPrediction(i).Visual(1:3,2).';
        %COP_impose(i,:)=ExternalForcesComputationResults.ExternalForcesPrediction_COP(i).Visual(1:3,2).';
        ExternalForcesComputationResults.ExternalForcesExperiments(i).Visual(6,2)
%     if (abs(ExternalForcesComputationResults.ExternalForcesExperiments(i).Visual(6,2))<50)
%         erreur_COP(i)=0;
%         erreur_COP_impo(i)=0;
%     else
         erreur_COP_2(i)=norm(COP_exp_2(i,:)-COP_pred_2(i,:));
%         erreur_COP_impo(i)=norm(COP_exp(i,:)-COP_impose(i,:));
%     end
end
end

plot(erreur_COP_2(:))
hold on
%plot(erreur_COP_impo(:))
xlabel('Frames')
ylabel('COP location error(m mais ça me parait énorme)')
title('COP location error for the left foot')
legend('erreur COP prédit', 'erreur COP poutant imposé')


%% Plot error of the results
max_forces=max(abs(F_Pred.(Solids{1})(:,:)-F_Xp.(Solids{1})(:,:)));

Directions={'X','Y','Z'};
for ii=1:3  
    subplot(2,3,ii)
    % LFoot
    plot(F_Pred.(Solids{1})(:,ii)-F_Xp.(Solids{1})(:,ii),'b-','LineWidth',2)
    hold on
%    plot(F_COP.(Solids{1})(:,ii)-F_Xp.(Solids{1})(:,ii),'r-','LineWidth',2)
%    plot(erreur_COP(:)*max_forces(ii)/max(erreur_COP(:)))
    
    xlim([0 203])
    xlabel('Frames')
    ylabel('Force (N)')
    title(['Force error applied on Left Foot on the ' Directions{ii} '-direction'])
end

max_moment=max(abs(M_Pred.(Solids{1})(:,:)-M_Xp.(Solids{1})(:,:)));

for ii=4:6
    subplot(2,3,ii)
    % LFoot
    plot(M_Pred.(Solids{1})(:,ii-3)-M_Xp.(Solids{1})(:,ii-3),'b-','LineWidth',2)
    hold on
%    plot(M_COP.(Solids{1})(:,ii-3)-M_Xp.(Solids{1})(:,ii-3),'r-','LineWidth',2)
 %   plot(erreur_COP(:)*max_moment(ii-3)/max(erreur_COP(:)))
    
    xlim([0 203])
    xlabel('Frames')
    ylabel('Moment (Nm)')
    title(['Moment error applied on Left Foot on the ' Directions{ii-3} '-direction'])
end

legend('Error prediction on external forces','Error prediction on external forces with COP location', 'Error COP location prediction')

%%
%% Plot error of the results, pied gauche
max_forces=max(abs(F_Pred.(Solids{2})(:,:)-F_Xp.(Solids{2})(:,:)));

Directions={'X','Y','Z'};
for ii=1:3
    subplot(2,3,ii)
    % LFoot
    plot(F_Pred.(Solids{2})(:,ii)-F_Xp.(Solids{2})(:,ii),'b-','LineWidth',2)
    hold on
%    plot(F_COP.(Solids{2})(:,ii)-F_Xp.(Solids{2})(:,ii),'r-','LineWidth',2)
%    plot(erreur_COP(:)*max_forces(ii)/max(erreur_COP(:)))
    
    xlim([0 203])
    xlabel('Frames')
    ylabel('Force (N)')
    title(['Force error applied on Left Foot on the ' Directions{ii} '-direction'])
end

max_moment=max(abs(M_Pred.(Solids{2})(:,:)-M_Xp.(Solids{2})(:,:)));

for ii=4:6
    subplot(2,3,ii)
    % LFoot
    plot(M_Pred.(Solids{2})(:,ii-3)-M_Xp.(Solids{2})(:,ii-3),'b-','LineWidth',2)
    hold on
%    plot(M_COP.(Solids{2})(:,ii-3)-M_Xp.(Solids{2})(:,ii-3),'r-','LineWidth',2)
%    plot(erreur_COP(:)*max_moment(ii-3)/max(erreur_COP(:)))
    
    xlim([0 203])
    xlabel('Frames')
    ylabel('Moment (Nm)')
    title(['Moment error applied on Left Foot on the ' Directions{ii-3} '-direction'])
end

legend('Error prediction on external forces','Error prediction on external forces with COP location', 'Error COP location prediction')


   figure
hold on
plot(COP_exp_1(:,2))
plot(COP_exp_2(:,2))
plot(COP_pred_1(:,2))
plot(COP_pred_2(:,2))
 legend('exp droit','exp gauche', 'pred droit','pred gauche')


figure
hold on
plot(COP_exp_1(:,1))
plot(COP_exp_2(:,1))
plot(COP_pred_1(:,1))
plot(COP_pred_2(:,1))
 legend('exp droit','exp gauche', 'pred droit','pred gauche')
%end