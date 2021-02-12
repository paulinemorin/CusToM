% Post Processing Pre manip BAHAMaS
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
% Authors : Antoine Muller, Charles Pontonnier, Pierre Puchaud,
% Georges Dumont and Louise Demestre
%_______________________________________________________

% Get the numbers of solids on which the forces are applied
num_s=[22 28]; %Pieds droit et gauche

% Name of the trial treated
filename_pred='SQA0004';
filename_exp='sqa_004_5000Hz';

% Loading data
load([filename_pred '/ExperimentalData.mat']); 
load(fullfile(filename_pred,'ExternalForcesComputationResults.mat'))
GRF_Pred=ExternalForcesComputationResults.ExternalForcesPrediction;
load(filename_exp)
e_plateform = 0.01;
center = SensorCenter(filename_pred,e_plateform);

% Experimental and prediction times
Temps_exp = Channel_1_Data;
Temps_pred = ExperimentalData.Time;
Nb_frames_pred = numel(ExternalForcesComputationResults.ExternalForcesPrediction);

% Get the forces applied on the solids (prediction)
F_Pred = zeros(Nb_frames_pred,3);
M_Pred = zeros(Nb_frames_pred,3);
for i=1:numel(num_s)
    cur_s=num_s(i); %LFoot and RFoot
    for jj_f=1:Nb_frames_pred
        F_Pred(jj_f,:) = F_Pred(jj_f,:) + GRF_Pred(jj_f).fext(cur_s).fext(:,1)';
        M_Pred(jj_f,:) = M_Pred(jj_f,:) + GRF_Pred(jj_f).fext(cur_s).fext(:,2)';
    end
end

%Départ au maximum d'amplitude (saut pour synchronisation)
% [max_exp,ind_exp] = max(Channel_4_Data);
% [max_pred,ind_pred] = max(F_Pred(:,3));
[min_exp,ind_exp] = min(Channel_4_Data);
[min_pred,ind_pred] = min(F_Pred(:,3));
Temps_exp = Temps_exp(ind_exp:end)-Temps_exp(ind_exp);
nb_frames_exp = length(Temps_exp);
FX_c = [Channel_2_Data(ind_exp:end) zeros(nb_frames_exp,1) zeros(nb_frames_exp,1) ones(nb_frames_exp,1)];
FY_c = [zeros(nb_frames_exp,1) Channel_3_Data(ind_exp:end) zeros(nb_frames_exp,1) ones(nb_frames_exp,1)];
FZ_c = [zeros(nb_frames_exp,1) zeros(nb_frames_exp,1) Channel_4_Data(ind_exp:end) ones(nb_frames_exp,1)];
MX_c = [Channel_5_Data(ind_exp:end) zeros(nb_frames_exp,1) zeros(nb_frames_exp,1) ones(nb_frames_exp,1)];
MY_c = [zeros(nb_frames_exp,1) Channel_6_Data(ind_exp:end) zeros(nb_frames_exp,1) ones(nb_frames_exp,1)];
MZ_c = [zeros(nb_frames_exp,1) zeros(nb_frames_exp,1) Channel_7_Data(ind_exp:end) ones(nb_frames_exp,1)];
center = center(ind_pred:end,:);
Temps_pred = Temps_pred(ind_pred:end)-Temps_pred(ind_pred);
F_Pred = F_Pred(ind_pred:end,:);
M_Pred = M_Pred(ind_pred:end,:);

%Exprimer les forces mesurées dans le repère monde
load([filename_pred '_T21.mat']); %Matrice de passage du repère monde au repère de la plateforme
T21 = T21(:,:,ind_pred:end);
[TC2] = ChgmtReperesCapteurPlateforme(T21,center); %Matrice de passage du repère du capteur au repère de la plateforme
FX_interp = zeros(Nb_frames_pred-ind_pred+1,4);
FY_interp = zeros(Nb_frames_pred-ind_pred+1,4);
FZ_interp = zeros(Nb_frames_pred-ind_pred+1,4);
MX_interp = zeros(Nb_frames_pred-ind_pred+1,4);
MY_interp = zeros(Nb_frames_pred-ind_pred+1,4);
MZ_interp = zeros(Nb_frames_pred-ind_pred+1,4);
for k = 1:size(FX_c,2)
    FX_interp(:,k) = interp1(Temps_exp,FX_c(:,k),Temps_pred)';
    FY_interp(:,k) = interp1(Temps_exp,FY_c(:,k),Temps_pred)';
    FZ_interp(:,k) = interp1(Temps_exp,FZ_c(:,k),Temps_pred)';
    MX_interp(:,k) = interp1(Temps_exp,MX_c(:,k),Temps_pred)';
    MY_interp(:,k) = interp1(Temps_exp,MY_c(:,k),Temps_pred)';
    MZ_interp(:,k) = interp1(Temps_exp,MZ_c(:,k),Temps_pred)';
end
FX1 = zeros(length(center),4);
FY1 = zeros(length(center),4);
FZ1 = zeros(length(center),4);
MX1 = zeros(length(center),4);
MY1 = zeros(length(center),4);
MZ1 = zeros(length(center),4);
for j = 1:length(center)
    FX1(j,:) = (inv(T21(:,:,j))*inv(TC2(:,:,j))*FX_interp(j,:)')';
    FY1(j,:) = (inv(T21(:,:,j))*inv(TC2(:,:,j))*FY_interp(j,:)')';
    FZ1(j,:) = (inv(T21(:,:,j))*inv(TC2(:,:,j))*FZ_interp(j,:)')';
    MX1(j,:) = (inv(T21(:,:,j))*inv(TC2(:,:,j))*MX_interp(j,:)')';
    MY1(j,:) = (inv(T21(:,:,j))*inv(TC2(:,:,j))*MY_interp(j,:)')';
    MZ1(j,:) = (inv(T21(:,:,j))*inv(TC2(:,:,j))*MZ_interp(j,:)')';
end
F_exp = [FX1(:,1)+FY1(:,1)+ FZ1(:,1) FX1(:,2)+FY1(:,2)+ FZ1(:,2) FX1(:,3)+FY1(:,3)+ FZ1(:,3)];
M_exp = [MX1(:,1)+MY1(:,1)+ MZ1(:,1) MX1(:,2)+MY1(:,2)+ MZ1(:,2) MX1(:,3)+MY1(:,3)+ MZ1(:,3)];

%Transport des moments au centre du capteur
M_Pred = M_Pred + cross(-center,F_Pred);

%% Plot the result
Directions={'X','Y','Z'};
for ii=1:3
    subplot(2,3,ii)
    plot(Temps_pred,F_Pred(:,ii),'r-','LineWidth',2)
    hold on  
    plot(Temps_pred,-F_exp(:,ii),'b-','LineWidth',2)
    
    xlim([0 str2double(Channel_1_Header.MaxLevel(1:2))])
    xlabel(strcat('Time (', Channel_1_Header.Unit, ')'))
    ylabel('Force (N)')
    title(['Force applied on Feets on the ' Directions{ii} '-direction'])
end

for ii=4:6
    subplot(2,3,ii)
    plot(Temps_pred,M_Pred(:,ii-3),'r-','LineWidth',2)
    hold on
    plot(Temps_pred,-M_exp(:,ii-3),'b-','LineWidth',2)
    
    xlim([0 str2double(Channel_1_Header.MaxLevel(1:2))])
    xlabel(strcat('Time (', Channel_1_Header.Unit, ')'))
    ylabel('Moment (Nm)')
    title(['Moment applied on Feets on the ' Directions{ii-3} '-direction'])
end

% %Tracé des résultats pour le résumé ISB2021
% subplot(1,2,1)
% plot(Temps_pred,F_Pred(:,3),'r-','LineWidth',2)
% hold on
% plot(Temps_pred,-F_exp(:,3),'b-','LineWidth',2)
% axis square
% xlim([0 13])
% xlabel(strcat('Time (', Channel_1_Header.Unit, ')'),'FontSize',30)
% ylabel('Force (N)','FontSize',30)
% set(gca, 'fontsize', 25);
% title('FZ applied on feet','FontSize',30)
% 
% subplot(1,2,2)
% plot(Temps_pred,M_Pred(:,2),'r-','LineWidth',2)
% hold on
% axis square
% plot(Temps_pred,-M_exp(:,2),'b-','LineWidth',2)
% xlim([0 13])
% xlabel(strcat('Time (', Channel_1_Header.Unit, ')'),'FontSize',30)
% ylabel('Moment (Nm)','FontSize',30)
% set(gca, 'fontsize', 25);
% title('MY applied on feet','FontSize',30)

legend('Predicted forces','Experimental forces','Location','southeast')

RMSE_FX = sqrt(immse(F_Pred(:,1),-F_exp(:,1)))
RMSE_FY = sqrt(immse(F_Pred(:,2),-F_exp(:,2)))
RMSE_FZ = sqrt(immse(F_Pred(:,3),-F_exp(:,3)))
RMSE_MX = sqrt(immse(M_Pred(:,1),-M_exp(:,1)))
RMSE_MY = sqrt(immse(M_Pred(:,2),-M_exp(:,2)))
RMSE_MZ = sqrt(immse(M_Pred(:,3),-M_exp(:,3)))