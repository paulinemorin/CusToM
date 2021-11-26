function [] = InsoleDataFile(filename)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

g=9.81;
f_cut=10;
freq=200;

%% Accélération marqueurs de pied
h = btkReadAcquisition([char(filename) '_all.c3d']);
[markers,markersInfo]=btkGetMarkers(h);

first_all=btkGetFirstFrame(h);

for k=1:3
    pos_R(:,k) = mean([markers.RTAR(:,k), markers.RTARI(:,k), markers.RTOE(:,k), markers.RHEE(:,k)],2);
    pos_R(:,k) = filt_data( pos_R(:,k),f_cut,freq);
end

for k=1:3
    pos_L(:,k) = mean([markers.LTAR(:,k), markers.LTARI(:,k), markers.LTOE(:,k), markers.LHEE(:,k)],2);
    pos_L(:,k) = filt_data( pos_L(:,k),f_cut,freq);
end


acc_R = derivee2(0.0050,(derivee2(0.0050,pos_R)));
acc_R_norm = vecnorm(acc_R')/1000;

acc_L = derivee2(0.0050,(derivee2(0.0050,pos_L)));
acc_L_norm = vecnorm(acc_L')/1000;


%% Insole data
data_full=importdata([filename '.txt']);

data = data_full.data;

for i = 2 :  size(data,1)-1
    for j=1: size(data,2)
        if isnan(data(i,j))
            data(i,j) = (data(i-1,j)+data(i+1,j))/2; 
        end
    end
end

LF_acc=data(:,18:20)'*(g);
RF_acc=data(:,40:42)'*(g);

% interpoler 200Hz
LF_acc_200 = double_freq(LF_acc)';
RF_acc_200 = double_freq(RF_acc)';

LF_acc_200(isnan(LF_acc_200))=0;
RF_acc_200(isnan(RF_acc_200))=0;

LF_acc_200_norm = vecnorm(LF_acc_200');
RF_acc_200_norm = vecnorm(RF_acc_200');

RF_acc_200_norm=filt_data(RF_acc_200_norm', f_cut, freq);
LF_acc_200_norm=filt_data(LF_acc_200_norm', f_cut, freq);

%% Synchronisation
delay_R=(finddelay(acc_R_norm,RF_acc_200_norm));
delay_L=(finddelay(acc_L_norm,LF_acc_200_norm));

% if norm(delay_R-delay_L)>50
%     delay_R=(finddelay(acc_R_norm,RF_acc_200_norm(1:length(RF_acc_200_norm)/2)));
%     delay_L=(finddelay(acc_L_norm,LF_acc_200_norm(1:length(LF_acc_200_norm)/2)));
% end

% figure 
% hold on
% plot(LF_acc_200_norm)
% plot(RF_acc_200_norm)
% plot(acc_R_norm)
% plot(acc_L_norm)
% legend('Semelle gauche','Semelle Droite','Mouvement Gauche','Mouvement Droit')


%% Croper selon le .c3d final
load(fullfile(filename,'ExperimentalData.mat'))

InsoleData.text_data=data_full.textdata;

for i = 1 : size(data,2)
    data_200(:,i)=double_freq(data(:,i)')';
end

if delay_L+ExperimentalData.LastFrame>length(data_200)
    delay_L=delay_R;
end

if delay_R+ExperimentalData.LastFrame>length(data_200)
    delay_R=delay_L;
end

InsoleData.data(:,2:23)=data_200(delay_L+ExperimentalData.FirstFrame-first_all:delay_L+ExperimentalData.LastFrame-first_all,2:23);
InsoleData.data(:,24:45)=data_200(delay_R+ExperimentalData.FirstFrame-first_all:delay_R+ExperimentalData.LastFrame-first_all,24:45);

delay_L_force=0;
delay_R_force=0;

% Loading external forces computation
load(fullfile(filename,'ExternalForcesComputationResults.mat'))


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

% Load experimental forces and forces from prediction algorithm.
GRF_Xp=ExternalForcesComputationResults.ExternalForcesExperiments;

% Get the forces applied on the solids
for ii=1:numel(num_s)
    cur_s=num_s(ii); %LFoot and RFoot
    for jj_f=1:ExperimentalData.LastFrame-ExperimentalData.FirstFrame
        % Experimental results
        F_Xp.(Solids{ii})(jj_f,:) = GRF_Xp(jj_f).fext(cur_s).fext(:,1)';
    end
end

%% External forces data
f_cut=5;

LF_force=InsoleData.data(:,21);
RF_force=InsoleData.data(:,43);

% InsoleData.data(end,21)=InsoleData.data(1,21);
% InsoleData.data(end,43)=InsoleData.data(1,43);

% LF_force=filt_data(InsoleData.data(:,21), f_cut, freq);
% RF_force=filt_data(InsoleData.data(:,43), f_cut, freq);


%% Synchronisation
delay_L_force=finddelay(F_Xp.LFoot(:,end),LF_force);
delay_R_force=finddelay(F_Xp.RFoot(:,end),RF_force);

if delay_L+delay_L_force+ExperimentalData.LastFrame>length(data_200) || delay_L_force>50
    delay_L_force=delay_R_force;
end

if delay_R+delay_L_force+ExperimentalData.LastFrame>length(data_200) || delay_R_force>50
    delay_R_force=delay_L_force;
end

% delay_R_force=0;
% delay_L_force=0;

InsoleData.data(:,2:23)=data_200(delay_L+delay_L_force+ExperimentalData.FirstFrame-first_all:delay_L+delay_L_force+ExperimentalData.LastFrame-first_all,2:23);
InsoleData.data(:,24:45)=data_200(delay_R+delay_R_force+ExperimentalData.FirstFrame-first_all:delay_R+delay_R_force+ExperimentalData.LastFrame-first_all,24:45);



% figure 
% hold on
% plot(InsoleData.data(:,21))
% plot(InsoleData.data(:,43))
% plot(F_Xp.LFoot(:,end))
% plot(F_Xp.RFoot(:,end))
% legend('Semelle gauche','Semelle Droite','Mouvement Gauche','Mouvement Droit')








save([filename '/InsoleData'],'InsoleData');
 
end

