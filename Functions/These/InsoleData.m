function [] = InsoleDataFile(filename)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

g=9.81;
f_cut=10;
freq=200;

%% Accélération marqueurs de pied
h = btkReadAcquisition([char(filename) '_all.c3d']);
[markers,markersInfo]=btkGetMarkers(h);

pos_RTAR =filt_data( markers.RTAR,f_cut,freq);
pos_LTAR =filt_data( markers.LTAR,f_cut,freq);

acc_RTAR = derivee2(0.0050,(derivee2(0.0050,pos_RTAR)));
acc_RTAR_norm = vecnorm(acc_RTAR');

acc_LTAR = derivee2(0.0050,(derivee2(0.0050,pos_RTAR)));
acc_LTAR_norm = vecnorm(acc_RTAR')/1000;

%% Insole data
data_full=importdata([filename '.txt']);

data = data_full.data;

LF_acc=data(:,18:20)'*(g);
RF_acc=data(:,40:42)'*(g);

% interpoler 200Hz
LF_acc_200 = double_freq(LF_acc)';
RF_acc_200 = double_freq(RF_acc)';

LF_acc_200(isnan(LF_acc_200))=0;
RF_acc_200(isnan(RF_acc_200))=0;

LF_acc_200_norm = vecnorm(LF_acc_200');
RF_acc_200_norm = vecnorm(LF_acc_200');

RF_acc_200_norm=filt_data(RF_acc_200_norm', f_cut, freq);
LF_acc_200_norm=filt_data(LF_acc_200_norm', f_cut, freq);

%% Synchronisation
delay_R=finddelay(acc_RTAR_norm,RF_acc_200_norm);
delay_L=finddelay(acc_LTAR_norm,LF_acc_200_norm);

%% Croper selon le .c3d final
load(fullfile(filename,'ExperimentalData.mat'))

InsoleData.text_data=data_full.textdata;

data_full.data_200=double_freq(data_full.data')';


InsoleData.data(:,2:23)=data_full.data_200(delay_L+ExperimentalData.FirstFrame:delay_L+ExperimentalData.LastFrame,2:23);
InsoleData.data(:,24:45)=data_full.data_200(delay_R+ExperimentalData.FirstFrame:delay_R+ExperimentalData.LastFrame,24:45);

save([filename '/InsoleData'],'InsoleData');
 
end

