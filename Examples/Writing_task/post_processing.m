%data post processing

% DPD = Deltoid Posterior : 5
% DAD = Deltoid Anterior : 3
% DAG = Deltoid Medial : 8
% BicD = Biceps : 2
% TriMD = Triceps medial : 7
% EPD = Wrist extensor muscle group (extensor digitorium, extensor carpi
% radialis, and ulnaris) : 1
% FPD = Wrist flexor muscle group  (flexor digitorium, flexor carpi
% radialis, and ulnaris)  : 6
% TraD = Trapezius : 4
%close all;

CusMuscles={{},{},{},{},{},{}};
nom_musclesEMG={'Trapeze',' Deltoid Posterior ',' Deltoid Anterior ','Deltoid Medial','Biceps ','Triceps Medial','Wrist extensor','Wrist flexor '};

cmap=colormap(colorcube(29));
set(groot, 'DefaultAxesColorOrder', cmap,'DefaultAxesFontSize',20,'DefaultLineLineWidth',2);
close all;

load('THI_trajectoire0008/MuscleForcesComputationResults.mat');
load('THI_trajectoire0008.mat');
load('THI_trajectoire0008/ExperimentalData.mat');
load('BiomechanicalModel.mat');

EMG_interest=THI_trajectoire0008.Analog.Data ;


MusCus_tot=[CusMuscles{:}];
MusCus_tot=[MusCus_tot{:}];
EMG_Cus=MuscleForcesComputationResults.MuscleActivations(MusCus_tot,:)';
nom_muscles={BiomechanicalModel.Muscles(MusCus_tot).name};

time=(0:size(EMG_interest,2)-1)*1/THI_trajectoire0008.Analog.Frequency;
time_capt=ExperimentalData.Time;


placement=1;
 for k= 1:length(nom_musclesEMG)
    RMS_EMG(:,k) = MYO_RMS(EMG_interest(k,:)',1000);
    RMS_EMG(:,k) =RMS_EMG(:,k) /max(RMS_EMG(:,k) );
%         for j=1:length(CusMuscles{k})
%            EMG_CusNorm(:,placement)=EMG_Cus(:,placement) /max(EMG_Cus(:,placement) );
%            placement=placement+1;
%         end
end

placement=1;
figure()
for k=1:length(nom_musclesEMG)
    subplot(4,2,k)
    plot(time,RMS_EMG(:,k))
    hold on
    leg={'Mesure'};
%     for j=1:length(CusMuscles{k})
%         plot(time_capt,EMG_CusNorm(:,placement))
%         R=corrcoef(RMS_EMG(1:10:end,k),EMG_CusNorm(:,placement));
%         leg=[leg,[nom_muscles{placement}, ', corr : ',num2str(R(1,2)),' max :', num2str(max(EMG_Cus(:,placement)))]]; 
%         placement=placement+1;
%     end
    title(nom_musclesEMG{k})
    legend(leg)
    xlabel("temps (s)")
end






% for k=1:length(CusMuscles)
%     R=corrcoef(RMS_EMG(1:10:end,k),EMG_Cus(:,k));
%     Rfin(k)=R(1,2);
% end
