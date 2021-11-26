function Foot_trial = Foot_trial_frame(filename,i,s)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
% i : frame
% s : solid

load(fullfile(filename,'ExperimentalData.mat'))

if s == 1 % pied gauche 
    Markers=[{'LHEE'},{'LTAR'},{'LTOE'},{'LTARI'}]';
else
     Markers=[{'RHEE'},{'RTAR'},{'RTOE'},{'RTARI'}]';
end

Name_list={ExperimentalData.MarkerPositions.name}';
Marker_liste=[];

for k= 1 : numel(Name_list)
    Marker_liste=[Marker_liste;Name_list{k}];
end

% Get the numbers of solids on which the forces are applied
[~,num_m]=intersect(Marker_liste,Markers,'stable');

position_HEE=ExperimentalData.MarkerPositions(num_m(1)).position_c3d;
position_TAR=ExperimentalData.MarkerPositions(num_m(2)).position_c3d;
position_TOE=ExperimentalData.MarkerPositions(num_m(3)).position_c3d;
position_TARI=ExperimentalData.MarkerPositions(num_m(4)).position_c3d;



y = (position_HEE(i,:)-position_TOE(i,:))/norm(position_HEE(i,:)-position_TOE(i,:));

if s==1 
    z_temp = (position_TARI(i,:)-position_TAR(i,:))/norm(position_TARI(i,:)-position_TAR(i,:));
else
    z_temp = (position_TAR(i,:)-position_TARI(i,:))/norm(position_TARI(i,:)-position_TAR(i,:));
end

x = cross(y, z_temp);
z = cross(x, y);

O = mean([position_TAR(i,:)',position_TARI(i,:)',position_HEE(i,:)',position_TOE(i,:)'],2);

Foot_trial = zeros(4);
Foot_trial(4,4)=1;
Foot_trial(1:3,4)= O;
Foot_trial(1:3,1:3)=[x',y',z'];

% figure
% hold on
% plot3(position_TAR(i,1),position_TAR(i,2),position_TAR(i,3),'*')
% plot3(position_TARI(i,1),position_TARI(i,2),position_TARI(i,3),'*')
% plot3(position_TOE(i,1),position_TOE(i,2),position_TOE(i,3),'*')
% plot3(position_HEE(i,1),position_HEE(i,2),position_HEE(i,3),'*')
% 
% legend('TAR','TARI','TOE','HEE')
% triad('matrix',Foot_trial,'scale',0.5,'Tag','Rmodel')

end

