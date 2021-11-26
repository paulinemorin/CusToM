function [Rot_R, R_pos, Rot_L, L_pos] = FootFrameMarker(filename, i)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
load([filename '/ExperimentalData.mat']); %#ok<LOAD>

load(fullfile(filename,'ExperimentalData.mat'))

Markers=[{'RHEE'},{'RTOE'}]';
Name_list={ExperimentalData.MarkerPositions.name}';
Marker_liste=[];

for k= 1 : numel(Name_list)
    Marker_liste=[Marker_liste;Name_list{k}];
end

% Get the numbers of solids on which the forces are applied
[~,num_m]=intersect(Marker_liste,Markers);

position_RHEE=ExperimentalData.MarkerPositions(num_m(1)).position_c3d(i,:);
position_RTOE=ExperimentalData.MarkerPositions(num_m(2)).position_c3d(i,:);

Markers=[{'RTAR'},{'RTARI'}]';
Name_list={ExperimentalData.MarkerPositions.name}';
Marker_liste=[];

for k= 1 : numel(Name_list)
    Marker_liste=[Marker_liste;Name_list{k}];
end

% Get the numbers of solids on which the forces are applied
[~,num_m]=intersect(Marker_liste,Markers);
position_RTAR=ExperimentalData.MarkerPositions(num_m(1)).position_c3d(i,:);
position_RTARI=ExperimentalData.MarkerPositions(num_m(2)).position_c3d(i,:);

Pts_R= [position_RTARI;position_RTAR;position_RHEE;position_RTOE];

R_pos = mean(Pts_R);

Z_R = position_RTAR - position_RTARI; 
Z_R = Z_R / norm( Z_R);
Y_R = position_RHEE - position_RTOE;
Y_R = Y_R / norm (Y_R);
X_R = cross (Y_R, Z_R);
Z_R = cross (X_R, Y_R);
Rot_R = [X_R', Y_R', Z_R'];





Markers=[{'LHEE'},{'LTOE'}]';
Name_list={ExperimentalData.MarkerPositions.name}';
Marker_liste=[];

for k= 1 : numel(Name_list)
    Marker_liste=[Marker_liste;Name_list{k}];
end

% Get the numbers of solids on which the forces are applied
[~,num_m]=intersect(Marker_liste,Markers);

position_LHEE=ExperimentalData.MarkerPositions(num_m(1)).position_c3d(i,:);
position_LTOE=ExperimentalData.MarkerPositions(num_m(2)).position_c3d(i,:);

Markers=[{'LTAR'},{'LTARI'}]';
Name_list={ExperimentalData.MarkerPositions.name}';
Marker_liste=[];

for k= 1 : numel(Name_list)
    Marker_liste=[Marker_liste;Name_list{k}];
end

% Get the numbers of solids on which the forces are applied
[~,num_m]=intersect(Marker_liste,Markers);
position_LTAR=ExperimentalData.MarkerPositions(num_m(1)).position_c3d(i,:);
position_LTARI=ExperimentalData.MarkerPositions(num_m(2)).position_c3d(i,:);

Pts_L= [position_LTARI;position_LTAR;position_LHEE;position_LTOE];

L_pos = mean(Pts_L);

Z_L = position_LTARI - position_LTAR; 
Z_L = Z_L / norm( Z_L);
Y_L = position_LHEE - position_LTOE;
Y_L = Y_L / norm (Y_L);
X_L = cross (Y_L, Z_L);

Rot_L = [X_L', Y_L', Z_L'];

end

