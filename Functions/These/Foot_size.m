function [Foot_length,Foot_large] = Foot_size(filename)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

load(fullfile(filename,'ExperimentalData.mat'))

Markers=[{'RHEE'},{'RTOE'}]';
Name_list={ExperimentalData.MarkerPositions.name}';
Marker_liste=[];

for k= 1 : numel(Name_list)
    Marker_liste=[Marker_liste;Name_list{k}];
end

% Get the numbers of solids on which the forces are applied
[~,num_m]=intersect(Marker_liste,Markers);

position_HEE=ExperimentalData.MarkerPositions(num_m(1)).position_c3d;
position_TOE=ExperimentalData.MarkerPositions(num_m(2)).position_c3d;

Foot_length=mean(vecnorm(position_HEE-position_TOE,2,2));

Markers=[{'RANI'},{'RANE'}]';
Name_list={ExperimentalData.MarkerPositions.name}';
Marker_liste=[];

for k= 1 : numel(Name_list)
    Marker_liste=[Marker_liste;Name_list{k}];
end

% Get the numbers of solids on which the forces are applied
[~,num_m]=intersect(Marker_liste,Markers);
position_ANI=ExperimentalData.MarkerPositions(num_m(1)).position_c3d;
position_ANE=ExperimentalData.MarkerPositions(num_m(2)).position_c3d;

Foot_large=mean(vecnorm(position_ANI-position_ANE,2,2));

end

