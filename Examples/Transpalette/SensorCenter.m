%Fonction permettant de repérer le centre du capteur 
%filename --> non du fichier traité
%e --> épaisseur de la plateforme
function center = SensorCenter(filename,e)
%% Loading structure data
load('AnalysisParameters');
acq = btkReadAcquisition(AnalysisParameters.Prediction.ReferenceFile);
markers = btkGetMarkers(acq);
markers_list = fieldnames(markers);
[real_markers_unsorted]=C3dProcessedData(filename, markers_list);
Time = real_markers_unsorted(1).time;
real_markers_unsorted = rmfield(real_markers_unsorted,'time');
freq=1/Time(2);

%Tri des marqueurs (Marqueurs dans l'ordre croissant selon l'axe priviliégié X2 puis marqueur pour créer Y2)
real_markers = real_markers_unsorted;
no_num = 0;
for index = 1:length(real_markers)
    name = char(real_markers(index).name);
    num_loc = (name > 'A');
    if all(num_loc==0)==1
       no_num = no_num + 1; 
    end
    if no_num > 1
        error('Labeling error (only 2 markers to form perpendicular axis with names without any number)')
    end
    loc_end = find(num_loc,1,'last');
    num = str2num(name(loc_end+1:end));
    real_markers(num) = real_markers_unsorted(index);
end

% Filtrage
f_cut = AnalysisParameters.Prediction.FilterCutOff;
if AnalysisParameters.General.FilterActive
    for i=1:numel(real_markers)
    	real_markers(i).position = filt_data(real_markers(i).position_c3d,f_cut,freq);
    end
else
    for i=1:numel(real_markers)
    	real_markers(i).position = real_markers(i).position_c3d;
    end
end
real_markers = rmfield(real_markers,'position_c3d');

%Repérage du centre du capteur
center = zeros(length(Time),3);
distance_lastMarker = norm(real_markers(end).position(1,:)-real_markers(end-1).position(1,:));
distance_firstMarker = norm(real_markers(end).position(1,:)-real_markers(1).position(1,:));
for i = 1:length(Time)
    Marker1 = real_markers(end).position(i,:);
    if distance_lastMarker > distance_firstMarker
        Marker2 = real_markers(end-1).position(i,:);
    else
        Marker2 = real_markers(1).position(i,:);
    end
    Vect = Marker2-Marker1;
    center(i,:) = Marker1+0.5*Vect-[0 0 e];
end
end
