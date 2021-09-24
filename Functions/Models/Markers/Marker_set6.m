function [Markers]=Marker_set6(varargin)
% Definition of the markers set used in XSENS
%
%   OUTPUT
%   - Markers: set of markers (see the Documentation for the structure) 
%________________________________________________________
%
% Licence
% Toolbox distributed under GPL 3.0 Licence
%________________________________________________________
%
% Authors : Antoine Muller, Charles Pontonnier, Pierre Puchaud and
% Georges Dumont
%________________________________________________________

s=cell(0);

% TRUNK & HEAD
s=[s;{'RASIS' 'RFWT' {'On';'Off';'On'};...
            'LASIS' 'LFWT' {'On';'Off';'On'};...
            'RPSIS' 'RBWT' {'On';'Off';'On'};...
            'LPSIS' 'LBWT' {'On';'Off';'On'};...
            'Lumbar' 'Pelvis_L5JointNode' {'Off';'On';'Off'}; ...
            'Thoracic' 'T10' {'Off';'On';'Off'};...
            'RightShoulder' 'RSHO' {'Off';'Off';'On'};...
            'LeftShoulder' 'LSHO' {'Off';'Off';'On'}; ...
            'Neck' 'C7' {'On';'On';'Off'}; ...
            'Sternum' 'STRN' {'Off';'On';'Off'};...
             'Xiphoid' 'CLAV'  {'Off';'Off';'Off'};  }];

Side1={'Right';'Left'};
Side2={'R';'L'};
% Arm
for i=1:2
    s=[s;{  ['p' Side2{i} 'Elbow'] [Side2{i} 'RAD'] {'On';'On';'Off'};...
            ['p' Side2{i} 'Ulna'] [Side2{i} 'WRB'] {'Off';'On';'Off'}; ...
            ['p' Side2{i} 'Radius'] [Side2{i} 'WRA'] {'Off';'Off';'Off'}; ...
            ['p' Side2{i} 'Wrist'] [Side2{i} 'Hand_WristJointNode'] {'Off';'Off';'Off'}; ...
        }];
end

% Leg
for i=1:2
    s=[s;{  ['p' Side1{i} 'ShankSuperior'] [Side2{i} 'KNE'] {'On';'On';'On'};...
            ['p' Side2{i} 'Patella'] [Side2{i} 'KNI'] {'Off';'On';'On'}; ...
            ['p' Side2{i} 'Hindfoot'] [Side2{i} 'ANE'] {'Off';'On';'Off'};...
            ['p' Side2{i} 'MidfootLateral'] [Side2{i} 'TAR'] {'On';'On';'On'}; ...
            ['p' Side2{i} 'MidfootMedial'] [Side2{i} 'TARI'] {'On';'On';'On'};...
            ['p' Side2{i} 'Heel'] [Side2{i} 'HEE'] {'Off';'On';'Off'}; ...
            ['p' Side2{i} 'Toe'] [Side2{i} 'TOE'] {'Off';'Off';'Off'}; ...
        }];
end

Markers=struct('name',{s{:,1}}','anat_position',{s{:,2}}','calib_dir',{s{:,3}}'); %#ok<CCAT1>

end