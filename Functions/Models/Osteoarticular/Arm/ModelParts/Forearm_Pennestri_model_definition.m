function [Human_model]= Forearm_Pennestri_model_definition(Human_model,k,Signe,Mass,varargin)
%% Fichier de description du mod�le de bras
% Human_model : partie du mod�le d�j� construite (si il existe)
% attachment_pt : nom du points d'attache (si il existe)
% k : coefficient multiplicateur pour le scaling lin�aire
% Signe : 'R' ou 'L' (Right ou Left)
% Mass : masse du mod�le complet
% scaling_choice  : choix de la m�thode de mise � l'�chelle des donn�es inertielles
% Density : densit� du corps

%% Variables de sortie :
% "enrichissement de la structure "Human_model""

%% Liste des solides
list_solid={'Elbow_J1' 'Elbow_J2' 'Radius' 'Ulna_J1' 'Ulna_J2' 'Ulna_J3' 'Ulna' 'UlnaHumerus'};

%% Choix bras droite ou gauche
if Signe == 'R'
    Mirror=[1 0 0; 0 1 0; 0 0 1];
else
    if Signe == 'L'
        Mirror=[1 0 0; 0 1 0; 0 0 -1];
    end
end

%% Incr�mentation du num�ro des groupes
% n_group=0;
% for i=1:numel(Human_model)
%     if size(Human_model(i).Group) ~= [0 0] %#ok<BDSCA>
%         n_group=max(n_group,Human_model(i).Group(1,1));
%     end
% end
% n_group=n_group+1;


%% Incr�mentation de la num�rotation des solides

s=size(Human_model,2)+1;  %#ok<NASGU> % num�ro du premier solide
for i=1:size(list_solid,2)      % num�rotation de chaque solide : s_"nom du solide"
    if i==1
        eval(strcat('s_',list_solid{i},'=s;'))
    else
        eval(strcat('s_',list_solid{i},'=s_',list_solid{i-1},'+1;'))
    end
end

% trouver le num�ro de la m�re � partir du nom du point d'attache : 'attachment_pt'
if numel(Human_model) == 0
    s_mother=0;
    pos_attachment_pt=[0 0 0]';
else
    attachment_pt=varargin{1};
    test=0;
    for i=1:numel(Human_model)
        for j=1:size(Human_model(i).anat_position,1)
            if strcmp(attachment_pt,Human_model(i).anat_position{j,1})
                s_mother=i;
                pos_attachment_pt=Human_model(i).anat_position{j,2}+Human_model(s_mother).c;
                test=1;
                break
            end
        end
        if i==numel(Human_model) && test==0
            error([attachment_pt ' is no existent'])
        end
    end
    if Human_model(s_mother).child == 0      % si la m�re n'a pas d'enfant
        Human_model(s_mother).child = eval(['s_' list_solid{1}]);    % l'enfant de cette m�re est ce solide
    else
        [Human_model]=sister_actualize(Human_model,Human_model(s_mother).child,eval(['s_' list_solid{1}]));   % recherche de la derni�re soeur
    end
end

%%                      D�finition des noeuds

% ------------------------- Radius ----------------------------------------

% Position des noeuds
Radius_ElbowJointNode = (k*[0 0.1741 0])*Mirror;
Radius_WristJointNode = (k*[0 -0.0887 0])*Mirror;
Radius_UlnaJointNode = (k*[0 -0.0777 -0.0382])*Mirror;

% ------------------------- Ulna ------------------------------------------

% Position des noeuds
Ulna_HumerusJointNode = (k*[0 0.1088 0])*Mirror;
Ulna_RadiusJointNode = (k*[0 -0.1430 0])*Mirror;

%%              D�finition des positions anatomiques

Radius_position_set = {...
    [Signe 'RAD'], k*Mirror*[0 0.15 0.023]'; ...
    [Signe 'WRA'], k*Mirror*[0 -0.09 0.03]'; ...
    [Signe 'Forearm_WristJointNode'], Radius_WristJointNode'; ...
    % BULLSHIT
    [Signe 'Radius_SupinatorBrevis_i'],[0 0 0]';...
    [Signe 'Radius_Brachialis'],[0 0 0]';...
     [Signe 'Radius_PronatorTeres_i'],[0 0 0]';...
     [Signe 'Radius_PronatorQuadrus_i'],[0  0 0]';...
    };

Ulna_position_set = {...
    [Signe 'WRB'], k*Mirror*[0 -0.1570 0]'; ...
    %BULLSHIT
    [Signe 'Ulna_SupinatorBrevis_o'],[0 0 0]';...
    [Signe 'Ulna_AbductorDigitiV_o'],[0 0 0]';...
    [Signe 'Ulna_Biceps_i'],[0 0 0]';...
  [Signe 'Ulna_Triceps_via5'],[0 0 0]';...
  [Signe 'Ulna_Triceps_i'] ,[0 0 0]';...
  [Signe 'Ulna_PronatorQuadrus_o'],[0 0 0]';...
    };

%Fin du bullshit

%%                     Mise � l'�chelle des inerties


%% ["Adjustments to McConville et al. and Young et al. body segment inertial parameters"] R. Dumas
% Forearm
Radius_Mass=Mass*(0.498)/(0.498+0.752);
Ulna_Mass=Mass*(0.752)/(0.498+0.752);
% ------------------------- Radius ----------------------------------------
Radius_radius_sagittal = k*0.033;
Radius_radius_transverse = k*0.079;
Radius_radius_longitudinal = k*0.079;
I_Radius=[Radius_radius_sagittal*Radius_radius_sagittal*Radius_Mass, Radius_radius_longitudinal*Radius_radius_longitudinal*Radius_Mass, Radius_radius_transverse*Radius_radius_transverse*Radius_Mass, 0, 0, 0];
% ------------------------- Ulna ------------------------------------------
Radius_ulna_sagittal = k*0.02;
Radius_ulna_transverse = k*0.0745;
Radius_ulna_longitudinal = k*0.0745;
I_Ulna=[Radius_ulna_sagittal*Radius_ulna_sagittal*Ulna_Mass, Radius_ulna_longitudinal*Radius_ulna_longitudinal*Ulna_Mass, Radius_ulna_transverse*Radius_ulna_transverse*Ulna_Mass, 0, 0, 0];



%% Cr�ation de la structure "Human_model"

num_solid=0;
%% Radius
% Elbow_J1
num_solid=num_solid+1;        % solide num�ro ...
name=list_solid{num_solid}; % nom du solide
eval(['incr_solid=s_' name ';'])  % num�ro du solide dans le mod�le
Human_model(incr_solid).name=[Signe name];
Human_model(incr_solid).sister=0;
Human_model(incr_solid).child=s_Elbow_J2;
Human_model(incr_solid).mother=s_mother;
Human_model(incr_solid).a=[0 0 1]';
Human_model(incr_solid).joint=1;
Human_model(incr_solid).limit_inf=0;
Human_model(incr_solid).limit_sup=pi;
Human_model(incr_solid).ActiveJoint=1;
Human_model(incr_solid).m=0;
Human_model(incr_solid).b=pos_attachment_pt;
Human_model(incr_solid).I=zeros(3,3);
Human_model(incr_solid).c=[0 0 0]';
Human_model(incr_solid).Visual=0;

% Elbow_J2
num_solid=num_solid+1;        % solide num�ro ...
name=list_solid{num_solid}; % nom du solide
eval(['incr_solid=s_' name ';'])  % num�ro du solide dans le mod�le
Human_model(incr_solid).name=[Signe name];
Human_model(incr_solid).sister=0;
Human_model(incr_solid).child=s_Radius;
Human_model(incr_solid).mother=s_Elbow_J1;
Human_model(incr_solid).a=[1 0 0]';
Human_model(incr_solid).joint=1;
Human_model(incr_solid).limit_inf=-pi/4;
Human_model(incr_solid).limit_sup=pi/4;
Human_model(incr_solid).ActiveJoint=0;
Human_model(incr_solid).m=0;
Human_model(incr_solid).b=[0 0 0]';
Human_model(incr_solid).I=zeros(3,3);
Human_model(incr_solid).c=[0 0 0]';
Human_model(incr_solid).Visual=0;

% Radius
num_solid=num_solid+1;        % solide num�ro ...
name=list_solid{num_solid}; % nom du solide
eval(['incr_solid=s_' name ';'])  % num�ro du solide dans le mod�le
Human_model(incr_solid).name=[Signe name];
Human_model(incr_solid).sister=0;
Human_model(incr_solid).child=s_Ulna_J1;
Human_model(incr_solid).mother=s_Elbow_J2;
Human_model(incr_solid).a=[0 1 0]';
Human_model(incr_solid).joint=1;
if Signe == 'R'
    Human_model(incr_solid).limit_inf=0;
    Human_model(incr_solid).limit_sup=pi;
else
    Human_model(incr_solid).limit_inf=-pi;
    Human_model(incr_solid).limit_sup=0;
end
Human_model(incr_solid).ActiveJoint=1;
Human_model(incr_solid).m=Radius_Mass;
% Human_model(incr_solid).Group=[n_group 2];
Human_model(incr_solid).b=[0 0 0]';
Human_model(incr_solid).I=[I_Radius(1) I_Radius(4) I_Radius(5); I_Radius(4) I_Radius(2) I_Radius(6); I_Radius(5) I_Radius(6) I_Radius(3)];
Human_model(incr_solid).c=-Radius_ElbowJointNode';
Human_model(incr_solid).anat_position=Radius_position_set;
Human_model(incr_solid).Visual=1;

%% Ulna
% Ulna_J1
num_solid=num_solid+1;        % solide num�ro ...
name=list_solid{num_solid}; % nom du solide
eval(['incr_solid=s_' name ';'])  % num�ro du solide dans le mod�le
Human_model(incr_solid).name=[Signe name];
Human_model(incr_solid).sister=0;
Human_model(incr_solid).child=s_Ulna_J2;
Human_model(incr_solid).mother=s_Radius;
Human_model(incr_solid).a=[0 1 0]';
Human_model(incr_solid).joint=2;
Human_model(incr_solid).limit_inf=-0.1;
Human_model(incr_solid).limit_sup=0.1;
Human_model(incr_solid).ActiveJoint=0;
Human_model(incr_solid).m=0;
Human_model(incr_solid).b=(Radius_UlnaJointNode-Radius_ElbowJointNode)';
Human_model(incr_solid).I=zeros(3,3);
Human_model(incr_solid).c=[0 0 0]';
Human_model(incr_solid).Visual=0;

% Ulna_J2
num_solid=num_solid+1;        % solide num�ro ...
name=list_solid{num_solid}; % nom du solide
eval(['incr_solid=s_' name ';'])  % num�ro du solide dans le mod�le
Human_model(incr_solid).name=[Signe name];
Human_model(incr_solid).sister=0;
Human_model(incr_solid).child=s_Ulna_J3;
Human_model(incr_solid).mother=s_Ulna_J1;
Human_model(incr_solid).a=[0 0 1]';
Human_model(incr_solid).joint=1;
Human_model(incr_solid).limit_inf=-pi;
Human_model(incr_solid).limit_sup=pi;
Human_model(incr_solid).ActiveJoint=0;
Human_model(incr_solid).m=0;
Human_model(incr_solid).b=[0 0 0]';
Human_model(incr_solid).I=zeros(3,3);
Human_model(incr_solid).c=[0 0 0]';
Human_model(incr_solid).Visual=0;

% Ulna_J3
num_solid=num_solid+1;        % solide num�ro ...
name=list_solid{num_solid}; % nom du solide
eval(['incr_solid=s_' name ';'])  % num�ro du solide dans le mod�le
Human_model(incr_solid).name=[Signe name];
Human_model(incr_solid).sister=0;
Human_model(incr_solid).child=s_Ulna;
Human_model(incr_solid).mother=s_Ulna_J2;
Human_model(incr_solid).a=[0 1 0]';
Human_model(incr_solid).joint=1;
Human_model(incr_solid).limit_inf=-pi;
Human_model(incr_solid).limit_sup=pi;
Human_model(incr_solid).ActiveJoint=0;
Human_model(incr_solid).m=0;
Human_model(incr_solid).b=[0 0 0]';
Human_model(incr_solid).I=zeros(3,3);
Human_model(incr_solid).c=[0 0 0]';
Human_model(incr_solid).Visual=0;

% Ulna
num_solid=num_solid+1;        % solide num�ro ...
name=list_solid{num_solid}; % nom du solide
eval(['incr_solid=s_' name ';'])  % num�ro du solide dans le mod�le
Human_model(incr_solid).name=[Signe name];
Human_model(incr_solid).sister=0;
Human_model(incr_solid).child=s_UlnaHumerus;
Human_model(incr_solid).mother=s_Ulna_J3;
Human_model(incr_solid).a=[1 0 0]';
Human_model(incr_solid).joint=1;
Human_model(incr_solid).limit_inf=-pi;
Human_model(incr_solid).limit_sup=pi;
Human_model(incr_solid).ActiveJoint=0;
Human_model(incr_solid).m=Ulna_Mass;
% Human_model(incr_solid).Group=[n_group 2];
Human_model(incr_solid).b=[0 0 0]';
Human_model(incr_solid).I=[I_Ulna(1) I_Ulna(4) I_Ulna(5); I_Ulna(4) I_Ulna(2) I_Ulna(6); I_Ulna(5) I_Ulna(6) I_Ulna(3)];
Human_model(incr_solid).c=-Ulna_RadiusJointNode';
Human_model(incr_solid).calib_k_constraint=s_Radius;
Human_model(incr_solid).anat_position=Ulna_position_set;
Human_model(incr_solid).Visual=1;

% UlnaHumerus
num_solid=num_solid+1;        % solide num�ro ...
name=list_solid{num_solid}; % nom du solide
eval(['incr_solid=s_' name ';'])  % num�ro du solide dans le mod�le
Human_model(incr_solid).name=[Signe name];
Human_model(incr_solid).sister=0;
Human_model(incr_solid).child=0;
Human_model(incr_solid).mother=s_Ulna;
Human_model(incr_solid).a=[0 0 1]';
Human_model(incr_solid).joint=1;
Human_model(incr_solid).limit_inf=-pi;
Human_model(incr_solid).limit_sup=pi;
Human_model(incr_solid).ClosedLoop = [Signe 'Humerus_UlnaJointNode'];
Human_model(incr_solid).ActiveJoint=0;
Human_model(incr_solid).m=0;
Human_model(incr_solid).b=(Ulna_HumerusJointNode-Ulna_RadiusJointNode)';
Human_model(incr_solid).I=zeros(3,3);
Human_model(incr_solid).c=[0 0 0]';
Human_model(incr_solid).Visual=0;

end