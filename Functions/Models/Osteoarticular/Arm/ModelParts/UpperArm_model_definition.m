function [Human_model]= UpperArm_model_definition(Human_model,k,Signe,Mass,varargin)
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
list_solid={'Scapula_J1' 'Scapula_J2' 'ScapulaThorax' 'Humerus' };

%% Choix bras droite ou gauche
if Signe == 'R'
    Mirror=[1 0 0; 0 1 0; 0 0 1];
else
    if Signe == 'L'
        Mirror=[1 0 0; 0 1 0; 0 0 -1];
    end
end

% %% Incr�mentation du num�ro des groupes
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

% ------------------------- Humerus ---------------------------------------

% Position des noeuds
Humerus_ghJointNode = (k*[0 0.1674 0])*Mirror;
Humerus_ElbowJointNode = (k*[0 -0.1674 0])*Mirror;
osim2antoine=[k (Humerus_ghJointNode(2)-Humerus_ElbowJointNode(2))/0.2904 k];
Humerus_RadiusJointNode = (k*[0 -0.1674 0.0191])*Mirror;
Humerus_UlnaJointNode = (k*[0 -0.1674 -0.0191])*Mirror;
Humerus_Brachioradialis = (k*[-0.006 -0.209 -0.007])*Mirror; %dans le repere gh Murray2001
% Humerus_Biceps = (k*[0.025 0.009 0.006])*Mirror; %dans le repere gh Murray2001
Humerus_BicepsL_via2 = (osim2antoine.*[0.02131 0.01793 0.01028])*Mirror;  %dans le repere OSIMarm26
Humerus_BicepsL_via3 = (osim2antoine.*[0.02378 -0.00511 0.01201])*Mirror;  %dans le repere OSIMarm26
Humerus_BicepsL_via4 = (osim2antoine.*[0.01345 -0.02827 0.00136])*Mirror;  %dans le repere OSIMarm26
Humerus_BicepsL_via5 = (osim2antoine.*[0.01068 -0.07736 -0.00165])*Mirror;  %dans le repere OSIMarm26
Humerus_BicepsL_via6 = (osim2antoine.*[0.01703 -0.12125 0.00024])*Mirror;  %dans le repere OSIMarm26
Humerus_BicepsS_via2 = (osim2antoine.*[0.01117 -0.07576 -0.01101])*Mirror;  %dans le repere OSIMarm26
Humerus_BicepsS_via3 = (osim2antoine.*[0.01703 -0.12125 -0.01079])*Mirror;  %dans le repere OSIMarm26
Humerus_Biceps_via7 = (osim2antoine.*[0.0228 -0.1754 -0.0063])*Mirror;  %dans le repere OSIMarm26

Humerus_ECRL = (k*[-0.005 -0.260 -0.002])*Mirror; %dans le repere gh Murray2001
% Humerus_Brachialis = (k*[0.008 -0.184 -0.013])*Mirror; %dans le repere gh Murray2001
Humerus_Brachialis = (k*[0.0068 -0.1739 -0.0036])*Mirror; %dans le repere OSIMarm26
Humerus_PronatorTeres = (k*[0.003 -0.270 -0.051])*Mirror; %dans le repere gh Murray2001

% Humerus_Triceps = (k*[-0.004 -0.039 -0.006])*Mirror; %dans le repere gh Murray2001
Humerus_TricepsLg_via1 = (osim2antoine.*[-0.02714 -0.11441 -0.00664])*Mirror;  %dans le repere OSIMarm26
Humerus_TricepsLat_o = (osim2antoine.*[-0.00599 -0.12646 0.00428])*Mirror;     %dans le repere OSIMarm26
Humerus_TricepsLat_via1 = (osim2antoine.*[-0.02344 -0.14528 0.00928])*Mirror;  %dans le repere OSIMarm26
Humerus_TricepsMed_o = (osim2antoine.*[-0.00838 -0.13695 -0.00906])*Mirror;    %dans le repere OSIMarm26
Humerus_TricepsMed_via1 = (osim2antoine.*[-0.02601 -0.15139 -0.0108])*Mirror;  %dans le repere OSIMarm26
Humerus_Triceps_via2 = (osim2antoine.*[-0.03184 -0.22637 -0.01217])*Mirror; %dans le repere OSIMarm26
Humerus_Triceps_via3 = (osim2antoine.*[-0.01743 -0.26757 -0.01208])*Mirror; %dans le repere OSIMarm26

%%              D�finition des positions anatomiques

Humerus_position_set = {...
    [Signe 'HUM'], k*Mirror*[0 -0.1674 -0.05]'; ...
    [Signe 'Humerus_RadiusJointNode'], Humerus_RadiusJointNode'; ...
    [Signe 'Humerus_UlnaJointNode'], Humerus_UlnaJointNode'; ...
    [Signe 'Humerus_ElbowJointNode'], Humerus_ElbowJointNode'; ...
    [Signe 'Humerus_ghJointNode'], Humerus_ghJointNode'; ...
    %     [Signe 'Humerus_Brachioradialis_o'], (Humerus_Brachioradialis+Humerus_ghJointNode)'; ...
    [Signe 'Humerus_Brachioradialis_o'], Humerus_RadiusJointNode'+[0 0.07 0]'; ...
    ...[Signe 'Humerus_Biceps'], (Humerus_Biceps+Humerus_ghJointNode)'; ...
    [Signe 'Humerus_BicepsL_via2'], (Humerus_BicepsL_via2+Humerus_ghJointNode)';
    [Signe 'Humerus_BicepsL_via3'], (Humerus_BicepsL_via3+Humerus_ghJointNode)';
    [Signe 'Humerus_BicepsL_via4'], (Humerus_BicepsL_via4+Humerus_ghJointNode)';
    [Signe 'Humerus_BicepsL_via5'], (Humerus_BicepsL_via5+Humerus_ghJointNode)';
    [Signe 'Humerus_BicepsL_via6'], (Humerus_BicepsL_via6+Humerus_ghJointNode)';
    [Signe 'Humerus_BicepsS_via2'], (Humerus_BicepsS_via2+Humerus_ghJointNode)';
    [Signe 'Humerus_BicepsS_via3'], (Humerus_BicepsS_via3+Humerus_ghJointNode)';
    [Signe 'Humerus_Biceps_via7'], (Humerus_Biceps_via7+Humerus_ghJointNode)';
    ...
    %     [Signe 'Humerus_ECRL_o'], (Humerus_ECRL+Humerus_ghJointNode)'; ...
    [Signe 'Humerus_ECRL_o'], Humerus_RadiusJointNode'+[0 0.03 0]'; ...
    [Signe 'Humerus_Brachialis_o'], (Humerus_Brachialis+Humerus_ghJointNode)'; ...
    %     [Signe 'Humerus_PronatorTeres_o'], (Humerus_PronatorTeres+Humerus_ghJointNode)'; ...
    [Signe 'Humerus_PronatorTeres_o'], Humerus_UlnaJointNode'+[0 0.02 0.01]'; ...
    ...
    ...[Signe 'Humerus_Triceps_o'], (Humerus_Triceps+Humerus_ghJointNode)'; ...
    [Signe 'Humerus_TricepsLg_via1'], (Humerus_TricepsLg_via1+Humerus_ghJointNode)';
    [Signe 'Humerus_TricepsLat_o'], (Humerus_TricepsLat_o+Humerus_ghJointNode)';
    [Signe 'Humerus_TricepsLat_via1'], (Humerus_TricepsLat_via1+Humerus_ghJointNode)';
    [Signe 'Humerus_TricepsMed_o'], (Humerus_TricepsMed_o+Humerus_ghJointNode)';
    [Signe 'Humerus_TricepsMed_via1'], (Humerus_TricepsMed_via1+Humerus_ghJointNode)';
    [Signe 'Humerus_Triceps_via2'], (Humerus_Triceps_via2+Humerus_ghJointNode)';
    [Signe 'Humerus_Triceps_via3'], (Humerus_Triceps_via3+Humerus_ghJointNode)';
    [Signe 'Humerus_Triceps_via4'], Humerus_ElbowJointNode' + k*[-0.028 0 0]';...
    % Partie placement bullshit A REFAIRE PROPREMENT
     [Signe 'Humerus_Coracobrachialis_i'],[ 0 0 0]';...
     [Signe 'Humerus_Deltoid_i'],[ 0 0 0]';...
    [Signe 'Humerus_LatissumusDorsi_i'],[ 0 0 0]';...
    [Signe 'Humerus_PectoralisMajor_i'],[ 0 0 0]';...
    [Signe 'Humerus_Supraspinatus_i'],[ 0 0 0]';...
    [Signe 'Humerus_Infraspinatus_i'],[ 0 0 0]';...
    [Signe 'Humerus_CubitalisAnterior_o'],[ 0 0 0]';...
    [Signe 'Humerus_FlexorCarpiUlnaris_o'],[ 0 0 0]';...
    [Signe 'Humerus_ExtensorCarpiUlnaris_o'],[ 0 0 0]';...
    [Signe 'Humerus_ExtensorDigitorum_o'],[ 0 0 0]';...
    [Signe 'Humerus_FlexorDigitorumSuperior_o'],[ 0 0 0]';...
    [Signe 'Humerus_FlexorCarpiRadialis_o'],[ 0 0 0]';...
    };



    % Partie placement bullshit A REFAIRE PROPREMENT
Scapula_position_set={...
    [Signe 'Scapula_Coracobrachialis_o'],[ 0 0 0]';...
   [Signe 'Scapula_Deltoid_o'],[ 0 0 0]';...
  [Signe 'Scapula_Supraspinatus_o'],[ 0 0 0]';...
   [Signe 'Scapula_Infraspinatus_o'],[ 0 0 0]';...    
   [Signe 'Scapula_BicepsL_o'],[ 0 0 0]';...
   [Signe 'Scapula_BicepsL_via1'],[ 0 0 0]';...
    [Signe 'Scapula_BicepsS_o'],[ 0 0 0]';...
    [Signe 'Scapula_BicepsS_via1'],[ 0 0 0]';...
    [Signe 'Scapula_Triceps_o'],[ 0 0 0]';...
    };

%Fin du bullshit

%%                     Mise � l'�chelle des inerties


%% ["Adjustments to McConville et al. and Young et al. body segment inertial parameters"] R. Dumas
% ------------------------- Humerus ---------------------------------------
Length_Humerus=norm(Humerus_ghJointNode-Humerus_ElbowJointNode);
[I_Humerus]=rgyration2inertia([31 14 32 6 5 2], Mass, [0 0 0], Length_Humerus, Signe);


%% Cr�ation de la structure "Human_model"

num_solid=0;
%% Humerus
% Shoulder_J1
num_solid=num_solid+1;        % solide num�ro ...
name=list_solid{num_solid}; % nom du solide
eval(['incr_solid=s_' name ';'])  % num�ro du solide dans le mod�le
Human_model(incr_solid).name=[Signe name];    % nom du solide ('R' ou 'L' en suffixe)
Human_model(incr_solid).sister=0;                       % sister : d�finit en entr�e de la fonction
Human_model(incr_solid).child=s_Scapula_J2;                 % child : Shank
Human_model(incr_solid).mother=s_mother;                       % mother : d�finit en entr�e de la fonction
Human_model(incr_solid).a=[0 1 0]';                          % rotation /x
Human_model(incr_solid).joint=1;
if Signe == 'R'
    Human_model(incr_solid).limit_inf=-pi/2;                       % but�e articulaire inf�rieure
    Human_model(incr_solid).limit_sup=pi;                      % but�e articulaire sup�rieure
else
    Human_model(incr_solid).limit_inf=-pi;                       % but�e articulaire inf�rieure
    Human_model(incr_solid).limit_sup=pi/2;                      % but�e articulaire sup�rieure
end
Human_model(incr_solid).ActiveJoint=1;
Human_model(incr_solid).m=0;                                 % masse de r�f�rence
Human_model(incr_solid).b=pos_attachment_pt;                    % position du point d'attache par rapport au rep�re parent
Human_model(incr_solid).I=zeros(3,3);                        % matrice d'inertie de r�f�rence
Human_model(incr_solid).c=[0 0 0]';                          % position du centre de masse dans le rep�re local
Human_model(incr_solid).calib_k_constraint=[];
Human_model(incr_solid).u=[];                          % rotation fixe selon l'axe u d'un angle theta (apr�s la rotation q)
Human_model(incr_solid).theta=[];
Human_model(incr_solid).KinematicsCut=[];              % coupure cin�matique
Human_model(incr_solid).ClosedLoop=[];
                 % si solide de fermeture de boucle : {num�ro du solide i sur lequel est attach� ce solide ; point d'attache (rep�re du solide i)}
Human_model(incr_solid).linear_constraint=[];
Human_model(incr_solid).Visual=0;

% Shoulder_J2
num_solid=num_solid+1;        % solide num�ro ...
name=list_solid{num_solid}; % nom du solide
eval(['incr_solid=s_' name ';'])  % num�ro du solide dans le mod�le
Human_model(incr_solid).name=[Signe name];
Human_model(incr_solid).sister=0;
Human_model(incr_solid).child=s_Humerus;
Human_model(incr_solid).mother=s_Scapula_J1;
Human_model(incr_solid).a=[1 0 0]';
Human_model(incr_solid).joint=1;
if Signe == 'R'
    Human_model(incr_solid).limit_inf=-pi;                       % but�e articulaire inf�rieure
    Human_model(incr_solid).limit_sup=pi/2;                      % but�e articulaire sup�rieure
else
    Human_model(incr_solid).limit_inf=-pi/2;                       % but�e articulaire inf�rieure
    Human_model(incr_solid).limit_sup=pi;                      % but�e articulaire sup�rieure
end
Human_model(incr_solid).ActiveJoint=1;
Human_model(incr_solid).m=0;
Human_model(incr_solid).b=[0 0 0]';
Human_model(incr_solid).I=zeros(3,3);
Human_model(incr_solid).c=[0 0 0]';
Human_model(incr_solid).Visual=0;
Human_model(incr_solid).anat_position=Scapula_position_set;



% ScapulaClavicule
num_solid=num_solid+1;        % solide num�ro ...
name=list_solid{num_solid}; % nom du solide
eval(['incr_solid=s_' name ';'])  % num�ro du solide dans le mod�le
Human_model(incr_solid).name=[Signe name];
Human_model(incr_solid).sister=0;
Human_model(incr_solid).child=0;
Human_model(incr_solid).mother=s_Scapula_J1;
Human_model(incr_solid).a=[0 0 1]';
Human_model(incr_solid).joint=1;
Human_model(incr_solid).limit_inf=-pi;
Human_model(incr_solid).limit_sup=pi;
Human_model(incr_solid).ClosedLoop = ['Thorax_Scapula_J1JointNode'];
Human_model(incr_solid).ActiveJoint=0;
Human_model(incr_solid).m=0;
Human_model(incr_solid).b=[0 0 0]';
Human_model(incr_solid).I=zeros(3,3);
Human_model(incr_solid).c=[0 0 0]';
Human_model(incr_solid).Visual=0;


% Humerus
num_solid=num_solid+1;        % solide num�ro ...
name=list_solid{num_solid}; % nom du solide
eval(['incr_solid=s_' name ';'])  % num�ro du solide dans le mod�le
Human_model(incr_solid).name=[Signe name];
Human_model(incr_solid).sister=0;
Human_model(incr_solid).child=0;
Human_model(incr_solid).mother=s_Scapula_J2;
Human_model(incr_solid).a=[0 1 0]';
Human_model(incr_solid).joint=1;
Human_model(incr_solid).limit_inf=-2*pi/3;
Human_model(incr_solid).limit_sup=2*pi/3;
Human_model(incr_solid).ActiveJoint=1;
% Human_model(incr_solid).Group=[n_group 1];
Human_model(incr_solid).m=Mass;
Human_model(incr_solid).b=[0 0 0]';
Human_model(incr_solid).I=[I_Humerus(1) I_Humerus(4) I_Humerus(5); I_Humerus(4) I_Humerus(2) I_Humerus(6); I_Humerus(5) I_Humerus(6) I_Humerus(3)];
Human_model(incr_solid).c=-Humerus_ghJointNode';
Human_model(incr_solid).anat_position=Humerus_position_set;
Human_model(incr_solid).Visual=1;
Human_model(incr_solid).L={[Signe 'Humerus_ghJointNode'];[Signe 'Humerus_ElbowJointNode']};





end