function [Muscles]=ArmMusclesPennestri(Muscles,Signe)
% Definition of an arm muscle model
%   This model contains 21 muscles
%
%	Based on:
%	- Pennestrì, E. , Stefanelli, R. , Valentini, P. P. , Vita, L.
%	Virtual musculo-skeletal model for the biomechanical
% analysis of the upper limb, Pennestrì2007

%   INPUT
%   - Muscles: set of muscles (see the Documentation for the structure)
%   - Signe: Signe of the arm model ('R' for right Signe or 'L' for left Signe)
%   OUTPUT
%   - Muscles: new set of muscles (see the Documentation for the structure)
%________________________________________________________
%
% Licence
% Toolbox distributed under GPL 3.0 Licence
%________________________________________________________
%
% Authors : Antoine Muller, Charles Pontonnier, Pierre Puchaud and
% Georges Dumont
%________________________________________________________


if strcmp(Signe,'Right')
    Signe = 'R';
else
    Signe = 'L';
end

s=cell(0);

s=[s;{
%    [Signe 'Trapezius'],240,0.1,[],[],[],{[Signe 'Humerus_Trapezius_o'];[Signe 'Thorax_Trapezius_i']},{};...   
   [Signe 'Anconeus'],350,0.027,[],[],[],{[Signe 'Humerus_AnconeusBF_o'];[Signe 'Humerus_AnconeusBF_VP1'];[Signe 'Ulna_Anconeus_VP2'];[Signe 'Ulna_Anconeus_i']},{};...
   [Signe 'FlexorCarpiUlnaris'],128.9,0.051,[],[],[],{[Signe 'Humerus_FlexorCarpiUlnarisBF_o'];[Signe 'Humerus_FlexorCarpiUlnarisBF_VP1'];[Signe 'Radius_FlexorCarpiUlnaris_VP2'];[Signe 'Radius_FlexorCarpiUlnaris_VP3'];[Signe 'Hand_FlexorCarpiUlnarisBF_VP4'];[Signe 'Hand_FlexorCarpiUlnarisBF_i']},{};...
   [Signe 'ExtensorCarpiUlnaris'],93.2,0.062,[],[],[],{[Signe 'Humerus_ExtensorCarpiUlnarisBF_o'];[Signe 'Humerus_ExtensorCarpiUlnarisBF_VP1'];[Signe 'Radius_ExtensorCarpiUlnaris_VP2'];[Signe 'Radius_ExtensorCarpiUlnaris_VP3'];[Signe 'Hand_ExtensorCarpiUlnarisBF_VP4'];[Signe 'Hand_ExtensorCarpiUlnarisBF_i']},{};...
   [Signe 'ExtensorCarpiRadialisLongus'],304.9,0.081,[],0.244,[],{[Signe 'Humerus_ExtensorCarpiRadialisLongusBF_o'];[Signe 'Humerus_ExtensorCarpiRadialisLongusBF_VP1'];[Signe 'Radius_ExtensorCarpiRadialisLongus_VP2'];[Signe 'Radius_ExtensorCarpiRadialisLongus_VP3'];[Signe 'Hand_ExtensorCarpiRadialisLongusBF_VP4'];[Signe 'Hand_ExtensorCarpiRadialisLongusBF_i']},{};...
   [Signe 'ExtensorCarpiRadialisBrevis'],100.5,0.059,[],0.2223,[],{[Signe 'Humerus_ExtensorCarpiRadialisBrevisBF_o'];[Signe 'Humerus_ExtensorCarpiRadialisBrevisBF_VP1'];[Signe 'Radius_ExtensorCarpiRadialisBrevis_VP2'];[Signe 'Radius_ExtensorCarpiRadialisBrevis_VP3'];[Signe 'Hand_ExtensorCarpiRadialisBrevisBF_VP4'];[Signe 'Hand_ExtensorCarpiRadialisBrevisBF_i']},{};...
   [Signe 'FlexorCarpiRadialis'],74,0.063,[],[],[],{[Signe 'Humerus_FlexorCarpiRadialisBF_o'];[Signe 'Humerus_FlexorCarpiRadialisBF_VP1'];[Signe 'Radius_FlexorCarpiRadialis_VP2'];[Signe 'Radius_FlexorCarpiRadialis_VP3'];[Signe 'Hand_FlexorCarpiRadialisBF_VP4'];[Signe 'Hand_FlexorCarpiRadialisBF_i']},{};...
   [Signe 'PronatorQuadratus'],75.5,0.028,[],[],[],{[Signe 'Radius_PronatorQuadratus_o'];[Signe 'Radius_PronatorQuadratus_VP1'];[Signe 'Ulna_PronatorQuadratus_VP2'];[Signe 'Ulna_PronatorQuadratus_i']},{};...
   [Signe 'SupinatorBrevis'],476,0.033,[],[],[],{[Signe 'Radius_SupinatorBrevis_o'];[Signe 'Radius_SupinatorBrevis_VP1'];[Signe 'Ulna_SupinatorBrevis_VP2'];[Signe 'Ulna_SupinatorBrevis_i']},{};...
   [Signe 'Brachialis'],987.3,0.086,4,0.0535,0,{[Signe 'Humerus_BrachialisBF_o'];[Signe 'Humerus_BrachialisBF_VP1'];[Signe 'Ulna_Brachialis_VP2'];[Signe 'Ulna_Brachialis_i']},{};...
   [Signe 'Brachioradialis'],261.3,0.173,[],[],[],{[Signe 'Humerus_BrachioradialisBF_o'];[Signe 'Humerus_BrachioradialisBF_VP1'];[Signe 'Radius_Brachioradialis_VP2'];[Signe 'Radius_Brachioradialis_i']},{};...
   [Signe 'PronatorTeres'],566.2,0.049,[],[],[],{[Signe 'Humerus_PronatorTeresBF_o'];[Signe 'Humerus_PronatorTeresBF_VP1'];[Signe 'Radius_PronatorTeres_VP2'];[Signe 'Radius_PronatorTeres_i']},{};...
   [Signe 'TricepsLat'],624.3,0.114,4,0.098,0.157,{[Signe 'Humerus_TricepsLatBF_o'];[Signe 'Humerus_TricepsLatBF_VP1'];[Signe 'Ulna_TricepsLat_VP2'];[Signe 'Ulna_TricepsLat_i']},{};... 
   [Signe 'TricepsMed'],624.3,0.114,4,0.098,0.157,{[Signe 'Humerus_TricepsMedBF_o'];[Signe 'Humerus_TricepsMedBF_VP1'];[Signe 'Ulna_TricepsMed_VP2'];[Signe 'Ulna_TricepsMed_i']},{};... arm26.osim
   [Signe 'PalmarisLongus'],26.7,0.064,4,0.098,0.157,{[Signe 'Humerus_PalmarisLongusBF_o'];[Signe 'Humerus_PalmarisLongusBF_VP1'];[Signe 'Radius_PalmarisLongus_VP2'];[Signe 'Radius_PalmarisLongus_VP3'];[Signe 'Hand_PalmarisLongusBF_VP4'];[Signe 'Hand_PalmarisLongusBF_i']},{};...

    
    % Fake muscles from (Seth et al., 2016)
    [Signe '_Deltoid_ant'],1218.9,0.0976,[],0.093,0.38397,{[Signe '_humerus_r_DELT1_r-P1'];[Signe '_humerus_r_DELT1_r-P1'];[Signe '_Scapula_DELT1-P3'];[Signe '_clavicle_r_DELT1_r-P4']},{['Wrap' Signe 'HumerusDelt']};... 
    [Signe '_Deltoid_mid'],1103.5,0.1078,[],0.1095,0.2618,{[Signe '_humerus_r_DELT2_r-P1'];[Signe '_Scapula_DELT2-P3'];[Signe '_Scapula_DELT2-P4']},{['Wrap' Signe 'ThoraxGH']};...   
    %     % Conservation du modèle de Holzbaur
% % on conserve les biceps du modèle de Holzbaur sauf qu'on part de la
%     % scapula pour le biceps short et glénoïde pour le biceps long
%     [Signe 'BicepsL'],624.3,0.1157,4,0.2723,0,{[Signe 'Scapula_BicepsL_o'];[Signe 'Scapula_BicepsL_via1'];[Signe 'Humerus_BicepsL_via2'];[Signe 'Humerus_BicepsL_via3'];[Signe 'Humerus_BicepsL_via4'];[Signe 'Humerus_BicepsL_via5'];[Signe 'Humerus_BicepsL_via6'];[Signe 'Humerus_Biceps_via7'];[Signe 'Ulna_Biceps_i']},{};... arm26.osim       
%     [Signe 'BicepsS'],435.56,0.1321,4,0.1923,0,{[Signe 'Scapula_BicepsS_o'];[Signe 'Scapula_BicepsS_via1'];[Signe 'Humerus_BicepsS_via2'];[Signe 'Humerus_BicepsS_via3'];[Signe 'Humerus_Biceps_via7'];[Signe 'Ulna_Biceps_i']},{};... arm26.osim    
%     % on conserve les biceps du modèle de Holzbaur sauf qu'on part de la
%     % scapula pour le triceps long
%     [Signe 'TricepsLg'],798.5,0.134,4,0.143,0.209,{[Signe 'Scapula_Triceps_o'];[Signe 'Humerus_TricepsLg_via1'];[Signe 'Humerus_Triceps_via2'];[Signe 'Humerus_Triceps_via3'];[Signe 'Humerus_Triceps_via4'];[Signe 'Ulna_TricepsMed_VP2'];[Signe 'Ulna_TricepsMed_i']},{};...       arm26.osim    
%      [Signe 'TricepsLat'],624.3,0.114,4,0.098,0.157,{[Signe 'Humerus_TricepsLat_o'];[Signe 'Humerus_TricepsLat_via1'];[Signe 'Humerus_Triceps_via2'];[Signe 'Humerus_Triceps_via3'];[Signe 'Humerus_Triceps_via4'];[Signe 'Ulna_Triceps_via5'];[Signe 'Ulna_Triceps_i']},{};... arm26.osim   
%      [Signe 'TricepsMed'],624.3,0.114,4,0.098,0.157,{[Signe 'Humerus_TricepsMed_o'];[Signe 'Humerus_TricepsMed_via1'];[Signe 'Humerus_Triceps_via2'];[Signe 'Humerus_Triceps_via3'];[Signe 'Humerus_Triceps_via4'];[Signe 'Ulna_Triceps_via5'];[Signe 'Ulna_Triceps_i']},{};... arm26.osim
  %         [Signe 'TricepsLg'],798.5,0.134,4,0.143,0.209,{[Signe 'Scapula_Triceps_o'];[Signe 'Humerus_TricepsLg_via1'];[Signe 'Humerus_Triceps_via2'];[Signe 'Humerus_Triceps_via3'];[Signe 'Humerus_Triceps_via4'];[Signe 'Ulna_Triceps_i']},{};...       arm26.osim    
%      [Signe 'TricepsLat'],624.3,0.114,4,0.098,0.157,{[Signe 'Humerus_TricepsLat_o'];[Signe 'Humerus_TricepsLat_via1'];[Signe 'Humerus_Triceps_via2'];[Signe 'Humerus_Triceps_via3'];[Signe 'Humerus_Triceps_via4'];[Signe 'Ulna_Triceps_i']},{};... arm26.osim   
%      [Signe 'TricepsMed'],624.3,0.114,4,0.098,0.157,{[Signe 'Humerus_TricepsMed_o'];[Signe 'Humerus_TricepsMed_via1'];[Signe 'Humerus_Triceps_via2'];[Signe 'Humerus_Triceps_via3'];[Signe 'Humerus_Triceps_via4'];[Signe 'Ulna_Triceps_i']},{};... arm26.osim
        }];


% Structure generation
Muscles=[Muscles;struct('name',{s{:,1}}','f0',{s{:,2}}','l0',{s{:,3}}',...
    'Kt',{s{:,4}}','ls',{s{:,5}}','alpha0',{s{:,6}}','path',{s{:,7}}','wrap',{s{:,8}}')]; %#ok<CCAT1>
end