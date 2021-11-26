function [AnalysisParameters]=AnalysisForcePrediction(AnalysisParameters,trialName)
if iscell(trialName)
    trialName=trialName{1};
end

%% AnalysisParameters
% AnalysisParameters=struct();
% AnalysisParameters.General.FilterActive=1;
% AnalysisParameters.General.FilterCutOff=10;
% AnalysisParameters.General.InputData=@C3dProcessedData;
% AnalysisParameters.General.InputDataOptions={};
%
% AnalysisParameters.CalibIK.filename=meta.TrialList;
% AnalysisParameters.CalibIK.Active=1;
% AnalysisParameters.CalibIK.Frames.Method=@UniformlyDistributed;
% AnalysisParameters.CalibIK.Frames.NbFrames=100;
% AnalysisParameters.CalibIK.LengthAdd={};
% AnalysisParameters.CalibIK.LengthDelete={};
% AnalysisParameters.CalibIK.AxisAdd={};
% AnalysisParameters.CalibIK.AxisDelete={};
% % AnalysisParameters.CalibIK.MarkersCalibModif={'RTAR',{'On';'On';'On'};'LTAR',{'On';'On';'On'}};
%
AnalysisParameters.filename={trialName};
% AnalysisParameters.IK.Active=1;
% AnalysisParameters.IK.Method=2; % LM
% AnalysisParameters.IK.FilterActive=1;
% AnalysisParameters.IK.FilterCutOff=5;

%AnalysisParameters.CalibID.Active=1;
AnalysisParameters.ID.InputData = 2;
AnalysisParameters.Prediction.FilterActive=true;
% AnalysisParameters.Prediction.FilterActive=false;
AnalysisParameters.Prediction.FilterCutOff=5;
AnalysisParameters.Prediction.PositionThreshold=0.05;
AnalysisParameters.Prediction.VelocityThreshold=0.80;
AnalysisParameters.Prediction.FrictionCoef=1.07; % d'après littérature
AnalysisParameters.Prediction.FrictionCoef=0.5;% parametre de base
AnalysisParameters.Prediction.ContactPoint=...
    {'LFootPrediction1',...
    'LFootPrediction10',...
    'LFootPrediction11',...
    'LFootPrediction12',...
    'LFootPrediction13',...
    'LFootPrediction14',...
    'LFootPrediction2',...
    'LFootPrediction3',...
    'LFootPrediction4',...
    'LFootPrediction5',...
    'LFootPrediction6',...
    'LFootPrediction7',...
    'LFootPrediction8',...
    'LFootPrediction9',...
    'RFootPrediction1',...
    'RFootPrediction10',...
    'RFootPrediction11',...
    'RFootPrediction12',...
    'RFootPrediction13',...
    'RFootPrediction14',...
    'RFootPrediction2',...
    'RFootPrediction3',...
    'RFootPrediction4',...
    'RFootPrediction5',...
    'RFootPrediction6',...
    'RFootPrediction7',...
    'RFootPrediction8',...
    'RFootPrediction9'};

% 
% % On charge les marqueurs des pieds
% [real_markers, nb_frame, Firstframe, Lastframe, f, missing_markers] =...
%     C3dProcessedData(trialName(1:end-4),{'RHEE','RTAR','RTARI','RTOE','LHEE','LTAR','LTARI','LTOE'});
% 
% % On vérifie qu'il y a bien le marqueur du talon
% % à droite
% if isempty(missing_markers) || isempty(intersect('RHEE',missing_markers))
%     RightFootMarker='RHEE';
% else
%     warning('RHEE is absent')
%     RightFootMarker='RTAR';
% end
% % à gauche
% if isempty(missing_markers) || isempty(intersect('LHEE',missing_markers))
%     LeftFootMarker='LHEE';
% else
%     warning('LHEE is absent')
%     LeftFootMarker='LTAR';
% end
% 
% [~,i_RFM]=intersect([real_markers.name]',RightFootMarker);
% [~,i_LFM]=intersect([real_markers.name]',LeftFootMarker);
% 
% I_foot=[i_RFM,i_LFM];
% 
% % On calcule la vitesse des marqueurs
% for ii = 1:numel(real_markers)
%     t=real_markers(ii).time(2);
%     real_markers(ii).vitesse=derivee2(t,real_markers(ii).position_c3d);
% end
% 
% % On charge les plateformes, elles vont nous dire laquelle est la plus
% % proche du marqueur à plus faible vitesse
% h = btkReadAcquisition(trialName);
% [ForceplatesData, ForceplatesInfo] = btkGetForcePlatforms(h);
% 
% if contains(trialName,'Marche')
%     
%     % Pour chaque plateforme : resample pour troncature
%     % for each platform change of sampling frequency
%     for i=1:numel(ForceplatesData)
%         for j=1:6
%             FieldsName = fieldnames(ForceplatesData(i).channels);
%             Data(i).FullData(:,j)=resample(ForceplatesData(i).channels.(FieldsName{j}),f,ForceplatesInfo(i).frequency);
%             Data(i).RawData(:,j) = Data(i).FullData(:,j);
%         end
%     end
%     num_plate=zeros(numel(ForceplatesData),1);
%     
%     %%%%%%%%%%%%%%%%%%
%     Threshold=40; % Newton
%     %%%%%%%%%%%%%%%%%%%%%%%%
%     % On cherche les instants quand les plateformes indiquent des valeurs
%     % supérieur à 40 N toutes les deux !
%     if numel(ForceplatesData)>2
%         warning('There is more than 2 forceplates')
%         keyboard
%     end
%     bool=[];
%     for ii=1:numel(ForceplatesData)
%         bool(:,ii)=abs(Data(ii).FullData(:,3))>Threshold;
%     end
%     bool = logical(prod(bool,2));
%     F = fieldnames(real_markers);
%     rmk=struct();
%     cmpt=1;
%     for ii=I_foot
%         for jj=1:length(F)
%             if ~iscell(real_markers(ii).(F{jj}))
%                 rmk(ii).(F{jj})=real_markers(ii).(F{jj})(bool,:);
%             else
%                 rmk(ii).(F{jj})=real_markers(ii).(F{jj});
%             end
%         end
%         [num_plate(cmpt)] = DetectFootOnForcePlate(rmk(ii),ForceplatesData,2);
%         cmpt=cmpt+1;
%     end
%     
%     AnalysisParameters.ExternalForces.Options{num_plate(1),:}='RFoot';
%     AnalysisParameters.ExternalForces.Options{num_plate(2),:}='LFoot';
%     if numel(ForceplatesData)==3
%         AnalysisParameters.ExternalForces.Options{3,:}='NoContact';
%     end
%     AnalysisParameters.ExternalForces.Options
%     
%     AnalysisParameters.ExternalForces.AnalysisParameters.ExternalForces.FilterActive =true;
%     AnalysisParameters.ExternalForces.FilterCutOff=5;
%     AnalysisParameters.ExternalForces.Method=@DataInC3D;
%     
% elseif contains(trialName,'Course')
%     
%     % TROUVER LES CONTACTS DES PIEDS
%     num_plate=zeros(numel(ForceplatesData),1);
%     cmpt=1;
%     for ii=I_foot
%         [num_plate(cmpt)] = DetectFootOnForcePlate(real_markers(ii),ForceplatesData,0.5);
%         cmpt=cmpt+1;
%     end
%    
%     AnalysisParameters.ExternalForces.Options{num_plate(1),:}='RFoot';
%     AnalysisParameters.ExternalForces.Options{num_plate(2),:}='LFoot';
%     if numel(ForceplatesData)==3
%         AnalysisParameters.ExternalForces.Options{3,:}='NoContact';
%     end
%     AnalysisParameters.ExternalForces.Options
%     
%     AnalysisParameters.ExternalForces.AnalysisParameters.ExternalForces.FilterActive =true;
%     AnalysisParameters.ExternalForces.FilterCutOff=5;
%     AnalysisParameters.ExternalForces.Method=@DataInC3D;
%     
%     
% elseif contains(trialName,'Direction') % Changement de direction
%     
%     % TROUVER LES CONTACTS DES PIEDS
%     num_plate=zeros(numel(ForceplatesData),1);
%     cmpt=1;
%     for ii=I_foot
%         [num_plate(cmpt)] = DetectFootOnForcePlate(real_markers(ii),ForceplatesData,0.5);
%         cmpt=cmpt+1;
%     end
%    
%     AnalysisParameters.ExternalForces.Options{num_plate(1),:}='RFoot';
%     AnalysisParameters.ExternalForces.Options{num_plate(2),:}='LFoot';
%     if numel(ForceplatesData)==3
%         AnalysisParameters.ExternalForces.Options{3,:}='NoContact';
%     end
%     AnalysisParameters.ExternalForces.Options
%     
%     AnalysisParameters.ExternalForces.AnalysisParameters.ExternalForces.FilterActive =true;
%     AnalysisParameters.ExternalForces.FilterCutOff=5;
%     AnalysisParameters.ExternalForces.Method=@DataInC3D;
%     
% elseif contains(trialName,'Saut')
%     
%     % TROUVER LES CONTACTS DES PIEDS
%     num_plate=zeros(numel(ForceplatesData),1);
%     cmpt=1;
%     for ii=I_foot
%         [num_plate(cmpt)] = DetectFootOnForcePlate(real_markers(ii),ForceplatesData,0.5);
%         cmpt=cmpt+1;
%     end
%    
%     AnalysisParameters.ExternalForces.Options{num_plate(1),:}='RFoot';
%     AnalysisParameters.ExternalForces.Options{num_plate(2),:}='LFoot';
%     if numel(ForceplatesData)==3
%         AnalysisParameters.ExternalForces.Options{3,:}='NoContact';
%     end
%     AnalysisParameters.ExternalForces.Options
%     
%     AnalysisParameters.ExternalForces.AnalysisParameters.ExternalForces.FilterActive =true;
%     AnalysisParameters.ExternalForces.FilterCutOff=5;
%     AnalysisParameters.ExternalForces.Method=@DataInC3D;
%     
% elseif contains(trialName,'ROM')
%     
%     % TROUVER LES CONTACTS DES PIEDS
%     num_plate=zeros(numel(ForceplatesData),1);
%     cmpt=1;
%     for ii=I_foot
%         [num_plate(cmpt)] = DetectFootOnForcePlate(real_markers(ii),ForceplatesData,0.5);
%         cmpt=cmpt+1;
%     end
%    
%     AnalysisParameters.ExternalForces.Options{num_plate(1),:}='RFoot';
%     AnalysisParameters.ExternalForces.Options{num_plate(2),:}='LFoot';
%     if numel(ForceplatesData)==3
%         AnalysisParameters.ExternalForces.Options{3,:}='NoContact';
%     end
%     AnalysisParameters.ExternalForces.Options
%     
%     AnalysisParameters.ExternalForces.AnalysisParameters.ExternalForces.FilterActive =true;
%     AnalysisParameters.ExternalForces.FilterCutOff=5;
%     AnalysisParameters.ExternalForces.Method=@DataInC3D;
%     
% elseif contains(trialName,'Statref') || contains(trialName,'statref')
%     
%     % TROUVER LES CONTACTS DES PIEDS
%     num_plate=zeros(numel(ForceplatesData),1);
%     cmpt=1;
%     for ii=I_foot
%         [num_plate(cmpt)] = DetectFootOnForcePlate(real_markers(ii),ForceplatesData,0.5);
%         cmpt=cmpt+1;
%     end
%    
%     AnalysisParameters.ExternalForces.Options{num_plate(1),:}='RFoot';
%     AnalysisParameters.ExternalForces.Options{num_plate(2),:}='LFoot';
%     if numel(ForceplatesData)==3
%         AnalysisParameters.ExternalForces.Options{3,:}='NoContact';
%     end
%     AnalysisParameters.ExternalForces.Options
%     
%     AnalysisParameters.ExternalForces.AnalysisParameters.ExternalForces.FilterActive =true;
%     AnalysisParameters.ExternalForces.FilterCutOff=5;
%     AnalysisParameters.ExternalForces.Method=@DataInC3D;
%     
% else
%     warning(['Type of task not detected for ' trialName])
%     keyboard
% end

end

% save('AnalysisParameters.mat','AnalysisParameters')
% save('ModelParameters','ModelParameters')