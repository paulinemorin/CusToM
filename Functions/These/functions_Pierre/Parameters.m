function [ModelParameters,AnalysisParameters]=Parameters(meta)

%% ModelParameters
ModelParameters=struct();

ModelParameters.Size=meta.taille; 
ModelParameters.Mass=meta.masse; 

ModelParameters.PelvisLowerTrunk=@PelvisLowerTrunk;
ModelParameters.UpperTrunk=@UpperTrunkClavicle;
ModelParameters.Head=@Head;
ModelParameters.LeftLeg=@Leg;
ModelParameters.RightLeg=@Leg;
ModelParameters.LeftArm=@Arm;
ModelParameters.RightArm=@Arm;

% ModelParameters.PelvisLowerTrunk=@Pelvis_TLEM;
% ModelParameters.UpperTrunk=@UpperTrunkClavicle;
% ModelParameters.Head=@Skull;
% ModelParameters.LeftLeg=@LegTLEM;
% ModelParameters.RightLeg=@LegTLEM;

ModelParameters.LeftArm=@Arm;
ModelParameters.RightArm=@Arm;

ModelParameters.Root='PelvisSacrum';

%ModelParameters.Markers=@Marker_set2;
ModelParameters.Markers=@Marker_set7;
ModelParameters.MarkersOptions=1;
ModelParameters.MarkersRemoved={};

ModelParameters.Muscles={@LegMuscles,@LegMuscles};
ModelParameters.MusclesOptions={'Right','Left'};
ModelParameters.Muscles={};
ModelParameters.MusclesOptions={};

%% AnalysisParameters
AnalysisParameters=struct();
AnalysisParameters.General.FilterActive=1;
AnalysisParameters.General.FilterCutOff=10;
AnalysisParameters.General.InputData=@C3dProcessedData;
AnalysisParameters.General.InputDataOptions={};
AnalysisParameters.General.Extension='*.c3d';

AnalysisParameters.CalibIK.filename=meta.TrialList;
AnalysisParameters.CalibIK.Active=1;
AnalysisParameters.CalibIK.Frames.Method=@UniformlyDistributed;
AnalysisParameters.CalibIK.Frames.NbFrames=50;
AnalysisParameters.CalibIK.LengthAdd={};
AnalysisParameters.CalibIK.LengthDelete={};
AnalysisParameters.CalibIK.AxisAdd={};
AnalysisParameters.CalibIK.AxisDelete={};

% AnalysisParameters.CalibIK.AxisDelete = {'RShank',[0.998807846858712,-0.00113107972216446;...
%     -0.00113107972216446,0.998926864935733;...
%     -0.0488016978406023,-0.0463016108610633];...
%     'RTalus',[0.931543714114037,-0.0246238845317481;...
%     -0.0246238845317481,0.991142731721629;...
%     0.362794670585433,0.130498082997241];...
%     'RFoot',[-0.463807411701659,-0.115201840940127;...
%     0.885479286838347,-0.0284449895563225;...
%     -0.0284449895563225,0.992934750330124];...
%     'LShank',[0.998807846858712,-0.00113107972216446;...
%     -0.00113107972216446,0.998926864935733;...
%     0.0488016978406023,0.0463016108610633];...
%     'LTalus',[0.931543714114037,-0.0246238845317481;...
%     -0.0246238845317481,0.991142731721629;...
%     -0.362794670585433,-0.130498082997241];...
%     'LFoot',[-0.463807411701659,-0.115201840940127;...
%     0.885479286838347,-0.0284449895563225;...
%     -0.0284449895563225,0.992934750330124]};

AnalysisParameters.CalibIK.MarkersCalibModif={};
% AnalysisParameters.CalibIK.MarkersCalibModif={'RTAR',{'On';'On';'On'};'LTAR',{'On';'On';'On'}};

AnalysisParameters.filename=meta.TrialList;
AnalysisParameters.IK.Active=1;
AnalysisParameters.IK.Method=2; % LM
AnalysisParameters.IK.FilterActive=1;
AnalysisParameters.IK.FilterCutOff=5;

AnalysisParameters.CalibID.Active=0;

%AnalysisParameters.ExternalForces.
% AnalysisParameters.ExternalForces.Options={'LFoot';'RFoot';'NoContact'};

AnalysisParameters.ID.Active=0;

AnalysisParameters.Muscles.Active=0;

AnalysisParameters.CalibMuscles=0;

% save('AnalysisParameters.mat','AnalysisParameters')
% save('ModelParameters','ModelParameters')