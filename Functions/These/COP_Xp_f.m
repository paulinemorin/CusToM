function [COP_Xp, Contact_detection_sum] = COP_Xp_f(filename,X)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
file=filename;
load('AnalysisParameters.mat')
filename=file;
AnalysisParameters.ExternalForces.Options{2,1}='RFoot';
AnalysisParameters.ExternalForces.Options{1,1}='LFoot';

% Get the solid names on which the forces are applied
Solids = AnalysisParameters.ExternalForces.Options;

% Loading the Biomechanicalmodel file
load('BiomechanicalModel.mat')
Human_model = BiomechanicalModel.OsteoArticularModel;


% Solid list extracted from the OsteoarticularModel
Solid_list = {BiomechanicalModel.OsteoArticularModel.name}';

% Get the numbers of solids on which the forces are applied
[~,num_s]=intersect(Solid_list,Solids);

% Loading external forces computation
load(fullfile(filename,'ExternalForcesComputationResults.mat'))


% Load experimental forces and forces from prediction algorithm.
GRF_Xp=ExternalForcesComputationResults.ExternalForcesExperiments;

load([filename '/InverseKinematicsResults.mat']); %#ok<LOAD>
q = InverseKinematicsResults.JointCoordinates';
if isfield(InverseKinematicsResults,'FreeJointCoordinates')
    q6dof = InverseKinematicsResults.FreeJointCoordinates';
else
    PelvisPosition = InverseKinematicsResults.PelvisPosition;
    PelvisOrientation = InverseKinematicsResults.PelvisOrientation;
end
load([filename '/ExperimentalData.mat']); %#ok<LOAD>


time = ExperimentalData.Time;
freq=1/time(2);
nbframe=size(q,1);

%% get rid of the 6DOF joint
if isfield(InverseKinematicsResults,'FreeJointCoordinates')
    Human_model(Human_model(end).child).mother = 0;
    Human_model=Human_model(1:(numel(Human_model)-6));
end

dt=1/freq;
dq=derivee2(dt,q);  % vitesses
ddq=derivee2(dt,dq);  % accélérations

%% Définition des données cinématiques du pelvis
% (position / vitesse / accélération / orientation / vitesse angulaire / accélération angulaire)
% Kinematical data for Pelvis (Position/speed/acceleration/angles/angular speed/angular acceleration)

if isfield(InverseKinematicsResults,'FreeJointCoordinates')
    p_pelvis=q6dof(:,1:3);  % frame i : p_pelvis(i,:)
    r_pelvis=cell(size(q6dof,1),1);
    for i=1:size(q6dof,1)
        r_pelvis{i}=Rodrigues([1 0 0]',q6dof(i,4))*Rodrigues([0 1 0]',q6dof(i,5))*Rodrigues([0 0 1]',q6dof(i,6)); % matrice de rotation en fonction des rotations successives (x,y,z) : frame i : r_pelvis{i}
    end
else
    p_pelvis = cell2mat(PelvisPosition)';
    r_pelvis  = PelvisOrientation';
end

%dR
dR=zeros(3,3,nbframe);
for ligne=1:3
    for colonne=1:3
        dR(ligne,colonne,:)=derivee2(dt,cell2mat(cellfun(@(b) b(ligne,colonne),r_pelvis,'UniformOutput',false)));
    end
end
w=zeros(nbframe,3);
for i=1:nbframe
    wmat=dR(:,:,i)*r_pelvis{i}';
    w(i,:)=[wmat(3,2),wmat(1,3),wmat(2,1)];
end

% v0
v=derivee2(dt,p_pelvis);
vw=zeros(nbframe,3);
for i=1:nbframe
    vw(i,:)=cross(p_pelvis(i,:),w(i,:));
end
v0=v+vw;

% dv0
dv0=derivee2(dt,v0);

% dw
dw=derivee2(dt,w);

%% Creation of a structure to add contact points
for i=1:numel(AnalysisParameters.Prediction.ContactPoint)
    Prediction(i).points_prediction_efforts = AnalysisParameters.Prediction.ContactPoint{i}; %#ok<AGROW>
end

Prediction=verif_Prediction_Humanmodel(Human_model,Prediction);
NbPointsPrediction = numel(Prediction);

%% Contact detection
Contact_detection = ContactDetectionInsole(filename);

Contact_detection_sum(1,:)=sum(Contact_detection(1:16,:),1);
Contact_detection_sum(2,:)=sum(Contact_detection(17:32,:),1);

%% Determination cOP Xp et detection contact sur plateforme

if GRF_Xp(1).Visual(4:6,1)==GRF_Xp(1).fext(num_s(1)).fext(:,1)
    Same_Xp=1;
else
    Same_Xp=-1;
end

for ii=1:numel(num_s)
    COP_Xp.(Solids{ii})=zeros(nbframe,3);
    for i=1:nbframe
        if Same_Xp==1
            if (GRF_Xp(i).Visual(6,ii)-GRF_Xp(i).Visual(3,ii))>75 && Contact_detection_sum(ii,i)>1
                COP_Xp.(Solids{ii})(i,:)=GRF_Xp(i).Visual(1:3,ii);
            else
                Contact_detection_sum(ii,i)=0;
            end
        else
            if (GRF_Xp(i).Visual(6,3-ii)-GRF_Xp(i).Visual(3,3-ii))>75 && Contact_detection_sum(ii,i)>1
                COP_Xp.(Solids{ii})(i,:)=GRF_Xp(i).Visual(1:3,3-ii);
            else
                Contact_detection_sum(ii,i)=0;
            end
        end
    end
end



% figure
% hold on
% 
% for ii=[1,2]
%     subplot(1,2,ii)
%     COP_Xp.(Solids{ii})(COP_Xp.(Solids{ii})==0) = NaN;
% %      x=CoP_Pos_G((ii-1)*3+2,:);
% %      y=CoP_Pos_G((ii-1)*3+1,:);
% %     plot(y,x,'*')
%     hold on
%     x=COP_Xp.(Solids{ii})(:,2);
%     y=COP_Xp.(Solids{ii})(:,1);
%     plot(y,x,'*')
%     
% end

end

