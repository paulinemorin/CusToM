

load('AnalysisParameters.mat')
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


% Name of the trial treated
%filename='ChgtDirection04';
filename='Step0001';
%filename='Trial001';

% Loading external forces computation
load(fullfile(filename,'ExternalForcesComputationResults.mat'))


% Load experimental forces and forces from prediction algorithm.
GRF_Xp=ExternalForcesComputationResults.ExternalForcesExperiments;
GRF_Pred=ExternalForcesComputationResults.ExternalForcesPrediction;

% Number of frames
Nb_frames = numel(GRF_Xp);
Contact_sum=zeros(numel(num_s),Nb_frames);

for ii=1:numel(num_s)
    COP_Pred.(Solids{ii})=zeros(Nb_frames,3);
    COP_Xp.(Solids{ii})=zeros(Nb_frames,3);
    COP_Insole.(Solids{ii})=zeros(Nb_frames,3);
end


% ordre solides pour visual identiques ou non exp et pred
if GRF_Pred(1).Visual(4:6,1)==GRF_Pred(1).fext(num_s(1)).fext(:,1)
    Same_Pred=1;
else
    Same_Pred=-1;
end

if GRF_Xp(1).Visual(4:6,1)==GRF_Xp(1).fext(num_s(1)).fext(:,1)
    Same_Xp=1;
else
    Same_Xp=-1;
end


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

%% Creation of a structure to add contact points
for i=1:numel(AnalysisParameters.Prediction.ContactPoint)
    Prediction(i).points_prediction_efforts = AnalysisParameters.Prediction.ContactPoint{i}; %#ok<AGROW>
end

Prediction=verif_Prediction_Humanmodel(Human_model,Prediction);
NbPointsPrediction = numel(Prediction);

%% Contact detection initialization
Contact_detection = zeros(NbPointsPrediction, nbframe);

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

PositionThreshold = AnalysisParameters.Prediction.PositionThreshold;
VelocityThreshold = AnalysisParameters.Prediction.VelocityThreshold;



%% taille de pied et semelles pour normaliser
[Foot_length,Foot_large] = Foot_size(filename);
[Foot_length_insole,Foot_large_insole] = Foot_size_insole(Foot_length,Foot_large);

%% COP Semelle
[Contact_detect, CoP_Pos] = ContactDetectionSole(filename, AnalysisParameters, Human_model,BiomechanicalModel, p_pelvis, r_pelvis, v0, w, dv0, dw, q, dq, ddq);

insole_foot_left=zeros(4);
a = [0 -1 0;
    0  0 -1;
    1 0  0];
insole_foot_left(1:3,1:3)=a';
insole_foot_left(1:3,4) =a'*[Foot_length_insole/2 0 0]';
insole_foot_left(4,4)=1;

insole_foot_right=zeros(4);
a = [0 -1 0;
    0  0 1;
    1 0  0];
insole_foot_right(1:3,1:3)=a';
insole_foot_right(1:3,4) =a'*[Foot_length_insole/2 0 0]';
insole_foot_right(4,4)=1;

insole_foot{1} = insole_foot_left;
insole_foot{2} = insole_foot_right;


for i=1:nbframe
    %attribution à chaque articulation de la position/vitesse/accélération (position/speed/acceleration for each joint)
    Human_model(1).p=p_pelvis(i,:)';
    Human_model(1).R=r_pelvis{i};
    Human_model(1).v0=v0(i,:)';
    Human_model(1).w=w(i,:)';
    Human_model(1).dv0=dv0(i,:)';
    Human_model(1).dw=dw(i,:)';
    for j=2:numel(Human_model)
        Human_model(j).q=q(i,j); %#ok<*SAGROW>
        Human_model(j).dq=dq(i,j);
        Human_model(j).ddq=ddq(i,j);
    end
    
    % Calcul positions / vitesses / accélération de chaque solide (computation of position/speed/acceleration for each solid)
    [Human_model,Prediction] = ForwardAllKinematicsPrediction(Human_model,Prediction,1);
    
    %% Threshold application
    for pred = 1:numel(Prediction)
        Prediction(pred).vitesse_temps(i)=sqrt(Prediction(pred).vitesse(1,:)^2+Prediction(pred).vitesse(2,:)^2+Prediction(pred).vitesse(3,:)^2); % Recuperation de la norme de la vitesse (repère monde)
        
        if Prediction(pred).pos_anim(3)<PositionThreshold && abs(Prediction(pred).vitesse_temps(i))<VelocityThreshold
            Contact_detection(pred,i)=1;
        end
    end
    
    
    for k = 1:numel(num_s)
        for l = 1:size(Contact_detection,1)
            if Prediction(l).num_solid==num_s(k)
                    Contact_sum(k, i)=Contact_sum(k, i)+Contact_detection(l,i);
            end
        end
    end
    
    R_g_to_f=zeros(4);
    % Get the forces applied on the solids
    for ii=1:numel(num_s)
       R = Human_model(num_s(ii)).R;
       R_g_to_f(1:3,1:3)=R';
       %R_g_to_f(1:3,1:3)=eye(3);
       Loc = R'* Human_model(num_s(ii)).p;
       %Loc = Human_model(num_s(ii)).p;
       %Loc=zeros(1,3);
       R_g_to_f(1:3,4) =-Loc;
       R_g_to_f(4,4)=1;
       cur_s=num_s(ii); %LFoot and RFoot

if Same_Xp==1
    if Same_Pred==1
        if norm(GRF_Xp(i).Visual(4:6,ii)-GRF_Xp(i).Visual(1:3,ii))>50 && Contact_sum(ii,i)>0
            a =  [R_g_to_f*[GRF_Xp(i).Visual(1:3,ii);1]]';
            COP_Xp.(Solids{ii})(i,:)= a(1:3);
            a =  [R_g_to_f*[GRF_Pred(i).Visual(1:3,ii);1]]';
            COP_Pred.(Solids{ii})(i,:) = a(1:3);
            a =  [R_g_to_f*insole_foot{ii}*[[Foot_length_insole,Foot_large_insole]*CoP_Pos(1:2,i);0;1]]';
            COP_Insole.(Solids{ii})(i,:) = a(1:3);
         end
    else
        if norm(GRF_Xp(i).Visual(4:6,ii)-GRF_Xp(i).Visual(1:3,ii))>50 && Contact_sum(ii,i)>0
            a =  [R_g_to_f*[GRF_Xp(i).Visual(1:3,ii);1]]';
            COP_Xp.(Solids{ii})(i,:)= a(1:3);
            a =  [R_g_to_f*[GRF_Pred(i).Visual(1:3,3-ii);1]]';
            COP_Pred.(Solids{ii})(i,:) = a(1:3);
            a =  [R_g_to_f*insole_foot{ii}*[[Foot_length_insole,Foot_large_insole]*CoP_Pos(1:2,i);0;1]]';
            COP_Insole.(Solids{ii})(i,:) = a(1:3);
        end
    end
else
    if Same_Pred==1
        if norm(GRF_Xp(i).Visual(4:6,3-ii)-GRF_Xp(i).Visual(1:3,3-ii))>50&&Contact_sum(ii,i)>0
            a =  [R_g_to_f*[GRF_Xp(i).Visual(1:3,3-ii);1]]';
            COP_Xp.(Solids{ii})(i,:)= a(1:3);
            a =  [R_g_to_f*[GRF_Pred(i).Visual(1:3,ii);1]]';
            COP_Pred.(Solids{ii})(i,:) = a(1:3);
            a =  [R_g_to_f*insole_foot{ii}*[[Foot_length_insole,Foot_large_insole]*CoP_Pos(1:2,i);0;1]]';
            COP_Insole.(Solids{ii})(i,:) = a(1:3);
        end
    else
        if norm(GRF_Xp(i).Visual(4:6,3-ii)-GRF_Xp(i).Visual(1:3,3-ii))>50&&Contact_sum(ii,i)>0
            a =  [R_g_to_f*[GRF_Xp(i).Visual(1:3,3-ii);1]]';
            COP_Xp.(Solids{ii})(i,:)= a(1:3);
            a =  [R_g_to_f*[GRF_Pred(i).Visual(1:3,3-ii);1]]';
            COP_Pred.(Solids{ii})(i,:) = a(1:3);
            a =  [insole_foot{ii}*[[[Foot_length_insole,Foot_large_insole].*CoP_Pos(1:2,i)']';0;1]]';
            COP_Insole.(Solids{ii})(i,:) = a(1:3);

        end
    end
end
    end
    
end


%% RMSE 

RMSE=zeros(numel(num_s),2);
% Get the forces applied on the solids
for ii=1:numel(num_s)
    cur_s=num_s(ii); %LFoot and RFoot
   
    RMSE(ii,1)=rms(COP_Pred.(Solids{ii})(:,2)-COP_Xp.(Solids{ii})(:,2));
    RMSE_norm(ii,1)=rms((COP_Pred.(Solids{ii})(:,2)-COP_Xp.(Solids{ii})(:,2))/Foot_length);
    RMSE(ii,2)=rms(COP_Pred.(Solids{ii})(:,3)-COP_Xp.(Solids{ii})(:,3));
    RMSE_norm(ii,2)=rms((COP_Pred.(Solids{ii})(:,3)-COP_Xp.(Solids{ii})(:,3))/Foot_large);
end

RMSE = mean(RMSE,1);
RMSE_norm = mean(RMSE_norm,1);


%% Plot the result
Directions={'X','Y','Z'};
figure();
sgtitle(['COP error position during ', filename])

set(gcf,'color','w')
for ii=1:3
%for ii=1:2
    subplot(1,3,ii)
    % LFoot
    plot(COP_Pred.(Solids{1})(:,ii),'b-','LineWidth',2)
    hold on
    plot(COP_Xp.(Solids{1})(:,ii),'b--','LineWidth',2)
    plot(COP_Insole.(Solids{1})(:,ii),'b.','LineWidth',2)
    % RFoot
    plot(COP_Pred.(Solids{2})(:,ii),'r-','LineWidth',2)
    plot(COP_Xp.(Solids{2})(:,ii),'r--','LineWidth',2)
    plot(COP_Insole.(Solids{2})(:,ii),'r.','LineWidth',2)
    
    
    %xlim([1 203])
    xlabel('Frames')
    ylabel('Position (m)')
    title(['COP position on the ' Directions{ii} '-direction'])
end

legend('Predicted on Left Foot','Experimental on Left Foot',...
    'Predicted on Right Foot','Experimental on Right Foot')

% On the graph, we can see that the predicted results are comparable to the
% experimental forces measured. However there is artefact when the foot is
% not on the plateforms. More work is needed on filters and algorithms to
% avoid it.