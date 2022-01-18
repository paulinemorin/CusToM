function X = CalibrationInsoleModel(filename)

options = optimoptions(@fmincon,'Algorithm','sqp','Display','final','GradObj','off','GradConstr','off','TolFun',1e-8,'TolX',1e-8,'MaxFunctionEvaluations', 10000);

X=zeros(12,1);
X0=zeros(12,1);

ubl = [-15;
-15;
-15;
-0.1;
-0.1;
-0.1];
ubr = [-15;
0;
-15;
-0.1;
-0.1;
-0.1];


ub =[ubl;ubr];

lbl = [15;
0;
0;
0;
0.1;
0.1];
lbr = [15;
15;
0;
0;
0.1;
0.1];

lb =[lbl;lbr];

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

% Loading external forces computation
load(fullfile(filename,'ExternalForcesComputationResults.mat'))

[start,stop] = croper(ExternalForcesComputationResults,AnalysisParameters, filename);

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

%% taille de pied et semelles pour normaliser
[Foot_length,Foot_large] = Foot_size(filename);
[Foot_length_insole,Foot_large_insole] = Foot_size_insole(Foot_length,Foot_large);
taille=[Foot_length_insole,Foot_large_insole];

%% COP Semelle
InsoleDataFile(filename)
load(fullfile(filename,'InsoleData.mat'))

CoP_Pos(1:2,:) = InsoleData.data(:,22:23)';
CoP_Pos(3:4,:) = InsoleData.data(:,44:45)';

%[Contact_detection, CoP_Pos] = ContactDetectionSoleSynchPlateforme(filename, AnalysisParameters, Human_model,BiomechanicalModel, p_pelvis, r_pelvis, v0, w, dv0, dw, q, dq, ddq);
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
            if (GRF_Xp(i).Visual(6,ii)-GRF_Xp(i).Visual(3,ii))>70 && Contact_detection_sum(ii,i)>0
                COP_Xp.(Solids{ii})(i,:)=GRF_Xp(i).Visual(1:3,ii);
            else
                Contact_detection_sum(ii,i)=0;
            end
        else
            if (GRF_Xp(i).Visual(6,3-ii)-GRF_Xp(i).Visual(3,3-ii))>70 && Contact_detection_sum(ii,i)>0
                COP_Xp.(Solids{ii})(i,:)=GRF_Xp(i).Visual(1:3,3-ii);
            else
                Contact_detection_sum(ii,i)=0;
            end
        end
    end
end

[time_g, time_d]= TimingFootFlat(Contact_detection_sum(:,start:stop));
time_g=time_g+start;
time_d=time_d+start;

%[X,~,~,output] = fmincon(@(X) RMSEInsoleModelCOP(X,CoP_Pos,PositionThreshold,VelocityThreshold,dw,dv0,dR,Human_model,Contact_detection,Prediction,nbframe,NbPointsPrediction,Same_Xp,Same_Pred,Foot_large,Foot_length,Foot_large_insole,Foot_length_insole,GRF_Pred,GRF_Xp,p_pelvis,r_pelvis,v0,w,Contact_sum,Solids,q,dq,ddq,num_s,filename,Contact_detect),X0,[],[],[],[],ub,lb,[],options);
%[X,~,~,output] = fmincon(@(X) CalibRMSEInsole(X,CoP_Pos,COP_Xp, dw,dv0,Human_model,Contact_detection_sum,Prediction,nbframe,Foot_large,Foot_length,taille, p_pelvis,r_pelvis,v0,w,Solids,q,dq,ddq,num_s,time_g,time_d),X0,[],[],[],[],ub,lb,[],options);

[X,~,~,output] = fmincon(@(X) CalibRMSEInsole(X,CoP_Pos,COP_Xp, dw,dv0,Human_model,Contact_detection_sum,Prediction,nbframe,Foot_large,Foot_length,taille, p_pelvis,r_pelvis,v0,w,Solids,q,dq,ddq,num_s,time_g,time_d, filename,start,stop),X0,[],[],[],[],[],[],[],options);

%[X,~,~,output] = fmincon(@(X) CalibRMSEInsole(X,CoP_Pos,COP_Xp, dw,dv0,Human_model,Contact_detection_sum,Prediction,nbframe,Foot_large,Foot_length,taille, p_pelvis,r_pelvis,v0,w,Solids,q,dq,ddq,num_s),X0,[],[],[],[],[],[],[],options);
%[X,~,~,output] = fminsearch(@(X) CalibRMSEInsole(X,CoP_Pos,COP_Xp, dw,dv0,Human_model,Contact_detection_sum,Prediction,nbframe,Foot_large,Foot_length,taille, p_pelvis,r_pelvis,v0,w,Solids,q,dq,ddq,num_s),X0,optimset);
%[X,~,~,output] = fminunc(@(X) CalibRMSEInsole(X,CoP_Pos,COP_Xp, dw,dv0,Human_model,Contact_detection_sum,Prediction,nbframe,Foot_large,Foot_length,taille, p_pelvis,r_pelvis,v0,w,Solids,q,dq,ddq,num_s),X0);

%% Définition passage modèle de pied à semelle


R1_l = rotx(X(1));
R2_l = roty(X(2));
R3_l = rotz(X(3));

Rx.(Solids{1}) =R3_l*R2_l*R1_l;
Locx.(Solids{1})=X(4:6);

R1_r = rotx(X(7));
R2_r = roty(X(8));
R3_r = rotz(X(9));

Rx.(Solids{2}) = R3_r*R2_r*R1_r;
Locx.(Solids{2})=X(10:12);

%% Transformation cop data semelles

for ii=[1,2]
    CoP_Pos_G((ii-1)*3+1,:)=[taille(1).*CoP_Pos((ii-1)*2+1,:)];
    % CoP_Pos_G((ii-1)*3+1,:)=[taille(1).*CoP_Pos((ii-1)*2+1,:)+taille(1)*0.5];
    CoP_Pos_G((ii-1)*3+2,:)=[taille(2).*CoP_Pos((ii-1)*2+2,:)];
    CoP_Pos_G((ii-1)*3+3,:)=zeros(1,length( CoP_Pos_G((ii-1)*3+2,:)));
end
 
% changement repère pou_r repère direct, pied droit semelle
CoP_Pos_G(2,:)=-CoP_Pos_G(2,:);

% changement repère de semelle à modèle de pied
a = [0 -1 0;
    0  0 -1;
    1 0  0];

% for i=1:nbframe  
%     % Get the forces applied on the solids
%     for ii=1:numel(num_s)
%         if Contact_detection_sum(ii,i)>0
%             CoP_Pos_G((ii-1)*3+1:(ii-1)*3+3,i)=a'*CoP_Pos_G((ii-1)*3+1:(ii-1)*3+3,i);
%         end
%     end
% end



for i=[time_g,time_d]
    [Rot_R, R_pos, Rot_L, L_pos] = FootFrameMarker(filename, i);

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
    
    % Get the forces applied on the solids
    if i == time_d
            %if i == 164
%         R.(Solids{2}) = Human_model(num_s(2)).R;
%         Loc.(Solids{2}) = Human_model(num_s(2)).p;
                R.(Solids{2}) = Rot_R;
       Loc.(Solids{2}) = R_pos';
%          Frame_Pied_Droit=Foot_trial_frame(filename,i,2);
%     R.(Solids{2})= Frame_Pied_Droit(1:3,1:3);
%     Loc.(Solids{2}) =Frame_Pied_Droit(1:3,4);
%     
    end
    if i == time_g
            %if i == 238
%         R.(Solids{1}) = Human_model(num_s(1)).R;
%         Loc.(Solids{1}) = Human_model(num_s(1)).p;
                        R.(Solids{1}) = Rot_L;
       Loc.(Solids{1}) = L_pos';
%         Frame_Pied_Gauche=Foot_trial_frame(filename,i,1);
%     R.(Solids{2})= Frame_Pied_Gauche(1:3,1:3);
%     Loc.(Solids{2}) =Frame_Pied_Gauche(1:3,4);
    
    end
end

for i=1:nbframe  
    % Get the forces applied on the solids
    for ii=1:numel(num_s)
        if Contact_detection_sum(ii,i)>0
            CoP_Pos_G((ii-1)*3+1:(ii-1)*3+3,i)=R.(Solids{ii})*(Rx.(Solids{ii})*(a'*CoP_Pos_G((ii-1)*3+1:(ii-1)*3+3,i)+Locx.(Solids{ii})))+Loc.(Solids{ii});
        end
    end
end

suivi_L=[];
suivi_R=[];

for i = start: stop
%for i = 1 : nbframe
    if Contact_detection_sum(1,i)>0
        suivi_L=[suivi_L,i];
    end
    if Contact_detection_sum(2,i)>0
        suivi_R=[suivi_R,i];
    end
end
suivi=[suivi_L,suivi_R];

for ii=1:numel(num_s)
    %cur_s=num_s(ii); %LFoot and RFoot
    COP_Insole.(Solids{ii})=CoP_Pos_G((ii-1)*3+1:(ii-1)*3+3,:);

    RMSE(ii,1)=rms(COP_Insole.(Solids{ii})(1,suivi((ii-1)*length(suivi_L)+1:end-(2-ii)*length(suivi_R)))'-COP_Xp.(Solids{ii})(suivi((ii-1)*length(suivi_L)+1:end-(2-ii)*length(suivi_R))',1));
    RMSE_norm(ii,1)=rms((COP_Insole.(Solids{ii})(1,suivi((ii-1)*length(suivi_L)+1:end-(2-ii)*length(suivi_R)))'-COP_Xp.(Solids{ii})(suivi((ii-1)*length(suivi_L)+1:end-(2-ii)*length(suivi_R))',1))/Foot_length);
    RMSE(ii,2)=rms(COP_Insole.(Solids{ii})(2,suivi((ii-1)*length(suivi_L)+1:end-(2-ii)*length(suivi_R)))'-COP_Xp.(Solids{ii})(suivi((ii-1)*length(suivi_L)+1:end-(2-ii)*length(suivi_R))',2));
    RMSE_norm(ii,2)=rms((COP_Insole.(Solids{ii})(2,suivi((ii-1)*length(suivi_L)+1:end-(2-ii)*length(suivi_R)))'-COP_Xp.(Solids{ii})(suivi((ii-1)*length(suivi_L)+1:end-(2-ii)*length(suivi_R))',2))/Foot_large);
end

for i = 1: nbframe
    if Contact_detection_sum(1,i)==0
        CoP_Pos_G(1:3,i)= [NaN; NaN; NaN];
    end
    if Contact_detection_sum(2,i)==0
        CoP_Pos_G(4:6,i)= [NaN; NaN; NaN];
    end
end

figure
hold on
set(gca,'color','w')
set(gcf,'color','w')
title_liste={'Left Foot','Right Foot'};

for ii=[1,2]
    subplot(1,2,ii)
    COP_Xp.(Solids{ii})(COP_Xp.(Solids{ii})==0) = NaN;
    x=CoP_Pos_G((ii-1)*3+2,suivi((ii-1)*length(suivi_L)+1:end-(2-ii)*length(suivi_R)));
    y=CoP_Pos_G((ii-1)*3+1,suivi((ii-1)*length(suivi_L)+1:end-(2-ii)*length(suivi_R)));
    plot(y,x,'*')
    hold on
    x=COP_Xp.(Solids{ii})(suivi((ii-1)*length(suivi_L)+1:end-(2-ii)*length(suivi_R)),2);
    y=COP_Xp.(Solids{ii})(suivi((ii-1)*length(suivi_L)+1:end-(2-ii)*length(suivi_R)),1);
    plot(y,x,'*')
        axis([nanmean(y)-0.5*taille(2) nanmean(y)+0.5*taille(2) nanmean(x)-0.5*taille(1) nanmean(x)+0.5*taille(1)])
    title(title_liste(ii))
    %axis([-0.5*taille(2) 0.5*taille(2) -0.5*taille(1) 0.5*taille(1)])
    ylabel({'y dans le repère global','longeur semelle (m)'})
    xlabel({'x dans le repère global','largeur semelle (m)'})
    
end

legend('CoP from insole','CoP from plateform')

figure
hold on

for ii=[1,2]
    %subplot(1,2,ii)
    COP_Xp.(Solids{ii})(COP_Xp.(Solids{ii})==0) = NaN;
    y=CoP_Pos_G((ii-1)*3+2,start:stop);
    z=CoP_Pos_G((ii-1)*3+3,start:stop);
    x=CoP_Pos_G((ii-1)*3+1,start:stop);
    plot3(x,y,z,'*')
    hold on
    y=COP_Xp.(Solids{ii})(start:stop,2);
    x=COP_Xp.(Solids{ii})(start:stop,1);
    z=COP_Xp.(Solids{ii})(start:stop,3);
    plot3(x,y,z,'*')
end


%figure
Rxframel=eye(4);
Rxframel(1:3,1:3)=Rx.(Solids{1});
Rxframel(1:3,4)=Locx.(Solids{1});

Rframel=eye(4);
Rframel(1:3,1:3)=R.(Solids{1});
Rframel(1:3,4)=Loc.(Solids{1});


triad('matrix',Rframel,'scale',0.5,'Tag','Rmodel')
hold on
triad('matrix',Rframel*Rxframel,'scale',0.5,'LineWidth',3,'Tag','Rsemelle')

Rxframer=eye(4);
Rxframer(1:3,1:3)=Rx.(Solids{2});
Rxframer(1:3,4)=Locx.(Solids{2});

Rframer=eye(4);
Rframer(1:3,1:3)=R.(Solids{2});
Rframer(1:3,4)=Loc.(Solids{2});


triad('matrix',Rframer,'scale',0.5,'Tag','Rmodel')
hold on
triad('matrix',Rframer*Rxframer,'scale',0.5,'LineWidth',3,'Tag','Rsemelle')

triad('matrix',eye(4))

end

