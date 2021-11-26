function [ExternalForcesComputationResults] = ExternalForcesPredictionCOPPlateform(filename, AnalysisParameters, BiomechanicalModel, ModelParameters)
% Prediction of ground reaction forces
%   Ground reaction forces are predicted from motion data.
%
%	Based on :
%	- Fluit, R., Andersen, M. S., Kolk, S., Verdonschot, N., & Koopman, H. F., 2014.
%	Prediction of ground reaction forces and moments during various activities of daily living. Journal of biomechanics, 47(10), 2321-2329.
%	- Skals, S., Jung, M. K., Damsgaard, M., & Andersen, M. S., 2017. 
%	Prediction of ground reaction forces and moments during sports-related movements. Multibody system dynamics, 39(3), 175-195.
%
%   INPUT
%   - filename: name of the file to process (character string)
%   - AnalysisParameters: parameters of the musculoskeletal analysis,
%   automatically generated by the graphic interface 'Analysis' 
%   - BiomechanicalModel: musculoskeletal model
%   - ModelParameters: parameters of the musculoskeletal model, automatically
%   generated by the graphic interface 'GenerateParameters' 
%   OUTPUT
%   - ExternalForcesComputationResults: results of the external forces
%   computation (see the Documentation for the structure)
%________________________________________________________
%
% Licence
% Toolbox distributed under GPL 3.0 Licence
%________________________________________________________
%
% Authors : Antoine Muller, Charles Pontonnier, Pierre Puchaud and
% Georges Dumont
%________________________________________________________

disp(['External Forces Prediction (' filename ') ...'])

% AnalysisParameters.ExternalForces.Options(2,1)={'LFoot'};
% AnalysisParameters.ExternalForces.Options(1,1)={'RFoot'};


%% Loading data
Human_model = BiomechanicalModel.OsteoArticularModel;
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


%% Contact detection

% if AnalysisParameters.Prediction.ContactDetection == 0
     Contact_detection = ContactDetectionThreshold(filename, AnalysisParameters, Human_model);

     % %elseif AnalysisParameters.Prediction.ContactDetection == 1
%    Contact_detection = ContactDetectionAutomaticThreshold(filename, AnalysisParameters, BiomechanicalModel);
% elseif AnalysisParameters.Prediction.ContactDetection == 2
%     Contact_detection = ContactDetectionOne(filename, AnalysisParameters);
% elseif AnalysisParameters.Prediction.ContactDetection == 3
%     [Contact_detection, Num] = ContactDetectionAddMarkers(filename, AnalysisParameters, Human_model);
% %elseif AnalysisParameters.Prediction.ContactDetection == 4
% %    Contact_detection = Perso(filename, AnalysisParameters);
% end

% X =[
% 
%     8.6620
%    -4.8632
%     4.1089
%    -0.0925
%     0.0429
%     0.0138
%    -9.9464
%     3.3038
%    -0.1245
%    -0.0790
%     0.0887
%    -0.0010];

% PositionCoP = CoPInsole(X,filename);
Contact_detection_sum(1,:)=sum(Contact_detection(1:16,:),1);
Contact_detection_sum(2,:)=sum(Contact_detection(17:32,:),1);
Contact_detection=Contact_detection_sum;

%[COP_Xp, Contact_detection] = COP_Xp_f(filename);
COP_Xp = COP_Xp_s(filename);


%% Gravity
g=[0 0 -9.81]';

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

% 
% [abc,PositionCoP] = ContactDetectionSole(filename, AnalysisParameters, Human_model,BiomechanicalModel, p_pelvis, r_pelvis, v0, w, dv0, dw,q, dq, ddq); 
% %[Contact_detection,PositionCoP] = ContactDetectionSole(filename, AnalysisParameters, Human_model,BiomechanicalModel, p_pelvis, r_pelvis, v0, w, dv0, dw,q, dq, ddq); 

 % Loading external forces computation
load(fullfile(filename,'ExternalForcesComputationResults.mat'))

% PositionCoP=zeros(4,nbframe);
% 
% for i=1:nbframe
%     % if norm(ExternalForcesComputationResults.ExternalForcesExperiments(i).Visual(4:6,1)-ExternalForcesComputationResults.ExternalForcesExperiments(i).Visual(1:3,1)) > 50
%         PositionCoP([1,2],i)=ExternalForcesComputationResults.ExternalForcesExperiments(i).Visual(1:2,2).';
%     % end
%      %if norm(ExternalForcesComputationResults.ExternalForcesExperiments(i).Visual(4:6,2)-ExternalForcesComputationResults.ExternalForcesExperiments(i).Visual(1:3,2)) > 50
%         PositionCoP([3,4],i)=ExternalForcesComputationResults.ExternalForcesExperiments(i).Visual(1:2,1).';
%      %end
%              
% end




% if ModelParameters.Size < 1.70
%     length_sole = 248.6 ;
%     width_sole = 90.2 ;
% else
%     length_sole = 274.2 ;
%     width_sole = 97.5 ;
% end
% 
% PositionCoP([1, 3],:) = PositionCoP([1, 3],:)*length_sole+0.5*length_sole*ones(2,nbframe);
% PositionCoP([2, 4],:) = PositionCoP([2, 4],:)*width_sole;

%% Creation of a structure to add contact points
AnalysisParameters.Prediction.ContactPoint={''};
AnalysisParameters.Prediction.ContactPoint(1)={'LFootMID'};
AnalysisParameters.Prediction.ContactPoint(2)={'RFootMID'};

for i=1:numel(AnalysisParameters.Prediction.ContactPoint)
    Prediction(i).points_prediction_efforts = AnalysisParameters.Prediction.ContactPoint{i}; %#ok<AGROW>
end
Prediction=verif_Prediction_Humanmodel(Human_model,Prediction);
NbPointsPrediction = numel(Prediction);

%% Initialisations pour la dynamique inverse

torques=zeros(nbframe,numel(Human_model));
f6dof=zeros(3,nbframe);
t6dof0=zeros(3,nbframe);
t6dof=t6dof0;
FContactDyn=struct('F',[],'T',[]);

%% Param�tres de l'optimisation fmincon pour probleme lineaire
if ~isempty(BiomechanicalModel.Muscles)
    Muscles = BiomechanicalModel.Muscles;
    idm = logical([Muscles.exist]);
    Nb_muscles=numel(Muscles(idm));
else
    Nb_muscles=0;
end





%% Initialisations des différents efforts et leur stockage
for f=1:nbframe
    for n=1:numel(Human_model)
        external_forces_pred(f).fext(n).fext=zeros(3,2); %#ok<AGROW>
    end
end

for pred = 1:NbPointsPrediction
    Prediction(pred).efforts_max=zeros(nbframe,3);
    Prediction(pred).efforts = zeros(nbframe,1);
end
Fx=zeros(NbPointsPrediction,nbframe);
Fy=zeros(NbPointsPrediction,nbframe);
Fz=zeros(NbPointsPrediction,nbframe);

%% Paramètres de l'optimisation fmincon pour probleme lineaire
X0= 1*zeros(3*NbPointsPrediction,1);
lb=-ones(3*NbPointsPrediction,1);
lb([3,6])=0;
ub=ones(3*NbPointsPrediction,1);
lb_init=lb; ub_init=ub;

options = optimoptions(@fmincon,'Algorithm','sqp','Display','off','GradObj','off','GradConstr','off','TolFun',1e-6,'TolX',1e-6);
optionsLM = optimset('Algorithm','Levenberg-Marquardt','Display','final','MaxIter',4e6,'MaxFunEval',5e6);


%% Calcul frame par frame
h = waitbar(0,['External Forces Prediction (' filename ')']);
Mass = ModelParameters.Mass;
coef_friction = AnalysisParameters.Prediction.FrictionCoef;

% Pos_gauche = zeros(nbframe,3);
% COP_gauche = zeros(nbframe,3);
% 
% Pos_droite = zeros(nbframe,3);
% COP_droite = zeros(nbframe,3);
% 
% COP_global = zeros(3,2);

for i=1:nbframe
    Activ_Contact_Point = sum(Contact_detection(:,i));
    
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
    
    
    %% Position COP
    
%     COP_left_insole = [PositionCoP(1,i) PositionCoP(2,i) 0]; % insole frame
%     COP_right_insole = [PositionCoP(4,i) PositionCoP(5,i) 0]; % insole frame
%     
%     
%     insole_foot_left = [0 -1 0;
%                      0  0 -1;
%                       1 0  0];
%                   
%     insole_foot_right = [0 -1 0;
%                      0  0 1;
%                       1 0  0];
%     
%     foot_global_right = Human_model(22).R ;
%     foot_global_left= Human_model(28).R ;
    
%     COP_global (:,1) = foot_global_left * insole_foot_left' * COP_left_insole';
%     COP_global (:,2) = foot_global_right * insole_foot_right' * COP_right_insole';
    
%      COP_global (:,1) = COP_left_insole';
%      COP_global (:,2) =  COP_right_insole';
     COP_global (:,1) = COP_Xp.LFoot(i,:)';
     COP_global (:,2) =  COP_Xp.RFoot(i,:)';
    
    
    %% Calcul des efforts maximaux disponibles (computation of maximum available effort)
    
    
%     Pos_gauche(i,:) = Prediction(1).pos_anim(1:3);
%     COP_gauche(i,:) = Pos_gauche(i,:);
%     COP_gauche(i,1)=COP_gauche(i,1)+PositionCoP(1,i);
%     COP_gauche(i,2)=COP_gauche(i,2)+PositionCoP(2,i);
% 
%     Pos_droite(i,:) = Prediction(2).pos_anim(1:3);
%     COP_droite(i,:) = Pos_droite(i,:);
%     COP_droite(i,1)=COP_droite(i,1)+PositionCoP(3,i);
%     COP_droite(i,2)=COP_droite(i,2)+PositionCoP(4,2);
    
    
    for pred = 1:numel(Prediction)
%         Prediction(pred).px(i)=Prediction(pred).pos_anim(1)+COP_global(1,pred);
%         Prediction(pred).py(i)=Prediction(pred).pos_anim(2)+COP_global(2,pred);
%         Prediction(pred).pz(i)=Prediction(pred).pos_anim(3);
%       


        Activ_Contact_Point_g = sum(Contact_detection(1,i));
        Activ_Contact_Point_d = sum(Contact_detection(2,i));
       



if  Activ_Contact_Point_d ==0
    zero_p=1;
    mmtom=0;
else
    zero_p=0;
    mmtom=1;
end

        Prediction(2).px(i)=zero_p*Human_model(22).p(1)+COP_global(1,2)*mmtom;
        Prediction(2).py(i)=zero_p*Human_model(22).p(2)+COP_global(2,2)*mmtom;
        Prediction(2).pz(i)=zero_p*Human_model(22).p(3)+COP_global(3,2)*mmtom;
 
        Prediction(2).pos_anim(1)=zero_p*Human_model(22).p(1)+COP_global(1,2)*mmtom;
        Prediction(2).pos_anim(2)=zero_p*Human_model(22).p(2)+COP_global(2,2)*mmtom;
        Prediction(2).pos_anim(3)=zero_p*Human_model(22).p(3)+COP_global(3,2)*mmtom; 
        
        
        
if i==50
    a=3;
end


        
if  Activ_Contact_Point_g ==0
    zero_p=1;
    mmtom=0;
else
    zero_p=0;
    mmtom=1;
end        
 
        Prediction(1).px(i)=zero_p*Human_model(28).p(1)+COP_global(1,1)*mmtom;
        Prediction(1).py(i)=zero_p*Human_model(28).p(2)+COP_global(2,1)*mmtom;
        Prediction(1).pz(i)=zero_p*Human_model(28).p(3)+COP_global(3,1)*mmtom;
        
        
%         Prediction(pred).pos_anim(1)=Prediction(pred).pos_anim(1)+COP_global(1,pred)*0.001;
%         Prediction(pred).pos_anim(2)=Prediction(pred).pos_anim(2)+COP_global(2,pred)*0.001; 
      


        Prediction(1).pos_anim(1)=zero_p*Human_model(28).p(1)+COP_global(1,1)*mmtom;
        Prediction(1).pos_anim(2)=zero_p*Human_model(28).p(2)+COP_global(2,1)*mmtom; 
        Prediction(1).pos_anim(3)=zero_p*Human_model(28).p(3)+COP_global(3,1)*mmtom;
        
   
        Activ_Contact_Point = Contact_detection(pred,i);
        
        if Activ_Contact_Point ~= 0
            % Cpi =( numel(Prediction)/Activ_Contact_Point )* Force_max_TOR(Contact_detection(pred,i),Mass);
            Cpi =14* 0.4*Mass*9.8;
        else
            Cpi = 0;
        end
        Fx(pred,i)=Cpi;
        Fy(pred,i)=Cpi;
        Fz(pred,i)=Cpi;
        Prediction(pred).efforts_max(i,1)=Cpi; %Fx
        Prediction(pred).efforts_max(i,2)=Cpi; %Fy
        Prediction(pred).efforts_max(i,3)=Cpi; %Fz
    end
    %Fmax=[Fx(:,i)' Fy(:,i)' Fz(:,i)'];
    Fmax=[Fx(1,i) Fy(1,i) Fz(1,i) Fx(2,i) Fy(2,i) Fz(2,i)];
    
    %% Direct optimisation by linearization of the dynamical condition.
    Activ_Contact_Point = sum(Contact_detection(:,i));
    if Activ_Contact_Point ~= 0
        
        A=zeros(6,3*numel(Prediction));
        b1=[0 0 0]';
        b2=[0 0 0]';
        
        [~,b1,b2]=InverseDynamicsSolid_lin(Human_model,g,1,b1,b2);
        bf=b1;
        bt=b2+cross(-p_pelvis(i,:)',b1); %on ramene les couples au niveau du pelvis (torques are expressed at pelvis point)
        b=[bf' bt']';
        
        for k = 1:numel(Prediction)
            % calcul des efforts
            A(1,1+(k-1)*3)=Prediction(k).efforts_max(i,1);
            A(2,2+(k-1)*3)=Prediction(k).efforts_max(i,2);
            A(3,3+(k-1)*3)=Prediction(k).efforts_max(i,3);
            % calcul des moments
            A(4,2+(k-1)*3)=-(Prediction(k).pz(i)-p_pelvis(i,3))*Prediction(k).efforts_max(i,2); %-pz*beta
            A(4,3+(k-1)*3)=(Prediction(k).py(i)-p_pelvis(i,2))*Prediction(k).efforts_max(i,3); %py*gamma
            A(5,1+(k-1)*3)=(Prediction(k).pz(i)-p_pelvis(i,3))*Prediction(k).efforts_max(i,1); %pz*alpha
            A(5,3+(k-1)*3)=-(Prediction(k).px(i)-p_pelvis(i,1))*Prediction(k).efforts_max(i,3); %-px*gamma
            A(6,1+(k-1)*3)=-(Prediction(k).py(i)-p_pelvis(i,2))*Prediction(k).efforts_max(i,1); %-py*alpha
            A(6,2+(k-1)*3)=(Prediction(k).px(i)-p_pelvis(i,1))*Prediction(k).efforts_max(i,2); %px*beta
        end
        
        
        %% Taking friction into account for every point to point link, |Fx|<coeff friction|Fz| et |Fy|<coeff friction|Fz|
        Afric=zeros(4*numel(Prediction),3*numel(Prediction));
        bfric=zeros(4*numel(Prediction),1);
        
        for k = 1:(numel(Prediction))
            Afric(1+(k-1)*4,1+(k-1)*3)=1*Prediction(k).efforts_max(i,1);
            Afric(2+(k-1)*4,2+(k-1)*3)=1*Prediction(k).efforts_max(i,2);
            Afric(1+(k-1)*4,3+(k-1)*3)=-coef_friction*Prediction(k).efforts_max(i,3);
            Afric(2+(k-1)*4,3+(k-1)*3)=-coef_friction*Prediction(k).efforts_max(i,2);
            Afric(3+(k-1)*4,1+(k-1)*3)=-1*Prediction(k).efforts_max(i,1);
            Afric(4+(k-1)*4,2+(k-1)*3)=-1*Prediction(k).efforts_max(i,2);
            Afric(3+(k-1)*4,3+(k-1)*3)=-coef_friction*Prediction(k).efforts_max(i,3);
            Afric(4+(k-1)*4,3+(k-1)*3)=-coef_friction*Prediction(k).efforts_max(i,3);
         end
        
        %% Minimizing sum of normalized efforts for each punctual joint while respecting dynamical equations and friction
        %alpha_beta_gamma = [0.5 0.5 0];
        alpha_beta_gamma = [0.9 0.1 0];
        
        Biomtemp=BiomechanicalModel;
    Biomtemp.OsteoArticularModel =Human_model;
    Aopti=A;
    bopti=b;
    Aopti(4:6,:)=0;
    bopti(4:6)=0;
    
    [X,~,~,output] = fmincon(@(X) MinimizationCost(X, Biomtemp, Prediction, i, external_forces_pred,Fmax, g, Nb_muscles , alpha_beta_gamma),X0,[],[],Aopti,bopti,lb,ub,[],options);

        
        
        
        
        
       % X = fmincon(@(X) sum(X.^2),X0,Afric,bfric,A,b,lb,ub,[],options);
        %X = fmincon(@(X) sum(X.^2),X0,[],[],A,b,lb,ub,[],options);
        
        %costfunctionprediction = @(X) CostFunctionPredictionLM(X,Afric,bfric,A,b,lb,ub);
        %X = lsqnonlin(costfunctionprediction,X0,[],[],optionsLM);
    
        
    else
        X= 1*zeros(3*NbPointsPrediction,1);
    end
    
    %% Optimisation de la prochaine minimisation
    lb=max(X-0.45,lb_init); %expérimentalement, les bornes ne varient pas de plus ou moins 0.45 (experimentaly, boundaries vary not more than 0.45)
    ub=min(X+0.45,ub_init);

    X0=X; %d'une frame à l'autre, on change très peu de position, donc de valeur d'effort (
    
    %% Récupération des forces finales, stockées d'abord dans Prediction (Final forces storage without prediction)
    for k = 1:numel(Prediction)
%         Prediction(k).efforts(i,1)=X(1+numel(Prediction)*(k-1))*Prediction(k).efforts_max(i,1);
%         Prediction(k).efforts(i,2)=X(2+numel(Prediction)*(k-1))*Prediction(k).efforts_max(i,2);
%         Prediction(k).efforts(i,3)=X(3+numel(Prediction)*(k-1))*Prediction(k).efforts_max(i,3);
        Prediction(k).efforts(i,1)=X(1+3*(k-1))*Prediction(k).efforts_max(i,1);
        Prediction(k).efforts(i,2)=X(2+3*(k-1))*Prediction(k).efforts_max(i,2);
        Prediction(k).efforts(i,3)=X(3+3*(k-1))*Prediction(k).efforts_max(i,3);
    end
    
    %% Calcul des efforts extérieurs tels qu’utilisés par la suite pour la dynamique
    %% Computation of external forces for use with dynamics
    external_forces_pred=addForces_Prediction_frame_par_frame_COP(X,external_forces_pred,Prediction,Fmax,i);
%     if AnalysisParameters.Prediction.ContactDetection == 3
%         %Coord = num2cell(good_X(:,i));
%         [external_forces_pred(i).Num_Markers] = Num{i}; % Abscisses des marqueurs en contact par rapport au marqueur origine de la structure
%         [external_forces_pred(i).p1] = Human_model(22).p;
%         [external_forces_pred(i).p2] = Human_model(28).p;
%     end
    
    waitbar(i/nbframe)
end

close(h)
disp(['... External Forces Prediction (' filename ') done'])

%% Filtrage des données

if AnalysisParameters.Prediction.FilterActive
    f_cut = AnalysisParameters.Prediction.FilterCutOff;
    % Conversion sous la forme d'une matrice (conversion into a matrix)
    for i=1:numel(external_forces_pred)
        for j=1:numel(external_forces_pred(i).fext)
            PredictionFx(i,j) = external_forces_pred(i).fext(j).fext(1,1); %#ok<AGROW>
            PredictionFy(i,j) = external_forces_pred(i).fext(j).fext(2,1); %#ok<AGROW>
            PredictionFz(i,j) = external_forces_pred(i).fext(j).fext(3,1); %#ok<AGROW>
            PredictionMx(i,j) = external_forces_pred(i).fext(j).fext(1,2); %#ok<AGROW>
            PredictionMy(i,j) = external_forces_pred(i).fext(j).fext(2,2); %#ok<AGROW> 
            PredictionMz(i,j) = external_forces_pred(i).fext(j).fext(3,2); %#ok<AGROW>
        end
    end
    % Filtrage
    PredictionFiltFx = filt_data(PredictionFx,f_cut,freq);
    PredictionFiltFy = filt_data(PredictionFy,f_cut,freq);
    PredictionFiltFz = filt_data(PredictionFz,f_cut,freq);
    PredictionFiltMx = filt_data(PredictionMx,f_cut,freq);
    PredictionFiltMy = filt_data(PredictionMy,f_cut,freq);
    PredictionFiltMz = filt_data(PredictionMz,f_cut,freq);
    % Remise sous la forme d'une structure (utilisée pour la dynamique inverse) (definition of a structure used for inverse dynamics)
    for i=1:numel(external_forces_pred)
        for j=1:numel(external_forces_pred(i).fext)
            external_forces_pred(i).fext(j).fext(1,1)=PredictionFiltFx(i,j);
            external_forces_pred(i).fext(j).fext(2,1)=PredictionFiltFy(i,j);
            external_forces_pred(i).fext(j).fext(3,1)=PredictionFiltFz(i,j);
            external_forces_pred(i).fext(j).fext(1,2)=PredictionFiltMx(i,j);
            external_forces_pred(i).fext(j).fext(2,2)=PredictionFiltMy(i,j);
            external_forces_pred(i).fext(j).fext(3,2)=PredictionFiltMz(i,j);
        end
    end
end

%% Pour animation (for animation purpose)

if ~any(strcmp('Visual',fieldnames(external_forces_pred)))
    external_forces_pred(1).Visual=[];
end
if ~isequal(AnalysisParameters.General.InputData, @MVNX_V3)
    for f=1:numel(external_forces_pred) % for every frame
%         % One global force
%             T = zeros(3,2);
%             for i=unique([Prediction.num_solid]) % for every solid
%                 T = T + external_forces_pred(f).fext(i).fext;
%             end
%             % CoP position
%             CoP = cross(T(:,1),T(:,2))/(norm(T(:,1))^2);
%             CoP = CoP - (CoP(3)/T(3,1))*T(:,1); % point on z=0
%             % external_forces structure
%             external_forces_pred(f).Visual = [external_forces_pred(f).Visual [CoP;T(:,1)]];
        % One force for each solid
            for i=unique([Prediction.num_solid]) % for every solid
                T = external_forces_pred(f).fext(i).fext;
                % CoP position
                pz_pieds = reshape([Prediction.pz], [numel([Prediction.pz])/length(Prediction) length(Prediction)]);
                CoP = cross(T(:,1),T(:,2))/(norm(T(:,1))^2);
                %CoP = CoP - (CoP(3)/T(3,1))*T(:,1); % point on z=0
                CoP = CoP - ((CoP(3)-mean(pz_pieds(f,:)))/T(3,1))*T(:,1); % point at average altitude of prediction points
                % external_forces structure
                external_forces_pred(f).Visual = [external_forces_pred(f).Visual [CoP;T(:,1)]];
            end
    end
else
    for f=1:numel(external_forces_pred) % for every frame
    % One force for each solid
        for i=unique([Prediction.num_solid]) % for every solid
            T = external_forces_pred(f).fext(i).fext;
            % CoP position
            CoP = cross(T(:,1),T(:,2))/(norm(T(:,1))^2);
            %CoP = CoP - (CoP(3)/T(3,1))*T(:,1); % point on z=0
            % external_forces structure
            external_forces_pred(f).Visual = [external_forces_pred(f).Visual [CoP;T(:,1)]];
        end
%     % One force for each foot
%         % Right foot (solids 52 and 55)
%             T = external_forces_pred(f).fext(52).fext + external_forces_pred(f).fext(55).fext;
%             CoP = cross(T(:,1),T(:,2))/(norm(T(:,1))^2);
%             CoP = CoP - (CoP(3)/T(3,1))*T(:,1); 
%             external_forces_pred(f).Visual = [external_forces_pred(f).Visual [CoP;T(:,1)]];
%         % Left foot (solids 64 and 67)
%             T = external_forces_pred(f).fext(64).fext + external_forces_pred(f).fext(67).fext;
%             CoP = cross(T(:,1),T(:,2))/(norm(T(:,1))^2);
%             CoP = CoP - (CoP(3)/T(3,1))*T(:,1); 
%             external_forces_pred(f).Visual = [external_forces_pred(f).Visual [CoP;T(:,1)]];
%     % One global force
%             T = external_forces_pred(f).fext(52).fext + external_forces_pred(f).fext(55).fext + ...
%                 external_forces_pred(f).fext(64).fext + external_forces_pred(f).fext(67).fext;
%             CoP = cross(T(:,1),T(:,2))/(norm(T(:,1))^2);
%             CoP = CoP - (CoP(3)/T(3,1))*T(:,1); 
%             external_forces_pred(f).Visual = [external_forces_pred(f).Visual [CoP;T(:,1)]];
    end
end

%% Sauvegarde des données (data saving)

if exist([filename '/ExternalForcesComputationResults.mat'],'file')
    load([filename '/ExternalForcesComputationResults.mat']); %#ok<LOAD>
end
ExternalForcesComputationResults.ExternalForcesPredictionCOP = external_forces_pred;

end