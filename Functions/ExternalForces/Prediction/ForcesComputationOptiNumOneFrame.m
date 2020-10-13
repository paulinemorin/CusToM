function [Fopt,Fmax] = ForcesComputationOptiNumOneFrame(BiomechanicalModel)
% Computation of the muscle forces estimation step by using an optimization method
%
%	Based on :
%	- Crowninshield, R. D., 1978.
%	Use of optimization techniques to predict muscle forces. Journal of Biomechanical Engineering, 100(2), 88-92.
%
%   INPUT

%   OUTPUT
%   - MuscleForcesComputationResults: results of the muscle forces
%   estimation step (see the Documentation for the structure)
%________________________________________________________
%
% Licence
% Toolbox distributed under GPL 3.0 Licence
%________________________________________________________
%
% Authors : Antoine Muller, Charles Pontonnier, Pierre Puchaud and
% Georges Dumont
%________________________________________________________

Muscles = BiomechanicalModel.Muscles;
HumanModel =  BiomechanicalModel.OsteoArticularModel;
q=[HumanModel .q];
torques =[HumanModel.torques];


Nb_q=size(q,1);

%existing muscles
idm = logical([Muscles.exist]);
Nb_muscles=numel(Muscles(idm));

%% computation of muscle moment arms from joint posture
% L0=zeros(Nb_muscles,1);
% Ls=zeros(Nb_muscles,1);
% for i=1:Nb_muscles
%     L0(i) = BiomechanicalModel.Muscles(i).l0;
%     Ls(i) = BiomechanicalModel.Muscles(i).ls;
% end

%Lmt(idm)   =   MuscleLengthComputationNum(BiomechanicalModel,q(:,i)); %dependant of every q (q_complete)
R   =   MomentArmsComputationNum(BiomechanicalModel,q,0.0001); %depend on reduced set of q (q_red)


%Lm = Lmt./(Ls./L0+1);
% Muscle length ratio to optimal length
%Lm_norm = Lm./L0;
% Muscle velocity
% Vm = gradient(Lm_norm)*freq;

idxj=find(sum(R(:,:,1),2)~=0)';

%% Computation of muscle forces (optimization)
% Optimisation parameters
Amin = zeros(Nb_muscles,1);
A0  = 0.5*ones(Nb_muscles,1);
Fmax = [Muscles(idm).f0]';
Amax = ones(Nb_muscles,1);
% Muscle Forces Matrices computation
%[Fa,Fp]=AnalysisParameters.Muscles.MuscleModel(Lm,Vm,Fmax);
% Solver parameters
options1 = optimoptions(@fmincon,'Algorithm','sqp','Display','off','GradObj','off','GradConstr','off','TolFun',1e-6,'MaxIterations',100000,'MaxFunEvals',100000);


if isfield(BiomechanicalModel.OsteoArticularModel,'ClosedLoop') && ~isempty([BiomechanicalModel.OsteoArticularModel.ClosedLoop])    
    k=BiomechanicalModel.GeometricalCalibration.k_calib;
        
    [solid_path1,solid_path2,num_solid,num_markers]=Data_ClosedLoop(BiomechanicalModel.OsteoArticularModel);

    dependancies=KinematicDependancy(BiomechanicalModel.OsteoArticularModel);
    % Closed-loop constraints
    KT=ConstraintsJacobian(BiomechanicalModel,q(:,1),solid_path1,solid_path2,num_solid,num_markers,k,0.0001,dependancies)';
    [idKT,~]=find(sum(KT(:,:,1),2)~=0);
    idq=intersect(idKT,idxj);
    % Adaptation of variables to closed-loop problem
    A0 = [A0 ; zeros(size(KT,2),1)];
    Amin = [Amin ;-inf*ones(size(KT,2),1)];
    Fmax = [Fmax ;inf*ones(size(KT,2),1)];
    Amax = [Amax ;inf*ones(size(KT,2),1)];
    % Moment arms and Active forces
    Aeq = [R(idq,:).*Fmax' KT(idq,:)];
    % Joint Torques
    beq = torques(idq);% - R(idq,:,1)*Fp(:,1);
    % First frame optimization
    [Aopt] = PolynomialFunction(A0, Aeq, beq, Amin, Amax, options1, 2, Fmax, Fmax);
    % Muscular activiy
%    Fopt = Fa.*Aopt(1:Nb_muscles,1)+Fp;
    Fopt = Fmax.*Aopt(1:Nb_muscles);
    
else
    % Moment arms and Active forces
    Aeq=R(idxj,:).*Fmax';
    % Joint Torques and Passive force
    beq=torques(idxj);% - R(idxj,:)*Fp;
    % First frame optimization
    [Aopt,exitflag] = PolynomialFunction(A0, Aeq, beq, Amin, Amax, options1,2, Fmax, Fmax);
    % Muscular activity
  %  Fopt= Fa.*Aopt+Fp;
       Fopt = Fmax.*Aopt(1:Nb_muscles);
       if exitflag==-2
           Fopt=Fopt*1000;
       end

    
end



end