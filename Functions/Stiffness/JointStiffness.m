function Kj = JointStiffness(BiomechanicalModel,MuscleConcerned,SolidConcerned,q,FMT,varargin)
%This functionreturns the task stiffness 
%Based on:
% Référence Stanev : Stanev, D., & Moustakas, K. (2019). Stiffness modulation of redundant musculoskeletal systems. Journal of Biomechanics, 85, 101–107. https://doi.org/10.1016/j.jbiomech.2019.01.017
%INPUT :
%- BiomechanicalModel : musculoskeletal model
%- MuscleConcerned : list of number of concerned muscles
%- SolidConcerned : list of number of concerned solids
%- q : vector of joint coordinates 
%- FMT : vector of muscle foces
%- varargin : supposed to take R when already calculated

km = zeros(length(FMT),1);
km(MuscleConcerned) = 23.4*FMT(MuscleConcerned)./[BiomechanicalModel.Muscles(MuscleConcerned).l0]';
Km = diag(km);
dq=0.0001;
if isempty(varargin{1})
    R  = MomentArmsComputationNum(BiomechanicalModel,q,dq)';
else
    R = varargin{1};
end

n = length(q);
m = length(BiomechanicalModel.Muscles);
dRdq = zeros(n,m,n);
for i = SolidConcerned
    qplus = q;
    qmoins = q;
    qplus(i) = q(i)+dq;
    qmoins(i) = q(i)-dq;
    Rplus = MomentArmsComputationNum(BiomechanicalModel,qplus,dq);
    Rmoins = MomentArmsComputationNum(BiomechanicalModel,qmoins,dq);
    dRdq(:,:,i) = (Rplus - Rmoins)/(2*dq);
end

Kj = zeros(n,n);
for i=1:n
    Kj(:,i) = -dRdq(:,:,i)*FMT;
end

Kj = Kj -R'*Km*R;
end