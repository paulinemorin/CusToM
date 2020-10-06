function Kt = TaskStiffness(BiomechanicalModel,MuscleConcerned, SolidConcerned,q,Fext,FMT,effector,varargin)
%This function returns the task stiffness 
%Based on:
% Référence Stanev : Stanev, D., & Moustakas, K. (2019). Stiffness modulation of redundant musculoskeletal systems. Journal of Biomechanics, 85, 101–107. https://doi.org/10.1016/j.jbiomech.2019.01.017
%INPUT :
%- BiomechanicalModel : musculoskeletal model
%- q : vector of joint coordinates 
%- Fext : vector of external forces
%- FMT : vector of muscle foces
%- effector : vector with the number of the solid effector and its number in anatomical position
%- varargin : supposed to take R when already calculate
%- MuscleConcerned : list of number of concerned muscles
%- SolidConcerned : list of number of concerned solids
%   Detailed explanation goes here

n=length(q);
Kj = JointStiffness(BiomechanicalModel,MuscleConcerned, SolidConcerned,q,FMT,varargin);
BiomechanicalModel.OsteoArticularModel(1).R = eye(3);
BiomechanicalModel.OsteoArticularModel(1).p = [0 0 0]';

for i=1:n
BiomechanicalModel.OsteoArticularModel(i).q = q(i);
end
dq=0.0001;

J = diffdXdq(effector, SolidConcerned, BiomechanicalModel, q, dq);
dJdq = zeros(n,3,n);
for i = SolidConcerned
    qplus = q;
    qmoins = q;
    qplus(i) = q(i)+dq;
    qmoins(i) = q(i)-dq;
    Jplus = diffdXdq(effector, SolidConcerned, BiomechanicalModel, qplus, dq);
    Jmoins = diffdXdq(effector, SolidConcerned, BiomechanicalModel, qmoins, dq);
    dJdq(:,:,i) = (Jplus -Jmoins)'/(2*dq);
end

for i=1:n
    Kt(:,i) = -dJdq(:,:,i)*Fext;
end
Kt = pinv(J')*(Kj-Kt)*pinv(J);
end