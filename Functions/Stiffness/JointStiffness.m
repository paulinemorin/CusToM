function Kj = JointStiffness(BiomechanicalModel,MuscleConcerned,FMT,R,dRdq)
%This functionreturns the task stiffness 
%Based on:
% R�f�rence Stanev : Stanev, D., & Moustakas, K. (2019). Stiffness modulation of redundant musculoskeletal systems. Journal of Biomechanics, 85, 101�107. https://doi.org/10.1016/j.jbiomech.2019.01.017
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

n = size(R,1);
Kj = zeros(n,n);
for i=1:n
    Kj(:,i) = -dRdq(:,:,i)*FMT;
end

Kj = Kj -R*Km*R';
end