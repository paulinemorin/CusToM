function Kt = TaskStiffness(BiomechanicalModel,MuscleConcerned,Fext,FMT,R,dRdq,J,dJdq)
%This function returns the task stiffness
%Based on:
% R�f�rence Stanev : Stanev, D., & Moustakas, K. (2019). Stiffness modulation of redundant musculoskeletal systems. Journal of Biomechanics, 85, 101�107. https://doi.org/10.1016/j.jbiomech.2019.01.017
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

n=size(R,1);
Kj = JointStiffness(BiomechanicalModel,MuscleConcerned,FMT,R,dRdq);

for i=1:n
    Kt(:,i) = -dJdq(:,:,i)*Fext;
end

Kt = pinv(J')*(Kj-Kt)*pinv(J);

Kt = (Kt+ Kt')/2;
end