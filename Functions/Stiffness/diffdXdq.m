function J = diffdXdq(effector, SolidConcerned, BiomechanicalModel, q, dq)
%This function returns the Jacobian of the SolidConcerned
%Based on:
% Référence Stanev : Stanev, D., & Moustakas, K. (2019). Stiffness modulation of redundant musculoskeletal systems. Journal of Biomechanics, 85, 101–107. https://doi.org/10.1016/j.jbiomech.2019.01.017
%INPUT :
%- effector : vector with the number of the solid effector and its number in anatomical position
%- SolidConcerned : list of number of concerned solids
%- BiomechanicalModel : musculoskeletal model
%- q : vector of joint coordinates 
%- dq : vector of joint speeds
%Detailed explanation goes here

J = zeros(3,length(q));

for i = SolidConcerned
    BiomechanicalModel.OsteoArticularModel(i).q = q(i)+dq;    
    Human_model = ForwardPositions(BiomechanicalModel.OsteoArticularModel,1);
    Pplus = Human_model(effector(1)).p+Human_model(effector(1)).R*Human_model(effector(1)).anat_position{effector(2),2};
    
    BiomechanicalModel.OsteoArticularModel(i).q = q(i)-dq;
    Human_model = ForwardPositions(BiomechanicalModel.OsteoArticularModel,1);
    Pmoins = Human_model(effector(1)).p+Human_model(effector(1)).R*Human_model(effector(1)).anat_position{effector(2),2};
    
    BiomechanicalModel.OsteoArticularModel(i).q = q(i);
    J(:,i) = (Pplus-Pmoins)'/(2*dq);
end