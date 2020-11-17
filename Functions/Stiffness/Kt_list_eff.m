function [Kt]=Kt_list_eff(BiomechanicalModel,MuscleConcerned,SolidConcerned,q,Fext,Fa,A,Fp,effector)
Kt = [];
for i_eff=1:numel(effector(:,1))
    Kt=[Kt TaskStiffness(BiomechanicalModel,MuscleConcerned(i_eff).list,SolidConcerned(i_eff).list,q,Fext,Fa.*A+Fp,effector(i_eff,:)];
end