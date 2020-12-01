function [Kt]=Kt_list_eff(BiomechanicalModel,MuscleConcerned,Fext,Fa,A,Fp,R,dRdq,J,dJdq)
Kt = [];
for i_eff=1:size(J,2)
    Kt=[Kt -funKtmax(A,BiomechanicalModel,MuscleConcerned(i_eff).list,Fext,Fa,Fp,R,dRdq{i_eff},J{i_eff},dJdq{i_eff})];
end



end