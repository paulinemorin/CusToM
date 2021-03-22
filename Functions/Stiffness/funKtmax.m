function fvalKt=funKtmax(A,BiomechanicalModel,MuscleConcerned,Fext,Fa,Fp,R,dRdq,J,dJdq)
FMT=Fa.*A(1:length(Fa))+Fp;
Kt=TaskStiffness(BiomechanicalModel,MuscleConcerned,Fext,FMT,R,dRdq,J,dJdq);
[V,D] = eig(Kt);
fvalKt=-norm(sum(V.*diag(D)));
end