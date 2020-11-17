function fvalKt=funKtmax(A,BiomechanicalModel,MuscleConcerned,SolidConcerned,q,Fext,Fa,Fp,effector)
FMT=Fa.*A+Fp;
Kt=TaskStiffness(BiomechanicalModel,MuscleConcerned,SolidConcerned,q,Fext,FMT,effector);
[V,D] = eig(Kt);
fvalKt=-norm(sum(V.*diag(D)));
end