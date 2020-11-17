function dJdq = diff2dXdq(effector, SolidConcerned, BiomechanicalModel, q, dp)

n=length(q);
dJdq = zeros(n,3,n);
for i = SolidConcerned
    qplus = q;
    qmoins = q;
    qplus(i) = q(i)+dp;
    qmoins(i) = q(i)-dp;
    Jplus = diffdXdq(effector, SolidConcerned, BiomechanicalModel, qplus, dp);
    Jmoins = diffdXdq(effector, SolidConcerned, BiomechanicalModel, qmoins, dp);
    dJdq(:,:,i) = (Jplus -Jmoins)'/(2*dp);
end

end