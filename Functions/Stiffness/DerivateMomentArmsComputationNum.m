function  dRdq = DerivateMomentArmsComputationNum(BiomechanicalModel,q,dp,SolidConcerned)

dRdq = zeros(length(q),numel(BiomechanicalModel.Muscles),length(q));
for j = SolidConcerned
    qplus = q;
    qmoins = q;
    qplus(j) = q(j)+dp;
    qmoins(j) = q(j)-dp;
    Rplus = MomentArmsComputationNum(BiomechanicalModel,qplus,dp);
    Rmoins = MomentArmsComputationNum(BiomechanicalModel,qmoins,dp);
    dRdq(:,:,j) = (Rplus - Rmoins)/(2*dp);
end


end