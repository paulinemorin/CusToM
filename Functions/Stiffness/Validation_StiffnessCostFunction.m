filename='THI_trajectoire0008';
load([filename '/InverseDynamicsResults']);
torques =InverseDynamicsResults.JointTorques;
Nb_frames=1/0.005; %size(torques,2);
Muscles = BiomechanicalModel.Muscles;
idm = logical([Muscles.exist]);
Nb_muscles=numel(Muscles(idm));
alpha_l=0.6:0.1:1;
A_avg=zeros(Nb_muscles,numel(alpha_l));
i=1;
for alpha=alpha_l
    AnalysisParameters.StiffnessPercent=alpha;
    ForcesComputationOptiNum(filename,BiomechanicalModel, AnalysisParameters);
    save(['MuscleForcesComputationResults_', num2str(alpha),'.mat'], 'MuscleForcesComputationResults')
    A_avg(:,i) = sum(MuscleForcesComputationResults.MuscleActivations,2)/Nb_frames;
    i=i+1;
end
figure
hold on
grid on
plot(0.6:0.1:1,A_avg)
xlabel('Pourcentage de raideur','FontSize',16);
ylabel('Moyenne temporelle de la raideur','FontSize',16);
legend(BiomechanicalModel.Osteoarticular);