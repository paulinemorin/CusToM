filename='THI_trajectoire0008';
load([filename '/InverseDynamicsResults']);
load('BiomechanicalModel.mat');
load('AnalysisParameters.mat');
torques =InverseDynamicsResults.JointTorques;
Muscles = BiomechanicalModel.Muscles;
idm = logical([Muscles.exist]);
Nb_muscles=numel(Muscles(idm));
alpha_l=0:0.1:1;
A_avg=zeros(Nb_muscles,numel(alpha_l));
i=1;
for alpha=alpha_l
    alpha
    AnalysisParameters.StiffnessPercent=alpha;
    MuscleForcesComputationResults = ForcesComputationOptiNum(filename,BiomechanicalModel, AnalysisParameters);
    save(['MuscleForcesComputationResults_', num2str(alpha),'.mat'], 'MuscleForcesComputationResults')
    A_avg(:,i) = mean(MuscleForcesComputationResults.MuscleActivations,2);
    i=i+1;
end
figure
hold on
grid on
plot(0:0.1:1,A_avg)
xlabel('Pourcentage de raideur','FontSize',16);
ylabel('Moyenne temporelle des activations','FontSize',16);
legend(BiomechanicalModel.Osteoarticular);