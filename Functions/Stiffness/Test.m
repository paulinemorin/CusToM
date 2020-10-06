MuscleConcerned = 1:41;
SolidConcerned = [];
for i=MuscleConcerned
   SolidConcerned = [SolidConcerned BiomechanicalModel.Muscles(i).num_solid'];
end
SolidConcerned = unique(SolidConcerned);
q = InverseKinematicsResults.JointCoordinates(:,1); %Joint coordinates of the first frame of cycling example
Fext = ExternalForcesComputationResults;
FMT = MuscleForcesComputationResults.MuscleForces(:,1); 
effector = [22, 3]; %RFOOT : solid_number 22 and marker anat_position : RTOE : 3
Fext = ExternalForcesComputationResults.ExternalForcesExperiments(1).fext(22);
Fext = Fext.fext(1:3,1); %external forces applied to the RFOOt at first frame
FMT = FMT*0;
Fext = [3;1;0];

Kt = TaskStiffness(BiomechanicalModel,MuscleConcerned,SolidConcerned,q,Fext,FMT,effector);