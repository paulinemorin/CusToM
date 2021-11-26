function [start,stop] = croper(ExternalForcesComputationResults,AnalysisParameters, filename)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

start=1;
stop=length(ExternalForcesComputationResults.ExternalForcesExperiments);

if contains(filename,'Run')||contains(filename,'Direction')
        
    AnalysisParameters.ExternalForces.Options(1,1)={'LFoot'};
    AnalysisParameters.ExternalForces.Options(2,1)={'RFoot'};
    
    FootOnPlates=AnalysisParameters.ExternalForces.Options;
    [Event,ID]=DetectRunPhases([char(filename) '.c3d'],FootOnPlates);
    
    start_1=ID(1);
    stop_1=ID(2);
    
    AnalysisParameters.ExternalForces.Options(1,1)={'RFoot'};
    AnalysisParameters.ExternalForces.Options(2,1)={'LFoot'};
    FootOnPlates=AnalysisParameters.ExternalForces.Options;
    [Event,ID]=DetectRunPhases([char(filename) '.c3d'],FootOnPlates);
    
    start_2=ID(1);
    stop_2=ID(2);
    
    start=min(start_1,start_2);
    stop=max(stop_1,stop_2);
    
elseif contains(filename,'Marche') || contains(filename,'Walk')% A faire pour les autres taches.
    AnalysisParameters.ExternalForces.Options(1,1)={'LFoot'};
    AnalysisParameters.ExternalForces.Options(2,1)={'RFoot'};

    FootOnPlates=AnalysisParameters.ExternalForces.Options;
    [Event,ID]=DetectOnlyOneFootPhase([char(filename) '.c3d'],FootOnPlates);
  
    start_1=ID(1);
    stop_1=ID(2);
    
    AnalysisParameters.ExternalForces.Options(1,1)={'RFoot'};
    AnalysisParameters.ExternalForces.Options(2,1)={'LFoot'};
    
    [Event,ID]=DetectOnlyOneFootPhase([char(filename) '.c3d'],FootOnPlates);
  
    start_2=ID(1);
    stop_2=ID(2);
    
    start=min(start_1,start_2);
    stop=max(stop_1,stop_2);
    
    if start == 0
        start=1;
    end
end
end

