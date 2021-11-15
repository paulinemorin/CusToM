function [] = CalibrateModelGenerationNum(ModelParameters,AnalysisParameters)
% Generation of the calibrated musculoskeletal model
%   A subject-specific model is generated. Additionnal computations used
%   for the analysis are also made.
%
%	Associated publication:
%	- Muller, A., Haering, D., Pontonnier, C., & Dumont, G., 2017. 
%	Non-invasive techniques for musculoskeletal model calibration. In 23�me Congr�s Fran�ais de M�canique.
%
%   INPUT
%   - ModelParameters: parameters of the musculoskeletal model,
%   automatically generated by the graphic interface 'GenerateParameters';
%   - AnalysisParameters: parameters of the musculoskeletal analysis,
%   automatically generated by the graphic interface 'Analysis'.
%   OUTPUT
%   The musculoskeletal model is automatically saved in the variable
%   'BiomechanicalModel'.
%   In this numerical version, the moment arms are no longer computed
%   analytically. The MusIC method is no longer activated
%________________________________________________________
%
% Licence
% Toolbox distributed under Licence GPL 3.0
%________________________________________________________
%
% Authors : Antoine Muller, Charles Pontonnier, Pierre Puchaud and
% Georges Dumont
%________________________________________________________

%% Model generation from a MVNX file
if isequal(AnalysisParameters.General.InputData, @MVNX_V3)
    disp('Model Generation ...')
    MVNXModelGeneration(ModelParameters, AnalysisParameters);
    disp('... Model Generation done')
    return;
end

%% Anthropometric model generation
disp('Anthropometric Model Generation ...')
[BiomechanicalModel.OsteoArticularModel, BiomechanicalModel.Markers, BiomechanicalModel.Muscles] = ModelGeneration(ModelParameters);
disp('... Anthropometric Model Generation done')

%% Geometrical calibration
if AnalysisParameters.CalibIK.Active
    disp('Geometrical Calibration ...')
    [BiomechanicalModel.OsteoArticularModel, BiomechanicalModel.GeometricalCalibration, BiomechanicalModel.Markers] = GeometricalCalibration(BiomechanicalModel.OsteoArticularModel, BiomechanicalModel.Markers, AnalysisParameters);
    disp('... Geometrical Calibration done')
end

%% Joint Calibration
if ~isempty(find(contains(...
        {BiomechanicalModel.OsteoArticularModel.name}','Patella'), 1))
    disp('Joint Calibration ...')
    [BiomechanicalModel]=CalibratePatellaJoint(BiomechanicalModel);
    disp('... Joint Calibration done')
end

%% Symbolic functions
disp('Preliminary Computations ...')
  [BiomechanicalModel.OsteoArticularModel] = Add6dof(BiomechanicalModel.OsteoArticularModel);
[BiomechanicalModel.OsteoArticularModel, ...
    BiomechanicalModel.Generalized_Coordinates] = SymbolicFunctionGenerationIK(BiomechanicalModel.OsteoArticularModel,BiomechanicalModel.Markers);
disp('... Preliminary Computations done')

%% Inertial calibration
if AnalysisParameters.CalibID.Active
    disp('Dynamic Calibration ...')
    [BiomechanicalModel] = DynamicCalibration(ModelParameters,AnalysisParameters, BiomechanicalModel);
    disp('... Dynamic Calibration done')
end


%% Moment arms matrix et muscular coupling
% not applicable here
% if numel(BiomechanicalModel.Muscles)
%     disp('Moment Arms Computation ...')
%     [BiomechanicalModel.MomentArms,BiomechanicalModel.MuscularCoupling] =...
%         MomentArmsComputation(BiomechanicalModel);
%     disp('... Moment Arms Computation done');
% end

%% Generation of the data base for computation of muscular forces by using MusIC method
% Not available here
% if numel(BiomechanicalModel.Muscles) && AnalysisParameters.Muscles.Method == 2
%     disp('MusIC Database Generation ...')
%     eval(['[BiomechanicalModel.MusICDatabase_' char(AnalysisParameters.Muscles.Costfunction) num2str(AnalysisParameters.Muscles.CostfunctionOptions) '_' num2str(AnalysisParameters.Muscles.DatabaseDensity(1)) '_' num2str(AnalysisParameters.Muscles.DatabaseDensity(2)) '] = MusICDatabaseGeneration(BiomechanicalModel, AnalysisParameters);']) 
%     disp('... MusIC Database Generation done')
% end

 %% Wrapping locations
 if ~isempty([BiomechanicalModel.Muscles])
    if ~isempty([BiomechanicalModel.Muscles.wrap]')
         [BiomechanicalModel] =...
             WrappingLocations(BiomechanicalModel);
    end
 end
 
 %% Muscular coupling computation
q=zeros(numel(BiomechanicalModel.OsteoArticularModel(:))-6,1)+0.01;
dp=0.001;


 if numel(BiomechanicalModel.Muscles)
     disp('Muscular Coupling Computation ...')
     [BiomechanicalModel.Coupling] =...
         MomentArmsComputationInit(BiomechanicalModel,q,dp);
     disp('... Muscular Coupling Computation done');
 end

 %% Muscular calibration
 if ~isempty([BiomechanicalModel.Muscles])
    [BiomechanicalModel]=MuscleCalibrationAnthropo(ModelParameters,AnalysisParameters,BiomechanicalModel);
 end
 
 
save('BiomechanicalModel','BiomechanicalModel');


%% Closedloop equations
% disp('Closed loop equations ...')
% BiomechanicalModel = AddClosedLoopEquations(BiomechanicalModel);
% disp('... Closed loop equations  done')

save('BiomechanicalModel','BiomechanicalModel');

end
