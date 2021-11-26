function []=EMG_traitementAll(meta)
addpath(genpath('C:\Users\ppuchaud\Documents\MATLAB\MOtoNMS'))

C3D2MAT()

h=btkReadAcquisition('ROM.c3d');
[hh,Anainfo]=btkGetAnalogs(h);
names=fieldnames(Anainfo.label);
Ind = cellfun(@(X)contains(X,'Voltage'),names);
f=names(Ind);
for ii=1:length(f)
    labels{ii}=Anainfo.label.(f{ii});
end

parameters.EMGsSelected.C3DLabels=labels;
parameters.EMGsSelected.OutputLabels=EMGlabelNames(labels,meta);

parameters.EMGOffset = 0;
parameters.OutputFileFormats.EMG='.mot';

tlist=cellfun(@(X)strrep(X,'.c3d',''),meta.TrialList_tasks,'UniformOutput',false);
parameters.MaxEmgTrialsList = tlist;
parameters.trialsList = tlist;

foldersPath.elaboration =fullfile(meta.path,'elaboration\');
mkdir(fullfile(meta.path,'elaboration'));
foldersPath.sessionData =fullfile(meta.path,'sessionData\');
mkdir(fullfile(meta.path,'sessionData'));

trialsList=parameters.trialsList;
foldersPath.trialOutput= mkOutputDir(foldersPath.elaboration,trialsList);

AnalogFrameRate=Anainfo.frequency;

[~,info]=btkGetMarkers(h);

for ii=1:length(meta.TrialList_tasks)
    AnalysisWindow{ii}.LabeledDataOffset= btkGetFirstFrame(h)-1;
    AnalysisWindow{ii}.startFrame=meta.FrameOfInterest(ii).ID(1);
    AnalysisWindow{ii}.endFrame=meta.FrameOfInterest(ii).ID(2);
    AnalysisWindow{ii}.rate=info.frequency;
end
EMGOFileFormat='.mat';

EMG_traitement();
end