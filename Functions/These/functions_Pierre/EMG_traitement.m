H = waitbar(0,'Elaborating data...Please wait!');
disp(' ')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('             EMG PROCESSING                    ')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
%Data needed:
%foldersPath,trialsName,trialsList,AnalogRawData,AnalogFrameRate,EMGLabels,
%AnalysisWindow,EMGOffset,MaxEmgTrialsList,EMGsSelected_C3DLabels,
%EMGsSelected_OutputLabels

%Ri-nomination from parameters

EMGsSelected_OutputLabels= parameters.EMGsSelected.OutputLabels;
EMGsSelected_C3DLabels= parameters.EMGsSelected.C3DLabels;
EMGOffset=parameters.EMGOffset;
MaxEmgTrialsList=parameters.MaxEmgTrialsList;

if isfield(parameters,'OutputFileFormats')
    EMGOFileFormat=parameters.OutputFileFormats.EMG;
else
    EMGOFileFormat='.mot';  %default EMG output file format
end

foldersPath.maxemg=[foldersPath.elaboration filesep 'maxemg'];
mkdir(foldersPath.maxemg)
%Loading Analog Raw Data from the choosen trials with the corresponding
%labels
[AnalogRawData, AnalogDataLabels, aFrames, aUnits]=loadMatData(foldersPath.sessionData, trialsList, 'AnalogData');

%Loading Analog Raw Data for EMG Max Computation from the trials list
if isequal(parameters.MaxEmgTrialsList,parameters.trialsList)
    AnalogRawForMax=AnalogRawData;
    AnalogLabelsForMax=AnalogDataLabels;
else
    [AnalogRawForMax, AnalogLabelsForMax]=loadMatData(foldersPath.sessionData, MaxEmgTrialsList, 'AnalogData');
end

%Loading Analog Data Labels
%load([foldersPath.sessionData 'AnalogDataLabels.mat'])
%NOTE: analog channels configuration may change according to the
%acquisition procedure (e.g. with Vicon), thus analog labels for each
%trial are loaded and used (as for markers)

%If there are EMGs --> processing
if (isempty(AnalogRawData)==0 && isempty(AnalogDataLabels)==0)
    %% --------------------------------------------------------------------
    %                   EMGs EXTRACTION and MUSCLES SELECTION
    %                   EMGs Arrangement for the Output file
    %----------------------------------------------------------------------
    
    for k=1:length(trialsList)
        
        EMGselectionIndexes{k}=findIndexes(AnalogDataLabels{k},EMGsSelected_C3DLabels);
        EMGsSelected{k}=AnalogRawData{k}(:,EMGselectionIndexes{k});
        EMGsUnits{k}=aUnits{k}(EMGselectionIndexes{k});
    end
    
    %The arrangement of EMG signals in the analog channels may change
    %among trials, thus EMGselectionIndexes may differ according to the
    %trials used for max computation
    
    for k=1:length(MaxEmgTrialsList)
        
        EMGselectionIndexesForMax{k}=findIndexes(AnalogLabelsForMax{k},EMGsSelected_C3DLabels);
        EMGsSelectedForMax{k}=AnalogRawForMax{k}(:,EMGselectionIndexesForMax{k});
    end
    
    %% ------------------------------------------------------------------------
    %                       EMG FILTERING: ENVELOPE
    %--------------------------------------------------------------------------
    %fcut for EMG assumed fixed (6Hz)
    EMGsEnvelope=EMGFiltering(EMGsSelected,AnalogFrameRate);
    EMGsEnvelopeForMax=EMGFiltering(EMGsSelectedForMax,AnalogFrameRate);
    
    %% ------------------------------------------------------------------------
    %                      EMG ANALYSIS WINDOW SELECTION
    %--------------------------------------------------------------------------
    
    [EMGsFiltered,EMGtime]=selectionData(EMGsEnvelope,AnalysisWindow,AnalogFrameRate,EMGOffset);
    
    %if trials for max computation are the same of those for elaboration, max
    %values are computed within the same analysis window, else all signals are
    %considered --> this is not the way!
    %if isequal(MaxEmgTrialsList,trialsList)
    %    EMGsForMax=selectionData(EMGsEnvelopeForMax,AnalysisWindow,AnalogFrameRate,EMGOffset);
    %else
    %TO DO: implement a way to select a different AnalysisWindow for
    %max EMG computation
    %The analysis window for max identification is the whole trial now:
    EMGsForMax=EMGsEnvelopeForMax;
    %end
    
    %% ------------------------------------------------------------------------
    %                        COMPUTE MAX EMG VALUES
    %--------------------------------------------------------------------------
    [MaxEMG_aframes, numMaxEMG_trials,MaxEMGvalues]=computeMaxEMGvalues(EMGsForMax);
    
    disp('Max values for selected emg signals have been computed')
    sMaxEMG_trials=MaxEmgTrialsList(numMaxEMG_trials);
    MaxEMG_time=MaxEMG_aframes/AnalogFrameRate;
    
    %print maxemg.txt
    printMaxEMGvalues(foldersPath.maxemg, EMGsSelected_C3DLabels, MaxEMGvalues, sMaxEMG_trials, MaxEMG_time);
    
    disp('Printed maxemg.txt')
    
    waitbar(6/7);
    
    %% ------------------------------------------------------------------------
    %                            NORMALIZE EMG
    %--------------------------------------------------------------------------
    NormEMG=normalizeEMG(EMGsFiltered,MaxEMGvalues);
    
    %% ------------------------------------------------------------------------
    %                          SAVING and PLOTTING
    %--------------------------------------------------------------------------
    
    %Manual method for the Windows Selection has no window offset
    EnvelopePlotting(EMGsFiltered,MaxEMGvalues,EMGsSelected_C3DLabels, EMGsUnits, foldersPath.trialOutput, AnalogFrameRate,EMGOffset)
    
    %storing all info related to max EMGs in a struct
    MaxEMGstruct.values=MaxEMGvalues;
    MaxEMGstruct.muscles=EMGsSelected_C3DLabels;
    MaxEMGstruct.aframes=MaxEMG_aframes;
    MaxEMGstruct.time=MaxEMG_time;
    MaxEMGstruct.trials=numMaxEMG_trials;
    MaxEMGstruct.trialNames=sMaxEMG_trials;
    
    maxEmgPlotting(EMGsSelectedForMax,EMGsForMax,EMGsUnits{1}{1},foldersPath.maxemg,AnalogFrameRate, MaxEMGstruct)
    
    % ------------------------------------------------------------------------
    %                            PRINT emg
    %--------------------------------------------------------------------------
    
    switch EMGOFileFormat
        
        case '.mat'
            
            for k=1:length(trialsList)
                
                printEMGtxt(foldersPath.trialOutput{k},EMGtime{k},NormEMG{k},EMGsSelected_OutputLabels);
            end
            
        case {'.sto','.mot'}
            
            for k=1:length(trialsList)
                
                printEMGmot(foldersPath.trialOutput{k},EMGtime{k},NormEMG{k},EMGsSelected_OutputLabels, EMGOFileFormat);
            end
            
        case '.mat'
            
            for k=1:length(trialsList)
                printEMGmat(foldersPath.trialOutput{k},EMGtime{k},NormEMG{k},EMGsSelected_OutputLabels, EMGOFileFormat);
            end
            
        otherwise
            error('ErrorTests:convertTest', ...
                ['----------------------------------------------------------------\nWARNING: EMG Output File Format not Available!\nChoose among: [' availableFileFormats ']. Please, check it in your elaboration.xml file'])
    end
    
    disp(['Printed emg' EMGOFileFormat])
    
    waitbar(7/7);
    close(H)
    
    % -------------------------------------------------------------------------
    %                           PLOTTING EMG
    %--------------------------------------------------------------------------
%     plotEMGChoice = questdlg('Do you want to plot EMGs Raw', ...
%         'Plotting EMGs', ...
%         'Yes','No','Yes');
    
%     if strcmp(plotEMGChoice,'Yes')
%         
        EMGsPlotting(EMGsSelected,EMGsEnvelope,AnalysisWindow,EMGsSelected_C3DLabels,EMGsUnits,foldersPath.trialOutput,AnalogFrameRate)
        disp('Plotted EMGs')
%     end
    
else
    waitbar(6/7);
    disp('Check your data and/or your configuration files: No EMG raw data to be processed')
    waitbar(7/7);
    close(h)
end
%% -------------------------------------------------------------------------

h = msgbox('Data Processing terminated successfully','Done!');
uiwait(h)

save_to_base(1)