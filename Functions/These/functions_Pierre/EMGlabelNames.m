function [labelsout]=EMGlabelNames(labels,meta)

% Initialisation
Muscles = struct('name',[],'f0',[],'l0',[],'Kt',[],'ls',[],'alpha0',[],'path',[],'wrap',[]);
Muscles(1) = [];

% add muscle sets
for i = 1:numel(meta.ModelParameters.Muscles)
    Muscles = meta.ModelParameters.Muscles{i}(Muscles,meta.ModelParameters.MusclesOptions{i});
end

for ii=1:length(labels)
    L=labels{ii};
    
    if contains(L,'.R') % Right
        S= 'R';
    elseif contains(L,'.L') % Left
        S='L';
    end
    
    if contains(L,'RF') % rectus femoris
        
        labelsout{ii}= [S 'RectusFemoris'];
        
    elseif contains(L,'BF') % biceps femoris
        
        labelsout{ii}= [S 'BicepsFemorisCaputLongum'];
        
    elseif contains(L,'TA') % tibialis anterior
        
        labelsout{ii}= [S 'TibialisAnterior'];
        
    elseif contains(L,'SOL') % Soleus
        
        labelsout{ii}= [S 'Soleus'];
        
    elseif contains(L,'GM') % Gluteus maximus
        
        labelsout{ii}= [S 'GluteusMaximus3'];
        
    elseif contains(L,'GasL') % Gastrocnemius Lateralis
        
        labelsout{ii}= [S 'Gastrocnemius'];
        
    elseif contains(L,'VL') % Vastus Lateralis
        
        labelsout{ii}= [S 'GluteusMaximus3'];
        
    elseif contains(L,'PL') % Peroneous Longus
        
        labelsout{ii}= [S 'PeroneusBrevis'];
    else
        warning('name of the emg not defined')
        keyboard
    end
end

end