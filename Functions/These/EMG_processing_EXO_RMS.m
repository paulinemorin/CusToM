function [RMS_EMG] = EMG_processing_EXO_RMS(EMG_raw,S_window)
    
m=numel(EMG_raw(:,1));
n=numel(EMG_raw(1,:));
RMS_EMG=zeros(m,n);

s=S_window; % s : nb de frame de la fenetre glissante -> 200ms à 200Hz 40 frames
    % .................................................STEP 1 REMOVE DC OFFSET 
    EMG_trend = EMG_raw;
    y=  EMG_trend ;    
    % .................................................STEP 2 WINDOWED RMS
    for i=1:m
        
        for j=1:n
            if i+s>m
                
                delta=i+s-m;
                RMS_EMG(i,j)=mean(y((i-delta):(i+s-delta),j));
                %RMS_EMG(i,j)=sqrt((1/s)*sum(y((i-delta):(i+s-delta),j).^2));
            else
                %RMS_EMG(i,j)=sqrt((1/s)*sum(y(i:i+s,j).^2));
                RMS_EMG(i,j)=mean((y(i:i+s,j)));
            end
        end
    end
    
    
    
end
