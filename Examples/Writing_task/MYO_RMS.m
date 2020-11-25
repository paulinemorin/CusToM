function [RMS_EMG] = MYO_RMS(EMG_raw,s)
%Pour plusieurs vecteurs colonnes d'EMG, regroupés dans une matrice
% s est donnée en nombre de frames, à adapter en fonction de la fréquence
% d'acquisition des EMG ~100ms pour membre sup, soit s=200 pour 2000Hz

m=numel(EMG_raw(1,:));
n=numel(EMG_raw(:,1));
RMS_EMG=zeros(n,m);




    % .................................................STEP 2 REMOVE DC OFFSET
    %EMG_trend = detrend(EMG_rawntch);
    EMG_trend = detrend(EMG_raw);
    
    %..................................................STEP 4 ECG FILTER  

%         y=filterECG(EMG_trend,fs_EMG,'ecg',ECG,...
%                     'approach','ica','order',0,...
%                     'rtop','off','yule','off','comb','off',...
%                     'verbose','on','debug','on');
%                 le{n}='ica(0)+ref';
                
for i=1:n
    
    if(i<=s/2)
        RMS_EMG(i,1)=rms(EMG_trend(1:i+s/2,1));
    elseif(i>n-s/2)
        RMS_EMG(i,1)=rms(EMG_trend(i-s/2:n,1));
    else
        RMS_EMG(i,1)=rms(EMG_trend(i-s/2:i+s/2,1));
        
    end
end
%   RMS_EMG(i,j)=sqrt((1/s)*sum(EMG_raw(i-s/2:i+s/2,j).^2));
% 
% n=10;
% s=6;
% a=ones(10,1);
% for i=1:n
%     i
%     if(i<=s/2)
%         length(a(1:i+s/2,1))
%     elseif((n-i)<s/2)
%         length(a(i-(s/2):n,1))
%     else
%         length(a(i-(s/2):i+(s/2),1))
%         
%     end
% end

          
%     end
%          

    
end
