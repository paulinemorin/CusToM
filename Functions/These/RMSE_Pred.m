function RMSE = RMSE_Pred(X,filename)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

[COP_Xp, Contact] = COP_Xp_f(filename,X);
[COP_Pred] = COP_Pred_f(filename,X);

Solids{2,1}='RFoot';
Solids{1,1}='LFoot';

nbframe = length(Contact);
for i = 1: nbframe
    for ii =  1:2
        if Contact(ii,i)==0
            COP_Pred.(Solids{ii})(i,:)= [NaN; NaN; NaN];
            COP_Xp.(Solids{ii})(COP_Xp.(Solids{ii})==0) = NaN;
        end
        COP_Xp.(Solids{ii})(COP_Xp.(Solids{ii})==0) = NaN;
        COP_Pred.(Solids{ii})(COP_Xp.(Solids{ii})==0) = NaN;
    end
end


suivi_L=[];
suivi_R=[];

for i = 1: nbframe
    if Contact(1,i)>0 && COP_Pred.(Solids{1})(i,1)>-0
        suivi_L=[suivi_L,i];
    end
    if Contact(2,i)>0  && COP_Pred.(Solids{2})(i,1)>-0
        suivi_R=[suivi_R,i];
    end
end
suivi=[suivi_L,suivi_R];


for ii=1:2
    RMSE(ii,1)=rms(COP_Pred.(Solids{ii})(suivi((ii-1)*length(suivi_L)+1:end-(2-ii)*length(suivi_R)),1)-COP_Xp.(Solids{ii})(suivi((ii-1)*length(suivi_L)+1:end-(2-ii)*length(suivi_R)),1));
    %RMSE_norm(ii,1)=rms((COP_Pred.(Solids{ii})(suivi((ii-1)*length(suivi_L)+1:end-(2-ii)*length(suivi_R)),1))/Foot_length-COP_Xp.(Solids{ii})(suivi((ii-1)*length(suivi_L)+1:end-(2-ii)*length(suivi_R)),1))/Foot_length);
    RMSE(ii,2)=rms(COP_Pred.(Solids{ii})(suivi((ii-1)*length(suivi_L)+1:end-(2-ii)*length(suivi_R)),2)-COP_Xp.(Solids{ii})(suivi((ii-1)*length(suivi_L)+1:end-(2-ii)*length(suivi_R)),2));
    %RMSE_norm(ii,2)=rms((COP_Pred.(Solids{ii})(suivi((ii-1)*length(suivi_L)+1:end-(2-ii)*length(suivi_R)),2))/Foot_large-COP_Xp.(Solids{ii})(suivi((ii-1)*length(suivi_L)+1:end-(2-ii)*length(suivi_R)),2))/Foot_large);
    RMSE(ii,3)=rms(COP_Pred.(Solids{ii})(suivi((ii-1)*length(suivi_L)+1:end-(2-ii)*length(suivi_R)),3)-COP_Xp.(Solids{ii})(suivi((ii-1)*length(suivi_L)+1:end-(2-ii)*length(suivi_R)),3));
    
end

RMSE = mean(RMSE,1);
% 
% 
% figure
% hold on
% 
% for ii=[1,2]
%     subplot(1,2,1)
%     COP_Xp.(Solids{ii})(COP_Xp.(Solids{ii})==0) = NaN;
%     COP_Pred.(Solids{ii})(COP_Pred.(Solids{ii})==0) = NaN;
%     
%     x=COP_Pred.(Solids{ii})(:,2);
%     y=COP_Pred.(Solids{ii})(:,1); 
%     plot(x,'*')
%     hold on
%     x=COP_Xp.(Solids{ii})(:,2);
%     plot(x,'*')
%     subplot(1,2,2)
%     plot(y,'*')
%     hold on
%     y=COP_Xp.(Solids{ii})(:,1); 
%     plot(y,'*')
% end
% 

% 
% figure
% hold on
% 
% for ii=[1,2]
%     subplot(1,2,ii)
%     COP_Xp.(Solids{ii})(COP_Xp.(Solids{ii})==0) = NaN;
%     x=COP_Pred.(Solids{ii})(:,2);
%     y=COP_Pred.(Solids{ii})(:,1);
%     plot(y,x,'*')
%     hold on
%     x=COP_Xp.(Solids{ii})(:,2);
%     y=COP_Xp.(Solids{ii})(:,1);
%     plot(y,x,'*')
%     
% end
% 
% 
% figure
% hold on
% 
% for ii=[1,2]
%     %subplot(1,2,ii)
%     COP_Xp.(Solids{ii})(COP_Xp.(Solids{ii})==0) = NaN;
%     y=COP_Pred.(Solids{ii})(:,2);
%     x=COP_Pred.(Solids{ii})(:,1);
%     z=COP_Pred.(Solids{ii})(:,3);
%     plot3(x,y,z,'*')
%     hold on
%     y=COP_Xp.(Solids{ii})(:,2);
%     x=COP_Xp.(Solids{ii})(:,1);
%     z=COP_Xp.(Solids{ii})(:,3);
%     plot3(x,y,z,'*')
% end



end


