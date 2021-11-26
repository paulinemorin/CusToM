
function error = error_COP_liste(X, filename)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here


[COP_Xp, Contact] = COP_Xp_f(filename,X);

CoP_Pos_G = CoPInsoleGlob(X,filename);


Solids{2,1}='RFoot';
Solids{1,1}='LFoot';

nbframe = length(Contact);
for i = 1: nbframe
    for ii =  1:2
        if Contact(ii,i)==0
            CoP_Pos_G(3*ii-2:3*ii-2+2,i)= [NaN; NaN; NaN];
        end
        COP_Xp.(Solids{ii})(COP_Xp.(Solids{ii})==0) = NaN;
    end
end

load('AnalysisParameters.mat')
load(fullfile(filename,'ExternalForcesComputationResults.mat'))
[start,stop] = croper(ExternalForcesComputationResults,AnalysisParameters, filename);


suivi_L=[];
suivi_R=[];

%for i = 1: nbframe
for i = start : stop
    if Contact(1,i)>0
        suivi_L=[suivi_L,i];
    end
    if Contact(2,i)>0
        suivi_R=[suivi_R,i];
    end
end
suivi=[suivi_L,suivi_R];

%% taille de pied et semelles pour normaliser
[Foot_length,Foot_large] = Foot_size(filename);
[Foot_length_insole,Foot_large_insole] = Foot_size_insole(Foot_length,Foot_large);
taille=[Foot_length_insole,Foot_large_insole];

error=[];
for ii=1:2
    %cur_s=num_s(ii); %LFoot and RFoot
    COP_Insole.(Solids{ii})=CoP_Pos_G((ii-1)*3+1:(ii-1)*3+3,:);

    RMSE(ii,1)=rms(COP_Insole.(Solids{ii})(1,suivi((ii-1)*length(suivi_L)+1:end-(2-ii)*length(suivi_R)))'-COP_Xp.(Solids{ii})(suivi((ii-1)*length(suivi_L)+1:end-(2-ii)*length(suivi_R)),1));
    RMSE_norm(ii,1)=rms((COP_Insole.(Solids{ii})(1,suivi((ii-1)*length(suivi_L)+1:end-(2-ii)*length(suivi_R)))'-COP_Xp.(Solids{ii})(suivi((ii-1)*length(suivi_L)+1:end-(2-ii)*length(suivi_R)),1))/Foot_length);
    RMSE(ii,2)=rms(COP_Insole.(Solids{ii})(2,suivi((ii-1)*length(suivi_L)+1:end-(2-ii)*length(suivi_R)))'-COP_Xp.(Solids{ii})(suivi((ii-1)*length(suivi_L)+1:end-(2-ii)*length(suivi_R)),2));
    RMSE_norm(ii,2)=rms((COP_Insole.(Solids{ii})(2,suivi((ii-1)*length(suivi_L)+1:end-(2-ii)*length(suivi_R)))'-COP_Xp.(Solids{ii})(suivi((ii-1)*length(suivi_L)+1:end-(2-ii)*length(suivi_R)),2))/Foot_large);
    RMSE(ii,3)=rms(COP_Insole.(Solids{ii})(3,suivi((ii-1)*length(suivi_L)+1:end-(2-ii)*length(suivi_R)))'-COP_Xp.(Solids{ii})(suivi((ii-1)*length(suivi_L)+1:end-(2-ii)*length(suivi_R)),3));
    

end

RMSE = mean(RMSE,1);

error_L(:,1)=(COP_Insole.(Solids{1})(1,suivi((1-1)*length(suivi_L)+1:end-(2-1)*length(suivi_R)))'-COP_Xp.(Solids{1})(suivi((1-1)*length(suivi_L)+1:end-(2-1)*length(suivi_R)),1));
error_L(:,2)=(COP_Insole.(Solids{1})(2,suivi((1-1)*length(suivi_L)+1:end-(2-1)*length(suivi_R)))'-COP_Xp.(Solids{1})(suivi((1-1)*length(suivi_L)+1:end-(2-1)*length(suivi_R)),2));
  
error_R(:,1)=(COP_Insole.(Solids{2})(1,suivi((2-1)*length(suivi_L)+1:end-(2-2)*length(suivi_R)))'-COP_Xp.(Solids{2})(suivi((2-1)*length(suivi_L)+1:end-(2-2)*length(suivi_R)),1));
error_R(:,2)=(COP_Insole.(Solids{2})(2,suivi((2-1)*length(suivi_L)+1:end-(2-2)*length(suivi_R)))'-COP_Xp.(Solids{2})(suivi((2-1)*length(suivi_L)+1:end-(2-2)*length(suivi_R)),2));

load(fullfile(filename,'InsoleData.mat'))

Force_LF=InsoleData.data(suivi((1-1)*length(suivi_L)+1:end-(2-1)*length(suivi_R)),21);
Force_RF=InsoleData.data(suivi((2-1)*length(suivi_L)+1:end-(2-2)*length(suivi_R)),43);

error=[error_L;error_R];
% 
% figure(5)
% hold on
% set(gca,'color','w')
% set(gcf,'color','w')
% plot(Force_LF,abs(error_L(:,1)),'*')
% plot(Force_RF,abs(error_R(:,1)),'*')
% title('Erreur sur l axe medio-lateral')
% xlabel('Force verticale totale captée par la semelle (N)')
% ylabel('Erreur (m)')
% 
% figure(6)
% hold on
% set(gca,'color','w')
% set(gcf,'color','w')
% plot(Force_LF,abs(error_L(:,2)),'*')
% plot(Force_RF,abs(error_R(:,2)),'*')
% title('Erreur sur l axe antero-posterieur')
% xlabel('Force verticale totale captée par la semelle (N)')
% ylabel('Erreur (m)')
% 
% 

% figure
% hold on
% 
% for ii=[1,2]
%     subplot(1,2,1)
%     COP_Xp.(Solids{ii})(COP_Xp.(Solids{ii})==0) = NaN;
%     x=CoP_Pos_G((ii-1)*3+2,:);
%     y=CoP_Pos_G((ii-1)*3+1,:);
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
% % 
% figure
% hold on
% set(gca,'color','w')
% set(gcf,'color','w')
% title_liste={'Pied gauche','Pied droit'};
% 
% %figure
% 
% for ii=[1,2]
%     subplot(1,2,ii)
%     COP_Xp.(Solids{ii})(COP_Xp.(Solids{ii})==0) = NaN;
%     x=CoP_Pos_G((ii-1)*3+2,start:stop);
%     y=CoP_Pos_G((ii-1)*3+1,start:stop);
%     plot(y,x,'*')
%     hold on
%     x=COP_Xp.(Solids{ii})(start:stop,2);
%     y=COP_Xp.(Solids{ii})(start:stop,1);
%     plot(y,x,'*')
%     
% end
% 
% for ii=[1,2]
%     subplot(1,2,ii)
%     COP_Xp.(Solids{ii})(COP_Xp.(Solids{ii})==0) = NaN;
%     x=CoP_Pos_G((ii-1)*3+2,start:stop);
%     y=CoP_Pos_G((ii-1)*3+1,start:stop);
%     plot(y,x,'*')
%     hold on
%     x=COP_Xp.(Solids{ii})(start:stop,2);
%     y=COP_Xp.(Solids{ii})(start:stop,1);
%     plot(y,x,'*')
%         axis([nanmean(y)-0.5*taille(2) nanmean(y)+0.5*taille(2) nanmean(x)-0.5*taille(1) nanmean(x)+0.5*taille(1)])
%     title(title_liste(ii))
%     %axis([-0.5*taille(2) 0.5*taille(2) -0.5*taille(1) 0.5*taille(1)])
%     ylabel({'y dans le repère global','longeur semelle (m)'})
%     xlabel({'x dans le repère global','largeur semelle (m)'})
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
%     y=CoP_Pos_G((ii-1)*3+2,start:stop);
%     z=CoP_Pos_G((ii-1)*3+3,start:stop);
%     x=CoP_Pos_G((ii-1)*3+1,start:stop);
%     plot3(x,y,z,'*')
%     hold on
%     y=COP_Xp.(Solids{ii})(start:stop,2);
%     x=COP_Xp.(Solids{ii})(start:stop,1);
%     z=COP_Xp.(Solids{ii})(start:stop,3);
%     plot3(x,y,z,'*')
% end



end

