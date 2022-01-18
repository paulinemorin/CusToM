function [RMSE] = CalibRMSEInsole(X,CoP_Pos,COP_Xp, dw,dv0,Human_model,Contact_detection_sum,Prediction,nbframe,Foot_large,Foot_length,taille, p_pelvis,r_pelvis,v0,w,Solids,q,dq,ddq,num_s,time_g,time_d, filename,start,stop)
%UNTITLED4 Summary of this function goes here
%   X = [theta_1; Theta_2; Theta_3; X; Y; Z]

%% Définition passage modèle de pied à semelle


R1_l = rotx(X(1));
R2_l = roty(X(2));
R3_l = rotz(X(3));

Rx.(Solids{1}) =R1_l*R2_l*R3_l;
Locx.(Solids{1})=X(4:6);

R1_r = rotx(X(7));
R2_r = roty(X(8));
R3_r = rotz(X(9));

Rx.(Solids{2}) = R1_r*R2_r*R3_r;
Locx.(Solids{2})=X(10:12);

%% Transformation cop data semelles

for ii=[1,2]
    CoP_Pos_G((ii-1)*3+1,:)=[taille(1).*CoP_Pos((ii-1)*2+1,:)];
       % CoP_Pos_G((ii-1)*3+1,:)=[taille(1).*CoP_Pos((ii-1)*2+1,:)+taille(1)*0.5];
    CoP_Pos_G((ii-1)*3+2,:)=[taille(2).*CoP_Pos((ii-1)*2+2,:)];
    CoP_Pos_G((ii-1)*3+3,:)=zeros(1,length( CoP_Pos_G((ii-1)*3+2,:)));
end
 
% changement repère pou_r repère direct, pied droit semelle
CoP_Pos_G(2,:)=-CoP_Pos_G(2,:);

% changement repère de semelle à modèle de pied
a = [0 -1 0;
    0  0 -1;
    1 0  0];

for i=[time_g,time_d]
    [Rot_R, R_pos, Rot_L, L_pos] = FootFrameMarker(filename, i);

    % Get the forces applied on the solids
    if i == time_d     
       R.(Solids{2}) = Rot_R;
       Loc.(Solids{2}) = R_pos';
    end
    if i == time_g
        R.(Solids{1}) = Rot_L;
        Loc.(Solids{1}) = L_pos';
    end
    
end

for i=1:nbframe  
    % Get the forces applied on the solids
    for ii=1:numel(num_s)
        if Contact_detection_sum(ii,i)>0
            CoP_Pos_G((ii-1)*3+1:(ii-1)*3+3,i)=R.(Solids{ii})* (Rx.(Solids{ii})*(a'*CoP_Pos_G((ii-1)*3+1:(ii-1)*3+3,i)+Locx.(Solids{ii})))+Loc.(Solids{ii}); 
        end
    end
end

suivi_L=[];
suivi_R=[];

for i = start: stop
%for i = 1 : nbframe
    if Contact_detection_sum(1,i)>0
        suivi_L=[suivi_L,i];
    end
    if Contact_detection_sum(2,i)>0
        suivi_R=[suivi_R,i];
    end
end
suivi=[suivi_L,suivi_R];

for ii=1:numel(num_s)
    %cur_s=num_s(ii); %LFoot and RFoot
    COP_Insole.(Solids{ii})=CoP_Pos_G((ii-1)*3+1:(ii-1)*3+3,:);

    RMSE(ii,1)=rms(COP_Insole.(Solids{ii})(1,suivi((ii-1)*length(suivi_L)+1:end-(2-ii)*length(suivi_R)))'-COP_Xp.(Solids{ii})(suivi((ii-1)*length(suivi_L)+1:end-(2-ii)*length(suivi_R)),1));
    RMSE_norm(ii,1)=rms((COP_Insole.(Solids{ii})(1,suivi((ii-1)*length(suivi_L)+1:end-(2-ii)*length(suivi_R)))'-COP_Xp.(Solids{ii})(suivi((ii-1)*length(suivi_L)+1:end-(2-ii)*length(suivi_R)),1))/Foot_length);
    RMSE(ii,2)=rms(COP_Insole.(Solids{ii})(2,suivi((ii-1)*length(suivi_L)+1:end-(2-ii)*length(suivi_R)))'-COP_Xp.(Solids{ii})(suivi((ii-1)*length(suivi_L)+1:end-(2-ii)*length(suivi_R)),2));
    RMSE_norm(ii,2)=rms((COP_Insole.(Solids{ii})(2,suivi((ii-1)*length(suivi_L)+1:end-(2-ii)*length(suivi_R)))'-COP_Xp.(Solids{ii})(suivi((ii-1)*length(suivi_L)+1:end-(2-ii)*length(suivi_R)),2))/Foot_large);
    RMSE(ii,3)=rms(COP_Insole.(Solids{ii})(3,suivi((ii-1)*length(suivi_L)+1:end-(2-ii)*length(suivi_R)))'-COP_Xp.(Solids{ii})(suivi((ii-1)*length(suivi_L)+1:end-(2-ii)*length(suivi_R)),3));
end

RMSE=sum(RMSE(:,1))+sum(RMSE(:,2))+sum(RMSE(:,3));

for i = 1: nbframe
    if Contact_detection_sum(1,i)==0
        CoP_Pos_G(1:3,i)= [NaN; NaN; NaN];
    end
    if Contact_detection_sum(2,i)==0
        CoP_Pos_G(4:6,i)= [NaN; NaN; NaN];
    end
end

% figure
% hold on
% 
% for ii=[1,2]
%     subplot(1,2,ii)
%     COP_Xp.(Solids{ii})(COP_Xp.(Solids{ii})==0) = NaN;
%     x=CoP_Pos_G((ii-1)*3+2,:);
%     y=CoP_Pos_G((ii-1)*3+1,:);
%     plot(y,x,'*')
%     hold on
%     x=COP_Xp.(Solids{ii})(:,2);
%     y=COP_Xp.(Solids{ii})(:,1);
%     plot(y,x,'*')
%     
% end
% 
% figure
% hold on
% set(gca,'color','w')
% set(gcf,'color','w')
% title_liste={'Left Foot','Right Foot'};
% 
% for ii=[1,2]
%     subplot(1,2,ii)
%     COP_Xp.(Solids{ii})(COP_Xp.(Solids{ii})==0) = NaN;
%     x=CoP_Pos_G((ii-1)*3+2,suivi((ii-1)*length(suivi_L)+1:end-(2-ii)*length(suivi_R)));
%     y=CoP_Pos_G((ii-1)*3+1,suivi((ii-1)*length(suivi_L)+1:end-(2-ii)*length(suivi_R)));
%     plot(y,x,'*')
%     hold on
%     x=COP_Xp.(Solids{ii})(suivi((ii-1)*length(suivi_L)+1:end-(2-ii)*length(suivi_R)),2);
%     y=COP_Xp.(Solids{ii})(suivi((ii-1)*length(suivi_L)+1:end-(2-ii)*length(suivi_R)),1);
%     plot(y,x,'*')
%         axis([nanmean(y)-0.5*taille(2) nanmean(y)+0.5*taille(2) nanmean(x)-0.5*taille(1) nanmean(x)+0.5*taille(1)])
%     title(title_liste(ii))
%     %axis([-0.5*taille(2) 0.5*taille(2) -0.5*taille(1) 0.5*taille(1)])
%     xlabel({'Medio-lateral position (m)'})
%     ylabel({'Antero-posterior position (m)'})
%     
% %     ylabel({'Medio-lateral position','longeur semelle (m)'})
% %     xlabel({'x dans le repère global','largeur semelle (m)'})
%     
% end
% 
% legend('CoP from insole','CoP from plateform')


end
