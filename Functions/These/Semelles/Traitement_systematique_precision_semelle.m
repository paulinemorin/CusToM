figure
hold on

for ii=1:3
    subplot(1,3,ii)
    plot(COP_Xp.LFoot(:,ii),'*')
    hold on
    plot(COP_Xp.RFoot(:,ii),'*')
end

figure
hold on

for ii=1:3
    subplot(1,3,ii)
    plot(COP_Insole.LFoot(:,ii))
    hold on
    plot(COP_Insole.RFoot(:,ii))
end


figure
hold on

for ii=[1,2]
    subplot(1,2,ii)
    COP_Xp.(Solids{ii})(COP_Xp.(Solids{ii})==0) = NaN;
    CoP_Pos(CoP_Pos==0) = NaN;
    taille=[Foot_length_insole,Foot_large_insole];
    x=[taille(1).*CoP_Pos((ii-1)*2+1,:).*Contact_bool(ii,:)];
    y=[taille(2).*CoP_Pos((ii-1)*2+2,:).*Contact_bool(ii,:)];
    plot(y,x,'*')
    hold on
    %plot(COP_Xp.(Solids{ii})(:,1),COP_Xp.Solids{ii}(:,2))

end

figure
hold on

for ii=[1,2]
    subplot(1,2,ii)
    COP_Xp.(Solids{ii})(COP_Xp.(Solids{ii})==0) = NaN;
    x=[COP_Xp.(Solids{ii})(:,2)];
    y=[COP_Xp.(Solids{ii})(:,1)];
    plot(y,x,'*')
    hold on
    %plot(COP_Xp.(Solids{ii})(:,1),COP_Xp.Solids{ii}(:,2))

end

figure
hold on

for i = 1: nbframe
    if Contact_detection_sum(1,i)==0
        CoP_Pos_G(1:3,i)= [NaN; NaN; NaN];
    end
    if Contact_detection_sum(2,i)==0
        CoP_Pos_G(4:6,i)= [NaN; NaN; NaN];
    end
end

for ii=[1,2]
    subplot(1,2,ii)
    COP_Xp.(Solids{ii})(COP_Xp.(Solids{ii})==0) = NaN;
    x=[CoP_Pos_G((ii-1)*3+2,:)];
    y=[CoP_Pos_G((ii-1)*3+1,:)];
    plot(y,x,'*')
    hold on
    %plot(COP_Xp.(Solids{ii})(:,1),COP_Xp.Solids{ii}(:,2))

end

Contact_bool=Contact_detection_sum;
Contact_bool(1,Contact_bool(1,:)>0) = 1;
Contact_bool(2,Contact_bool(2,:)>0) = 1;

for i=1:nbframe
    for ii=1:2
    if (GRF_Xp(i).Visual(6,ii)-GRF_Xp(i).Visual(3,ii))<50%&&Contact_detection_sum(ii,i)>0
        Contact_bool(ii,i) = 0;
        %Contact_bool(2,i) = 0;
    end
    end
end

figure
hold on

for ii=[1,2]
    subplot(1,2,ii)
    taille=[Foot_length_insole,Foot_large_insole];
    plot(taille(ii).*CoP_Pos(ii,:).*Contact_bool(1,:),'*')
    hold on
    plot(taille(ii).*CoP_Pos(ii+2,:).*Contact_bool(2,:),'*')
    %plot(Contact_detection_sum(1,:)/10)
    %plot(Contact_detection_sum(2,:)/10)
end

figure
hold on

for ii=1:3
    subplot(1,3,ii)
    plot(Contact_bool(1,:)'.*COP_Xp.LFoot(:,ii))
    hold on
    plot(Contact_bool(2,:)'.*COP_Xp.RFoot(:,ii))
end

figure
hold on

for ii=1:3
    subplot(1,3,ii)
    plot(CoP_Pos_G(ii,:),'*')
    hold on
    plot(Contact_bool(1,:)'.*COP_Xp.LFoot(:,ii),'*')
    plot(CoP_Pos_G(ii+3,:),'*')
    plot(Contact_bool(2,:)'.*COP_Xp.RFoot(:,ii),'*')
end

figure
hold on

for ii=1:3
    subplot(1,3,ii)
    plot(COP_Xp.LFoot(:,ii))
    hold on
    plot(COP_Xp.RFoot(:,ii))
end


figure
hold on
plot(Contact_detection_sum(1,:)*100)
plot(Contact_detection_sum(2,:)*100)
legend('Gauche','Droite')

figure
hold on

for ii=1:3
    subplot(1,3,ii)
    plot(GRF.LFoot(:,ii))
    hold on
    plot(GRF.RFoot(:,ii))
end



figure
Rxframel=eye(4);
Rxframel(1:3,1:3)=Rx.(Solids{1});
Rxframel(1:3,4)=Locx.(Solids{1});

Rframel=eye(4);
Rframel(1:3,1:3)=R.(Solids{1});
Rframel(1:3,4)=Loc.(Solids{1});


triad('matrix',Rframel)
hold on
triad('matrix',Rframel*Rxframel)

Rxframer=eye(4);
Rxframer(1:3,1:3)=Rx.(Solids{2});
Rxframer(1:3,4)=Locx.(Solids{2});

Rframer=eye(4);
Rframer(1:3,1:3)=R.(Solids{2});
Rframer(1:3,4)=Loc.(Solids{2});


triad('matrix',Rframer,'scale',0.5)
hold on
triad('matrix',Rframer*Rxframer,'scale',0.5)

triad('matrix',eye(4))



Rbase=eye(4);
Rbase(1:3,1:3)=R.(Solids{1});
triad('matrix',Rbase)
%triad('matrix',eye(4))


clear all
close all
clc
cheminglobal=pwd;
cd(cheminglobal);


doss=dir(cheminglobal);
doss={doss([doss.isdir]').name}';
doss(1:2)=[];

SUJETS=doss;
nb_s=length(SUJETS);


for ii=6:nb_s
    cd(fullfile(cheminglobal,SUJETS{ii}));
    Results = struct();
    load AnalysisParameters
    load('ModelParameters.mat')
    load('BiomechanicalModel.mat')
    filenameliste=AnalysisParameters.filename;
    
    
    L=length(filenameliste);
    for j = 1:length(filenameliste)
        filename_calib=filenameliste{j};
        filename_calib=strrep(filename_calib, '.c3d','');
        if contains(filename_calib,'Direction')
            filename_calib_type='ChangmntDirection';
        elseif contains(filename_calib,'Run')
            filename_calib_type='Run';
        elseif contains(filename_calib,'Walk')
            filename_calib_type='Walk';
        end
        
        % %filename='ChangmntDirectionVinyl0002';
        % %filename='WalkVinyl0001';
        X = CalibrationInsoleModel(filename_calib);
        
        
        for i = 1:length(filenameliste)
            filename=filenameliste{i};
            filename=strrep(filename, '.c3d','');
            if contains(filename,'Direction')
                filename_type='ChangmntDirection';
            elseif contains(filename,'Run')
                filename_type='Run';
            elseif contains(filename,'Walk')
                filename_type='Walk';
            end
            
            Results(i+(j-1)*L).filename_calib=filename_calib;
            Results(i+(j-1)*L).filename_calib_type=filename_calib_type;
            Results(i+(j-1)*L).filename=filename;
            Results(i+(j-1)*L).filename_type=filename_type;
            Results(i+(j-1)*L).RMSE_Insole = RMSE_COP(X, filename);
            Results(i+(j-1)*L).RMSE_Insole_liste = error_COP_liste(X, filename);
            Results(i+(j-1)*L).Fz_Insole_liste = Fz_COP(X, filename);
            % Results(i++(j-1)*15).RMSE_Predi = RMSE_Pred(X, filename);
            %Results(i).RMSE_Pied_Frame = RMSE_COP_Pied(X, filename);
            Results(i+(j-1)*L).RMSE_Insole_x = Results(i+(j-1)*L).RMSE_Insole(1);
            Results(i+(j-1)*L).RMSE_Insole_y = Results(i+(j-1)*L).RMSE_Insole(2);
            %         [ExternalForcesComputationResults] = ExternalForcesPredictionCOP(filename, AnalysisParameters, BiomechanicalModel, ModelParameters);
            %         SaveDataExternalForces(filename,ExternalForcesComputationResults);
            %         [RMSE_Pred, RMSE_Pred_COP] = RMSE_PostProcess(filename);
            %         Results(i+(j-1)*L).RMSE_Pred=RMSE_Pred;
            %         Results(i+(j-1)*L).RMSE_Pred_COP=RMSE_Pred_COP;
            
        end
    end
save('matlab_Results_liste.mat','Results')
end


% Results = struct();
% load AnalysisParameters
% load('ModelParameters.mat')
% load('BiomechanicalModel.mat')
% filenameliste=AnalysisParameters.filename;
% 
% 
% L=length(filenameliste);
% for j = 1:length(filenameliste)
%     filename_calib=filenameliste{j};
%     filename_calib=strrep(filename_calib, '.c3d','');
%     if contains(filename_calib,'Direction')
%         filename_calib_type='ChangmntDirection';
%     elseif contains(filename_calib,'Run')
%         filename_calib_type='Run';
%     elseif contains(filename_calib,'Walk')
%         filename_calib_type='Walk';
%     end
%     
%     % %filename='ChangmntDirectionVinyl0002';
%     % %filename='WalkVinyl0001';
%     X = CalibrationInsoleModel(filename_calib);
%     
%     
%     for i = 1:length(filenameliste)
%         filename=filenameliste{i};
%         filename=strrep(filename, '.c3d','');
%         if contains(filename,'Direction')
%             filename_type='ChangmntDirection';
%         elseif contains(filename,'Run')
%             filename_type='Run';
%         elseif contains(filename,'Walk')
%             filename_type='Walk';
%         end
%         
%         Results(i+(j-1)*L).filename_calib=filename_calib;
%         Results(i+(j-1)*L).filename_calib_type=filename_calib_type;
%         Results(i+(j-1)*L).filename=filename;
%         Results(i+(j-1)*L).filename_type=filename_type;
%         Results(i+(j-1)*L).RMSE_Insole = RMSE_COP(X, filename);
%         Results(i+(j-1)*L).RMSE_Insole_liste = error_COP_liste(X, filename);
%         Results(i+(j-1)*L).Fz_Insole_liste = Fz_COP(X, filename);
%            % Results(i++(j-1)*15).RMSE_Predi = RMSE_Pred(X, filename);
%            %Results(i).RMSE_Pied_Frame = RMSE_COP_Pied(X, filename);
%         Results(i+(j-1)*L).RMSE_Insole_x = Results(i+(j-1)*L).RMSE_Insole(1);
%         Results(i+(j-1)*L).RMSE_Insole_y = Results(i+(j-1)*L).RMSE_Insole(2);
% %         [ExternalForcesComputationResults] = ExternalForcesPredictionCOP(filename, AnalysisParameters, BiomechanicalModel, ModelParameters);
% %         SaveDataExternalForces(filename,ExternalForcesComputationResults);
% %         [RMSE_Pred, RMSE_Pred_COP] = RMSE_PostProcess(filename);
% %         Results(i+(j-1)*L).RMSE_Pred=RMSE_Pred;
% %         Results(i+(j-1)*L).RMSE_Pred_COP=RMSE_Pred_COP;
% 
%     end
% end


T = struct2table(Results);
Calib_type_order = {'ChangmntDirection','Run','Walk'};
T.filename_calib_type = categorical(T.filename_calib_type,Calib_type_order);
T.filename_type = categorical(T.filename_type);
T.filename = categorical(T.filename);

% figure
% boxchart(T.filename_type,T.RMSE_Insole_x,'GroupByColor',T.filename_calib_type)
% title('Erreur sur la position du centre de pression')
% %boxchart(T.filename_calib,T.RMSE_Insole_x)
% ylabel('RMSE selon l axe x')
% xlabel('Type de l essai pour la calibration')
% legend


figure
group = T.filename_type.*T.filename_calib;
Calib_type_order = {'ChangmntDirection','','','','','Run','','','','','Walk','','','','','ChangmntDirection','','','','','Run','','','','','Walk','','','','','ChangmntDirection','','','','','Run','','','','','Walk','','','',''};

boxchart(group,T.RMSE_Insole_x,'GroupByColor',T.filename_calib_type)
title('Erreur sur la position du centre de pression')
%boxchart(T.filename_calib,T.RMSE_Insole_x)
ylabel('RMSE selon l axe x')
xlabel('Type de l essai traité')
legend



figure
boxchart(group,T.RMSE_Insole_y,'GroupByColor',T.filename_calib_type)
title('Erreur sur la position du centre de pression')
%boxchart(T.filename_calib,T.RMSE_Insole_x)
ylabel('RMSE selon l axe y')
xlabel('Type de l essai traité')
legend
% 
% boxchart(T.filename,T.RMSE_Insole_y,'GroupByColor',T.filename_calib_type)
% 
% title('Erreur sur la position du centre de pression')
% %boxchart(T.filename_calib,T.RMSE_Insole_x)
% ylabel('RMSE selon l axe y')
% xlabel('Type de l essai traité')
% legend



