clear all
close all
clc
cheminglobal=pwd;
cd(cheminglobal);

%Chargement des tailles
excel='BAHAMaS_sujets.xlsx';
% excel='Tableau_sujets_Militaire.xlsx';
%excel='Tableau_sujets_Civil.xlsx';
tab=importdata(excel);

names = tab.textdata(2:end,1);
taille = tab.data(:,2);  %adapté
masse = tab.data(:,3);
% taille = tab.data(:,7);
% masse = tab.data(:,8);

doss=dir(cheminglobal);
doss={doss([doss.isdir]').name}';
doss(1:2)=[];

SUJETS=doss;
nb_s=length(SUJETS);

if exist(fullfile(pwd, 'meta.mat'),'file')==2
    load('meta')
else
    %% METADATA TREATEMENT
    meta=struct();
    for ii=1:nb_s
        meta(ii).name=SUJETS{ii};
        [~,ind]=intersect(names,SUJETS{ii});
        meta(ii).taille = taille(ind,1);
        meta(ii).masse = masse(ind,1);
        meta(ii).path = fullfile(cheminglobal,SUJETS{ii});
    end
    %% TRIAL LIST
    for ii=1:nb_s
        cd(meta(ii).path)
        c3dFiles = dir ([meta(ii).path filesep '*.c3d']);
        meta(ii).TrialList={c3dFiles.name}';
        meta(ii).TrialList_tasks=setdiff(meta(ii).TrialList,{'ROM.c3d','Statref.c3d'});
    end
    
    %% GENERATE PARAMTERS
    for ii=1:nb_s
        cd(meta(ii).path)
        
        meta(ii).ModelParameters=struct();
        meta(ii).AnalysisParameters=struct();
        
        [meta(ii).ModelParameters,...
            meta(ii).AnalysisParameters] = Parameters(meta(ii));
        
         for jj=1:length(meta(ii).TrialList_tasks)
            meta(ii).ExternalForces(jj).AnalysisParameters = ...
                AnalysisExternalForces( meta(ii).AnalysisParameters, ...
                meta(ii).TrialList(jj));
%             
%             %if contains(meta(ii).TrialList(jj),'Marche')
%             %if contains(meta(ii).TrialList(jj),'Course')
%             if contains(meta(ii).TrialList_tasks(jj),'Course') || contains(meta(ii).TrialList_tasks(jj),'Marche')...
%                     || contains(meta(ii).TrialList_tasks(jj),'ChgtDirection') || contains(meta(ii).TrialList_tasks(jj),'Saut')
%                 
%                 
                 meta(ii).ExternalForcesPrediction(jj).AnalysisParameters =...
                     AnalysisForcePrediction(meta(ii).AnalysisParameters, ...
                     meta(ii).TrialList(jj));

%                 meta(ii).ExternalForcesPrediction_plus(jj).AnalysisParameters =...
%                     AnalysisForcePrediction(meta(ii).AnalysisParameters, ...
%                     meta(ii).TrialList(jj));
%                 
%             end
%             
         end
    end
    save(fullfile(cheminglobal,'meta'),'meta')
end
%% Frames of interest - event
for ii=1:nb_s
    cd(meta(ii).path)
    for jj=1:length(meta(ii).TrialList_tasks)
        
        t= meta(ii).TrialList(jj);
        
        if contains(t,'Marche') || contains(t,'Walk')% A faire pour les autres taches.
            FootOnPlates = FootOnPlates_ref;
            %FootOnPlates = meta(ii).ExternalForces(jj).AnalysisParameters.ExternalForces.Options;
            [meta(ii).FrameOfInterest(jj).Event,...
                meta(ii).FrameOfInterest(jj).ID]=...
                DetectOnlyOneFootPhase(t,FootOnPlates);
            
        elseif contains(t,'Course') || contains(t,'Direction')|| contains(t,'Run')
            FootOnPlates = FootOnPlates_ref;
            %FootOnPlates = meta(ii).ExternalForces(jj).AnalysisParameters.ExternalForces.Options;
            [meta(ii).FrameOfInterest(jj).Event,...
                meta(ii).FrameOfInterest(jj).ID]=...
                DetectRunPhases(t,FootOnPlates);
            
        elseif contains(t,'Saut')
        %if contains(t,'Saut')
            
        
            FootOnPlates = meta(ii).ExternalForces(jj).AnalysisParameters.ExternalForces.Options;
            
            [meta(ii).FrameOfInterest(jj).Event,...
                meta(ii).FrameOfInterest(jj).ID]=...
                DetectJumpPhases(t,FootOnPlates);
             
             
             cd(meta(ii).path)
             meta(ii).FrameOfInterest(jj).ID(1)=DetectJumpOffPrediction(meta(ii).TrialList{jj}(1:end-4),meta(ii).ExternalForcesPrediction_plus(jj).AnalysisParameters,...
                    meta(ii).ModelParameters,meta(ii).FrameOfInterest(jj).ID);
             
                          
            
        end
    end
end
save(fullfile(cheminglobal,'meta'),'meta')

%% TRAITEMENT EMG
% for ii=1:nb_s % Ajouter pour ne pas tenir compte des EMGs mauvais pour les max !!
%     cd(meta(ii).path)
%     EMG_traitementAll(meta(ii));
% end
%% CALIBRATION GEOMETRIQUE + INVERSE KINEMATICS
for ii=1:nb_s
    cd(meta(ii).path)
    
    CalibrateModelGenerationNum(meta(ii).ModelParameters,meta(ii).AnalysisParameters);
    InverseKinematics(meta(ii).AnalysisParameters);
    
end

%% EXTERNAL FORCE COMPUTATION + INVERSE DYNAMICS
for ii=10:nb_s
    cd(meta(ii).path)
%     for jj=1:length(meta(ii).TrialList)
%         if ~contains(meta(ii).TrialList(jj),'ROM')
             for jj=1:length(meta(ii).TrialList_tasks)
        if ~contains(meta(ii).TrialList_tasks(jj),'ROM')
             ExternalForcesComputation(meta(ii).ExternalForces(jj).AnalysisParameters);
             InverseDynamics(meta(ii).ExternalForces(jj).AnalysisParameters);
%             ExternalForcesComputation(meta(ii).ExternalForcesPrediction(jj).AnalysisParameters,...
%                     meta(ii).ModelParameters);
%             
            %Prediction si marche
            %if contains(meta(ii).TrialList_tasks(jj),'Marche')
            %if contains(meta(ii).TrialList_tasks(jj),'Course') || contains(meta(ii).TrialList_tasks(jj),'Marche')...
                   % || contains(meta(ii).TrialList_tasks(jj),'ChgtDirection') 
            %if contains(meta(ii).TrialList_tasks(jj),'Saut') 
                %ExternalForcesComputation(meta(ii).ExternalForcesPrediction(jj).AnalysisParameters,...
                  %  meta(ii).ModelParameters);
%             ExternalForcesComputation_plus(meta(ii).ExternalForcesPrediction_plus(jj).AnalysisParameters,...
%                 meta(ii).ModelParameters);
%             ExternalForcesComputation_prop_adaptatif(meta(ii).ExternalForcesPrediction_plus(jj).AnalysisParameters,...
%                 meta(ii).ModelParameters);
%             ExternalForcesComputation_prop_k_fixe(meta(ii).ExternalForcesPrediction_plus(jj).AnalysisParameters,...
%                 meta(ii).ModelParameters);


            ExternalForcesComputation(meta(ii).ExternalForcesPrediction(jj).AnalysisParameters,...
                 meta(ii).ModelParameters);

%             load([meta(ii).TrialList_tasks{jj}(1:end-4) '/InverseDynamicsResults']);
%                 ResID(jj).Ext=InverseDynamicsResults;
%                 
%                 InverseDynamics(meta(ii).ExternalForcesPrediction(jj).AnalysisParameters);
%                 load([meta(ii).TrialList_tasks{jj}(1:end-4) '/InverseDynamicsResults']);
%                 ResID(jj).Pred=InverseDynamicsResults;
%             %end
        end
    end
  %  save('ResID','ResID')
end

%% POSTPROCESS PREDICTION MARCHE

for ii=1:nb_s
    cd(meta(ii).path)
    for jj=1:length(meta(ii).TrialList_tasks)
        if ~contains(meta(ii).TrialList(jj),'ROM')&& ~contains(meta(ii).TrialList(jj),'Step')
        filename = meta(ii).TrialList(jj);
        match='.c3d';
        filename =string( erase(filename,match) );
        %[meta(ii).RMS_Threshold(jj,:),meta(ii).RMS_InsoleDetection(jj,:)]=PostProcessISB(jj,filename,meta(ii).FrameOfInterest(jj).ID);
        [meta(ii).Max(jj,:),meta(ii).Moy(jj,:)]=MaxMoy(jj,filename,meta(ii).FrameOfInterest(jj).ID);
        
        end
    end
end

RMS.RMS_Threshold_Direction=[];
RMS.RMS_Insole_Direction=[];

RMS.RMS_Threshold_Run=[];
RMS.RMS_Insole_Run=[];

RMS.RMS_Threshold_Walk=[];
RMS.RMS_Insole_Walk=[];


RMS.RMS_Threshold=[];
RMS.RMS_Insole=[];

%RMS=[];

for ii=1:nb_s
    for jj=1:length(meta(ii).TrialList_tasks)
        L_T=[(meta(ii).RMS_Threshold(jj,:)),"Threshold"];
        L_I=[(meta(ii).RMS_InsoleDetection(jj,:)),"Insole"];
        if ~contains(meta(ii).TrialList(jj),'ROM')&& ~contains(meta(ii).TrialList(jj),'Step')
            if contains(meta(ii).TrialList(jj),'Direction')
                L_T=[L_T,"Direction"];
                L_I=[L_I,"Direction"];
                RMS.RMS_Threshold_Direction=[RMS.RMS_Threshold_Direction;meta(ii).RMS_Threshold(jj,:)];
                RMS.RMS_Insole_Direction=[RMS.RMS_Insole_Direction;meta(ii).RMS_InsoleDetection(jj,:)];
            elseif contains(meta(ii).TrialList(jj),'Run')
                L_T=[L_T,"Run"];
                L_I=[L_I,"Run"];
                RMS.RMS_Threshold_Run=[RMS.RMS_Threshold_Run;meta(ii).RMS_Threshold(jj,:)];
                RMS.RMS_Insole_Run=[RMS.RMS_Insole_Run;meta(ii).RMS_InsoleDetection(jj,:)];
            elseif contains(meta(ii).TrialList(jj),'Walk')
                L_T=[L_T,"Walk"];
                L_I=[L_I,"Walk"];
                RMS.RMS_Threshold_Walk=[RMS.RMS_Threshold_Walk;meta(ii).RMS_Threshold(jj,:)];
                RMS.RMS_Insole_Walk=[RMS.RMS_Insole_Walk;meta(ii).RMS_InsoleDetection(jj,:)];
            end 
            RMS.RMS_Threshold=[RMS.RMS_Threshold;meta(ii).RMS_Threshold(jj,:)];
            RMS.RMS_Insole=[RMS.RMS_Insole;meta(ii).RMS_InsoleDetection(jj,:)];
        end
        
    end
end


filename ='RMS.xlsx';
writematrix(RMS,filename)

filename ='RMS_Threshold_Direction_2.xlsx';
writematrix(RMS.RMS_Threshold_Direction,filename)
filename ='RMS_Insole_Direction_2.xlsx';
writematrix(RMS.RMS_Insole_Direction,filename)

filename ='RMS_Threshold_Run_2.xlsx';
writematrix(RMS.RMS_Threshold_Run,filename)
filename ='RMS_Insole_Run_2.xlsx';
writematrix(RMS.RMS_Insole_Run,filename)

filename ='RMS_Threshold_Walk_2.xlsx';
writematrix(RMS.RMS_Threshold_Walk,filename)
filename ='RMS_Insole_Walk_2.xlsx';
writematrix(RMS.RMS_Insole_Walk,filename)


filename ='RMS_Threshold_2.xlsx';
writematrix(RMS.RMS_Threshold,filename)
filename ='RMS_Insole_2.xlsx';
writematrix(RMS.RMS_Insole,filename)


%%

liste=[1 2 3 4 5 6 7 8 9 10];

RMS_tot_1_0=mean(meta(ii).RMS_total((liste-ones(1,length(meta(ii).TrialList)))*5+ones(1,length(meta(ii).TrialList)),:),1);
RMS_tot_75_25=mean(meta(ii).RMS_total((liste-ones(1,length(meta(ii).TrialList)))*5+2*ones(1,length(meta(ii).TrialList)),:),1);
RMS_tot_50_50=mean(meta(ii).RMS_total((liste-ones(1,length(meta(ii).TrialList)))*5+3*ones(1,length(meta(ii).TrialList)),:),1);
RMS_tot_25_75=mean(meta(ii).RMS_total((liste-ones(1,length(meta(ii).TrialList)))*5+4*ones(1,length(meta(ii).TrialList)),:),1);
RMS_tot_0_1=mean(meta(ii).RMS_total((liste-ones(1,length(meta(ii).TrialList)))*5+5*ones(1,length(meta(ii).TrialList)),:),1);

RMS_tot=[RMS_tot_1_0;
    RMS_tot_75_25;
    RMS_tot_50_50;
    RMS_tot_25_75;
    RMS_tot_0_1];


RMS_dynamic_1_0_liste=meta(ii).RMS_dynamic((liste-ones(1,length(meta(ii).TrialList)))*5+ones(1,length(meta(ii).TrialList)),:);
RMS_dynamic_1_0=mean(meta(ii).RMS_dynamic((liste-ones(1,length(meta(ii).TrialList)))*5+ones(1,length(meta(ii).TrialList)),:),1);

RMS_dynamic_75_25_liste=meta(ii).RMS_dynamic((liste-ones(1,length(meta(ii).TrialList)))*5+2*ones(1,length(meta(ii).TrialList)),:);
RMS_dynamic_75_25=mean(meta(ii).RMS_dynamic((liste-ones(1,length(meta(ii).TrialList)))*5+2*ones(1,length(meta(ii).TrialList)),:),1);

RMS_dynamic_50_50_liste=meta(ii).RMS_dynamic((liste-ones(1,length(meta(ii).TrialList)))*5+3*ones(1,length(meta(ii).TrialList)),:);
RMS_dynamic_50_50=mean(meta(ii).RMS_dynamic((liste-ones(1,length(meta(ii).TrialList)))*5+3*ones(1,length(meta(ii).TrialList)),:),1);

RMS_dynamic_25_75_liste=meta(ii).RMS_dynamic((liste-ones(1,length(meta(ii).TrialList)))*5+4*ones(1,length(meta(ii).TrialList)),:);
RMS_dynamic_25_75=mean(meta(ii).RMS_dynamic((liste-ones(1,length(meta(ii).TrialList)))*5+4*ones(1,length(meta(ii).TrialList)),:),1);

RMS_dynamic_0_1_liste=meta(ii).RMS_dynamic((liste-ones(1,length(meta(ii).TrialList)))*5+5*ones(1,length(meta(ii).TrialList)),:);
RMS_dynamic_0_1=mean(meta(ii).RMS_dynamic((liste-ones(1,length(meta(ii).TrialList)))*5+5*ones(1,length(meta(ii).TrialList)),:),1);

RMS_dynamic=[RMS_dynamic_1_0;
    RMS_dynamic_75_25;
    RMS_dynamic_50_50;
    RMS_dynamic_25_75;
    RMS_dynamic_0_1];

RMS_static_1_0_liste=[meta(ii).RMS_static((liste-ones(1,length(meta(ii).TrialList)))*5+ones(1,length(meta(ii).TrialList)),:)];
RMS_static_1_0=mean(meta(ii).RMS_static((liste-ones(1,length(meta(ii).TrialList)))*5+ones(1,length(meta(ii).TrialList)),:),1);

RMS_static_75_25_liste=meta(ii).RMS_static((liste-ones(1,length(meta(ii).TrialList)))*5+2*ones(1,length(meta(ii).TrialList)),:);
RMS_static_75_25=mean(meta(ii).RMS_static((liste-ones(1,length(meta(ii).TrialList)))*5+2*ones(1,length(meta(ii).TrialList)),:),1);

RMS_static_50_50_liste=meta(ii).RMS_static((liste-ones(1,length(meta(ii).TrialList)))*5+3*ones(1,length(meta(ii).TrialList)),:);
RMS_static_50_50=mean(meta(ii).RMS_static((liste-ones(1,length(meta(ii).TrialList)))*5+3*ones(1,length(meta(ii).TrialList)),:),1);

RMS_static_25_75_liste=meta(ii).RMS_static((liste-ones(1,length(meta(ii).TrialList)))*5+4*ones(1,length(meta(ii).TrialList)),:);
RMS_static_25_75=mean(meta(ii).RMS_static((liste-ones(1,length(meta(ii).TrialList)))*5+4*ones(1,length(meta(ii).TrialList)),:),1);

RMS_static_0_1_liste=meta(ii).RMS_static((liste-ones(1,length(meta(ii).TrialList)))*5+5*ones(1,length(meta(ii).TrialList)),:);
RMS_static_0_1=mean(meta(ii).RMS_static((liste-ones(1,length(meta(ii).TrialList)))*5+5*ones(1,length(meta(ii).TrialList)),:),1);

Directions={'Medio-lateral','Antero-posterior','Vertical','Frontal','Sagittal','Transversal'};

   
figure
set(gcf,'color','w')
for k = 1:6
    subplot(1,6,k)
%     data={[RMS_static_1_0_liste(:,k),RMS_static_75_25_liste(:,k),RMS_static_50_50_liste(:,k),RMS_static_25_75_liste(:,k),RMS_static_0_1_liste(:,k)],...
%         [RMS_dynamic_1_0_liste(:,k), RMS_dynamic_75_25_liste(:,k),RMS_dynamic_50_50_liste(:,k),RMS_dynamic_25_75_liste(:,k), RMS_dynamic_0_1_liste(:,k)]};
%     
    %data={[RMS_Insole_Direction(:,k),RMS_Threshold_Direction(:,k)]};
    data={[RMS_Insole_Run(:,k),RMS_Threshold_Run(:,k)]};
    %data={[RMS_Insole_Walk(:,k),RMS_Threshold_Walk(:,k)]};
    h=boxplotGroup(data, 'SecondaryLabels', {'II','KT'}, ...
    'PrimaryLabels',{'', ''}, 'InterGroupSpace', 1);
%     set(h.axis.Children(1).Children,'Color', 'b')
%     set(h.axis.Children(2).Children,'Color', 'r')
    title({[''],[ Directions{k}],[ ' component']})
    %xlabel(['Minimization paramters',newline])
    if k==2
        ylabel('RMSE (N)')
        
    end
        if k==4
        ylabel('RMSE (N.m)')
        
    end
end
box_vars = findall(gca,'Tag','Box');
hLegend = legend(box_vars([1,25]), {'Dynamic phase','Static phase'});

rRMSE=RMS;
for i= 1:length(rRMSE.RMS_Threshold_Walk)
    poids=meta(rRMSE.RMS_Threshold_Walk(i,7)).masse;
    rRMSE.RMS_Threshold_Walk(i,1:6)= rRMSE.RMS_Threshold_Walk(i,1:6)/poids;
end
for i= 1:length(rRMSE.RMS_Insole_Walk)
    poids=meta(rRMSE.RMS_Insole_Walk(i,7)).masse;
    rRMSE.RMS_Insole_Walk(i,1:6)= rRMSE.RMS_Insole_Walk(i,1:6)/poids;
end
for i= 1:length(rRMSE.RMS_Threshold_Run)
    poids=meta(rRMSE.RMS_Threshold_Run(i,7)).masse;
    rRMSE.RMS_Threshold_Run(i,1:6)= rRMSE.RMS_Threshold_Run(i,1:6)/poids;
end
for i= 1:length(rRMSE.RMS_Insole_Run)
    poids=meta(rRMSE.RMS_Insole_Run(i,7)).masse;
    rRMSE.RMS_Insole_Run(i,1:6)= rRMSE.RMS_Insole_Run(i,1:6)/poids;
end

  
figure
set(gcf,'color','w')
for k = 1:2
    subplot(3,2,k+2)
    data={[rRMSE.RMS_Insole_Run(:,1+(k-1)*3:3+(k-1)*3)],[rRMSE.RMS_Threshold_Run(:,1+(k-1)*3:3+(k-1)*3)]};
    
    if k==1
        h=boxplotGroup(data, 'SecondaryLabels', {'Fx','Fy','Fz'}, ...
            'PrimaryLabels',{'II','IK'}, 'InterGroupSpace', 1,'whisker',30);
        ylabel('rRMSE (N/kg)')
        %title({'Estimated Force RMSE','for run trials'})
        set(h.axis.Children(1).Children,'Color', 'b')
        set(h.axis.Children(2).Children,'Color', 'r')
        
axis([0 9 0 5])
       % yt = get(gca, 'YTick');
       yt = 310;
        yt = 3.5;
%axis([xlim    0  ceil(max(yt)*1.2)])
xt = get(gca, 'XTick');
hold on
%plot(xt([1 2]), [1 1]*max(yt)*1.1, '-k',  mean(xt([1 2])), max(yt)*1.15, '*k')
plot(xt([3 4]), [1 1]*max(yt)*1.1, '-k',  mean(xt([3 4])), max(yt)*1.15, '*k')
plot(xt([5 6]), [1 1]*max(yt)*1.1, '-k',  mean(xt([5 6])), max(yt)*1.15, '*k')
hold off
    else
        h=boxplotGroup(data, 'SecondaryLabels', {'Mx','My','Mz'}, ...
            'PrimaryLabels',{'II','IK'}, 'InterGroupSpace', 1,'whisker',30);
        ylabel('rRMSE (N.m/kg)')
        %title({'Estimated Moment RMSE','for run trials'})
        set(h.axis.Children(1).Children,'Color', 'b')
        set(h.axis.Children(2).Children,'Color', 'r')
        yt = 340;
         yt = 4.5;
        % yt = get(gca, 'YTick');
%axis([xlim    0  ceil(max(yt)*1.2)])
axis([xlim    0  6])
xt = get(gca, 'XTick');
hold on
plot(xt([1 2]), [1 1]*max(yt)*1.1, '-k',  mean(xt([1 2])), max(yt)*1.15, '*k')
plot(xt([3 4]), [1 1]*max(yt)*1.1, '-k',  mean(xt([3 4])), max(yt)*1.15, '*k')
plot(xt([5 6]), [1 1]*max(yt)*1.1, '-k',  mean(xt([5 6])), max(yt)*1.15, '*k')
hold off

    end
    
end

for k = 1:2
    subplot(3,2,k+4)
    data={[rRMSE.RMS_Insole_Walk(:,1+(k-1)*3:3+(k-1)*3)],[rRMSE.RMS_Threshold_Walk(:,1+(k-1)*3:3+(k-1)*3)]};
    
    if k==1
        h=boxplotGroup(data, 'SecondaryLabels', {'Fx','Fy','Fz'}, ...
            'PrimaryLabels',{'II','IK'}, 'InterGroupSpace', 1,'whisker',30);
        ylabel('rRMSE (N/kg)')
        %title({'Estimated Force RMSE','for walk trials'})
        set(h.axis.Children(1).Children,'Color', 'b')
        set(h.axis.Children(2).Children,'Color', 'r')
        
%axis([0 9 0 350])
        yt = get(gca, 'YTick');
         yt =3.5;
%axis([xlim    0  ceil(max(yt)*1.2)])
axis([0 9 0 5])
xt = get(gca, 'XTick');
hold on
plot(xt([1 2]), [1 1]*max(yt)*1.1, '-k',  mean(xt([1 2])), max(yt)*1.15, '*k')
plot(xt([3 4]), [1 1]*max(yt)*1.1, '-k',  mean(xt([3 4])), max(yt)*1.15, '*k')
plot(xt([5 6]), [1 1]*max(yt)*1.1, '-k',  mean(xt([5 6])), max(yt)*1.15, '*k')
hold off
    else
        h=boxplotGroup(data, 'SecondaryLabels', {'Mx','My','Mz'}, ...
            'PrimaryLabels',{'II','IK'}, 'InterGroupSpace', 1,'whisker',30);
        ylabel('rRMSE (N.m/kg)')
        %title({'Estimated Moment RMSE','for walk trials'})
        set(h.axis.Children(1).Children,'Color', 'b')
        set(h.axis.Children(2).Children,'Color', 'r')

        yt = 4.5;
axis([0 9 0 6])
  
xt = get(gca, 'XTick');
hold on
plot(xt([1 2]), [1 1]*max(yt)*1.1, '-k',  mean(xt([1 2])), max(yt)*1.15, '*k')
plot(xt([3 4]), [1 1]*max(yt)*1.1, '-k',  mean(xt([3 4])), max(yt)*1.15, '*k')
plot(xt([5 6]), [1 1]*max(yt)*1.1, '-k',  mean(xt([5 6])), max(yt)*1.15, '*k')
hold off

    end 
end
  
figure
set(gcf,'color','w')
for k = 1:2
    subplot(3,2,k)
    data={[RMS.RMS_Insole_Direction(:,1+(k-1)*3:3+(k-1)*3)],[RMS.RMS_Threshold_Direction(:,1+(k-1)*3:3+(k-1)*3)]};
    
    if k==1
        h=boxplotGroup(data, 'SecondaryLabels', {'Fx','Fy','Fz'}, ...
            'PrimaryLabels',{'II','IK'}, 'InterGroupSpace', 1,'whisker',30);
        ylabel('RMSE (N)')
        %title({'Estimated Force RMSE','for run with direction', 'change trials'})
        set(h.axis.Children(1).Children,'Color', 'b')
        set(h.axis.Children(2).Children,'Color', 'r')
        
%axis([0 9 0 300])
       % yt = get(gca, 'YTick');
         yt = 230;
          yt = 325;
axis([0 9 0 400])
         %axis([xlim    0  ceil(max(yt)*1.2)])
xt = get(gca, 'XTick');
hold on
plot(xt([1 2]), [1 1]*max(yt)*1.1, '-k',  mean(xt([1 2])), max(yt)*1.15, '*k')
plot(xt([3 4]), [1 1]*max(yt)*1.1, '-k',  mean(xt([3 4])), max(yt)*1.15, '*k')
plot(xt([5 6]), [1 1]*max(yt)*1.1, '-k',  mean(xt([5 6])), max(yt)*1.15, '*k')
hold off
    else
        h=boxplotGroup(data, 'SecondaryLabels', {'Mx','My','Mz'}, ...
            'PrimaryLabels',{'II','IK'}, 'InterGroupSpace', 1,'whisker',30);
        ylabel('RMSE (N.m)')
        %title({'Estimated Moment RMSE','for run with direction', 'change trials'})
        set(h.axis.Children(1).Children,'Color', 'b')
        set(h.axis.Children(2).Children,'Color', 'r')
      %  axis([0 9 0 350])
        yt = get(gca, 'YTick');
         yt = 325;
axis([0 9 0 400])
        %axis([xlim    0  ceil(max(yt)*1.2)])
xt = get(gca, 'XTick');
hold on
plot(xt([1 2]), [1 1]*max(yt)*1.1, '-k',  mean(xt([1 2])), max(yt)*1.15, '*k')
%plot(xt([3 4]), [1 1]*max(yt)*1.1, '-k',  mean(xt([3 4])), max(yt)*1.15, '*k')
plot(xt([5 6]), [1 1]*max(yt)*1.1, '-k',  mean(xt([5 6])), max(yt)*1.15, '*k') % p=0.024
hold off

    end
%     set(h.axis.Children(1).Children,'Color', 'b')
%     ho = findobj(h.axis.Children(1),'tag','Outliers');
%     set(ho,'MarkerEdgeColor','b')
%     set(h.axis.Children(2).Children,'Color', 'r')
    
    %title({[''],[ Directions{k}],[ ' component']})
    %xlabel(['Minimization paramters',newline])
    
end

%figure
%set(gcf,'color','w')
for k = 1:2
    subplot(3,2,k+2)
    data={[RMS.RMS_Insole_Run(:,1+(k-1)*3:3+(k-1)*3)],[RMS.RMS_Threshold_Run(:,1+(k-1)*3:3+(k-1)*3)]};
    
    if k==1
        h=boxplotGroup(data, 'SecondaryLabels', {'Fx','Fy','Fz'}, ...
            'PrimaryLabels',{'II','IK'}, 'InterGroupSpace', 1,'whisker',30);
        ylabel('RMSE (N)')
        %title({'Estimated Force RMSE','for run trials'})
        set(h.axis.Children(1).Children,'Color', 'b')
        set(h.axis.Children(2).Children,'Color', 'r')
        
axis([0 9 0 400])
       % yt = get(gca, 'YTick');
       yt = 310;
        yt = 300;
%axis([xlim    0  ceil(max(yt)*1.2)])
xt = get(gca, 'XTick');
hold on
%plot(xt([1 2]), [1 1]*max(yt)*1.1, '-k',  mean(xt([1 2])), max(yt)*1.15, '*k')
plot(xt([3 4]), [1 1]*max(yt)*1.1, '-k',  mean(xt([3 4])), max(yt)*1.15, '*k')
plot(xt([5 6]), [1 1]*max(yt)*1.1, '-k',  mean(xt([5 6])), max(yt)*1.15, '*k')
hold off
    else
        h=boxplotGroup(data, 'SecondaryLabels', {'Mx','My','Mz'}, ...
            'PrimaryLabels',{'II','IK'}, 'InterGroupSpace', 1,'whisker',30);
        ylabel('RMSE (N.m)')
        %title({'Estimated Moment RMSE','for run trials'})
        set(h.axis.Children(1).Children,'Color', 'b')
        set(h.axis.Children(2).Children,'Color', 'r')
        yt = 340;
         yt = 325;
        % yt = get(gca, 'YTick');
%axis([xlim    0  ceil(max(yt)*1.2)])
axis([xlim    0  400])
xt = get(gca, 'XTick');
hold on
plot(xt([1 2]), [1 1]*max(yt)*1.1, '-k',  mean(xt([1 2])), max(yt)*1.15, '*k')
plot(xt([3 4]), [1 1]*max(yt)*1.1, '-k',  mean(xt([3 4])), max(yt)*1.15, '*k')
plot(xt([5 6]), [1 1]*max(yt)*1.1, '-k',  mean(xt([5 6])), max(yt)*1.15, '*k')
hold off

    end
%     set(h.axis.Children(1).Children,'Color', 'b')
%     ho = findobj(h.axis.Children(1),'tag','Outliers');
%     set(ho,'MarkerEdgeColor','b')
%     set(h.axis.Children(2).Children,'Color', 'r')
    
    %title({[''],[ Directions{k}],[ ' component']})
    %xlabel(['Minimization paramters',newline])
    
end

%figure
%set(gcf,'color','w')
for k = 1:2
    subplot(3,2,k+4)
    data={[RMS.RMS_Insole_Walk(:,1+(k-1)*3:3+(k-1)*3)],[RMS.RMS_Threshold_Walk(:,1+(k-1)*3:3+(k-1)*3)]};
    
    if k==1
        h=boxplotGroup(data, 'SecondaryLabels', {'Fx','Fy','Fz'}, ...
            'PrimaryLabels',{'II','IK'}, 'InterGroupSpace', 1,'whisker',30);
        ylabel('RMSE (N)')
        %title({'Estimated Force RMSE','for walk trials'})
        set(h.axis.Children(1).Children,'Color', 'b')
        set(h.axis.Children(2).Children,'Color', 'r')
        
%axis([0 9 0 350])
        yt = get(gca, 'YTick');
         yt = 325;
%axis([xlim    0  ceil(max(yt)*1.2)])
axis([0 9 0 400])
xt = get(gca, 'XTick');
hold on
plot(xt([1 2]), [1 1]*max(yt)*1.1, '-k',  mean(xt([1 2])), max(yt)*1.15, '*k')
plot(xt([3 4]), [1 1]*max(yt)*1.1, '-k',  mean(xt([3 4])), max(yt)*1.15, '*k')
plot(xt([5 6]), [1 1]*max(yt)*1.1, '-k',  mean(xt([5 6])), max(yt)*1.15, '*k')
hold off
    else
        h=boxplotGroup(data, 'SecondaryLabels', {'Mx','My','Mz'}, ...
            'PrimaryLabels',{'II','IK'}, 'InterGroupSpace', 1,'whisker',30);
        ylabel('RMSE (N.m)')
        %title({'Estimated Moment RMSE','for walk trials'})
        set(h.axis.Children(1).Children,'Color', 'b')
        set(h.axis.Children(2).Children,'Color', 'r')
        %yt = get(gca, 'YTick');
        yt = 325;
axis([0 9 0 400])
        %axis([xlim    0  ceil(max(yt)*1.2)])
xt = get(gca, 'XTick');
hold on
plot(xt([1 2]), [1 1]*max(yt)*1.1, '-k',  mean(xt([1 2])), max(yt)*1.15, '*k')
plot(xt([3 4]), [1 1]*max(yt)*1.1, '-k',  mean(xt([3 4])), max(yt)*1.15, '*k')
plot(xt([5 6]), [1 1]*max(yt)*1.1, '-k',  mean(xt([5 6])), max(yt)*1.15, '*k')
hold off

    end
%     set(h.axis.Children(1).Children,'Color', 'b')
%     ho = findobj(h.axis.Children(1),'tag','Outliers');
%     set(ho,'MarkerEdgeColor','b')
%     set(h.axis.Children(2).Children,'Color', 'r')
    
    %title({[''],[ Directions{k}],[ ' component']})
    %xlabel(['Minimization paramters',newline])
    
end

plot(xt([3 4]), [1 1]*max(yt)*1.1, '-k',  mean(xt([3 4])), max(yt)*1.15, '*k')


yt = get(gca, 'YTick');
axis([xlim    0  ceil(max(yt)*1.2)])
xt = get(gca, 'XTick');
hold on
plot(xt([2 3]), [1 1]*max(yt)*1.1, '-k',  mean(xt([2 3])), max(yt)*1.15, '*k')
hold off

ho = findobj(h.axis.Children(1).Children,'tag','Outliers');
set(ho,'MarkerSize',34,'MarkerEdgeColor','b')
set(ho,'Color','b')




RMS_static=[RMS_static_1_0;
    RMS_static_75_25;
    RMS_static_50_50;
    RMS_static_25_75;
    RMS_static_0_1];



%% Calcul RMS


for ii=1:nb_s
    cd(meta(ii).path)
    name=meta(ii).name;
    RMS.name.Prediction=zeros(2,6);
    RMS.name.Prediction_plus=zeros(2,6);
    GFR_M_L=zeros(length(meta(ii).TrialList),6);
    GFR_M_L_plus=zeros(length(meta(ii).TrialList),6);
    GFR_M_L_prop_adaptatif=zeros(length(meta(ii).TrialList),6);
    GFR_M_L_prop_k_fixe=zeros(length(meta(ii).TrialList),6);
    
    GFR_M_R=zeros(length(meta(ii).TrialList),6);
    GFR_M_R_plus=zeros(length(meta(ii).TrialList),6);
    GFR_M_R_prop_adaptatif=zeros(length(meta(ii).TrialList),6);
    GFR_M_R_prop_k_fixe=zeros(length(meta(ii).TrialList),6);
    
    load('BiomechanicalModel.mat')
    for jj=1:length(meta(ii).TrialList)
        %if contains(meta(ii).TrialList(jj),'Marche') || contains(meta(ii).TrialList(jj),'Course')...
         %       || contains(meta(ii).TrialList_tasks(jj),'ChgtDirection') 
        %if contains(meta(ii).TrialList(jj),'Course')
           [rms,rms_plus, rms_prop_adaptatif, rms_prop_k_fixe] = RMSBetweenEventsOfInterest(BiomechanicalModel,...
                meta(ii).ExternalForces(jj).AnalysisParameters,...
                meta(ii).TrialList{jj}(1:end-4),meta(ii).FrameOfInterest(jj));
           GFR_M_L(jj,:)=rms(1,:);
           GFR_M_L_plus(jj,:)=rms_plus(1,:);
           GFR_M_L_prop_adaptatif(jj,:)=rms_prop_adaptatif(1,:);
           GFR_M_L_prop_k_fixe(jj,:)=rms_prop_k_fixe(1,:);
           GFR_M_R(jj,:)=rms(2,:);
           GFR_M_R_plus(jj,:)=rms_plus(2,:);
           GFR_M_R_prop_adaptatif(jj,:)=rms_prop_adaptatif(2,:);
           GFR_M_R_prop_k_fixe(jj,:)=rms_prop_k_fixe(2,:);
           
        %end
    end
    
    
    Course = [];
    ChgtDirection = [];
    Marche=[];
    Saut=[];
    for jj=1:length(meta(ii).TrialList)
        if contains(meta(ii).TrialList(jj),'Course')
            Course = [Course,jj];
        elseif contains(meta(ii).TrialList(jj),'Direction')
            ChgtDirection = [ChgtDirection,jj];
        elseif contains(meta(ii).TrialList(jj),'Marche')
            Marche = [Marche,jj];
        elseif contains(meta(ii).TrialList(jj),'Saut')
            Saut = [Saut,jj];
        end
    end
    
    
    meta(ii).RMS.Course=[mean(GFR_M_L(Course,:));mean(GFR_M_R(Course,:))];
    meta(ii).RMS_plus.Course=[mean(GFR_M_L_plus(Course,:));mean(GFR_M_R_plus(Course,:))];
    meta(ii).RMS_prop_adaptatif.Course=[mean(GFR_M_L_prop_adaptatif(Course,:));mean(GFR_M_R_prop_adaptatif(Course,:))];
    meta(ii).RMS_prop_k_fixe.Course=[mean(GFR_M_L_plus(Course,:));mean(GFR_M_R_plus(Course,:))];
    
    meta(ii).RMS.ChgtDirection=[mean(GFR_M_L(ChgtDirection,:));mean(GFR_M_R(ChgtDirection,:))];
    meta(ii).RMS_plus.ChgtDirection=[mean(GFR_M_L_plus(ChgtDirection,:));mean(GFR_M_R_plus(ChgtDirection,:))];
    meta(ii).RMS_prop_adaptatif.ChgtDirection=[mean(GFR_M_L_prop_adaptatif(ChgtDirection,:));mean(GFR_M_R_prop_adaptatif(ChgtDirection,:))];
    meta(ii).RMS_prop_k_fixe.ChgtDirection=[mean(GFR_M_L_plus(ChgtDirection,:));mean(GFR_M_R_plus(ChgtDirection,:))];
    
    meta(ii).RMS.Marche=[mean(GFR_M_L(Marche,:));mean(GFR_M_R(Marche,:))];
    meta(ii).RMS_plus.Marche=[mean(GFR_M_L_plus(Marche,:));mean(GFR_M_R_plus(Marche,:))];
    meta(ii).RMS_prop_adaptatif.Marche=[mean(GFR_M_L_prop_adaptatif(Marche,:));mean(GFR_M_R_prop_adaptatif(Marche,:))];
    meta(ii).RMS_prop_k_fixe.Marche=[mean(GFR_M_L_plus(Marche,:));mean(GFR_M_R_plus(Marche,:))];
    
    meta(ii).RMS.Saut=[mean(GFR_M_L(Saut,:));mean(GFR_M_R(Saut,:))];
    meta(ii).RMS_plus.Saut=[mean(GFR_M_L_plus(Saut,:));mean(GFR_M_R_plus(Saut,:))];
    meta(ii).RMS_prop_adaptatif.Saut=[mean(GFR_M_L_prop_adaptatif(Saut,:));mean(GFR_M_R_prop_adaptatif(Saut,:))];
    meta(ii).RMS_prop_k_fixe.Saut=[mean(GFR_M_L_plus(Saut,:));mean(GFR_M_R_plus(Saut,:))];
    
    
    meta(ii).STD.Course=[std([GFR_M_L(Course,:);GFR_M_R(Course,:)])].*(1/meta(ii).masse*ones(1,6));
    meta(ii).STD_plus.Course=[std([GFR_M_L_plus(Course,:);GFR_M_R_plus(Course,:)])].*(1/meta(ii).masse*ones(1,6));
    meta(ii).STD_prop_adaptatif.Course=[std([GFR_M_L_prop_adaptatif(Course,:);GFR_M_R_prop_adaptatif(Course,:)])].*(1/meta(ii).masse*ones(1,6));
    meta(ii).STD_prop_k_fixe.Course=[std([GFR_M_L_plus(Course,:);GFR_M_R_plus(Course,:)])].*(1/meta(ii).masse*ones(1,6));
    
    meta(ii).STD.ChgtDirection=[std([GFR_M_L(ChgtDirection,:);GFR_M_R(ChgtDirection,:)])].*(1/meta(ii).masse*ones(1,6));
    meta(ii).STD_plus.ChgtDirection=[std([GFR_M_L_plus(ChgtDirection,:);GFR_M_R_plus(ChgtDirection,:)])].*(1/meta(ii).masse*ones(1,6));
    meta(ii).STD_prop_adaptatif.ChgtDirection=[std([GFR_M_L_prop_adaptatif(ChgtDirection,:);GFR_M_R_prop_adaptatif(ChgtDirection,:)])].*(1/meta(ii).masse*ones(1,6));
    meta(ii).STD_prop_k_fixe.ChgtDirection=[std([GFR_M_L_plus(ChgtDirection,:);GFR_M_R_plus(ChgtDirection,:)])].*(1/meta(ii).masse*ones(1,6));
    
    meta(ii).STD.Marche=[std([GFR_M_L(Marche,:);GFR_M_R(Marche,:)])].*(1/meta(ii).masse*ones(1,6));
    meta(ii).STD_plus.Marche=[std([GFR_M_L_plus(Marche,:);GFR_M_R_plus(Marche,:)])].*(1/meta(ii).masse*ones(1,6));
    meta(ii).STD_prop_adaptatif.Marche=[std([GFR_M_L_prop_adaptatif(Marche,:);GFR_M_R_prop_adaptatif(Marche,:)])].*(1/meta(ii).masse*ones(1,6));
    meta(ii).STD_prop_k_fixe.Marche=[std([GFR_M_L_plus(Marche,:);GFR_M_R_plus(Marche,:)])].*(1/meta(ii).masse*ones(1,6));
    
    meta(ii).STD.Saut=[std([GFR_M_L(Saut,:);GFR_M_R(Saut,:)])].*(1/meta(ii).masse*ones(1,6));
    meta(ii).STD_plus.Saut=[std([GFR_M_L_plus(Saut,:);GFR_M_R_plus(Saut,:)])].*(1/meta(ii).masse*ones(1,6));
    meta(ii).STD_prop_adaptatif.Saut=[std([GFR_M_L_prop_adaptatif(Saut,:);GFR_M_R_prop_adaptatif(Saut,:)])].*(1/meta(ii).masse*ones(1,6));
    meta(ii).STD_prop_k_fixe.Saut=[std([GFR_M_L_plus(Saut,:);(GFR_M_R_plus(Saut,:))])].*(1/meta(ii).masse*ones(1,6));
    
    
    
    
    
%             RMS.name.Prediction_plus(1,:)=mean(GFR_M_L_plus);
%     RMS.name.Prediction_prop_adaptatif(1,:)=mean(GFR_M_L_prop_adaptatif);
%     RMS.name.Prediction_prop_k_fixe(1,:)=mean(GFR_M_L_prop_k_fixe);
%     
%     RMS.name.Prediction(2,:)=mean(GFR_M_R);
%     RMS.name.Prediction_plus(2,:)=mean(GFR_M_R_plus);
%     RMS.name.Prediction_prop_adaptatif(2,:)=mean(GFR_M_R_prop_adaptatif);
%     RMS.name.Prediction_prop_k_fixe(2,:)=mean(GFR_M_R_prop_k_fixe);
end


Course_RMS=[];
Course_RMS_plus=[];
Course_RMS_prop_adaptatif=[];
Course_RMS_prop_k_fixe=[];

for ii=1:nb_s
    Course_RMS=[Course_RMS;meta(ii).RMS.Course/meta(ii).masse];
    Course_RMS_plus=[Course_RMS_plus;meta(ii).RMS_plus.Course/meta(ii).masse];
    Course_RMS_prop_adaptatif=[Course_RMS_prop_adaptatif;meta(ii).RMS_prop_adaptatif.Course/meta(ii).masse];
    Course_RMS_prop_k_fixe=[Course_RMS_prop_k_fixe;meta(ii).RMS_prop_k_fixe.Course/meta(ii).masse];
end

% for ii=1:nb_s
%     Course_RMS=[Course_RMS;meta(ii).RMS.Course];
%     Course_RMS_plus=[Course_RMS_plus;meta(ii).RMS_plus.Course];
%     Course_RMS_prop_adaptatif=[Course_RMS_prop_adaptatif;meta(ii).RMS_prop_adaptatif.Course];
%     Course_RMS_prop_k_fixe=[Course_RMS_prop_k_fixe;meta(ii).RMS_prop_k_fixe.Course];
% end

Course_RMS=mean(Course_RMS)
Course_RMS_plus=mean(Course_RMS_plus)
Course_RMS_prop_adaptatif=mean(Course_RMS_prop_adaptatif)
Course_RMS_prop_k_fixe=mean(Course_RMS_prop_k_fixe)


ChgtDirection_RMS=[];
ChgtDirection_RMS_plus=[];
ChgtDirection_RMS_prop_adaptatif=[];
ChgtDirection_RMS_prop_k_fixe=[];

for ii=1:nb_s
    ChgtDirection_RMS=[ChgtDirection_RMS;meta(ii).RMS.ChgtDirection/meta(ii).masse];
    ChgtDirection_RMS_plus=[ChgtDirection_RMS_plus;meta(ii).RMS_plus.ChgtDirection/meta(ii).masse];
    ChgtDirection_RMS_prop_adaptatif=[ChgtDirection_RMS_prop_adaptatif;meta(ii).RMS_prop_adaptatif.ChgtDirection/meta(ii).masse];
    ChgtDirection_RMS_prop_k_fixe=[ChgtDirection_RMS_prop_k_fixe;meta(ii).RMS_prop_k_fixe.ChgtDirection/meta(ii).masse];
end

ChgtDirection_RMS=mean(ChgtDirection_RMS)
ChgtDirection_RMS_plus=mean(ChgtDirection_RMS_plus)
ChgtDirection_RMS_prop_adaptatif=mean(ChgtDirection_RMS_prop_adaptatif)
ChgtDirection_RMS_prop_k_fixe=mean(ChgtDirection_RMS_prop_k_fixe)

Marche_RMS=[];
Marche_RMS_plus=[];
Marche_RMS_prop_adaptatif=[];
Marche_RMS_prop_k_fixe=[];

for ii=1:nb_s
    Marche_RMS=[Marche_RMS;meta(ii).RMS.Marche/meta(ii).masse];
    Marche_RMS_plus=[Marche_RMS_plus;meta(ii).RMS_plus.Marche/meta(ii).masse];
    Marche_RMS_prop_adaptatif=[Marche_RMS_prop_adaptatif;meta(ii).RMS_prop_adaptatif.Marche/meta(ii).masse];
    Marche_RMS_prop_k_fixe=[Marche_RMS_prop_k_fixe;meta(ii).RMS_prop_k_fixe.Marche/meta(ii).masse];
end

Marche_RMS=mean(Marche_RMS)
Marche_RMS_plus=mean(Marche_RMS_plus)
Marche_RMS_prop_adaptatif=mean(Marche_RMS_prop_adaptatif)
Marche_RMS_prop_k_fixe=mean(Marche_RMS_prop_k_fixe)



Saut_RMS=[];
Saut_RMS_plus=[];
Saut_RMS_prop_adaptatif=[];
Saut_RMS_prop_k_fixe=[];

for ii=1:nb_s
    Saut_RMS=[Saut_RMS;meta(ii).RMS.Saut/meta(ii).masse];
    Saut_RMS_plus=[Saut_RMS_plus;meta(ii).RMS_plus.Saut/meta(ii).masse];
    Saut_RMS_prop_adaptatif=[Saut_RMS_prop_adaptatif;meta(ii).RMS_prop_adaptatif.Saut/meta(ii).masse];
    Saut_RMS_prop_k_fixe=[Saut_RMS_prop_k_fixe;meta(ii).RMS_prop_k_fixe.Saut/meta(ii).masse];
end

Saut_RMS=mean(Saut_RMS)
Saut_RMS_plus=mean(Saut_RMS_plus)
Saut_RMS_prop_adaptatif=mean(Saut_RMS_prop_adaptatif)
Saut_RMS_prop_k_fixe=mean(Saut_RMS_prop_k_fixe)




%% Calcul rRMS


for ii=1:nb_s
    cd(meta(ii).path)
    
    GFR_M_L=zeros(length(meta(ii).TrialList),6);
    GFR_M_L_plus=zeros(length(meta(ii).TrialList),6);
    GFR_M_L_prop_adaptatif=zeros(length(meta(ii).TrialList),6);
    GFR_M_L_prop_k_fixe=zeros(length(meta(ii).TrialList),6);
    
    GFR_M_R=zeros(length(meta(ii).TrialList),6);
    GFR_M_R_plus=zeros(length(meta(ii).TrialList),6);
    GFR_M_R_prop_adaptatif=zeros(length(meta(ii).TrialList),6);
    GFR_M_R_prop_k_fixe=zeros(length(meta(ii).TrialList),6);
    
    load('BiomechanicalModel.mat')
    for jj=1:length(meta(ii).TrialList)
        %if contains(meta(ii).TrialList(jj),'Marche') || contains(meta(ii).TrialList(jj),'Course')...
         %       || contains(meta(ii).TrialList_tasks(jj),'ChgtDirection') 
        %if contains(meta(ii).TrialList(jj),'Course')
           [rrms,rrms_plus, rrms_prop_adaptatif, rrms_prop_k_fixe] = rRMSBetweenEventsOfInterestPlus(BiomechanicalModel,...
                meta(ii).ExternalForces(jj).AnalysisParameters,...
                meta(ii).TrialList{jj}(1:end-4),meta(ii).FrameOfInterest(jj));
           GFR_M_L(jj,:)=rrms(1,:);
           GFR_M_L_plus(jj,:)=rrms_plus(1,:);
           GFR_M_L_prop_adaptatif(jj,:)=rrms_prop_adaptatif(1,:);
           GFR_M_L_prop_k_fixe(jj,:)=rrms_prop_k_fixe(1,:);
           GFR_M_R(jj,:)=rrms(2,:);
           GFR_M_R_plus(jj,:)=rrms_plus(2,:);
           GFR_M_R_prop_adaptatif(jj,:)=rrms_prop_adaptatif(2,:);
           GFR_M_R_prop_k_fixe(jj,:)=rrms_prop_k_fixe(2,:);
           
        %end
    end
    
    
    Course = [];
    ChgtDirection = [];
    Marche=[];
    Saut=[];
    for jj=1:length(meta(ii).TrialList)
        if contains(meta(ii).TrialList(jj),'Course')
            Course = [Course,jj];
        elseif contains(meta(ii).TrialList(jj),'Direction')
            ChgtDirection = [ChgtDirection,jj];
        elseif contains(meta(ii).TrialList(jj),'Marche')
            Marche = [Marche,jj];
        elseif contains(meta(ii).TrialList(jj),'Saut')
            Saut = [Saut,jj];
        end
    end
    
    
    meta(ii).rRMS.Course=[mean(GFR_M_L(Course,:));mean(GFR_M_R(Course,:))];
    meta(ii).rRMS_plus.Course=[mean(GFR_M_L_plus(Course,:));mean(GFR_M_R_plus(Course,:))];
    meta(ii).rRMS_prop_adaptatif.Course=[mean(GFR_M_L_prop_adaptatif(Course,:));mean(GFR_M_R_prop_adaptatif(Course,:))];
    meta(ii).rRMS_prop_k_fixe.Course=[mean(GFR_M_L_plus(Course,:));mean(GFR_M_R_plus(Course,:))];
    
    meta(ii).rRMS.ChgtDirection=[mean(GFR_M_L(ChgtDirection,:));mean(GFR_M_R(ChgtDirection,:))];
    meta(ii).rRMS_plus.ChgtDirection=[mean(GFR_M_L_plus(ChgtDirection,:));mean(GFR_M_R_plus(ChgtDirection,:))];
    meta(ii).rRMS_prop_adaptatif.ChgtDirection=[mean(GFR_M_L_prop_adaptatif(ChgtDirection,:));mean(GFR_M_R_prop_adaptatif(ChgtDirection,:))];
    meta(ii).rRMS_prop_k_fixe.ChgtDirection=[mean(GFR_M_L_plus(ChgtDirection,:));mean(GFR_M_R_plus(ChgtDirection,:))];
    
    meta(ii).rRMS.Marche=[mean(GFR_M_L(Marche,:));mean(GFR_M_R(Marche,:))];
    meta(ii).rRMS_plus.Marche=[mean(GFR_M_L_plus(Marche,:));mean(GFR_M_R_plus(Marche,:))];
    meta(ii).rRMS_prop_adaptatif.Marche=[mean(GFR_M_L_prop_adaptatif(Marche,:));mean(GFR_M_R_prop_adaptatif(Marche,:))];
    meta(ii).rRMS_prop_k_fixe.Marche=[mean(GFR_M_L_plus(Marche,:));mean(GFR_M_R_plus(Marche,:))];
    
    meta(ii).rRMS.Saut=[mean(GFR_M_L(Saut,:));mean(GFR_M_R(Saut,:))];
    meta(ii).rRMS_plus.Saut=[mean(GFR_M_L_plus(Saut,:));mean(GFR_M_R_plus(Saut,:))];
    meta(ii).rRMS_prop_adaptatif.Saut=[mean(GFR_M_L_prop_adaptatif(Saut,:));mean(GFR_M_R_prop_adaptatif(Saut,:))];
    meta(ii).rRMS_prop_k_fixe.Saut=[mean(GFR_M_L_plus(Saut,:));mean(GFR_M_R_plus(Saut,:))];
    
    
    
    
    
%             RMS.name.Prediction_plus(1,:)=mean(GFR_M_L_plus);
%     RMS.name.Prediction_prop_adaptatif(1,:)=mean(GFR_M_L_prop_adaptatif);
%     RMS.name.Prediction_prop_k_fixe(1,:)=mean(GFR_M_L_prop_k_fixe);
%     
%     RMS.name.Prediction(2,:)=mean(GFR_M_R);
%     RMS.name.Prediction_plus(2,:)=mean(GFR_M_R_plus);
%     RMS.name.Prediction_prop_adaptatif(2,:)=mean(GFR_M_R_prop_adaptatif);
%     RMS.name.Prediction_prop_k_fixe(2,:)=mean(GFR_M_R_prop_k_fixe);
end


Course_rRMS=[];
Course_rRMS_plus=[];
Course_rRMS_prop_adaptatif=[];
Course_rRMS_prop_k_fixe=[];

for ii=1:nb_s
    Course_rRMS=[Course_rRMS;meta(ii).rRMS.Course];
    Course_rRMS_plus=[Course_rRMS_plus;meta(ii).rRMS_plus.Course];
    Course_rRMS_prop_adaptatif=[Course_rRMS_prop_adaptatif;meta(ii).rRMS_prop_adaptatif.Course];
    Course_rRMS_prop_k_fixe=[Course_rRMS_prop_k_fixe;meta(ii).rRMS_prop_k_fixe.Course];
end

Course_rRMS=mean(Course_rRMS)*100
Course_rRMS_plus=mean(Course_rRMS_plus)*100
Course_rRMS_prop_adaptatif=mean(Course_rRMS_prop_adaptatif)*100
Course_rRMS_prop_k_fixe=mean(Course_rRMS_prop_k_fixe)*100



ChgtDirection_rRMS=[];
ChgtDirection_rRMS_plus=[];
ChgtDirection_rRMS_prop_adaptatif=[];
ChgtDirection_rRMS_prop_k_fixe=[];

for ii=1:nb_s
    ChgtDirection_rRMS=[ChgtDirection_rRMS;meta(ii).rRMS.ChgtDirection];
    ChgtDirection_rRMS_plus=[ChgtDirection_rRMS_plus;meta(ii).rRMS_plus.ChgtDirection];
    ChgtDirection_rRMS_prop_adaptatif=[ChgtDirection_rRMS_prop_adaptatif;meta(ii).rRMS_prop_adaptatif.ChgtDirection];
    ChgtDirection_rRMS_prop_k_fixe=[ChgtDirection_rRMS_prop_k_fixe;meta(ii).rRMS_prop_k_fixe.ChgtDirection];
end

ChgtDirection_rRMS=mean(ChgtDirection_rRMS)*100
ChgtDirection_rRMS_plus=mean(ChgtDirection_rRMS_plus)*100
ChgtDirection_rRMS_prop_adaptatif=mean(ChgtDirection_rRMS_prop_adaptatif)*100
ChgtDirection_rRMS_prop_k_fixe=mean(ChgtDirection_rRMS_prop_k_fixe)*100

Marche_rRMS=[];
Marche_rRMS_plus=[];
Marche_rRMS_prop_adaptatif=[];
Marche_rRMS_prop_k_fixe=[];

for ii=1:nb_s
    Marche_rRMS=[Marche_rRMS;meta(ii).rRMS.Marche];
    Marche_rRMS_plus=[Marche_rRMS_plus;meta(ii).rRMS_plus.Marche];
    Marche_rRMS_prop_adaptatif=[Marche_rRMS_prop_adaptatif;meta(ii).rRMS_prop_adaptatif.Marche];
    Marche_rRMS_prop_k_fixe=[Marche_rRMS_prop_k_fixe;meta(ii).rRMS_prop_k_fixe.Marche];
end

Marche_rRMS=mean(Marche_rRMS)*100
Marche_rRMS_plus=mean(Marche_rRMS_plus)*100
Marche_rRMS_prop_adaptatif=mean(Marche_rRMS_prop_adaptatif)*100
Marche_rRMS_prop_k_fixe=mean(Marche_rRMS_prop_k_fixe)*100



Saut_rRMS=[];
Saut_rRMS_plus=[];
Saut_rRMS_prop_adaptatif=[];
Saut_rRMS_prop_k_fixe=[];

for ii=1:nb_s
    Saut_rRMS=[Saut_rRMS;meta(ii).rRMS.Saut];
    Saut_rRMS_plus=[Saut_rRMS_plus;meta(ii).rRMS_plus.Saut];
    Saut_rRMS_prop_adaptatif=[Saut_rRMS_prop_adaptatif;meta(ii).rRMS_prop_adaptatif.Saut];
    Saut_rRMS_prop_k_fixe=[Saut_rRMS_prop_k_fixe;meta(ii).rRMS_prop_k_fixe.Saut];
end

Saut_rRMS=mean(Saut_rRMS)*100
Saut_rRMS_plus=mean(Saut_rRMS_plus)*100
Saut_rRMS_prop_adaptatif=mean(Saut_rRMS_prop_adaptatif)*100
Saut_rRMS_prop_k_fixe=mean(Saut_rRMS_prop_k_fixe)*100



%% Pearson correlation

for ii=1:nb_s
    cd(meta(ii).path)
    
    GFR_M_L=zeros(length(meta(ii).TrialList),6);
    GFR_M_L_plus=zeros(length(meta(ii).TrialList),6);
    GFR_M_L_prop_adaptatif=zeros(length(meta(ii).TrialList),6);
    GFR_M_L_prop_k_fixe=zeros(length(meta(ii).TrialList),6);
    
    GFR_M_R=zeros(length(meta(ii).TrialList),6);
    GFR_M_R_plus=zeros(length(meta(ii).TrialList),6);
    GFR_M_R_prop_adaptatif=zeros(length(meta(ii).TrialList),6);
    GFR_M_R_prop_k_fixe=zeros(length(meta(ii).TrialList),6);
    
    load('BiomechanicalModel.mat')
    for jj=1:length(meta(ii).TrialList)
        %if contains(meta(ii).TrialList(jj),'Marche') || contains(meta(ii).TrialList(jj),'Course')...
         %       || contains(meta(ii).TrialList_tasks(jj),'ChgtDirection') 
        %if contains(meta(ii).TrialList(jj),'Course')
           [pearson,pearson_plus, pearson_prop_adaptatif, pearson_prop_k_fixe] = PearsonBetweenEventsOfInterestPlus(BiomechanicalModel,...
                meta(ii).ExternalForces(jj).AnalysisParameters,...
                meta(ii).TrialList{jj}(1:end-4),meta(ii).FrameOfInterest(jj));
           GFR_M_L(jj,:)=pearson(1,:);
           GFR_M_L_plus(jj,:)=pearson(1,:);
           GFR_M_L_prop_adaptatif(jj,:)=pearson_prop_adaptatif(1,:);
           GFR_M_L_prop_k_fixe(jj,:)=pearson_prop_k_fixe(1,:);
           GFR_M_R(jj,:)=pearson(2,:);
           GFR_M_R_plus(jj,:)=pearson_plus(2,:);
           GFR_M_R_prop_adaptatif(jj,:)=pearson_prop_adaptatif(2,:);
           GFR_M_R_prop_k_fixe(jj,:)=pearson_prop_k_fixe(2,:);
           
        %end
    end
    
    
    Course = [];
    ChgtDirection = [];
    Marche=[];
    Saut=[];
    for jj=1:length(meta(ii).TrialList)
        if contains(meta(ii).TrialList(jj),'Course')
            Course = [Course,jj];
        elseif contains(meta(ii).TrialList(jj),'Direction')
            ChgtDirection = [ChgtDirection,jj];
        elseif contains(meta(ii).TrialList(jj),'Marche')
            Marche = [Marche,jj];
        elseif contains(meta(ii).TrialList(jj),'Saut')
            Saut = [Saut,jj];
        end
    end
    
    
    meta(ii).Pearson.Course=[mean(GFR_M_L(Course,:));mean(GFR_M_R(Course,:))];
    meta(ii).Pearson_plus.Course=[mean(GFR_M_L_plus(Course,:));mean(GFR_M_R_plus(Course,:))];
    meta(ii).Pearson_prop_adaptatif.Course=[mean(GFR_M_L_prop_adaptatif(Course,:));mean(GFR_M_R_prop_adaptatif(Course,:))];
    meta(ii).Pearson_prop_k_fixe.Course=[mean(GFR_M_L_plus(Course,:));mean(GFR_M_R_plus(Course,:))];
    
    meta(ii).Pearson.ChgtDirection=[mean(GFR_M_L(ChgtDirection,:));mean(GFR_M_R(ChgtDirection,:))];
    meta(ii).Pearson_plus.ChgtDirection=[mean(GFR_M_L_plus(ChgtDirection,:));mean(GFR_M_R_plus(ChgtDirection,:))];
    meta(ii).Pearson_prop_adaptatif.ChgtDirection=[mean(GFR_M_L_prop_adaptatif(ChgtDirection,:));mean(GFR_M_R_prop_adaptatif(ChgtDirection,:))];
    meta(ii).Pearson_prop_k_fixe.ChgtDirection=[mean(GFR_M_L_plus(ChgtDirection,:));mean(GFR_M_R_plus(ChgtDirection,:))];
    
    meta(ii).Pearson.Marche=[mean(GFR_M_L(Marche,:));mean(GFR_M_R(Marche,:))];
    meta(ii).Pearson_plus.Marche=[mean(GFR_M_L_plus(Marche,:));mean(GFR_M_R_plus(Marche,:))];
    meta(ii).Pearson_prop_adaptatif.Marche=[mean(GFR_M_L_prop_adaptatif(Marche,:));mean(GFR_M_R_prop_adaptatif(Marche,:))];
    meta(ii).Pearson_prop_k_fixe.Marche=[mean(GFR_M_L_plus(Marche,:));mean(GFR_M_R_plus(Marche,:))];
    
    meta(ii).Pearson.Saut=[mean(GFR_M_L(Saut,:));mean(GFR_M_R(Saut,:))];
    meta(ii).Pearson_plus.Saut=[mean(GFR_M_L_plus(Saut,:));mean(GFR_M_R_plus(Saut,:))];
    meta(ii).Pearson_prop_adaptatif.Saut=[mean(GFR_M_L_prop_adaptatif(Saut,:));mean(GFR_M_R_prop_adaptatif(Saut,:))];
    meta(ii).Pearson_prop_k_fixe.Saut=[mean(GFR_M_L_plus(Saut,:));mean(GFR_M_R_plus(Saut,:))];
    
    
    
    
    meta(ii).STD_Pearson.Course=[std([GFR_M_L(Course,:);GFR_M_R(Course,:)])]
    meta(ii).STD_Pearson_plus.Course=[std([GFR_M_L_plus(Course,:);GFR_M_R_plus(Course,:)])]
    meta(ii).STD_Pearson_prop_adaptatif.Course=[std([GFR_M_L_prop_adaptatif(Course,:);GFR_M_R_prop_adaptatif(Course,:)])];
    meta(ii).STD_Pearson_prop_k_fixe.Course=[std([GFR_M_L_plus(Course,:);GFR_M_R_plus(Course,:)])];
    
    meta(ii).STD_Pearson.ChgtDirection=[std([GFR_M_L(ChgtDirection,:);GFR_M_R(ChgtDirection,:)])];
    meta(ii).STD_Pearson_plus.ChgtDirection=[std([GFR_M_L_plus(ChgtDirection,:);GFR_M_R_plus(ChgtDirection,:)])];
    meta(ii).STD_Pearson_prop_adaptatif.ChgtDirection=[std([GFR_M_L_prop_adaptatif(ChgtDirection,:);GFR_M_R_prop_adaptatif(ChgtDirection,:)])];
    meta(ii).STD_Pearson_prop_k_fixe.ChgtDirection=[std([GFR_M_L_plus(ChgtDirection,:);GFR_M_R_plus(ChgtDirection,:)])];
    
    meta(ii).STD_Pearson_.Marche=[std([GFR_M_L(Marche,:);GFR_M_R(Marche,:)])];
    meta(ii).STD_Pearson_plus.Marche=[std([GFR_M_L_plus(Marche,:);GFR_M_R_plus(Marche,:)])];
    meta(ii).STD_Pearson_prop_adaptatif.Marche=[std([GFR_M_L_prop_adaptatif(Marche,:);GFR_M_R_prop_adaptatif(Marche,:)])];
    meta(ii).STD_Pearson_prop_k_fixe.Marche=[std([GFR_M_L_plus(Marche,:);GFR_M_R_plus(Marche,:)])];
    
    meta(ii).STD_Pearson.Saut=[std([GFR_M_L(Saut,:);GFR_M_R(Saut,:)])];
    meta(ii).STD_Pearson_plus.Saut=[std([GFR_M_L_plus(Saut,:);GFR_M_R_plus(Saut,:)])];
    meta(ii).STD_Pearson_prop_adaptatif.Saut=[std([GFR_M_L_prop_adaptatif(Saut,:);GFR_M_R_prop_adaptatif(Saut,:)])];
    meta(ii).STD_Pearson_prop_k_fixe.Saut=[std([GFR_M_L_plus(Saut,:);(GFR_M_R_plus(Saut,:))])];
    
    
    
    
    
%             RMS.name.Prediction_plus(1,:)=mean(GFR_M_L_plus);
%     RMS.name.Prediction_prop_adaptatif(1,:)=mean(GFR_M_L_prop_adaptatif);
%     RMS.name.Prediction_prop_k_fixe(1,:)=mean(GFR_M_L_prop_k_fixe);
%     
%     RMS.name.Prediction(2,:)=mean(GFR_M_R);
%     RMS.name.Prediction_plus(2,:)=mean(GFR_M_R_plus);
%     RMS.name.Prediction_prop_adaptatif(2,:)=mean(GFR_M_R_prop_adaptatif);
%     RMS.name.Prediction_prop_k_fixe(2,:)=mean(GFR_M_R_prop_k_fixe);
end


Course_Pearson=[];
Course_Pearson_plus=[];
Course_Pearson_prop_adaptatif=[];
Course_Pearson_prop_k_fixe=[];

for ii=1:nb_s
    Course_Pearson=[Course_Pearson;meta(ii).Pearson.Course];
    Course_Pearson_plus=[Course_Pearson_plus;meta(ii).Pearson_plus.Course];
    Course_Pearson_prop_adaptatif=[Course_Pearson_prop_adaptatif;meta(ii).Pearson_prop_adaptatif.Course];
    Course_Pearson_prop_k_fixe=[Course_Pearson_prop_k_fixe;meta(ii).Pearson_prop_k_fixe.Course];
end

% for ii=1:nb_s
%     Course_RMS=[Course_RMS;meta(ii).RMS.Course];
%     Course_RMS_plus=[Course_RMS_plus;meta(ii).RMS_plus.Course];
%     Course_RMS_prop_adaptatif=[Course_RMS_prop_adaptatif;meta(ii).RMS_prop_adaptatif.Course];
%     Course_RMS_prop_k_fixe=[Course_RMS_prop_k_fixe;meta(ii).RMS_prop_k_fixe.Course];
% end

Course_Pearson=mean(Course_Pearson)
Course_Pearson_plus=mean(Course_Pearson_plus)
Course_Pearson_prop_adaptatif=mean(Course_Pearson_prop_adaptatif)
Course_Pearson_prop_k_fixe=mean(Course_Pearson_prop_k_fixe)


ChgtDirection_Pearson=[];
ChgtDirection_Pearson_plus=[];
ChgtDirection_Pearson_prop_adaptatif=[];
ChgtDirection_Pearson_prop_k_fixe=[];

for ii=1:nb_s
    ChgtDirection_Pearson=[ChgtDirection_Pearson;meta(ii).Pearson.ChgtDirection];
    ChgtDirection_Pearson_plus=[ChgtDirection_Pearson_plus;meta(ii).Pearson_plus.ChgtDirection];
    ChgtDirection_Pearson_prop_adaptatif=[ChgtDirection_Pearson_prop_adaptatif;meta(ii).Pearson_prop_adaptatif.ChgtDirection];
    ChgtDirection_Pearson_prop_k_fixe=[ChgtDirection_Pearson_prop_k_fixe;meta(ii).Pearson_prop_k_fixe.ChgtDirection];
end

ChgtDirection_Pearson=mean(ChgtDirection_Pearson)
ChgtDirection_Pearson_plus=mean(ChgtDirection_Pearson_plus)
ChgtDirection_Pearson_prop_adaptatif=mean(ChgtDirection_Pearson_prop_adaptatif)
ChgtDirection_Pearson_prop_k_fixe=mean(ChgtDirection_Pearson_prop_k_fixe)

Marche_Pearson=[];
Marche_Pearson_plus=[];
Marche_Pearson_prop_adaptatif=[];
Marche_Pearson_prop_k_fixe=[];

for ii=1:nb_s
    Marche_Pearson=[Marche_Pearson;meta(ii).Pearson.Marche];
    Marche_Pearson_plus=[Marche_Pearson_plus;meta(ii).Pearson_plus.Marche];
    Marche_Pearson_prop_adaptatif=[Marche_Pearson_prop_adaptatif;meta(ii).Pearson_prop_adaptatif.Marche];
    Marche_Pearson_prop_k_fixe=[Marche_Pearson_prop_k_fixe;meta(ii).Pearson_prop_k_fixe.Marche];
end

Marche_Pearson=mean(Marche_Pearson)
Marche_Pearson_plus=mean(Marche_Pearson_plus)
Marche_Pearson_prop_adaptatif=mean(Marche_Pearson_prop_adaptatif)
Marche_Pearson_prop_k_fixe=mean(Marche_Pearson_prop_k_fixe)



Saut_Pearson=[];
Saut_Pearson_plus=[];
Saut_Pearson_prop_adaptatif=[];
Saut_Pearson_prop_k_fixe=[];

for ii=1:nb_s
    Saut_Pearson=[Saut_Pearson;meta(ii).Pearson.Saut];
    Saut_Pearson_plus=[Saut_Pearson_plus;meta(ii).Pearson_plus.Saut];
    Saut_Pearson_prop_adaptatif=[Saut_Pearson_prop_adaptatif;meta(ii).Pearson_prop_adaptatif.Saut];
    Saut_Pearson_prop_k_fixe=[Saut_Pearson_prop_k_fixe;meta(ii).Pearson_prop_k_fixe.Saut];
end

Saut_Pearson=mean(Saut_Pearson)
Saut_Pearson_plus=mean(Saut_Pearson_plus)
Saut_Pearson_prop_adaptatif=mean(Saut_Pearson_prop_adaptatif)
Saut_Pearson_prop_k_fixe=mean(Saut_Pearson_prop_k_fixe)







for ii=1:nb_s
    cd(meta(ii).path)
    for jj=1:length(meta(ii).TrialList)
        load([meta(ii).TrialList{jj}(1:end-4) '/InverseKinematicsResults']);
        disp([meta(ii).TrialList{jj}(1:end-4) ' '...
            num2str(mean(mean(InverseKinematicsResults.ReconstructionError)))])
    end
end