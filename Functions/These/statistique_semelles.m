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

Precision_Semelles=[];
Precision_Semelles_reduit=[];

for ii=1:nb_s
    cd(fullfile(cheminglobal,SUJETS{ii}));
    load('matlab_Results.mat')
    Suj = [];
    L=[];
    for j=1:size(Results,2)
        Results(j).Sujet = SUJETS{ii};
        if contains(Results(j).filename_type, Results(j).filename_calib_type) && size(Results(j).filename_type,2)==size(Results(j).filename_calib_type,2)
            L=[L;Results(j)];
        end
    end
    T = struct2table(Results);
    Precision_Semelles=[Precision_Semelles; T];
    
    L = struct2table(L); 
    Precision_Semelles_reduit=[Precision_Semelles_reduit; L];
end

%writetable(Precision_Semelles,'Precision_Semelles_4.xlsx')
%writetable(Precision_Semelles_reduit,'Precision_Semelles_calib_par_cat.xlsx')


T = Precision_Semelles;
%T = struct2table(Precision_Semelles);
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


figure(1)
set(gca,'color','w')
set(gcf,'color','w')
%group = T.filename_calib_type.*T.filename_type;
grouptext = {'Sidecut Maneuver','Run','Walk'};
group = T.filename_calib_type;
Calib_type_order = {'Sidecut Maneuver','','','','','Run','','','','','Walk','','','','','Sidecut Maneuver','','','','','Run','','','','','Walk','','','','','Sidecut Maneuver','','','','','Run','','','','','Walk','','','',''};

boxchart(group,T.RMSE_Insole_x,'GroupByColor',T.filename_type)

%
%boxchart(group,T.RMSE_Insole_x,'GroupByColor',T.filename_calib_type)
%title('Erreur sur la position du centre de pression')
%boxchart(T.filename_calib,T.RMSE_Insole_x)
ylabel('RMSE AM (m)')
xlabel('Calibration Trial Type')
xticklabels(grouptext)
legend('Sidecut Maneuver','Run','Walk')
% axe=gca;
% Unique_name_groupe=unique(group,'stable');
% x.Xaxis.Categories={Unique_name_groupe(1),'1',Unique_name_groupe(2),'2',Unique_name_groupe(3)};
% arrayfun(@(x)xline(x,'r-','LineWidth',1.5),categorical({'1','2'}))
% ax.XTickLabel([2,4]) = {''};
% 

figure(2)
set(gca,'color','w')
set(gcf,'color','w')
A=boxchart(group,T.RMSE_Insole_y,'GroupByColor',T.filename_type);

%boxchart(group,T.RMSE_Insole_y,'GroupByColor',T.filename_calib_type)
%title('Erreur sur la position du centre de pression')
%boxchart(T.filename_calib,T.RMSE_Insole_x)
ylabel('RMSE AP (m)')
xlabel('Calibration Trial Type')
legend('Sidecut Maneuver','Run','Walk')
xticklabels(grouptext)

%% Etude selon la force
% clear all
% close all
clc
cheminglobal=pwd;
cd(cheminglobal);


doss=dir(cheminglobal);
doss={doss([doss.isdir]').name}';
doss(1:2)=[];

SUJETS=doss;
nb_s=length(SUJETS);

Precision_Semelles=[];
Precision_Semelles_reduit=[];

for ii=1:nb_s
    cd(fullfile(cheminglobal,SUJETS{ii}));
    load('matlab_Results_liste.mat')
    Suj = [];
    L=[];
    for j=1:size(Results,2)
        Results(j).Sujet = SUJETS{ii};
        if contains(Results(j).filename, Results(j).filename_calib) && size(Results(j).filename,2)==size(Results(j).filename_calib,2)
        %if contains(Results(j).filename_type, Results(j).filename_calib_type) && size(Results(j).filename_type,2)==size(Results(j).filename_calib_type,2)
            L=[L;Results(j)];
        end
    end
    T = struct2table(Results);
    Precision_Semelles=[Precision_Semelles; T];
    
    L = struct2table(L); 
    Precision_Semelles_reduit=[Precision_Semelles_reduit; L];
end

Precision_force=struct();
Nb_reduit = size(Precision_Semelles_reduit,1);
Fz=[];
Fz_round=[];
RMSE_x=[];
RMSE_y=[];
Contact=[];
for i=1:Nb_reduit
    Fz=[Fz;Precision_Semelles_reduit.Fz_Insole_liste{i,1}];
    Fz_round=[Fz_round;round(Precision_Semelles_reduit.Fz_Insole_liste{i,1}*2,-2)/2];
    RMSE_x=[RMSE_x;abs(Precision_Semelles_reduit.RMSE_Insole_liste{i,1}(:,1))];
    RMSE_y=[RMSE_y;abs(Precision_Semelles_reduit.RMSE_Insole_liste{i,1}(:,2))];
    Contact=[Contact;Precision_Semelles_reduit.Contact_liste{i,1}];
end
Precision_force.Fz=Fz;
Precision_force.Fz_round=Fz_round;
Precision_force.RMSE_x=RMSE_x;
Precision_force.RMSE_y=RMSE_y;
Precision_force.Contact=Contact;

Precision_force = struct2table(Precision_force); 


%writetable(Precision_force,'Precision_force_6.xlsx')

figure(3)
set(gca,'color','w')
set(gcf,'color','w')

