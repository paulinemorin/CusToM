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
        %if contains(Results(j).filename_type ,{'Walk'})
        %if contains(Results(j).filename_type ,{'Run'})
        %if contains(Results(j).filename_type ,{'ChangmntDirection'})
        if contains(Results(j).filename, Results(j).filename_calib) && size(Results(j).filename,2)==size(Results(j).filename_calib,2)
        else%if %contains(Results(j).filename_type ,{'Walk'})
            L=[L;Results(j)];
        end
    end
    T = struct2table(Results);
    Precision_Semelles=[Precision_Semelles; T];
    
    L = struct2table(L); 
    Precision_Semelles_reduit=[Precision_Semelles_reduit; L];
end
Run_x=[];
Dir_x=[];
Walk_x=[];

Run_y=[];
Dir_y=[];
Walk_y=[];
for i=1:size(Precision_Semelles_reduit,1)
    if contains(Precision_Semelles_reduit.filename_calib_type(i),'ChangmntDirection')
        Dir_x=[Dir_x;Precision_Semelles_reduit.RMSE_Insole_x(i)];
        Dir_y=[Dir_y;Precision_Semelles_reduit.RMSE_Insole_y(i)];
    elseif contains(Precision_Semelles_reduit.filename_calib_type(i),'Run')
        Run_x=[Run_x;Precision_Semelles_reduit.RMSE_Insole_x(i)];
        Run_y=[Run_y;Precision_Semelles_reduit.RMSE_Insole_y(i)];
    elseif contains(Precision_Semelles_reduit.filename_calib_type(i),'Walk')
        Walk_x=[Walk_x;Precision_Semelles_reduit.RMSE_Insole_x(i)];
        Walk_y=[Walk_y;Precision_Semelles_reduit.RMSE_Insole_y(i)];
    end
end
    
RMSE_x_tot_Friedman={Dir_x,Walk_x,Run_x};
RMSE_y_tot_Friedman={Dir_y,Walk_y,Run_y};

%writetable(Precision_Semelles_reduit,'Precision_Semelles_calib_self_2.xlsx')
%writetable(Precision_Semelles,'Precision_Semelles_4.xlsx')
%writetable(Precision_Semelles_reduit,'Precision_Semelles_dirdction_sans_self.xlsx')
%% Claire
Calib_type_order = {'ChangmntDirection','1','Run','2','Walk'};
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






figure(3)
set(gca,'color','w')
set(gcf,'color','w')
%group = T.filename_calib_type.*T.filename_type;
grouptext = {'Sidecut Maneuver','1','Run','2','Walk'};
group =T.filename_calib_type;



Calib_type_order = {'Sidecut Maneuver','','','','','Run','','','','','Walk','','','','','Sidecut Maneuver','','','','','Run','','','','','Walk','','','','','Sidecut Maneuver','','','','','Run','','','','','Walk','','','',''};





boxchart(group,T.RMSE_Insole_x,'GroupByColor',T.filename_type)





%
%boxchart(group,T.RMSE_Insole_x,'GroupByColor',T.filename_calib_type)
%title('Erreur sur la position du centre de pression')
%boxchart(T.filename_calib,T.RMSE_Insole_x)
ylabel('RMSE AM (m)')
xlabel('Calibration Trial Type')
xticklabels(grouptext)
axe=gca;
% Unique_name_groupe=unique(group,'stable');
% x.Xaxis.Categories={Unique_name_groupe(1),'1',Unique_name_groupe(2),'2',Unique_name_groupe(3)};
arrayfun(@(x)xline(x,'r-','LineWidth',1.5),categorical({'1','2 '}))
axe.XTickLabel([2,4]) = {''};
%



legend('Sidecut Maneuver','Run','Walk')

%% debut
%T = Precision_Semelles;
T = Precision_Semelles_reduit;
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


%boxchart(group,T.RMSE_Insole_x,'GroupByColor',T.filename_calib_type)
%title('Erreur sur la position du centre de pression')
%boxchart(T.filename_calib,T.RMSE_Insole_x)
ylabel('RMSE AM (m)')
xlabel('Calibration Trial Type')
xticklabels(grouptext)
legend('Sidecut Maneuver','Run','Walk')


figure(4)
set(gca,'color','w')
set(gcf,'color','w')

Y1 = NaN(size(T,1),3);
Y2=NaN(size(T,1),3);
Y3=NaN(size(T,1),3);
% for i=1:size(T,1)
%     if T.filename_calib_type(i) =={'ChangmntDirection'}
%         j=1;
%     elseif T.filename_calib_type(i) =={'Run'}
%         j=2;
%     elseif T.filename_calib_type(i) =={'Walk'}
%         j=3;
%     end
%     if T.filename_type(i) =={'ChangmntDirection'}
%         Y1(i,j) = T.RMSE_Insole_x(i);
%     elseif T.filename_type(i) =={'Run'}
%         Y2(i,j)=T.RMSE_Insole_x(i);
%     elseif T.filename_type(i) =={'Walk'}
%         Y3(i,j)=T.RMSE_Insole_x(i);
%     end
% end

for i=1:size(T,1)
    if T.filename_type(i) =={'ChangmntDirection'}
        j=1;
    elseif T.filename_type(i) =={'Run'}
        j=2;
    elseif T.filename_type(i) =={'Walk'}
        j=3;
    end
    if T.filename_calib_type(i) =={'ChangmntDirection'}
        Y1(i,j) = T.RMSE_Insole_x(i);
    elseif T.filename_calib_type(i) =={'Run'}
        Y2(i,j)=T.RMSE_Insole_x(i);
    elseif T.filename_calib_type(i) =={'Walk'}
        Y3(i,j)=T.RMSE_Insole_x(i);
    end
end

A= [0.171 0.195 0.201]/norm([0.171 0.195 0.201]);
color={[0.171 0.195 0.201],[0.244 0.212 0.212],[0.204 0.190 0.159]};
color={'#601A4A','#EE442F','#63ACBE'};
color={[0.171 0.195 0.201]*4,[0.244*2 0.212*2 0.212*2]*1.5,[0.204*2 0.190*2 0.159*2]*2};
boxplot_custom(Y1,Y2,Y3,grouptext,'mode',3,'list_legends',grouptext,'outlier_multiplier',Inf,'fillcolor',color)
%boxchart(group,T.RMSE_Insole_x,'GroupByColor',T.filename_type) 'fillcolor','b'

ylabel('RMSE ML (m)')
xlabel('Trial Type')
xticklabels(grouptext)
%legend('Sidecut Maneuver','Run','Walk')


yt = get(gca, 'YTick');
axis([xlim    0  ceil(max(yt)*1.2)])
xt = get(gca, 'XTick');
hold on
%plot(xt([1 2]), [1 1]*max(yt)*1.1, '-k',  mean(xt([1 2])), max(yt)*1.15, '*k')
plot([0.785 1], [1 1]*max(yt)*1, '-k',  [0.8625 0.8925 0.9275], [1 1 1]*max(yt)*1.05, '*k') % mean([0.785 1]), max(yt)*0.85, '*k')
plot([0.785 1.215], [1 1]*max(yt)*0.9, '-k', [0.965 1 1.035], [1 1 1]*max(yt)*0.95, '*k')% mean([0.785 1.215]), max(yt)*0.95, '*k')

plot([1.785 2], [1 1]*max(yt)*0.9, '-k',  [1.8625 1.8925 1.9275], [1 1 1]*max(yt)*0.95, '*k') %mean([1.785 2]), max(yt)*1.05, '*k')
plot([1.785 2.215], [1 1]*max(yt)*0.8, '-k', [1.965 2 2.035], [1 1 1]*max(yt)*0.85, '*k')
plot([2 2.215], [1 1]*max(yt)*0.7, '-k',  [2.1075], [1]*max(yt)*0.75, '*k') %mean([2 2.215]), max(yt)*0.95, '*k')

plot([2.785 3], [1 1]*max(yt)*1.0, '-k', [2.8625 2.8925 2.9275], [1 1 1]*max(yt)*1.05, '*k') % mean([2.785 3]), max(yt)*0.85, '*k')
plot([2.785 3.215], [1 1]*max(yt)*0.8, '-k', [2.965 3 3.035], [1 1 1]*max(yt)*0.85, '*k')% mean([2.785 3.215]), max(yt)*0.95, '*k')
plot([3 3.215], [1 1]*max(yt)*0.9, '-k', [3.0775 3.1075 3.1375], [1 1 1]*max(yt)*0.95, '*k') % mean([3 3.215]), max(yt)*0.75, '*k')

% plot([0.785 1], [1 1]*max(yt)*1, '-k',  [0.8625 0.8925 0.9275], [1 1 1]*max(yt)*1.05, '*k') % mean([0.785 1]), max(yt)*0.85, '*k')
% plot([0.785 1.215], [1 1]*max(yt)*0.9, '-k', [0.965 1 1.035], [1 1 1]*max(yt)*0.95, '*k')% mean([0.785 1.215]), max(yt)*0.95, '*k')
% 
% plot([1.785 2], [1 1]*max(yt)*1, '-k',  [1.8625 1.8925 1.9275], [1 1 1]*max(yt)*1.05, '*k') %mean([1.785 2]), max(yt)*1.05, '*k')
% plot([1.785 2.215], [1 1]*max(yt)*1.1, '-k', [1.965 2 2.035], [1 1 1]*max(yt)*1.15, '*k')
% plot([2 2.215], [1 1]*max(yt)*0.9, '-k',  [2.0775 2.1075 2.1375], [1 1 1]*max(yt)*0.95, '*k') %mean([2 2.215]), max(yt)*0.95, '*k')
% 
% plot([2.785 3], [1 1]*max(yt)*0.8, '-k', [2.8625 2.8925 2.9275], [1 1 1]*max(yt)*0.85, '*k') % mean([2.785 3]), max(yt)*0.85, '*k')
% plot([2.785 3.215], [1 1]*max(yt)*0.9, '-k', [2.965 3 3.035], [1 1 1]*max(yt)*0.95, '*k')% mean([2.785 3.215]), max(yt)*0.95, '*k')
% plot([3 3.215], [1 1]*max(yt)*0.7, '-k', [3.0775 3.1075 3.1375], [1 1 1]*max(yt)*0.75, '*k') % mean([3 3.215]), max(yt)*0.75, '*k')


figure(5)
set(gca,'color','w')
set(gcf,'color','w')

Y1 = NaN(size(T,1),3);
Y2=NaN(size(T,1),3);
Y3=NaN(size(T,1),3);
for i=1:size(T,1)
    if T.filename_type(i) =={'ChangmntDirection'}
        j=1;
    elseif T.filename_type(i) =={'Run'}
        j=2;
    elseif T.filename_type(i) =={'Walk'}
        j=3;
    end
    if T.filename_calib_type(i) =={'ChangmntDirection'}
        Y1(i,j) = T.RMSE_Insole_y(i);
    elseif T.filename_calib_type(i) =={'Run'}
        Y2(i,j)=T.RMSE_Insole_y(i);
    elseif T.filename_calib_type(i) =={'Walk'}
        Y3(i,j)=T.RMSE_Insole_y(i);
    end
end
% 
% for i=1:size(T,1)
%     if T.filename_calib_type(i) =={'ChangmntDirection'}
%         j=1;
%     elseif T.filename_calib_type(i) =={'Run'}
%         j=2;
%     elseif T.filename_calib_type(i) =={'Walk'}
%         j=3;
%     end
%     if T.filename_type(i) =={'ChangmntDirection'}
%         Y1(i,j) = T.RMSE_Insole_y(i);
%     elseif T.filename_type(i) =={'Run'}
%         Y2(i,j)=T.RMSE_Insole_y(i);
%     elseif T.filename_type(i) =={'Walk'}
%         Y3(i,j)=T.RMSE_Insole_y(i);
%     end
% end

A= [0.171 0.195 0.201]/norm([0.171 0.195 0.201]);
color={[0.171 0.195 0.201],[0.244 0.212 0.212],[0.204 0.190 0.159]};
color={'#601A4A','#EE442F','#63ACBE'};
color={[0.171 0.195 0.201]*4,[0.244*2 0.212*2 0.212*2]*1.5,[0.204*2 0.190*2 0.159*2]*2};
boxplot_custom(Y1,Y2,Y3,grouptext,'mode',3,'list_legends',grouptext,'outlier_multiplier',Inf,'fillcolor',color)
yt = get(gca, 'YTick');
axis([xlim    0  ceil(max(yt)*1.2)])
xt = get(gca, 'XTick');
hold on

plot([1 1.215], [1 1]*max(yt)*0.8, '-k',  [1.0775 1.1075 1.1375], [1 1 1]*max(yt)*0.85, '*k') % mean([0.785 1]), max(yt)*0.85, '*k')
plot([0.785 1.215], [1 1]*max(yt)*0.9, '-k', [0.965 1 1.035], [1 1 1]*max(yt)*0.95, '*k')% mean([0.785 1.215]), max(yt)*0.95, '*k')

plot([1.785 2], [1 1]*max(yt)*1, '-k',  [1.8925], [1]*max(yt)*1.05, '*k') %mean([1.785 2]), max(yt)*1.05, '*k')
plot([1.785 2.215], [1 1]*max(yt)*1.1, '-k', [ 1.9825 2.0125 ], [1 1]*max(yt)*1.15, '*k')
plot([2 2.215], [1 1]*max(yt)*0.9, '-k',  [2.0775 2.1075 2.1375], [1 1 1]*max(yt)*0.95, '*k') %mean([2 2.215]), max(yt)*0.95, '*k')

% plot([1 1.215], [1 1]*max(yt)*1.1, '-k',  [1.0775 1.1075 1.1375], [1 1 1]*max(yt)*1.15, '*k') % mean([0.785 1]), max(yt)*0.85, '*k')
% plot([0.785 1.215], [1 1]*max(yt)*1.2, '-k', [0.965 1 1.035], [1 1 1]*max(yt)*1.25, '*k')% mean([0.785 1.215]), max(yt)*0.95, '*k')
% 
% plot([1.785 2], [1 1]*max(yt)*1, '-k',  [1.8625 1.8925 1.9275], [1 1 1]*max(yt)*1.05, '*k') %mean([1.785 2]), max(yt)*1.05, '*k')
% plot([1.785 2.215], [1 1]*max(yt)*1.1, '-k', [ 2 ], [ 1]*max(yt)*1.15, '*k')
% plot([2 2.215], [1 1]*max(yt)*0.9, '-k',  [2.0775 2.1075 2.1375], [1 1 1]*max(yt)*0.95, '*k') %mean([2 2.215]), max(yt)*0.95, '*k')
% 
% plot([2.785 3], [1 1]*max(yt)*1, '-k')%, [2.8625 2.8925 2.9275], [1 1 1]*max(yt)*1.05, '*k') % mean([2.785 3]), max(yt)*0.85, '*k')
% plot([2.785 3.215], [1 1]*max(yt)*1.1, '-k', [2.965 3 3.035], [1 1 1]*max(yt)*1.15, '*k')% mean([2.785 3.215]), max(yt)*0.95, '*k')
% plot([3 3.215], [1 1]*max(yt)*0.9, '-k', [ 3.1075], [1 ]*max(yt)*0.95, '*k') % mean([3 3.215]), max(yt)*0.75, '*k')


%boxchart(group,T.RMSE_Insole_x,'GroupByColor',T.filename_calib_type)
%title('Erreur sur la position du centre de pression')
%boxchart(T.filename_calib,T.RMSE_Insole_x)
ylabel('RMSE AP (m)')
xlabel('Trial Type')
xticklabels(grouptext)
%legend('Sidecut Maneuver','Run','Walk')
% axe=gca;
% Unique_name_groupe=unique(group,'stable');
% x.Xaxis.Categories={Unique_name_groupe(1),'1',Unique_name_groupe(2),'2',Unique_name_groupe(3)};
% arrayfun(@(x)xline(x,'r-','LineWidth',1.5),categorical({'1','2'}))
% ax.XTickLabel([2,4]) = {''};
% 
% 
% figure(2)
% set(gca,'color','w')
% set(gcf,'color','w')
% A=boxchart(group,T.RMSE_Insole_y,'GroupByColor',T.filename_type);

%boxchart(group,T.RMSE_Insole_y,'GroupByColor',T.filename_calib_type)
%title('Erreur sur la position du centre de pression')
%boxchart(T.filename_calib,T.RMSE_Insole_x)
% ylabel('RMSE AP (m)')
% xlabel('Trial Type')
% legend('Sidecut Maneuver','Run','Walk')
% xticklabels(grouptext)

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
%             Fz=[Results(j).Fz_Insole_liste];
%             RMSE_x=[abs(Results(j).RMSE_Insole_liste(:,1))];
%             RMSE_y=[abs(Results(j).RMSE_Insole_liste(:,2))];
%             for lm=1:length(RMSE_x)
%                 if Fz(lm)>1500 && RMSE_x(lm)>0.04
%                     A=ii
%                     B=j
%                 end
%                     
%             end
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
Fz_round_2=[];
RMSE_x=[];
RMSE_y=[];
Contact=[];
for i=1:Nb_reduit
    Fz=[Fz;Precision_Semelles_reduit.Fz_Insole_liste{i,1}];
    %Fz_round=[Fz_round;round(Precision_Semelles_reduit.Fz_Insole_liste{i,1}*2,-2)/2];
    Fz_round=[Fz_round;round(Precision_Semelles_reduit.Fz_Insole_liste{i,1},-2)];
    %Fz_round_2=[Fz_round_2;round(Precision_Semelles_reduit.Fz_Insole_liste{i,1}*2,-2)/2];
        Fz_round_2=[Fz_round_2;round(Precision_Semelles_reduit.Fz_Insole_liste{i,1}+50,-2)-50];
    RMSE_x=[RMSE_x;abs(Precision_Semelles_reduit.RMSE_Insole_liste{i,1}(:,1))];
    RMSE_y=[RMSE_y;abs(Precision_Semelles_reduit.RMSE_Insole_liste{i,1}(:,2))];
    Contact=[Contact;Precision_Semelles_reduit.Contact_liste{i,1}];
end
Precision_force.Fz=Fz;
Precision_force.Fz_round=Fz_round;
Precision_force.Fz_round_2=Fz_round_2;
Precision_force.RMSE_x=RMSE_x;
Precision_force.RMSE_y=RMSE_y;
Precision_force.Contact=Contact;

Precision_force = struct2table(Precision_force); 


%writetable(Precision_force,'Precision_force_6.xlsx')

% Y=NaN(length(Fz_round),length(unique(Fz_round)));
% X=NaN(length(Fz_round),length(unique(Fz_round)));

% for i=1:length(Fz_round)
%     j=find(Fz_round(i)==unique(Fz_round));
%     Y(i,j)=RMSE_y(i);
%     X(i,j)=RMSE_x(i);
% end

% Cate_liste=unique(Fz_round);
% %Cate={num2str(unique(Fz_round))};
% Cate=cell(1,length(unique(Fz_round)));
% for j=1:length(unique(Fz_round))
%     Cate{1,j}=num2str(Cate_liste(j));
% end

Y=NaN(length(Fz_round_2),length(unique(Fz_round_2)));
X=NaN(length(Fz_round_2),length(unique(Fz_round_2)));

for i=1:length(Fz_round_2)
    j=find(Fz_round_2(i)==unique(Fz_round_2));
    Y(i,j)=RMSE_y(i);
    X(i,j)=RMSE_x(i);
end

Cate_liste=unique(Fz_round_2);
%Cate={num2str(unique(Fz_round))};
Cate=cell(1,length(unique(Fz_round_2)));
for j=1:length(unique(Fz_round_2))
    Cate{1,j}=num2str(Cate_liste(j));
end
figure(6)
subplot(2,1,1);


% y1 = linspace(0,2400,24);
% x1=zeros(1,24);
% 
% plot(y1,x1)
set(gca,'color','w')
set(gcf,'color','w')
boxplot_custom(Y,Cate,'outlier_multiplier',Inf)
ylabel('Error AP (m)')
xlabel('Vertical force (N)')

Number=[];
for k=1:24
    [a,b]=find(Y(:,k)>0);
Number = [Number,length(a) ];
end

%axis([0 2400 0 0.015])

subplot(2,1,2);
%figure(7)
set(gca,'color','w')
set(gcf,'color','w')
boxplot_custom(X,Cate,'outlier_multiplier',Inf)
ylabel('Error ML (m)')
xlabel('Vertical force (N)')

ylim([0 0.15])

%% selon contact
Y=NaN(length(Precision_force.Contact),length(unique(Precision_force.Contact)));
X=NaN(length(Precision_force.Contact),length(unique(Precision_force.Contact)));

for i=1:length(Fz_round)
    j=find(Precision_force.Contact(i)==unique(Precision_force.Contact));
    Y(i,j)=RMSE_y(i);
    X(i,j)=RMSE_x(i);
end

Cate_liste=unique(Precision_force.Contact);
%Cate={num2str(unique(Fz_round))};
Cate=cell(1,length(unique(Precision_force.Contact)));
for j=1:length(unique(Precision_force.Contact))
    Cate{1,j}=num2str(Cate_liste(j));
end
figure(6)
subplot(2,1,1);
set(gca,'color','w')
set(gcf,'color','w')
boxplot_custom(Y,Cate,'outlier_multiplier',Inf)
ylabel('Error AP (m)')
xlabel('Number of activated pressure cells')

subplot(2,1,2);
%figure(7)
set(gca,'color','w')
set(gcf,'color','w')
boxplot_custom(X,Cate,'outlier_multiplier',Inf)
ylabel('Error ML (m)')
xlabel('Number of activated pressure cells')



Number=[];
for k=1:24
    [a,b]=find(Y(:,k)>0);
Number = [Number,length(a) ];
end

% figure(5)
% set(gca,'color','w')
% set(gcf,'color','w')
% Precision_force.Fz_round = categorical(Precision_force.Fz_round);
% group =Precision_force.Fz_round ;
% A=boxchart(group,Precision_force.RMSE_x,'GroupByColor',Precision_force.Fz_round);
% 
% %boxchart(group,T.RMSE_Insole_y,'GroupByColor',T.filename_calib_type)
% %title('Erreur sur la position du centre de pression')
% %boxchart(T.filename_calib,T.RMSE_Insole_x)
% % ylabel('RMSE AP (m)')
% % xlabel('Calibration Trial Type')
% % legend('Sidecut Maneuver','Run','Walk')
% % xticklabels(grouptext)

