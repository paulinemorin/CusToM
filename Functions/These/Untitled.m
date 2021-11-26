
clear
filename='StraightRunVinyl0003';
%filename='WalkVinyl0001';
X = CalibrationInsoleModel(filename);
%
%X=zeros(12,1);
[COP_Xp, Contact] = COP_Xp_f(filename,X);

CoP_size_feet = CoPInsole(X,filename);
CoP_size_feet_G = CoPInsoleGlob(X,filename);

[Foot_length,Foot_large] = Foot_size(filename);
[Foot_length_insole,Foot_large_insole] = Foot_size_insole(Foot_length,Foot_large);
taille=[Foot_length_insole,Foot_large_insole];

%% Transformation cop data semelles

for ii=[1,2]
    CoP_size_feet((ii-1)*2+1,:)=[taille(1).*CoP_size_feet((ii-1)*2+1,:)];
    CoP_size_feet((ii-1)*2+2,:)=[taille(2).*CoP_size_feet((ii-1)*2+2,:)];
end

% 

start = find(sum(Contact,1)>0,1)-10;
stop = find(sum(Contact,1)>0);
stop = stop(length(stop))+10;

CoP_size_feet(:,1:start)=NaN;
CoP_size_feet(:,stop:end)=NaN;

CoP_size_feet_G(:,1:start)=NaN;
CoP_size_feet_G(:,stop:end)=NaN;




Solids{2,1}='RFoot';
Solids{1,1}='LFoot';

figure
hold on

for j = 1:2
    for i=1:3
        subplot(2,3,i+(j-1)*3)
        COP_Xp.(Solids{j})(COP_Xp.(Solids{j})==0) = NaN;
        x=CoP_size_feet_G(i+(j-1)*3,:);        
        plot(x,'*')
        hold on
        x=COP_Xp.(Solids{j})(:,i);
        plot(x,'*')
    end
end

figure
hold on
set(gca,'color','w')
set(gcf,'color','w')
title_liste={'Pied gauche','Pied droit'};

for i=1:2
    subplot(1,2,i)
%     x=COP_Xp.(Solids{j})(:,(i-1)*2+1);
%     y=COP_Xp.(Solids{j})(:,(i-1)*2+2);
%     plot(x,y,'*')
%     hold on
    x=CoP_size_feet((i-1)*2+1,:);
    y=CoP_size_feet((i-1)*2+2,:);
    plot(y,x,'*')
    title(title_liste(i))
    %axis([-0.5*taille(2) 0.5*taille(2) -0.5*taille(1) 0.5*taille(1)])
    axis([nanmean(y)-0.5*taille(2) nanmean(y)+0.5*taille(2) nanmean(x)-0.5*taille(1) nanmean(x)+0.5*taille(1)])
    
    ylabel({'composante antero-posterieure','longeur semelle (m)'})
    xlabel({'coomposante medio-latérale','largeur semelle (m)'})
end



figure
hold on
set(gca,'color','w')
set(gcf,'color','w')
title_liste={'Pied gauche','Pied droit'};

for i=1:2
    subplot(1,2,i)
    x=CoP_size_feet_G((i-1)*3+1,:);
    y=CoP_size_feet_G((i-1)*3+2,:);
    plot(x,y,'*')
    hold on
    COP_Xp.(Solids{i})(COP_Xp.(Solids{i})==0) = NaN;
    x=COP_Xp.(Solids{i})(:,2);
    y=COP_Xp.(Solids{i})(:,1);


    plot(y,x,'*')
    hold on
    axis([nanmean(y)-0.5*taille(2) nanmean(y)+0.5*taille(2) nanmean(x)-0.5*taille(1) nanmean(x)+0.5*taille(1)])
    title(title_liste(i))
    %axis([-0.5*taille(2) 0.5*taille(2) -0.5*taille(1) 0.5*taille(1)])
    ylabel({'y dans le repère global','longeur semelle (m)'})
    xlabel({'x dans le repère global','largeur semelle (m)'})
end


