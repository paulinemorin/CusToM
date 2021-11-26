function Contact_detection = ContactDetectionInsole(filename)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

g=9.81;
f_cut=5;

InsoleDataFile(filename)
load(fullfile(filename,'InsoleData.mat'))

nbframe=length(InsoleData.data);

Pressure_data(1:16,:)=InsoleData.data(:,2:17)';
Pressure_data(17:32,:)=InsoleData.data(:,24:39)';

%% Position CoP
liste=[22,23,44,45];
CoP=[InsoleData.data(:,liste)]';

%% NaN et bruit, les semelles ont un code qui détecte le contact et 
%renseigne la position du CoP à l'origine s'il n'y a pas de contact
for i = 2 : nbframe
    if CoP(1,i)==0 && CoP(2,i)==0
        Pressure_data(1:16,i)=zeros(16,1);
    end
    if CoP(3,i)==0 && CoP(4,i)==0
        Pressure_data(17:32,i)=zeros(16,1);
    end
    if isnan(CoP(3,i))
        CoP(3:4,i) = CoP(3:4,i-1);
    end
        if isnan(CoP(1,i))
        CoP(1:2,i) = CoP(1:2,i-1);
    end
end

Pressure_data(Pressure_data<=2)=0;
%Pressure_data(Pressure_data<=1.5)=0;
Pressure_data(isnan(Pressure_data))=1;

for i = 1:32
    for j=3:nbframe-2
        if (Pressure_data(i,j-1)==0||Pressure_data(i,j-2)==0) && (Pressure_data(i,j+1)==0||Pressure_data(i,j+2)==0)
            Pressure_data(i,j)=0;
        end
    end
end

%% Tableau contact from pressur synch
Contact_detection = isinf(Pressure_data*inf);


%% Retirer point illogiques

%Tableau de groupement de points
Orteil_1=[9,10,11,12,13,14,15,16];
Orteil_2=Orteil_1+ones(1,length(Orteil_1))*16;

Talon_1=[1,2,3,4];
Talon_2=Talon_1+ones(1,length(Talon_1))*16;

Milieu_1=[5,6,7,8];
Milieu_2=Milieu_1+ones(1,length(Milieu_1))*16;

for i = 1 : nbframe
    if sum(Contact_detection(Orteil_1,i))<3
        Contact_detection(Orteil_1,i)=zeros(length(Orteil_1),1);
    end
    if sum(Contact_detection(Orteil_2,i))<3
        Contact_detection(Orteil_2,i)=zeros(length(Orteil_2),1);
    end  
    if sum(Contact_detection(Talon_1,i))<2
        Contact_detection(Talon_1,i)=zeros(length(Talon_1),1);
    end
    if sum(Contact_detection(Talon_2,i))<2
        Contact_detection(Talon_2,i)=zeros(length(Talon_2),1);
    end
    
    if sum(Contact_detection(Milieu_1,i))<2
        Contact_detection(Milieu_1,i)=zeros(length(Milieu_1),1);
    end
    if sum(Contact_detection(Milieu_2,i))<2
        Contact_detection(Milieu_2,i)=zeros(length(Milieu_2),1);
    end
end

end

