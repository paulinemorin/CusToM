function [event,ID]=DetectJumpPhases(trialName,FootOnPlates)
if iscell(trialName)
    trialName=trialName{1};
end

SideFP{1,1} = FootOnPlates{1}(1);
SideFP{2,1} = FootOnPlates{2}(1);

h = btkReadAcquisition(trialName);
[~,Mk_info]=btkGetMarkers(h);
Mkstruct = C3dProcessedData(trialName(1:end-4),{'RBWT','RKNE','RANE'});

[~,info] = btkGetForcePlatforms(h);


r = Mk_info.frequency/info(1).frequency;
event=struct();

%% JUMP OFF
[btkEvent]=btkGetEvents(h);
Firstframe = btkGetFirstFrame(h);
if isfield(btkEvent,'General_Jump_off')
    [t]=btkGetEventsValues(h);
    
    ID(1)=round(t*Mk_info.frequency,0) - Firstframe;
    
    event(1).frame = ID(1);
    event(1).t= ID(1)*1/Mk_info.frequency;
    % HeelStrike (HS) ou ToeOff (TO)
    event(1).type='JO';
    % Cot� de la jambe 'R' ou 'L'
    event(1).leg_side={'R','L'};
    % Num�re de la plateforme
    event(1).FP={'1','2'};
else
    warning('no jump off labelled')
    %keyboard
end


%% LANDING
% On charge torseur des efforts des deux plateformes
% threshold � 40N
% Tous les efforts sont mis � zeros lorsque l'effort en z est inf�rieur �
% 40 N
Threshold = 40;
[GroundReactionWrenches] = btkGetGroundReactionWrenches(h,Threshold);

for ii=1:numel(GroundReactionWrenches)
    Ind(ii) = find(GroundReactionWrenches(ii).F(:,3)~=0,1,'first');
end
Ind = min(Ind);
event(2).frame =round(Ind*r,0);
event(2).t= round(Ind*r,0)*1/Mk_info.frequency;
% Landing (L)
event(2).type='L';
% Cot� de la jambe 'R' ou 'L'
event(2).leg_side={'R','L'};
% Num�re de la plateforme
event(2).FP={'1','2'};

%% STANDING
% apr�s le max
% mass = (Norm2(GroundReactionWrenches(1).F)+Norm2(GroundReactionWrenches(2).F))/9.81;
% mass_fin = mass(end);
% [~,ind_max]=max(mass,[],1);
n=0.05;

% figure
% % plot(mass(ind_max:end));
% subplot(2,1,1)
% plot((1:length(mass))*r,mass(1:end)); 
% hold on
% plotHline(mass_fin,'k')
% plotHline(mass_fin*(1-n),'b')
% plotHline(mass_fin*(1+n),'b')
% 
% ind_stabilisation = find(mass(ind_max:end)> mass_fin*(1-n) & mass(ind_max:end)<mass_fin*(1+n),1,'first');

KP = Mkstruct(1).position_c3d-Mkstruct(2).position_c3d; % Knee-pelvis
KA = Mkstruct(1).position_c3d-Mkstruct(3).position_c3d; % Knee-Ankle
Angle = acosd(dot(KP',KA')'./(Norm2(KA).*Norm2(KP)));

[~,ind_max]=max(Angle,[],1);
A=Angle-Angle(end);
AA=A(ind_max:end);
figure
plot(AA); hold on
plotHline(0,'k')
plotHline(1,'b')
plotHline(-1,'b')
ind_standing = find(AA> -1 & AA< 1 ,1,'first');
plot(ind_standing,AA(ind_standing),'o')

ind_standing = ind_standing + ind_max;

event(3).frame=round((ind_max+ind_standing),0);
event(3).t=round((ind_max+ind_standing),0)*1/Mk_info.frequency;
% Standing
event(3).type='S';
% Cot� de la jambe 'R' ou 'L'
event(3).leg_side={'R','L'};
% Num�re de la plateforme
event(3).FP={'1','2'};

ID(2)=length(A);

end