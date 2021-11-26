function [Event,ID] = DetectRunPhases(trialName,FootOnPlates)
if iscell(trialName)
    trialName=trialName{1};
end
SideFP{1,1} = FootOnPlates{1}(1);
SideFP{2,1} = FootOnPlates{2}(1);

h = btkReadAcquisition(trialName);
[~,Mk_info]=btkGetMarkers(h);

[~,info] = btkGetForcePlatforms(h);
r = Mk_info.frequency / info(1).frequency;

% On charge torseur des efforts des deux plateformes
% threshold à 40N
% Tous les efforts sont mis à zeros lorsque l'effort en z est inférieur à
% 40 N
Threshold = 40;
[GroundReactionWrenches] = btkGetGroundReactionWrenches(h,Threshold);

hs_to={'HS','TO'};
for ii=1:size(SideFP,1)
    ind(1)= find(GroundReactionWrenches(ii).F(:,3)~=0,1,'first'); % Colonnes 1 2 3, x,y,z
    ind(2)= find(GroundReactionWrenches(ii).F(:,3)~=0,1,'last');
    
    for jj=1:length(ind)
        ij = sub2ind([size(SideFP,1) length(ind)],jj,ii);
 
        event(ij).frame=round(ind(jj)*r,0);
        % HeelStrike (HS) ou ToeOff (TO)
        event(ij).type=char(hs_to{jj});
        % Coté de la jambe 'R' ou 'L'
        event(ij).leg_side=SideFP{ii,1};
        % Numére de la plateforme
        event(ij).FP=ii;
    end
    
end

%% On réordonne les évenements
Event=struct('frame',[],'type',char(),'leg_side',char(),'FP',[]);
[~,z]=sort([event.frame]);
Event=event(z);

%% Compatibilité dynamique inverse
%phases d'envol
n = 15;
if Event(1).frame-n <=0
    deb=1;
else
    deb = Event(1).frame-n ;
end
if length(GroundReactionWrenches(ii).F(:,3))-(Event(end).frame+n )<0
    fin=round(length(GroundReactionWrenches(ii).F(:,3))*r,0);
else
    fin = (Event(end).frame+n );
end

ID=[deb fin];
end

%A stocker dans meta.FrameOfInterest(ii).Event

%Et enfin mettre les frames compatible avec la dynamique inverse,
% tu peux prendre quelques frames avant après le contact des plateformes
% pour les phases d'envols pour la course et le changment de direction
%meta.FrameOfInterest.ID = [ num de la frame 1, num de la frame 2];




