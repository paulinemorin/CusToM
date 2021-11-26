function [event,ID]=DetectOnlyOneFootPhase(trialName,FootOnPlates)
if iscell(trialName)
    trialName=trialName{1};
end

SideFP{1,1} = FootOnPlates{1}(1);
SideFP{2,1} = FootOnPlates{2}(1);
h = btkReadAcquisition(trialName);
[~,Mk_info]=btkGetMarkers(h);
[~,info] = btkGetForcePlatforms(h);
[GroundReactionWrenches] = btkGetGroundReactionWrenches(h,40);%threshold of 10N

for ii =1:numel(GroundReactionWrenches)
    %    dCOP(:,ii)= filt_data(Norm2(GroundReactionWrenches(ii).P),50,info(ii).frequency)  ;
    dCOP(:,ii)= Norm2(GroundReactionWrenches(ii).P);
    %     dCOP2(:,ii)= resample(Norm2(GroundReactionWrenches(ii).P),Mk_info.frequency,info.frequency);
    vCOP(:,ii)= round(derivee2(1/info(ii).frequency,dCOP(:,ii)),5);
    dCOPcut{ii}=[find(vCOP(:,ii)~=0) dCOP(vCOP(:,ii)~=0,ii)];
    vCOPcut{ii}=[find(vCOP(:,ii)~=0) vCOP(vCOP(:,ii)~=0,ii)];
    %     aCOP(:,ii)= round(derivee2(1/info(ii).frequency,vCOP(:,ii)),5);
end

% Sens de la marche
if dCOPcut{1}(1,1)<dCOPcut{2}(1,1) % Plateforme 1 puis 2
    
    toeoff = find(vCOP(:,1)~=0, 1, 'last' );
    heelstrike = find(vCOP(:,2)~=0, 1, 'first' );
    
    frame(1) = heelstrike;
    typ{1}=char('HS');
    leg_side{1}=char(SideFP{2,1});
    PF(1) = 2;
    
    frame(2) = toeoff;
    typ{2}=char('TO');
    leg_side{2}=char(SideFP{1,1});
    PF(2) = 1;
    
elseif dCOPcut{1}(1,1)>dCOPcut{2}(1,1) % Plateforme 2 puis 1
    
    heelstrike = find(vCOP(:,1)~=0, 1, 'first' );
    toeoff = find(vCOP(:,2)~=0, 1, 'last' );
    
    frame(1) = heelstrike;
    typ{1}=char('HS');
    leg_side{1}=char(SideFP{1,1});
    PF(1) = 1;
    
    frame(2) = toeoff;
    typ{2}=char('TO');
    leg_side{2}=char(SideFP{2,1});
    PF(2) = 2;
end

% Couper les parties superposées
rmv_frames = intersect(dCOPcut{1}(:,1),dCOPcut{2}(:,1));

for ii=1:numel(GroundReactionWrenches)
    [~,ind]=intersect(dCOPcut{ii}(:,1),rmv_frames);
    dCOPcut{ii}(ind,:)=[];
end
% Couper la fin si la variation de distance du COP est trop importante
% Variation de distance du CoP sur la plateforme inferieur a un seuil
Threshold = 400; %mm
if dCOPcut{1}(1,1)<dCOPcut{2}(1,1)
    if abs(max(dCOPcut{1}(:,2))-min(dCOPcut{1}(:,2)))>Threshold
        dCOPcut{1}(1:20,:)=[];
    end
    if abs(max(dCOPcut{2}(:,2))-min(dCOPcut{2}(:,2)))>Threshold
        dCOPcut{2}(end-19:end,:)=[];
    end
end
if dCOPcut{1}(1,1)>dCOPcut{2}(1,1)
    if abs(max(dCOPcut{2}(:,2))-min(dCOPcut{2}(:,2)))>Threshold
        dCOPcut{2}(1:20,:)=[];
    end
    if abs(max(dCOPcut{1}(:,2))-min(dCOPcut{1}(:,2)))>Threshold
        dCOPcut{1}(end-19:end,:)=[];
    end
end
figure;
plot(dCOPcut{1}(:,1),dCOPcut{1}(:,2)); hold on
plot(dCOPcut{2}(:,1),dCOPcut{2}(:,2))
% plot(dCOPcut{1}(:,1),fct_moyenne_glissante(dCOPcut{1}(:,2),21))
% plot(dCOPcut{2}(:,1),fct_moyenne_glissante(dCOPcut{2}(:,2),21))
% plot(vCOPcut{1}(:,1),filt_data(vCOPcut{1}(:,2),2,info(ii).frequency)); hold on
% plot(vCOPcut{2}(:,1),vCOPcut{2}(:,2))
hs_to={'HS','TO'};
for ii=1:numel(GroundReactionWrenches)
    x=dCOPcut{ii}(:,1);
    y=fct_moyenne_glissante(dCOPcut{ii}(:,2),25);
    r(ii)=info(ii).frequency/Mk_info.frequency;
    if abs(max(y)-min(y))>Threshold
        
        X=x(1:r(ii):end);
        Y=y(1:r(ii):end);
        [~,ind] = Umethod(X,Y,'on');
        ind = ind*r(ii);
        for jj=1:length(ind)
            frame(end+1) = x(ind(jj));
            typ{end+1}=char(hs_to{jj});
            PF(end+1) = ii;
            if dCOPcut{1}(1,1)<dCOPcut{2}(1,1)% PF 1 to PF2
                leg_side{end+1}=char(SideFP{jj,1});
            else % PF 2 to PF 1
                leg_side{end+1}=char(SideFP{length(ind)+1-jj,1});
            end
        end
    else
        %keyboard
        frame(end+1) = x(1);
        typ{end+1}=char(hs_to{1});
        PF(end+1) = ii;
        leg_side{end+1} = SideFP{ii,1} ;
    end
    % figure
    % plot(x,y,'k-'); hold on
    % plot(x(ind),y(ind),'bo')
end

event=struct('frame',[],'type',char(),'leg_side',char(),'FP',[]);
[~,z]=sort(frame);

for ii=1:length(frame)
    event(ii).frame=round(frame(z(ii))/r(PF(z(ii))),0); % to the sample frequency of MoCap
    event(ii).type=typ(z(ii));
    event(ii).leg_side=leg_side(z(ii));
    event(ii).FP=PF(z(ii));
end

% Compatibilité Inverse dynamics

nb_pieds(1)=length(unique([event([event.FP]==1).leg_side]));
nb_pieds(2)=length(unique([event([event.FP]==2).leg_side]));
if nb_pieds(1)==1 && nb_pieds(2)==1
    ID = [event(1).frame,event(end).frame];
elseif nb_pieds(1)==2 && nb_pieds(2)==1
    frames=[event([event.FP]==2).frame];
    event_fp2=event([event.FP]==1);
    if frames(1)<event_fp2(1).frame
        frames=[frames,event_fp2(1:2).frame];
    else
        frames=[frames,event_fp2(end-1:end).frame];
    end
    ID = [min(frames),max(frames)];
elseif nb_pieds(1)==1 && nb_pieds(2)==2
    frames=[event([event.FP]==1).frame];
    event_fp1=event([event.FP]==2);
    if frames(1)<event_fp1(1).frame
        frames=[frames,event_fp1(1:2).frame];
    else
        frames=[frames,event_fp1(end-1:end).frame];
    end
    ID = [min(frames),max(frames)];
else
    [~,d]=min(abs(frame(3:end)-frame(1)));
    [~,f]=min(abs(frame(3:end)-frame(2)));
    ID = [frame(d+2),frame(f+2)];
end
  ID = round(ID/r(1),0);
end

%Threshold
% Threshold2=1000;
% for ii=1:numel(GroundReactionWrenches)
%     x=dCOPcut{ii}(:,1);
%     y=fct_moyenne_glissante(dCOPcut{ii}(:,2),25);
%     dy=derivee2(1/info(ii).frequency,y);
%     find(dy>Threshold2);
%     % Il reste a dire quand le seuil de vitesse est trop important.
%     if dCOPcut{1}(1,1)<dCOPcut{2}(1,1)
%         if abs(max(dCOPcut{1}(:,2))-min(dCOPcut{1}(:,2)))>Threshold
%             dCOPcut{1}(1:100,:)=[];
%         end
%         if abs(max(dCOPcut{2}(:,2))-min(dCOPcut{2}(:,2)))>Threshold
%             dCOPcut{2}(end-99:end,:)=[];
%         end
%
%     elseif dCOPcut{1}(1,1)>dCOPcut{2}(1,1)
%         if abs(max(dCOPcut{2}(:,2))-min(dCOPcut{2}(:,2)))>Threshold
%             dCOPcut{2}(1:100,:)=[];
%         end
%         if abs(max(dCOPcut{1}(:,2))-min(dCOPcut{1}(:,2)))>Threshold
%             dCOPcut{1}(end-99:end,:)=[];
%         end
%     end
%     plot(x,dy)
% end

% pts=findchangepts(dCOPcut{ii}(:,2),'Statistic','linear','MinThreshold',20);
% pts=findchangepts(y,'MaxNumChanges',2,'Statistic','rms')
% pts=findchangepts(y,'Statistic','linear','MinThreshold',0.5*10^5)
% plot(dCOPcut{ii}(:,1),dCOPcut{ii}(:,2),'-b',...
%     x(pts),y(pts),'rx')
% plot(aCOP);
% findpeaks(vCOP(:,1))


% end