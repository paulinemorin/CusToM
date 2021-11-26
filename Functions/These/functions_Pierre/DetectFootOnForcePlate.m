function [plate,frame]=DetectFootOnForcePlate(real_markers,ForceplatesData,Threshold)
if nargin<3
    Threshold=0.5;
end
% Vitesse du marqueur
V=Norm2(real_markers.vitesse);
% On filtre pour eviter les chocs
try % si on peut;
    data = filt_data(V,5,1/(real_markers.time(2)-real_markers.time(1)));
catch
   data =V; 
end
%         plot(data); hold on
% on detect les changement de signe
ind = find(diff(data(1:end-1)).*diff(data(2:end))<0);
if isempty(ind)
    warning(['No sign variation on velocity of the marker ' char(real_markers.name)])
    [~,ind]=min(data);
else
end
%         plot(ind,data(ind),'o')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% threshold %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Vitesse limite

Ind_Vpossible = find(V(ind)<Threshold);
%         plot(ind(V(ind)<Threshold),data(ind(V(ind)<Threshold)),'.')
if isempty(Ind_Vpossible)
    warning('No local minimum found for the velocity of the marker')
    keyboard
end

% Position des marqueurs pour des vitesses de marqueurs très faible
for jj=1:length(Ind_Vpossible)
    pos(jj,:)=real_markers.position_c3d(ind(Ind_Vpossible(jj)),:);
end

% figure; ShowPlates(ForceplatesData)
% hold on; plot3(pos(:,1),pos(:,2),pos(:,3),'o');

% On cherche de quelle plateforme le marqueurs est le plus proche
for jj=1:length(ForceplatesData)
    d_pos(:,jj) = Norm2(pos - mean([ForceplatesData(jj).corners/1000]'));
end



[val,I]=min(d_pos,[],1);
[~,II]=min(val);
instant=I(II); plate=II;
frame = ind(Ind_Vpossible(instant));

end

function ShowPlates(ForceplatesData)

x_fp = []; y_fp = []; z_fp = [];
for i=1:numel(ForceplatesData)
        x_fp = [x_fp ForceplatesData(i).corners(1,:)'/1000]; %#ok<AGROW> % mm -> m
        y_fp = [y_fp ForceplatesData(i).corners(2,:)'/1000]; %#ok<AGROW> % mm -> m
        z_fp = [z_fp ForceplatesData(i).corners(3,:)'/1000]; %#ok<AGROW> % mm -> m
        XX=mean([ForceplatesData(i).corners]'/1000);
        plot3(XX(:,1),XX(:,2),XX(:,3)+0.05,'o','MarkerFaceColor','k'); hold on
end
patch(x_fp,y_fp,z_fp,[1 1 1]); axisGeom
end