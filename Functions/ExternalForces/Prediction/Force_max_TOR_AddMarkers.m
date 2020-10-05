function Cpi = Force_max_TOR_AddMarkers(py2,pz2,vp2,Mass,zcrit,vcrit,offset,width)
% Maximal force available at a contact point for the prediction of the
% ground reaction forces for unconventional floors requiring additional
% markers
%   if the vertical position and the normal of the velocity of the point is
%   lower than thresholds, the maximal force available is equal to 40% of
%   the subject weight
%
%   INPUT
%   - py2: relative position of the point according to Y2
%   - pz2: relative position of the point according to Z2
%   - vp2: norm of the relative velocity of the point
%   - Mass: subject mass
%   - zcrit: vectical position threshold
%   - vcrit: norm of the velocity threshold
%   - offset: offset between structure markers and contact surface
%   - width: structure width
%________________________________________________________
%
% Licence
% Toolbox distributed under GPL 3.0 Licence
%________________________________________________________
%
% Authors : Antoine Muller, Charles Pontonnier, Pierre Puchaud and
% Georges Dumont
%________________________________________________________

prop=0.4;
if 0<=abs(py2) && abs(py2)<=width && abs(pz2)<=zcrit+offset && vp2<=vcrit
    Cpi=prop* Mass*9.81; 
else
    Cpi=0;
end

end