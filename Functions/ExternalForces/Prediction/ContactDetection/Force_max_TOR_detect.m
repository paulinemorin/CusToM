function Cpi = Force_max_TOR_detect(contact_boolean,Mass)
% Maximal force available at a contact point for the prediction of the ground reaction forces
%   if the vertical position and the normal of the velocity of the point is
%   lower than thresholds, the maximal force available is equal to 40% of
%   the subject weight
%
%   INPUT
%   - pz: vectical position of the point
%   - vp: norm of the velocity of the point
%   - Mass: subject mass
%   - zcrit: vectical position threshold
%   - vcrit: norm of the velocity threshold
%________________________________________________________
%
% Licence
% Toolbox distributed under GPL 3.0 Licence
%________________________________________________________
%
% Authors : Antoine Muller, Charles Pontonnier, Pierre Puchaud,
% Georges Dumont and Pauline Morin
%________________________________________________________

prop=0.4;

if contact_boolean
    Cpi=prop* Mass*9.81; 
else
    Cpi=0;
end

end