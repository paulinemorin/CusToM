function Cpi = Force_max_TOR(contact_boolean,Mass)
% Maximal force available at a contact point for the prediction of the ground reaction forces
%   if the point has been detected as in contact,
%   the maximal force available is equal to 40% of
%   the subject weight
%
%   INPUT
%   - contact_boolean: activation of the considered contact point
%   - Mass: subject mass
%________________________________________________________
%
% Licence
% Toolbox distributed under GPL 3.0 Licence
%________________________________________________________
%
% Authors : Antoine Muller, Charles Pontonnier, Pierre Puchaud, Pauline
% Morin and Georges Dumont 
%________________________________________________________

prop=0.4;

if contact_boolean
    Cpi=prop* Mass*9.81; 
else
    Cpi=0;
end

end