function [L,Typ] = distance_point(Point1,Bone1,Point2,Bone2,Human_model,q)
% Computation of the distance between two points
%
%   INPUT
%   - Point1: position of the anatomical position of the first studied
%   point on its solid 
%   - Bone1: number of the solid containing the first studied point
%   - Point2: position of the anatomical position of the second studied
%   point on its solid  
%   - Bone2: number of the solid containing the second studied point
%   - Human_model: osteo-articular model (see the Documentation for the structure)
%   - q: vector of joint coordinates at a given instant
%   OUTPUT
%   - L: distance between the two studied points
%   - Typ: relative orientation of the two points (-1 ou 1)
%________________________________________________________
%
% Licence
% Toolbox distributed under GPL 3.0 Licence
%________________________________________________________
%
% Authors : Antoine Muller, Charles Pontonnier, Pierre Puchaud and
% Georges Dumont
%________________________________________________________

if Bone1 == Bone2
    A=Human_model(Bone1).anat_position{Point1,2};
    B=Human_model(Bone2).anat_position{Point2,2};
    L = norm(B-A);
    Typ = sign(B(2)-A(2));
    return
end

% finding the index of common antecedent (ppac:plus proche antecedent commun)
t = 1;
t1 = find_solid_path(Human_model,Bone1,1);
t2 = find_solid_path(Human_model,Bone2,1);
while t<=min(numel(t1),numel(t2)) && t1(t)==t2(t)
    ppac = t1(t);
    ppac_t = t;
    t = t + 1;
end

% minimal paths
solid1 = t1(ppac_t:end);
solid2 = t2(ppac_t:end);

% Direct kinematics initialization (Initialisation de la cinématique directe)
Human_model(ppac).p = zeros(3,1);
Human_model(ppac).R = eye(3);

% Computation of coordinates of point A (calcul des coordonnées du point A)
if numel(solid1)==1
    A = Human_model(Bone1).c+Human_model(Bone1).anat_position{Point1,2};
else
    for n=solid1(2:end)
        m = Human_model(n).mother;
        Human_model(n).p = Human_model(m).R * Human_model(n).b + Human_model(m).p;
        Human_model(n).R = Human_model(m).R * Rodrigues(Human_model(n).a,q(n)) * Rodrigues(Human_model(n).u,Human_model(n).theta);
    end
    l = solid1(end);
    Human_model(l).c_global = Human_model(l).p + Human_model(l).R * Human_model(l).c;
    A = Human_model(l).c_global + Human_model(l).R * Human_model(Bone1).anat_position{Point1,2};
end

% Computation of coordinates of point B (calcul des coordonnées du point B)
if numel(solid2)==1
    B = Human_model(Bone2).c+Human_model(Bone2).anat_position{Point2,2};
else
    for n=solid2(2:end)
        m = Human_model(n).mother;
        Human_model(n).p = Human_model(m).R * Human_model(n).b + Human_model(m).p;
        Human_model(n).R = Human_model(m).R * Rodrigues(Human_model(n).a,q(n)) * Rodrigues(Human_model(n).u,Human_model(n).theta);
    end
    l = solid2(end);
    Human_model(l).c_global = Human_model(l).p + Human_model(l).R * Human_model(l).c;
    B = Human_model(l).c_global + Human_model(l).R * Human_model(Bone2).anat_position{Point2,2};
end

%Distance between A and B
L = norm(B-A);
Typ = sign(B(2)-A(2));

end