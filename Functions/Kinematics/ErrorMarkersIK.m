function [error] = ErrorMarkersIK(q,real_markers,f,list_markers,Rcut,pcut)
% Computation of reconstruction error for the inverse kinematics step
%   Computation of the distance between the position of each experimental 
%   marker and the position of the corresponded model marker
%
%   INPUT
%   - q: vector of joint coordinates at a given instant
%   - nb_cut: number of geometrical cut done in the osteo-articular model
%   - real_markers: 3D position of experimental markers
%   - f: current frame
%   - list_markers: list of the marker names
%   - Rcut: pre-initialization of Rcut
%   - pcut: pre-initialization of pcut
%   OUTPUT
%   - error: distance between the position of each experimental marker and
%   the position of the corresponded model marker 
%________________________________________________________
%
% Licence
% Toolbox distributed under GPL 3.0 Licence
%________________________________________________________
%
% Authors : Antoine Muller, Charles Pontonnier, Pierre Puchaud and
% Georges Dumont
%________________________________________________________
[Rcut,pcut]=fcut(q,pcut,Rcut);

error=zeros(numel(list_markers),1);
for m=1:numel(list_markers)
    error(m,:) = norm(feval([list_markers{m} '_Position'],q,pcut,Rcut) - real_markers(m).position(f,:)');
end

end