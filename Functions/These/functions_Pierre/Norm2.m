function B = Norm2(A, dim)
%
% Calcul la norme d'un tableau de vecteurs en ligne
%

if nargin == 1
    dim = 2 ;
end

B = sqrt(sum(A.^2,dim)) ;
%
% fin de la fonction