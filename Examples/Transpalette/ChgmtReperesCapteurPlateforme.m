%Changement de repère entre les repères du capteur et de la plateforme
%Pré manip BAHAMaS
%center --> origine 

function [TC2] = ChgmtReperesCapteurPlateforme(T21,center)
phi = (45*pi)/180; %D'après les données du constructeur
nbframes = size(center,1);
TC2 = zeros(4,4,nbframes);
for i = 1:nbframes
    X_capteur = cos(phi)*[1 0 0]+sin(phi)*[0 1 0];
    X_capteur = X_capteur/norm(X_capteur);
    Y_capteur = -sin(phi)*[1 0 0]+cos(phi)*[0 1 0];
    Y_capteur = Y_capteur/norm(Y_capteur);
    Z_capteur = [0 0 1];
    R2C = [X_capteur' Y_capteur' Z_capteur'];
    T12 = inv(T21(:,:,i));
    %Center et TC2 exprimées dans le repère monde
    %Différence des deux donne l'origine de capteur dans 2 mais exprimée dans le repère 2
    %On la multiplie par R12 pour l'avoir dans le repère 2
    %T2C est donc bien exprimée dans le repère 2
    T2C = [R2C T21(1:3,1:3,i)*(center(i,:)'-T12(1:3,end)); 0 0 0 1];
    TC2(:,:,i) = inv(T2C);
end