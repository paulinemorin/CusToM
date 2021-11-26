function [b]=fct_moyenne_glissante(a,n) % a : vecteur colonne ? lisser, n : valeur pour la moyenne glissante (doit ?tre impair)
 n2=(n-1)/2;
 a2=zeros(size(a,1)+2*n2,1); a2(n2+1:size(a2,1)-n2,1)=a;
   for i =1:n2 % extrapolation sur n2 valeur avant et apr?s la matrice a (permet de ne pas perdre de valeurs et eviter les effets de bords)
     a2(n2+1-i,1)=a2(n2+1,1)-(a2(n2+1+i,1)-a2(n2+1));
     a2(size(a2,1)-n2+i,1)=a2(size(a2,1)-n2,1)-(a2(size(a2,1)-n2-i,1)-a2(size(a2,1)-n2,1));
   end % end for
   % lissage - sens horaire
   a3=a2;
   for i=n2+1:(size(a2,1)-n2)
     a3(i,1)=mean(a2(i-n2:i+n2,1));
   end % end for 
   % lissage - sens anti-horaire
   a4=a3;
   for i=n2:size(a2,1)-(n2+1)
     a4(size(a2,1)-i,1)=mean(a3(size(a2,1)-i-n2:size(a2,1)-i+n2,1));
   end % end for 
  b=a4(n2+1:size(a2,1)-n2);
end