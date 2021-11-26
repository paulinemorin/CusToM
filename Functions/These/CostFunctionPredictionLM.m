function cost = CostFunctionPrediction(X,Afric,bfric,A,b,lb,ub)

boundspen = zeros(length(X),1);

idxsup = X>ub;
idxinf = X<lb;
    
boundspen(idxsup)=(X(idxsup)-ub(idxsup)).^2;
boundspen(idxinf)=(X(idxinf) - lb(idxinf)).^2;


inequaliy_constraint = Afric*X - bfric;
idx = inequaliy_constraint > 0;

frictionpen = zeros(length(inequaliy_constraint),1);
frictionpen(idx) =  inequaliy_constraint(idx).^2;





%cost = [ 100*X ; A*X - b ; boundspen ];
cost = [ 100*X ; 100*(A*X - b) ; boundspen ; frictionpen]';

end