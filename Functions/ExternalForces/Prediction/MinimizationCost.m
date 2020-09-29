function cost = MinimizationCost(X, Human_model, Prediction, i,  external_forces_pred, f6dof, t6dof0,  Fmax, g)
    alpha =1;
    beta = 1;
    
    for k = 1:numel(Prediction)
        Prediction(k).efforts(i,1)=X(k)*Prediction(k).efforts_max(i,1);
        Prediction(k).efforts(i,2)=X(k+numel(Prediction))*Prediction(k).efforts_max(i,2);
        Prediction(k).efforts(i,3)=X(k+2*numel(Prediction))*Prediction(k).efforts_max(i,3);
    end
    
    external_forces_pred=addForces_Prediction_frame_par_frame(X,external_forces_pred,Prediction,Fmax,i);
    
    [Human_model,f6dof(:,i),t6dof0(:,i)]=InverseDynamicsSolid(Human_model,external_forces_pred(i).fext,g,1);
    
    cost = alpha * sum(X.^2) + beta * sum([Human_model.torques]);
    
    
return
