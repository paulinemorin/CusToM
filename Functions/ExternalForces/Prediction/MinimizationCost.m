function cost = MinimizationCost(X, BiomechanicalModel, Prediction, i,  external_forces_pred, f6dof, t6dof0,  Fmax, g, Nb_muscles , alpha_beta_gamma)
alpha = alpha_beta_gamma(1);
beta =  alpha_beta_gamma(2);
gamma =  alpha_beta_gamma(3);

Pred_length=3*numel(Prediction);


if beta
    Human_model = BiomechanicalModel.OsteoArticularModel;
    
    for k = 1:numel(Prediction)
        Prediction(k).efforts(i,1)=X(k)*Prediction(k).efforts_max(i,1);
        Prediction(k).efforts(i,2)=X(k+numel(Prediction))*Prediction(k).efforts_max(i,2);
        Prediction(k).efforts(i,3)=X(k+2*numel(Prediction))*Prediction(k).efforts_max(i,3);
    end
    
    external_forces_pred=addForces_Prediction_frame_par_frame(X,external_forces_pred,Prediction,Fmax,i);
    
    [Human_model,f6dof(:,i),t6dof0(:,i)]=InverseDynamicsSolid(Human_model,external_forces_pred(i).fext,g,1);
    
    cost = alpha * sum(X(1:Pred_length).^2) + beta * sum([Human_model.torques].^2) + gamma * sum(X(1+Pred_length:Nb_muscles).^2); 
    
else
    cost = alpha * sum(X(1:Pred_length).^2) + gamma * sum(X(1+Pred_length:Nb_muscles).^2); 
end




return
