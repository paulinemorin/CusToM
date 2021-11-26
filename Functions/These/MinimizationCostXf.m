function cost = MinimizationCostXf(X, X_f, BiomechanicalModel, Prediction, i,  external_forces_pred, Fmax, g, Nb_muscles , alpha_beta_gamma)
alpha = alpha_beta_gamma(1);
beta =  alpha_beta_gamma(2);
gamma =  alpha_beta_gamma(3);

Pred_length=2*numel(Prediction);
X_f=[X;X_f(2*numel(Prediction)+1:3*numel(Prediction))];

if beta
    Human_model = BiomechanicalModel.OsteoArticularModel;
    
    for k = 1:numel(Prediction)
        Prediction(k).efforts(i,1)=X_f(k)*Prediction(k).efforts_max(i,1);
        Prediction(k).efforts(i,2)=X_f(k+numel(Prediction))*Prediction(k).efforts_max(i,2);
        Prediction(k).efforts(i,3)=X_f(k+2*numel(Prediction))*Prediction(k).efforts_max(i,3);
    end
    
    external_forces_pred=addForces_Prediction_frame_par_frame(X_f,external_forces_pred,Prediction,Fmax,i);
    
    [Human_model,~,~]=InverseDynamicsSolid(Human_model,external_forces_pred(i).fext,g,1);
    
    cost = alpha * sum(X(1:Pred_length).^2) + beta * sum([Human_model.torques].^2) + gamma * sum(X(1+Pred_length:Nb_muscles).^2); 
    
else
    cost = alpha * sum(X(1:Pred_length).^2) + gamma * sum(X(1+Pred_length:Nb_muscles).^2); 
end




return
