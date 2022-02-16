function [c,ceq]=Prediction_activation_dependancy(X, Human_model, Prediction, i,  external_forces_pred,  Fmax, g, R, G,Fmmax)

c=[];



for k = 1:numel(Prediction)
    Prediction(k).efforts(i,1)=X(k)*Prediction(k).efforts_max(i,1);
    Prediction(k).efforts(i,2)=X(k+numel(Prediction))*Prediction(k).efforts_max(i,2);
    Prediction(k).efforts(i,3)=X(k+2*numel(Prediction))*Prediction(k).efforts_max(i,3);
end

external_forces_pred=addForces_Prediction_frame_par_frame(X,external_forces_pred,Prediction,Fmax,i);

[Human_model, ~ ,~]=InverseDynamicsSolid(Human_model,external_forces_pred(i).fext,g,1);


torques =[0 ,  Human_model.torques];

a = X(3*numel(Prediction) + 1:end);

ceq = G*(R.*Fmmax*a - torques');


end