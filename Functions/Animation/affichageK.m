Kt = rand(3,3);
[V,D] = eig(Kt);

F_ef=[];V_ef=[];
for i_for=1:3
        X_array = V(1,i_for);
        Y_array = V(2,i_for);
        Z_array = V(3,i_for);
        F_ef = [F_ef; [1 2]+size(V_ef,1)]; %#ok<AGROW>
        V_ef = [V_ef; [X_array' Y_array' Z_array']]; %#ok<AGROW>
end
if isfield(AnimateParameters,'Mode')  && (isequal(AnimateParameters.Mode, 'Figure') ...
        || isequal(AnimateParameters.Mode, 'GenerateParameters') ...
        || isequal(AnimateParameters.Mode, 'GenerateAnimate'))
    finv = figure('visible','off');
    Ext = gpatch(F_ef,V_ef,[],Colors.color_vect_force,1,4);
    copyobj(Ext,ax);
    close(finv);
elseif f==f_affich(1)
    Ext = gpatch(F_ef,V_ef,[],Colors.color_vect_force,1,4);
end
animStruct.Handles{f} = [animStruct.Handles{f} Ext];
animStruct.Props{f} = {animStruct.Props{f}{:},'Vertices'};
animStruct.Set{f} = {animStruct.Set{f}{:},V_ef};
    