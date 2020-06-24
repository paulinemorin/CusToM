function [varargout] = PlotAnimation(ModelParameters, AnimateParameters)
% Generation of an animation
%
%   INPUT
%   - ModelParameters: parameters of the musculoskeletal model,
%   automatically generated by the graphic interface 'GenerateParameters' 
%   - AnimateParameters: parameters of the animation, automatically
%   generated by the graphic interface 'GenerateAnimate'
%________________________________________________________
%
% Licence
% Toolbox distributed under GPL 3.0 Licence
%________________________________________________________
%
% Authors : Antoine Muller, Charles Pontonnier, Pierre Puchaud and
% Georges Dumont
%________________________________________________________

[filename,DataXSens,q,q6dof,Markers_set,Muscles,real_markers,PelvisPosition,PelvisOrientation,EnableModel,Human_model,AnalysisParameters,InverseKinematicsResults,ExperimentalData,BiomechanicalModel]=DataExtractionForAnimation(AnimateParameters,ModelParameters);

% AnimateParameters
options=OptionsChoices(BiomechanicalModel,AnimateParameters);

% Preliminary computations
if (options.seg_anim || options.AnatLandmark || isfield(options,'num_mus')) && ~DataXSens %anatomical position where other segments are attached
    [Human_model] = anat_position_solid_repere(Human_model,find(~[Human_model.mother]));
end


%% Mesh of bones
if options.bone_anim
    [Human_model]=BoneView(DataXSens,Human_model,BiomechanicalModel,ModelParameters,options.Segment);
end


%% Figure

if isfield(AnimateParameters,'Mode')  && (isequal(AnimateParameters.Mode, 'Figure') ...
        || isequal(AnimateParameters.Mode, 'Picture'))
    fig=figure('outerposition',[483,60,456*1.5,466*1.5]);
    ax=gca;
elseif (isfield(AnimateParameters,'Mode')  &&  isequal(AnimateParameters.Mode, 'cFigure')) 
    fig=cFigure; % from GIBBON
    view(3); axis equal; axis tight; axis vis3d; grid on; box on;
    camlight headlight; axis off; axis manual; lighting gouraud;
    ax=gca;
    ax.Clipping = 'off';
    drawnow;
elseif isfield(AnimateParameters,'Mode')  &&  (isequal(AnimateParameters.Mode, 'GenerateAnimate') || isequal(AnimateParameters.Mode, 'GenerateParameters'))
    ax = AnimateParameters.ax; 
    fig=ax.Parent;
    camlight(ax, 'headlight'); lighting(ax,'gouraud');
%     material(ax, 'metal');
end

% Frames to display
if isfield(AnimateParameters,'Mode') && (isequal(AnimateParameters.Mode, 'Picture') ...
        || isequal(AnimateParameters.Mode, 'GenerateAnimate') ...
        || isequal(AnimateParameters.Mode, 'GenerateParameters'))
    f_affich = AnimateParameters.PictureFrame;
else
    f_affich = 1:1:size(q,2);
end

%Initialization animStruct
animStruct=struct();
if (isfield(AnimateParameters,'Mode')  && ~isequal(AnimateParameters.Mode, 'GenerateParameters') &&...
        isfield(AnimateParameters,'Noc3d') &&  ~AnimateParameters.Noc3d)
    animStruct.Time=ExperimentalData.Time;
end

animStruct.Handles=cell(1,size(q,2));
animStruct.Props=cell(1,size(q,2));
animStruct.Set=cell(1,size(q,2));


%% Animation frame by frame

[animStruct]=AnimationFramebyFrame(ax,fig,filename,AnalysisParameters,ModelParameters,AnimateParameters,DataXSens,q,q6dof,PelvisPosition,PelvisOrientation,Markers_set,f_affich,Muscles,animStruct,real_markers,BiomechanicalModel,Human_model);












if isfield(AnimateParameters,'Mode')  && isequal(AnimateParameters.Mode, 'Figure')
    close all
    v=VideoWriter([filename '.avi']);
    v.FrameRate=1/(3*ExperimentalData.Time(2));
    open(v)
    writeVideo(v,M);
    close(v)
elseif (isfield(AnimateParameters,'Mode')  && isequal(AnimateParameters.Mode, 'cFigure') ) 
    anim8(fig,animStruct);
end

% varargout
if isfield(AnimateParameters,'Mode')  && isequal(AnimateParameters.Mode, 'GenerateParameters')
    varargout{1} = Human_model;
    varargout{2} = Markers_set;
    varargout{3} = EnableModel;
end
% save('L','L')




end

