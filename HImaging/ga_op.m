function [x,fval,exitflag,output,population,score] = ga_op(nvars,lb,ub,PopulationSize_Data,MaxGenerations_Data)
%% This is an auto generated MATLAB file from Optimization Tool.

%% Start with the default options
options = optimoptions('ga');
%% Modify options setting
options = optimoptions(options,'PopulationSize', PopulationSize_Data);
options = optimoptions(options,'MaxGenerations', MaxGenerations_Data);
options = optimoptions(options,'FitnessScalingFcn', {  @fitscalingshiftlinear [] });
options = optimoptions(options,'SelectionFcn', @selectionroulette);
options = optimoptions(options,'CrossoverFcn', @crossoverscattered);
options = optimoptions(options,'MutationFcn', @mutationadaptfeasible);
options = optimoptions(options,'Display', 'off');
options = optimoptions(options,'PlotFcn', { @gaplotbestf });
options = optimoptions(options,'UseVectorized', false);
options = optimoptions(options,'UseParallel', true);
% options = optimoptions('ga', 'Display', 'off', 'OutputFcn', @saveB);

% options = optimoptions('ga', 'OutputFcn', @myoutputfcn);

[x,fval,exitflag,output,population,score] = ...
ga(@optimizeA1,nvars,[],[],[],[],lb,ub,[],[],options);

% 


end
