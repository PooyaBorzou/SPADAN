function [x,fval,exitflag,output] = optimizer_simplex5(x0)
%% This is an auto generated MATLAB file from Optimization Tool.

%% Start with the default options
global figure_show
options = optimset;
%% Modify options setting
if(figure_show)
options = optimset(options,'PlotFcns', {  @optimplotx @optimplotfval });
end
options = optimset(options,'Display', 'iter');
options = optimset(options,'MaxIter', 80);

[x,fval,exitflag,output] = ...
fminsearch(@eval_sqr_fmin_v20,x0,options);
