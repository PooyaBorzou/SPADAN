function [x,fval,exitflag,output] = optimizer_Ifmin5(x0)
%% This is an auto generated MATLAB file from Optimization Tool.
global figure_show
MaxFunEvals_Data = 10000;
MaxIter_Data = 10000;
%% Start with the default options
options = optimoptions('fminunc');
%% Modify options setting
options = optimoptions(options,'MaxFunEvals', MaxFunEvals_Data);
options = optimoptions(options,'MaxIter', MaxIter_Data);
if(figure_show)
options = optimoptions(options,'PlotFcns', {  @optimplotx @optimplotfunccount @optimplotfval @optimplotstepsize @optimplotfirstorderopt });
end
options = optimoptions(options,'Algorithm', 'quasi-newton');
[x,fval,exitflag,output] = ...
fminunc(@eval_sqr_fmin5,x0,options);
