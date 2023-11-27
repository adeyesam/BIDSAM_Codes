%% File for simulating model for any set of values of inputs
% Please ensure you have the following in the same directory as this file:
% 1) "simulate_model_dynamic" function file, 2) "BasisInterpretation"
% function file and 3) "solution_workspace" mat file.

% Step 1: declare the set of inputs (u) for which you desire to use the model
% to predict outputs. This can be for just a time instant or multiple
% instances e.g load('u.mat'), u=[2 5 5]'....

% Step 2: declare the rank (slnrank) of the model you desire to run from the set of
% hierarchically ranked models the algorithm has identified. e.g 
slnrank = 4;

% Step 3: declare initial values of output variables e.g
y_initial = data(:,1);

% Step 3: run simulation
y = simulate_model_dynamic(u, y_initial, slnrank);
