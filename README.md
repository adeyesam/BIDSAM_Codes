# BML_Codes
This repository contains MATLAB codes of Bayesian ML algorithm developed by Adeyemo and Bhattacharyya (2023) for identification and simulation of sparse steady state and dynamic data-driven models. 

These codes will come in versions with continuos efforts to enhance user friendliness.

# Sample data
Sample steady state and dynamic data are uploaded in the folders "Sample Steady State Data" and "Sample Dynamic Data".
In each MATLAB data file named 'dataset' are the measured output variables named as 'data' and input variables named as 'u'.
In each of these, the rows of 'u' are the input variables while the columns are the time instants. 
![image](https://github.com/adeyesam/BML_Codes/assets/148823677/aaf28c87-3d13-41ed-bb01-2583af4be74b)


The measured outputs come in similar fashion.


# BML_V0
## Model training
The primary file for training a model (steady state/dynamic) is the Run_main file. At this point the user may wish to specify
1. "est_sel": indicies of available data to be used for model training. The default setting for steady state modeling is 70% randomly selected points from available data while for dynamic modeling, it is set to 1-800. This can be adjusted depending on the amount of data considered.

2. "Nworker": the number or workers to be employed for parallel computation in the Branch and bound algorithm. This depends on the choice of the user and the number of available processors on the machine used for running the codes. The default here is set to 8.

3. "Ndes": the number of desired top rank models the user desires the algorithm to return after searching by the Branch and bound algorithm.

4. The user may rarely need to consider changing the value of "R" in the "EstimateModel" function file. This determines the relative weighting of penalty due to model error and model complexity. For varying cases (both synthetic and real plant data) for which this algorithm has been tested, values in the range of 0.04 - 0.06 have yielded satisfactory results.

## Model simulation
The function files "simulate_model_dynamic" and "simulate_model_steady" will facilitate the interpretation and simulation of the hierarchichally ranked models obtained from the algorithm automatically saved as "solution_workspace" in the current directory.

For the steady state model, the user needs to specify the values of the input values for which prediction of the output variables are desired and the rank of the model that is desired to be tested (specified as "solnrank"). The file can give both short-term (one time instant) and long-term (multiple time instances) predictions as may be needed by the user and specified by the values of "inputs' fed in while calling the function file as follows:

y = simulate_model_steady(inputs, solnrank)

Simulation of the dynamic model follows similar manner except for the need to specify the initial conditions as "y_initial" while calling the function file as follows:

y = simulate_model_dynamic(inputs, y_initial, solnrank)

### Important!
1. The top rank model may not necessarily give the best fit as the criteria for ranking models seek a balance between model fitness (as measured by log likelihood function), model size and complexity (measured as a function of the number of parameters and their estimated covariance). Hence the user may need to examine a few of the top rank models to select the most suitable model based on the model accuracy and desired interpretability as reported by the chosen basis functions.

2. The model training is done with normalized values. Hence, during simulation, the input variables are first scaled while the normalized predictions are converted back to the normal scale before reporting back. The normalization is done based on the "buffered" range of each input and output variables fed into the algorithm.


## Reference
Adeyemo S, Bhattacharyya D “Optimal Dynamic Model Selection and Bayesian Parameter Estimation for Nonlinear Systems”, Comput. Chem. Eng., (Reviewed)
