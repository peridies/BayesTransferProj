Read Me
To run the code for data analysis 
Please first download the following Matlab file exchange packages on Mathworks
1. Bayes factor https://www.mathworks.com/matlabcentral/fileexchange/69794-bayesfactor
2. violin plot https://www.mathworks.com/matlabcentral/fileexchange/45134-violin-plot
3. point_to_line_distance https://www.mathworks.com/matlabcentral/fileexchange/64396-point-to-line-distance
4. (elective) subtightplot  https://www.mathworks.com/matlabcentral/fileexchange/39664-subtightplot
5. (elective) hex2rgb https://www.mathworks.com/matlabcentral/fileexchange/46289-rgb2hex-and-hex2rgb

Data include raw and processed data
4 folders [exp1_discovery, exp1_validation, exp2_dsicovery, and exp2_validation] 
Each folder contains 
.csv : unprocessed raw data including likelihood-only task, learing and transfer phases 
.mat : processed experimental and simulated [linear-model and exemplar-model] data.  

Code
experiment 1 
0. input4statistic_PSPartLearn.m and input4statistic_PSTransfer.m: compute slope for each prior/likelihood combination & each participant  
1. SlopeAnalyses_PS.m: plot figure 3.a & 3.b and statistics
2. TransferScore_PS.m: plot figure 3.c & 3.d and statistics
3. Model_PS: compare Bayesian, linear, likelihood-only linear and exemplar (number of samples = 5) models. 

experiment 2 
0. input4statistic_IEPartLearn.m and input4statistic_IETransfer.m: compute slope for each prior/likelihood combination & each participant  
1. SlopeAnalyses_IE.m: plot figure 4.a & 4.b and statistics
2. TransferScore_IE.m: plot figure 4.c & 4.d and statistics
3. Model_IE: compare Bayesian, linear, likelihood-only linear and exemplar (number of samples = 5) models. 

%%other codes
1. SubjectiveVarGamma_PS.m & SubjectiveVarGamma_IE.m show how we calculated subject-specific prior variances  
2. Modelled exp1 data ('OldNewTable_Exemplar_P.mat') already provided. 
To create your own simulated exemplar model data. Use scripts and functions in the exemplar model folder    
 
 

