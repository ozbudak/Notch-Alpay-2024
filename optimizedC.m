% Generate random perturbations within ±5%
pb1 = 1+(0.05 * (2*rand(1,64) - 1)); 
pb2 = 1+(0.05 * (2*rand(1,64) - 1));

% Define the objective function to be optimized using the genetic algorithm
getScore = @(parameters)getscoredC(parameters,pb1, pb2);

% Set options for the genetic algorithm (GA) to use parallel processing
options=optimoptions('ga','UseParallel',true);

% Define the lower bounds (lb) for the parameters
lb = [15 0.1 15 0.1 5  0.1 5  0.1 0.1 0.1 0.1 15 0.1 5  0.1 15 0.1 5  0.1 15 0.1 5  0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.003 0.0003 0.003 0.0003 0.003 0.0003 0.003 0.0003 0.003 0.0003 0.003 0.0003 0.003 0.0003 0.003 0.0003 0.003 0.0003 30   30   30   30   30   5  9  9  5  1 5  1 1 0  0  0  0];
% Define the upper bounds (ub) for the parameters
ub = [65 0.4 65 0.4 60 0.4 60 0.4 0.4 0.4 0.4 65 0.4 60 0.4 65 0.4 60 0.4 65 0.4 60 0.4 0.4 0.4 0.4 0.4 0.4 0.4 0.3   0.003  0.3   0.003  0.3   0.003  0.3   0.003  0.3   0.003  0.3   0.003  0.3   0.003  0.3   0.003  0.3   0.003  1500 1500 1500 1500 1500 12 27 27 12 2 12 2 2 5  5  5  5];

% Run the genetic algorithm (GA) to find the optimal parameters
% 'x' contains the optimal parameters found by the GA
% 'fval' is the value of the objective function at the optimal parameters
[x,fval]=ga(getScore,64,[],[],[],[],lb,ub,[],[],options);

%% The following part is not required; it is used to get multiple parameter sets.
% count=1;
% while count<=200 %set the number of parameter sets
%     [x,fval]=ga(getScore,64,[],[],[],[],lb,ub,[],[],options);
%     if fval<-10
%         name="setbest"+count+".mat";
%         save(name,'x','fval','pb1','pb2');
%         count=count+1;
%     end
% end
