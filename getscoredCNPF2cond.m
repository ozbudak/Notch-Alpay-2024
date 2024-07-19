function score=getscoredCNPF2cond(parameters, pb1, pb2)
% Function to evaluate the score of a parameter set for No Positive Feedback model
% Inputs:
%   parameters - The parameters to be evaluated
%   pb1, pb2 - Perturbation factors for the model
% Output:
%   score - The score of the parameter set

score=0;% Initialize the score to zero

% Simulate the model with the given parameters
[t,mh1Matrix]=dCmodelNPF(parameters,pb1,pb2);
% Check if her1 mRNA has sustained oscillations
if checkSusOsc(t,mh1Matrix(1,:))==true
    score=score+1;
end

% Check if her1 mRNA oscillation period is around 30 min
period=checkPeriod(t,mh1Matrix(1,:));
score=score+1/(abs(period-30)+1);
score=-score;
end