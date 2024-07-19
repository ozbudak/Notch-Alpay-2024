function score=getscoredCNPF3cond(parameters, pb1, pb2)
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

% Check if her1 mRNA synchronization score is above 0.8
cell1=mh1Matrix(1,:);
cell2=mh1Matrix(2,:);
sync = corr(cell1', cell2');
if sync>=0.8
    score=score+1;
end
score=-score;
end