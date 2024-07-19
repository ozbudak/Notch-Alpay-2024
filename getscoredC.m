function score=getscoredC(parameters, pb1, pb2)
% Function to evaluate the score of a parameter set for a model
% Inputs:
%   parameters - The parameters to be evaluated
%   pb1, pb2 - Perturbation factors for the model
% Output:
%   score - The score of the parameter set

score=0; % Initialize the score to zero

% Simulate the model with the given parameters
[t,mh1Matrix]=dCmodelPF(parameters,pb1,pb2);
% Check if her1 mRNA has sustained oscillations
if checkSusOsc(t,mh1Matrix(1,:))==true
    score=score+1;
end

% Check if her1 mRNA oscillation period is around 30 min
period=checkPeriod(t,mh1Matrix(1,:));
score=score+1/(abs(period-30)+1);

if score>1
    % Check the effect of removing deltaC
    paramC=parameters;
    paramC(1)=0;
    [t,mh1Matrix]=dCmodelPF(paramC,pb1,pb2);
    if checkSusOsc(t,mh1Matrix(1,:))==true
        score=score+1;
    end
    time1=syncBreak(mh1Matrix(1,:),mh1Matrix(2,:));

    % Check the effect of removing deltaD
    paramD=parameters;
    paramD(3)=0;
    [t,mh1Matrix]=dCmodelPF(paramD,pb1,pb2);
    if checkSusOsc(t,mh1Matrix(1,:))==true
        score=score+1;
    end
    time2=syncBreak(mh1Matrix(1,:),mh1Matrix(2,:));

    % Check if the synchrony break time in the deltaC mutant 
    % is earlier than that in the deltaD mutant
    if ~isnan(time1) && ~isnan(time2)
        if time1 <= time2-90
            score = score + 1;
        end
    end

    % Check the effect of removing deltaC, deltaD
    paramCD=paramC;
    paramCD(3)=0;
    [t,mh1Matrix]=dCmodelPF(paramCD,pb1,pb2);
    if checkSusOsc(t,mh1Matrix(1,:))==true
        score=score+1;
    end
    time3=syncBreak(mh1Matrix(1,:),mh1Matrix(2,:));

    % Check if the synchrony break time in the double mutant 
    % is close to that in the deltaC mutant
    if ~isnan(time1) && ~isnan(time3)
        if abs(time1-time3)<=5
            score = score + 1;
        end
    end
end

if score>6
    % Check the effect of removing Her1, Her7
    paramph1ph7=parameters;
    paramph1ph7(14)=0;
    paramph1ph7(18)=0;
    [t,mh1Matrix,mh7Matrix]=dCmodelPF(paramph1ph7,pb1,pb2);
    % Calculate the average combined mRNA level of her1 and her7
    mh17=mh1Matrix(1,:)+mh7Matrix(1,:);
    level0=(sum(mh17(1,101:end)))/200;
    score=score+1/(abs(level0-97)+1);

    % Check the effect of removing Her1, Her7, and deltaC
    paramph1ph7C=paramph1ph7;
    paramph1ph7C(1)=0;
    [t,mh1Matrix,mh7Matrix]=dCmodelPF(paramph1ph7C,pb1,pb2);
    mh17=mh1Matrix(1,:)+mh7Matrix(1,:);
    level1=(sum(mh17(1,101:end)))/200;
    if level0/level1>1.3 && level0/level1<1.4
        score=score+1;
    end

    % Check the effect of removing Her1, Her7, and deltaD
    paramph1ph7D=paramph1ph7;
    paramph1ph7D(3)=0;
    [t,mh1Matrix,mh7Matrix]=dCmodelPF(paramph1ph7D,pb1,pb2);
    mh17=mh1Matrix(1,:)+mh7Matrix(1,:);
    level2=(sum(mh17(1,101:end)))/200;
    if level0/level2>1.74 && level0/level2<1.84
        score=score+1;
    end

    % Check the effect of removing Her1, Her7, deltaC, and deltaD
    paramph1ph7CD=paramph1ph7;
    paramph1ph7CD(1)=0;
    paramph1ph7CD(3)=0;
    [t,mh1Matrix,mh7Matrix]=dCmodelPF(paramph1ph7CD,pb1,pb2);
    mh17=mh1Matrix(1,:)+mh7Matrix(1,:);
    level3=(sum(mh17(1,101:end)))/200;
    if level2/level3>0.95 && level2/level3<1.05
        score=score+1;
    end
end
score=-score;
end