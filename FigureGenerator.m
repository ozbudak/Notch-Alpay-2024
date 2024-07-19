% Initialize arrays
time1_array = zeros(1, 10);
time2_array = zeros(1, 10);
time3_array = zeros(1, 10);
level1_array = zeros(1, 10);
level2_array = zeros(1, 10);
level3_array = zeros(1, 10);

N=10; % number of parameter sets you have

for i = 1:N
    % Load the .mat file
    filename = sprintf('setbest%d.mat', i);
    data = load(filename);
    
    % Extract x from the loaded data (assuming x is the variable name)
    parameters = data.x;
    pb1 = data.pb1;
    pb2 = data.pb2;
    
    % Find deltaC mutant synchrony break time 
    paramC = parameters;
    paramC(1) = 0;
    [t, mh1Matrix] = dCmodel(paramC, pb1, pb2);
    time1 = syncBreak(mh1Matrix(1,:), mh1Matrix(2,:));
    time1_array(i) = time1;
  
    % Find deltaD mutant synchrony break time 
    paramD = parameters;
    paramD(3) = 0;
    [t, mh1Matrix] = dCmodel(paramD, pb1, pb2);
    time2 = syncBreak(mh1Matrix(1,:), mh1Matrix(2,:));
    time2_array(i) = time2;

    % Find deltaC, deltaD double mutant synchrony break time 
    paramCD=paramC;
    paramCD(3)=0;
    [t,mh1Matrix]=dCmodel(paramCD,pb1,pb2);
    time3=syncBreak(mh1Matrix(1,:),mh1Matrix(2,:));
    time3_array(i) = time3;

    % Find deltaC-triple mutant average mRNA level
    paramph1ph7=parameters;
    paramph1ph7(14)=0;
    paramph1ph7(18)=0;
    paramph1ph7C=paramph1ph7;
    paramph1ph7C(1)=0;
    [t,mh1Matrix,mh7Matrix]=dCmodel(paramph1ph7C,pb1,pb2);
    mh17=mh1Matrix(1,:)+mh7Matrix(1,:);
    level1=(sum(mh17(1,101:end)))/200;
    level1_array(i) = level1;

    % Find deltaD-triple mutant average mRNA level
    paramph1ph7D=paramph1ph7;
    paramph1ph7D(3)=0;
    [t,mh1Matrix,mh7Matrix]=dCmodel(paramph1ph7D,pb1,pb2);
    mh17=mh1Matrix(1,:)+mh7Matrix(1,:);
    level2=(sum(mh17(1,101:end)))/200;
    level2_array(i) = level2;

    % Find quadruple mutant average mRNA level
    paramph1ph7CD=paramph1ph7;
    paramph1ph7CD(1)=0;
    paramph1ph7CD(3)=0;
    [t,mh1Matrix,mh7Matrix]=dCmodel(paramph1ph7CD,pb1,pb2);
    mh17=mh1Matrix(1,:)+mh7Matrix(1,:);
    level3=(sum(mh17(1,101:end)))/200;
    level3_array(i) = level3;
end

% Calculate the mean and standard error of time1_array
mean_time1 = mean(time1_array);
se_time1 = std(time1_array) / sqrt(length(time1_array));

% Calculate the mean and standard error of time2_array
mean_time2 = mean(time2_array);
se_time2 = std(time2_array) / sqrt(length(time2_array));

% Calculate the mean and standard error of time3_array
mean_time3 = mean(time3_array);
se_time3 = std(time3_array) / sqrt(length(time3_array));

% Calculate the mean and standard error of level1_array
mean_level1 = mean(level1_array);
se_level1 = std(level1_array) / sqrt(length(level1_array));

% Calculate the mean and standard error of level2_array
mean_level2 = mean(level2_array);
se_level2 = std(level2_array) / sqrt(length(level2_array));

% Calculate the mean and standard error of level3_array
mean_level3 = mean(level3_array);
se_level3 = std(level3_array) / sqrt(length(level3_array));

% Create a bar plot with error bars
figure(1);
hold on;
% Create bar objects
b1 = bar(1, mean_time1);
b2 = bar(2, mean_time2);
b3 = bar(3, mean_time3);

% Set bar colors
b1.FaceColor = 'w'; % white
b2.FaceColor = [0.8 0.8 0.8]; % light grey
b3.FaceColor = [0.5 0.5 0.5]; % dark grey

errorbar(1, mean_time1, 2*se_time1, 'k', 'LineStyle', 'none');
errorbar(2, mean_time2, 2*se_time2, 'k', 'LineStyle', 'none');
errorbar(3, mean_time3, 2*se_time3, 'k', 'LineStyle', 'none');
hold off;

% Customize the plot
set(gca, 'XTick', [1 2 3], 'XTickLabel', {'deltaC^{-/-}', 'deltaD^{-/-}','deltaC^{-/-};deltaD^{-/-}'});
yticks([0 150 300]);
ax=gca;
ax.TickLength = [0 0];
ax.FontSize = 12;
xlabel('Mutant','FontSize', 14);
ylabel('Synchrony Break Time','FontSize', 14)
box on;


% Create a bar plot with error bars
figure(2);
hold on;
% Create bar objects
b1 = bar(1, mean_level1);
b2 = bar(2, mean_level2);
b3 = bar(3, mean_level3);

% Set bar colors
b1.FaceColor = 'w'; % white
b2.FaceColor = [0.8 0.8 0.8]; % light grey
b3.FaceColor = [0.5 0.5 0.5]; % dark grey

errorbar(1, mean_level1, 2*se_level1, 'k', 'LineStyle', 'none');
errorbar(2, mean_level2, 2*se_level2, 'k', 'LineStyle', 'none');
errorbar(3, mean_level3, 2*se_level3, 'k', 'LineStyle', 'none');
hold off;

% Customize the plot
set(gca, 'XTick', [1 2 3], 'XTickLabel', {'deltaC^{-/-}', 'deltaD^{-/-}','deltaC^{-/-};deltaD^{-/-}'});
yticks([0 40 80]);
ax=gca;
ax.TickLength = [0 0];
ax.FontSize = 12;
xlabel('Mutant','FontSize', 14);
ylabel('her1+her7 mRNA','FontSize', 14)
box on;

saveas(figure(1), 'Average_Synchrony Break.png');
saveas(figure(2), 'Average_mRNA Level.png');