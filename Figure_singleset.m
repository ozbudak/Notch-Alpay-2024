% Load data from a single set file
data = load('setbest5.mat');

% Extract parameters and pb1, pb2 from the loaded data
parameters = data.x;
pb1 = data.pb1;
pb2 = data.pb2;

paramC = parameters;
paramC(1) = 0;
[t, mh1Matrix] = dCmodel(paramC, pb1, pb2);
time1 = syncBreak(mh1Matrix(1,:), mh1Matrix(2,:));

% Plot the her1 mRNA levels over time for Cell_1 and Cell_2
figure(1);
plot(t, mh1Matrix(1,:), 'b', 'DisplayName', 'Cell#1','LineWidth', 2);
hold on;
plot(t, mh1Matrix(2,:), 'r', 'DisplayName', 'Cell#2','LineWidth', 2);

% Add a vertical line to indicate the synchrony break time
xline(time1, '--k', 'DisplayName', 'Synchrony Break Time','LineWidth', 2);

% Add labels, title, and legend
xlabel('Time','FontSize', 14);
ylabel('her1 mRNA','FontSize', 14);

% Specify the position of the legend
legend('show', 'Location', 'northwest','FontSize', 10);

% Adjust the axis range
ylim([0 45]);
yticks([20 40]);
xticks([0 150 300]);

ax=gca;
ax.TickLength = [0 0];
ax.FontSize = 12;

hold off;

% Repeat for deltaD mutant and double mutant
paramD=parameters;
paramD(3)=0;
[t,mh1Matrix]=dCmodel(paramD,pb1,pb2);
time2=syncBreak(mh1Matrix(1,:),mh1Matrix(2,:));

figure(2);
plot(t, mh1Matrix(1,:), 'b', 'DisplayName', 'Cell#1','LineWidth', 2);
hold on;
plot(t, mh1Matrix(2,:), 'r', 'DisplayName', 'Cell#2','LineWidth', 2);

xline(time2, '--k', 'DisplayName', 'Synchrony Break Time','LineWidth', 2);

xlabel('Time','FontSize', 14);
ylabel('her1 mRNA','FontSize', 14);

legend('show', 'Location', 'northwest','FontSize', 10);

ylim([0 45]);
yticks([20 40]);
xticks([0 150 300]);

ax=gca;
ax.TickLength = [0 0];
ax.FontSize = 12;

hold off;

paramCD=paramC;
paramCD(3)=0;
[t,mh1Matrix]=dCmodel(paramCD,pb1,pb2);
time3=syncBreak(mh1Matrix(1,:),mh1Matrix(2,:));

figure(3);
plot(t, mh1Matrix(1,:), 'b', 'DisplayName', 'Cell#1','LineWidth', 2);
hold on;
plot(t, mh1Matrix(2,:), 'r', 'DisplayName', 'Cell#2','LineWidth', 2);

xline(time3, '--k', 'DisplayName', 'Synchrony Break Time','LineWidth', 2);

xlabel('Time','FontSize', 14);
ylabel('her1 mRNA','FontSize', 14);

legend('show', 'Location', 'northwest','FontSize', 10);


ylim([0 45]);
yticks([20 40]);
xticks([0 150 300]);

ax=gca;
ax.TickLength = [0 0];
ax.FontSize = 12;

hold off;

saveas(figure(1), 'Synchrony_deltaC.png');
saveas(figure(2), 'Synchrony_deltaD.png');
saveas(figure(3), 'Synchrony_double.png');


