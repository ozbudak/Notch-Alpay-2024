%% Optimization
% Generate random perturbations within ±5%
pb1 = 1+(0.05 * (2*rand(1,64) - 1));  
pb2 = 1+(0.05 * (2*rand(1,64) - 1));

% Set options for the genetic algorithm (GA) to use parallel processing
options=optimoptions('ga','UseParallel',true);

% Define the lower bounds (lb) for the parameters
lb=[16.7284566166532	0.290979002611650	15	0.189097440707135	5	0.100000000000000	5	0.187693641091572	0.184270218704426	0.264939739969374	0.252538385744830	19.1740729120170	0.215671347967585	37.4720925024376	0.100000000000000	15	0.294716105834435	5	0.153611820664325	15	0.141566996608334	6.46188494579985	0.100000000000000	0.270579197851299	0.100000000000000	0.359328256234682	0.261287121721257	0.100000000000000	0.377782142213363	0.00300000000000000	0.000300000000000000	0.0941343904350536	0.00264326449705189	0.0300814734561588	0.00230296687621812	0.00300000000000000	0.00151774014244727	0.00300000000000000	0.000300000000000000	0.00300000000000000	0.000423731838358995	0.167808199755981	0.00117616966402214	0.0258463289929095	0.00146443171284367	0.00300000000000000	0.000300000000000000	44.4185745210511	30	700.772612086524	477.922505191913	1216.06119120049	6.75986054361099	9.39395706817813	23.7951217078933	5	1.15394390544709	9.87547785266050	1.09845141483228	1	0.360570302872564	0	1.42715340200494	0];
% Define the upper bounds (ub) for the parameters
ub=[26.0796959139297	0.400000000000000	18.7265674491828	0.400000000000000	60	0.400000000000000	31.4416119980394	0.400000000000000	0.400000000000000	0.275642816377563	0.372011287101638	65	0.400000000000000	60	0.400000000000000	47.0899848335048	0.400000000000000	60	0.400000000000000	65	0.383141472929458	48.7406513073387	0.282158344778565	0.313344108600953	0.215440683727082	0.384725030217659	0.386818211088618	0.400000000000000	0.400000000000000	0.300000000000000	0.00300000000000000	0.300000000000000	0.00286881311818320	0.0482729824177737	0.00284412434377824	0.300000000000000	0.00300000000000000	0.132595919888828	0.00300000000000000	0.300000000000000	0.00203121768136093	0.300000000000000	0.00224702057593868	0.270383365066799	0.00261752083851082	0.300000000000000	0.00104545292651286	55.7091773895632	1500	1130.43536995080	817.695088946271	1480.51249163070	11.2231060728182	16.6063229307699	27	7.47606881518977	1.15849416565971	12	2	2	1.80992464076270	0.364396538435602	2.97870968416463	4.69782667331410];

% Define the objective function to be optimized using the genetic algorithm
getScore = @(parameters)getscoredC2cond(parameters,pb1, pb2);

% Submit 50 optimization rounds for Model with Positive Feedback
count=1;
i=0;
sync=zeros(1,2);
while count<=50
   [x,fval,exitflag,output,population,scores]=ga(getScore,64,[],[],[],[],lb,ub,[],[],options);
   if fval==-2
       i=i+1;
       [t,mh1Matrix]=dCmodelPF(x,pb1,pb2);
       cell1=mh1Matrix(1,:);
       cell2=mh1Matrix(2,:);
       sync(1,i) = corr(cell1', cell2'); % Store Synchronization Score
   end
   count=count+1;
end

% Repeat for Model without Positive Feedback
getScore1 = @(parameters)getscoredCNPF2cond(parameters,pb1, pb2);
count1=1;
j=0;
sync1=zeros(1,2);
while count1<=50
    [x,fval,exitflag,output,population,scores]=ga(getScore1,64,[],[],[],[],lb,ub,[],[],options);
    if fval==-2
       j=j+1;
       [t,mh1Matrix]=dCmodelNPF(x,pb1,pb2);
       cell1=mh1Matrix(1,:);
       cell2=mh1Matrix(2,:);
       sync1(1,j) = corr(cell1', cell2');
   end
    count1=count1+1;
end
%% Generate Figures & T-test
% Calculate the mean and standard error
mean_sync = mean(sync);
se_sync = std(sync) / sqrt(length(sync));

mean_sync1 = mean(sync1);
se_sync1 = std(sync1) / sqrt(length(sync1));

% Conduct t-test
[h,p,ci,stats] = ttest2(sync, sync1);

% Display t-test result
fprintf('T-test result:\n');
fprintf('h = %d\n', h);
fprintf('p = %.4f\n', p);
fprintf('ci = [%.4f, %.4f]\n', ci(1), ci(2));
fprintf('tstat = %.4f\n', stats.tstat);
fprintf('df = %d\n', stats.df);

% Create a bar plot with error bars
figure(1);
hold on;
% Create bar objects
b1 = bar(1, mean_sync);
b2 = bar(2, mean_sync1);

% Set bar colors
b1.FaceColor = 'w'; % white
b2.FaceColor = [0.8 0.8 0.8]; % light grey

errorbar(1, mean_sync, 2*se_sync, 'k', 'LineStyle', 'none');
errorbar(2, mean_sync1, 2*se_sync1, 'k', 'LineStyle', 'none');
hold off;

% Customize the plot
set(gca, 'XTick', [1 2], 'XTickLabel', {'Positive Feedback', 'No Positive Feedback'});
yticks([0 0.25 0.5 0.75 1]);
ax=gca;
ax.TickLength = [0 0];
ax.FontSize = 12;
xlabel('Models','FontSize', 14);
ylabel('Synchrony Score','FontSize', 14)
box on;

saveas(figure(1), '2 condition.png');
name="Pass_ttest_2.mat";
save(name,'i','j','h','p','ci','stats');