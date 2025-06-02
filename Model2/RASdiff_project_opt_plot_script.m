%% Parameters optimization and simulations with optimized parameters
clc; clear all; close all

%% Load the project into SBPD toolbox
RASdiff_project = IQMprojectSB('RASdiff_project')

% Keep the original project unchanged
RASdiff_opt = RASdiff_project;

%% Select parameters to estimate and corresponding bounds
%  Global parameters (define only those which will be optimized)
%  Name  Lower bounds Upper bounds
boundFactor= 1 %3
boundFactorNarrow= 1.5
paramdata = {
   'mP'     0.58/boundFactorNarrow      0.58*boundFactorNarrow
   'mPTT'    0     1*boundFactor %remove?
   'mPcP'     0      1*boundFactor %remove?
   'mPcT'     0      1*boundFactor
   'mT'     0.68/boundFactorNarrow      0.68*boundFactorNarrow
   'mTDD'    0     1*boundFactor
   'mTcT'     0      1*boundFactor %remove?
   'mTcD'     0      1*boundFactor
   'd'     0       1*boundFactor
%    'fRasCil'    0     2 
   };

%  Local parameters
paramdatalocal = {
%    'mTD'     0       1*boundFactor %%% for Fig3 local param mTD
    };

%  Initial conditions
icdata = {
       'PROG'        50       150
       'TRANS'     4000     12000
       'DIFF'         0       150
       'PROGc'      250       750
       'TRANSc'    1000      3000
    };

% Choose number of iterations
repetition=100 %1/3/10/100

indExpToOpt=[1:2] %1:2
figureTitles={...
    'high serum only, Fig2F',...
    'diff with serum switch, Fig.2A',...
    'diff with serum switch, Fig.3A-I, CTRL'...
    'diff with serum switch, Fig.3A-C, siKRAS'...
    'diff with serum switch, Fig.3D-F, siNRAS'...
    'diff with serum switch, Fig.3G-I, siHRAS'...
    }

%% Define parameters to collect results
estim = []; % collect result from parameter estimations
optimalcost = []; % collect only optimal cost of each parameter set
optParams = []; % collect all parameter values of each set
optICs=nan(numel(indExpToOpt),size(icdata,1),repetition);

% Model and experiment settings
estim.modelindex=1;
estim.experiments.indices=[indExpToOpt];
estim.experiments.weight=ones(1,numel(indExpToOpt));

% Perform optimization (Global [pswarmSB] then Local [simplexSB])
for i=1:repetition
    
    disp(['iteration number = ' num2str((i))]) % counter number
    
    % Always start with default model setting
    RASdiff_opt = RASdiff_project;
    
    % Optimization parameters
    estim.integrator.options.abstol = 1e-006;
    estim.integrator.options.reltol = 1e-006;
    estim.integrator.options.minstep = 0;
    estim.integrator.options.maxstep = Inf;
    estim.integrator.options.maxnumsteps = 1000;
    
    estim.displayFlag = 2;
    estim.scalingFlag = 2; %0/2 %doc IQMparameterestimation 
    estim.timescalingFlag = 0;
    estim.initialconditionsFlag = 0; %1
    
    estim.parameters = paramdata;
    estim.parameterslocal = paramdatalocal;
    estim.initialconditions = icdata;
    
    % Method of optimization and number of maxfunevals [Global]
    estim.optimization.method = 'pswarmIQM';
    estim.optimization.options.maxfunevals = 10000*1;
    
    % Run estimation and get optimized project [pswarmSB] & get result
    output = IQMparameterestimation(RASdiff_opt,estim);
    RASdiff_opt = output.projectopt;

    % Method of optimization and number of maxfunevals [Local]
    estim.optimization.method = 'simplexIQM';
    estim.optimization.options.maxfunevals = 3000*1;

    % Run estimation and get optimized project [simplexSB] & get result
    output = IQMparameterestimation(RASdiff_opt,estim);
    RASdiff_opt = output.projectopt;
    
    % Collect results just only FVALopt and optimized parameter values
    optimalcost=[optimalcost; output.FVALopt];
    optParams=[optParams; output.Popt];
    optICs(:,:,i)=output.ICopt;
end

%% Plot selected optimal solutions
% Set desirable threshold
threshold=0.10 %inf %0.05

for counter2=1:numel(indExpToOpt)
% Get model structure and experiment then merge them together
model_structure = IQMgetmodel(RASdiff_opt);
experiment = IQMgetexperiment(RASdiff_opt,indExpToOpt(counter2)); % Select measurement (only 1)
model = IQMmergemodexp(model_structure, experiment);

% Define states to plot 
plotstates=[1:numel(IQMstates(model))]; % plot all 3 states
% plotvariables=[2,4:numel(IQMvariables(model))]; % plot all variables except serumSwitchTime + sumPTD
[temp,plotvariables]=setdiff(IQMvariables(model),{'serumSwitch','rasDiff','sumPTD'},'stable');

% Plot figures using all parameters which pass the threshold

figure(counter2)

for j=1:length(optimalcost)
    
    disp(['Plotting parameters set = ' num2str((j))]); % counter number
    
    if (optimalcost(j)<threshold) % if optimalcost pass threshold
        
        % plug in parameters and simulate with time > 120 min (300 min in this setting)
        model2 = IQMparameters(model,output.parameters,optParams(j,:));
        model2 = IQMinitialconditions(model2,optICs(counter2,:,j));
        simulation = IQMsimulate(model2,5);
        
        % start plotting each state for each parameter set
        for k=1:length(plotstates);
            plotstate=plotstates(k);
            subplot(3,4,k), hold on, plot(simulation.time,simulation.statevalues(:,plotstate),'LineWidth',2),
            legend(simulation.states(plotstate)), hold off
            hold on
        end
        for k2=1:length(plotvariables);
            plotstate=plotvariables(k2);
            subplot(3,4,k+k2), hold on, plot(simulation.time,simulation.variablevalues(:,plotstate),'LineWidth',2),
            legend(simulation.variables(plotstate)), hold off
            hold on
        end
        
    elseif (optimalcost(j)>=threshold) % if optimalcost didn't pass threshold
        disp(['Not included'])
    end
    
end

% Add the marker and error from experimental data to compare with the simulation

% Extract measurement data from datasheet inside the project
[time,componentNames,values,minvalues,maxvalues]=IQMmeasurementdata(IQMgetmeasurement(RASdiff_project,indExpToOpt(counter2)));

% Plot measurement data with range of S.D. as errorbars on the plot
for counter=1:numel(componentNames)
    if isempty(find(ismember(simulation.variables,componentNames(counter))))
       temp=find(ismember(simulation.states,componentNames(counter)));
    else
%         temp=find(ismember(simulation.variables,componentNames(counter)))+k-2; %-1
        temp=find(ismember(simulation.variables,componentNames(counter)))+k-numel(IQMvariables(model))+numel(plotvariables); %-1
    end
    subplot(3,4,temp), hold on, errorbar(time,values(:,counter),values(:,counter)-minvalues(:,counter),maxvalues(:,counter)-values(:,counter),'r.','MarkerSize',20), hold off
hold on 

sgtitle(figureTitles(indExpToOpt(counter2)))
end

end

%% param & ICs stats
figure
boxplot(optParams(optimalcost<=threshold,:))
ylabel({'parameter value'},'fontweight','bold','fontsize',12)
set(gca,'xtick',1:size(optParams,2),'xticklabel',paramdata(:,1),'fontweight','bold','fontsize',12);
xtickangle(45)

for counter=1:numel(indExpToOpt)
    figure
    temp=optICs(:,:,optimalcost<=threshold);
    temp=squeeze(temp(counter,:,:))';
    boxplot(temp)
    ylabel({'initial condition value'},'fontweight','bold','fontsize',12)
    set(gca,'xtick',1:size(optParams,2),'xticklabel',icdata(:,1),'fontweight','bold','fontsize',12);
    xtickangle(45)
    title(['Experiment: ' num2str(indExpToOpt(counter))])
end

%%
optimalcost_sel=optimalcost(optimalcost<=threshold,:)
optParams_sel=optParams(optimalcost<=threshold,:)
% disp(' ')
% disp('Symmetrical self-renewal: PROG -> PROG + PROG')
% disp('Asymmetrical dividion: PROG -> PROG + TRANS')
% disp(' ')
% disp('Fraction of asymmetrical divisions: xPT')
% disp(' ')
% disp('Parameter name:')
% disp(paramdata(:,1)')
% disp('Mean:')
% disp(mean(optParams_sel))
% 
nrDatapoints=0;
for counter=1:numel(indExpToOpt)
   temp2=IQMmeasurementdata(IQMgetmeasurement(RASdiff_opt,indExpToOpt(counter)));
   nrDatapoints=nrDatapoints+numel(temp2);
end
% disp('Standard error of the mean:')
% disp(std(optParams_sel)/sqrt(nrDatapoints))

% AIC
p=size(optParams_sel,2) %number of parameters

% AIC example
% err=[0.1 0.15 0.2 0.01 0.001]
% p=numel(optParams) %number of parameters
% n=numel(err)
% AIC = n*log(sum(err.*err)/n) + 2*p
% AICc = AIC + 2*p*(p+1)/(n-p-1)

AIC = nrDatapoints*log(optimalcost_sel) + 2*p
% correction for small sample sizes
AICc = AIC + 2*p*(p+1)/(nrDatapoints-p-1)
figure
hist(AIC)
title('AIC')

% Then the quantity exp((AICmin ? AICi)/2) can be interpreted as being
% proportional to the probability that the ith model minimizes the
% (estimated) information loss.
% As an example, suppose that there are three candidate models, whose AIC
% values are 100, 102, and 110. Then the second model is 
% exp((100 - 102)/2) = 0.368 times as probable as the first model to 
% minimize the information loss.
AICcs=AICc;
% AICcs=[-15.8 -12.6]
temp = exp((min(AICcs) - AICcs)/2);
resProb = temp ./ sum(temp)

%% reporting
disp('min(optimalcost)')
min(optimalcost)
disp('min(AIC)')
min(AIC)
disp('min(AICc)')
min(AICc)
