%% Plotting of final fits with optimized parameters
clc; clear all; close all

% simple model
load('WORKSPACE260924_simpleModel.mat')
threshold=0.041
dataSet=1

% full model: estimate fAsym
load('WORKSPACE_exp12_fAsym_100.mat')
threshold=min(optimalcost)*4/3 %threshold=0.06
dataSet=2

% full model: exp3-6
load('WORKSPACE_exp3456_fAsym_ICs_100.mat')
threshold=0.05
dataSet=3

simulationTime=[0:0.05:5];

%% Plot selected optimal solutions
% Set desirable threshold
% threshold=0.10 %inf %0.05

% Define states to plot
model_structure = IQMgetmodel(RASdiff_opt);
experiment = IQMgetexperiment(RASdiff_opt,indExpToOpt(counter2)); % Select measurement (only 1)
model = IQMmergemodexp(model_structure, experiment);
plotstates=[1:numel(IQMstates(model))]; % plot all 3 states
% plotvariables=[2,4:numel(IQMvariables(model))]; % plot all variables except serumSwitchTime + sumPTD
[plotVariableNames,plotvariables]=setdiff(IQMvariables(model),{'serumSwitch','rasDiff','sumPTD'},'stable');
% IQMvariables(model)r

plotDataAll=zeros(numel(simulationTime),numel(plotstates)+numel(plotvariables),sum(optimalcost<threshold),numel(indExpToOpt));
% plotData2=plotData1;
for counter2=1:numel(indExpToOpt)
    counter3=0;
    % Get model structure and experiment then merge them together
    model_structure = IQMgetmodel(RASdiff_opt);
    experiment = IQMgetexperiment(RASdiff_opt,indExpToOpt(counter2)); % Select measurement (only 1)
    model = IQMmergemodexp(model_structure, experiment);
%     IQMvariables(model)
    [~,~,plotvariables]=intersect(plotVariableNames,IQMvariables(model),'stable');

    % Plot figures using all parameters which pass the threshold
    
    figure(counter2)
    
    for j=1:length(optimalcost)
        
        disp(['Plotting parameters set = ' num2str((j))]); % counter number
        
        if (optimalcost(j)<threshold) % if optimalcost pass threshold
            counter3=counter3+1;
            
            % plug in parameters and simulate with time > 120 min (300 min in this setting)
            model2 = IQMparameters(model,output.parameters,optParams(j,:));
            if exist('optICs')
                model2 = IQMinitialconditions(model2,optICs(counter2,:,j));
            end
            %             simulation = IQMsimulate(model2,5);
            simulation = IQMsimulate(model2,simulationTime); %simulation time
            [C,IA,IB]=intersect(simulationTime,simulation.time);
            plotDataAll(:,:,counter3,counter2)=[simulation.statevalues(IB,plotstates), simulation.variablevalues(IB,plotvariables)];
            
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

size(plotDataAll)

%% median + 68% confidence interval
% close all

fontSize=14
lineWidth=1
printAs='-djpeg'  %'-depsc' %'-djpeg'
titleOn=1

for counter=1:numel(indExpToOpt)
    plotData=plotDataAll(:,:,:,counter);
    mPlotData=median(plotData,3);
    
    figure('color','white')
    for counter6=[10:1:68] %conf levels
        %             fadingColor=[1 (counter6+20)/100 (counter6+20)/100];
        fadingColor=[1 (counter6-10)/58*0.8+0.1 (counter6-10)/58*0.8+0.1]; %for range [10, 68]
        lbPlotData=prctile(plotData,50-counter6/2,3); %68% range of fits / predictions (similar standard deviation)
        ubPlotData=prctile(plotData,50+counter6/2,3);
        for counter2=1:length(plotstates)
            plotState=plotstates(counter2);
            subplot(3,4,counter2), hold on, plot(simulationTime,lbPlotData(:,counter2),'-','LineWidth',lineWidth,'Color',fadingColor), hold off
            subplot(3,4,counter2), hold on, plot(simulationTime,ubPlotData(:,counter2),'-','LineWidth',lineWidth,'Color',fadingColor), hold off
        end
        for counter3=1:length(plotvariables)
            plotState=plotvariables(counter3);
            subplot(3,4,counter2+counter3), hold on, plot(simulationTime,lbPlotData(:,counter2+counter3),'-','LineWidth',lineWidth,'Color',fadingColor), hold off
            subplot(3,4,counter2+counter3), hold on, plot(simulationTime,ubPlotData(:,counter2+counter3),'-','LineWidth',lineWidth,'Color',fadingColor), hold off
        end
    end
    for counter2=1:length(plotstates)
        plotState=plotstates(counter2);
        subplot(3,4,counter2), hold on, plot(simulationTime,mPlotData(:,plotState),'r--','LineWidth',lineWidth*3)
        if titleOn, title(simulation.states(plotState)), end
%         set(gca,'FontSize',fontSize,'XTickLabel',{'-2','0','2','4','6','8'}), hold off
        set(gca,'FontSize',fontSize), hold off
        temp=axis; temp(2)=5; axis(temp)
    end
    for counter3=1:length(plotvariables)
        plotState=plotvariables(counter3);
        subplot(3,4,counter2+counter3), hold on, plot(simulationTime,mPlotData(:,counter2+counter3),'r--','LineWidth',lineWidth*3)
        if titleOn, title(simulation.variables(plotState)), end
        set(gca,'FontSize',fontSize), hold off
        temp=axis; temp(2)=5;
        if dataSet==1
            if counter==2
                if plotState==6
                    temp(3)=70;
                elseif plotState==7
                    temp(4)=30;
                end
            end
        end
        if dataSet==2
            if counter==2
                if plotState==6
                    temp(4)=2.5;
                elseif plotState==8
                    temp(4)=30;
                elseif plotState==9
                    temp(4)=0.8;
                end
            end
        end
        if dataSet==3
            if counter==1
                if plotState==7
                    temp(3)=0.7;
                end
            elseif counter==2
                if plotState==7
                    temp(3)=0.45;
                elseif plotState==8
                    temp(3)=75;
                elseif plotState==9
                    temp(4)=22;
                end
            elseif counter==4
                if plotState==7
                    temp(3)=0.65;
                elseif plotState==8
                    temp(3)=85;
                end
            end
        end
        axis(temp)
    end
    
    % Plot measurement data with range of S.D. as errorbars on the plot
    [time,componentNames,values,minvalues,maxvalues]=IQMmeasurementdata(IQMgetmeasurement(RASdiff_project,indExpToOpt(counter)));
    for counter4=1:numel(componentNames)
        if isempty(find(ismember(simulation.variables,componentNames(counter4))))
            temp=find(ismember(simulation.states,componentNames(counter4)));
        else
            %         temp=find(ismember(simulation.variables,componentNames(counter)))+k-2; %-1
            temp=find(ismember(simulation.variables,componentNames(counter4)))+k-numel(IQMvariables(model))+numel(plotvariables); %-1
        end
        
        temp=find(ismember([simulation.states(plotstates),simulation.variables(plotvariables)],componentNames(counter4)));
        
        subplot(3,4,temp), hold on, errorbar(time,values(:,counter4),values(:,counter4)-minvalues(:,counter4),maxvalues(:,counter4)-values(:,counter4),'k.','MarkerSize',30), hold off
        hold on
    end
    sgtitle(figureTitles(indExpToOpt(counter)))
end

%% keep only Fig1 and print
h=gcf
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,'-dpdf','-r0','Fig_Exp1.pdf')
