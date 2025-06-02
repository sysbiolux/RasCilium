%% Plotting of final fits with optimized parameters
clc; clear all; close all

% simple model
% load('C:\Users\thomas.sauter\OneDrive - University of Luxembourg\work_other\Projects\Daniel_2024\model_230924_simpleModel\WORKSPACE260924_simpleModel.mat')
% threshold=0.041

% full model: estimate fAsym
% load('WORKSPACE_exp12_fAsym_100.mat')
% threshold=min(optimalcost)*4/3 %threshold=0.06

% full model: exp3-6
load('WORKSPACE_exp3456_fAsym_ICs_100.mat')
threshold=0.05

simulationTime=[0:0.05:5];

%% Plot selected optimal solutions
% Set desirable threshold
% threshold=0.10 %inf %0.05

% Define states to plot
model_structure = IQMgetmodel(RASdiff_opt);
experiment = IQMgetexperiment(RASdiff_opt,indExpToOpt(counter2)); % Select measurement (only 1)
model = IQMmergemodexp(model_structure, experiment);
plotstates=[1:numel(IQMstates(model))]; % plot all 3 states
plotstates=[]
% plotvariables=[2,4:numel(IQMvariables(model))]; % plot all variables except serumSwitchTime + sumPTD
[plotVariableNames,plotvariables]=setdiff(IQMvariables(model),{'serumSwitch','rasDiff','sumPTD'},'stable');
% IQMvariables(model)r
[plotVariableNames,plotvariables]=intersect(IQMvariables(model),{'pPROG','pTRANS','pDIFF','pCil'},'stable');

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
    
% % %     figure(counter2)
    
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
% % %             for k=1:length(plotstates);
% % %                 plotstate=plotstates(k);
% % %                 subplot(3,4,k), hold on, plot(simulation.time,simulation.statevalues(:,plotstate),'LineWidth',2),
% % %                 legend(simulation.states(plotstate)), hold off
% % %                 hold on
% % %             end
% % %             for k2=1:length(plotvariables);
% % %                 plotstate=plotvariables(k2);
% % %                 subplot(3,4,k+k2), hold on, plot(simulation.time,simulation.variablevalues(:,plotstate),'LineWidth',2),
% % %                 legend(simulation.variables(plotstate)), hold off
% % %                 hold on
% % %             end
            
        elseif (optimalcost(j)>=threshold) % if optimalcost didn't pass threshold
            disp(['Not included'])
        end
        
    end
    
    % Add the marker and error from experimental data to compare with the simulation
    
    % Extract measurement data from datasheet inside the project
% % %     [time,componentNames,values,minvalues,maxvalues]=IQMmeasurementdata(IQMgetmeasurement(RASdiff_project,indExpToOpt(counter2)));
    
    % Plot measurement data with range of S.D. as errorbars on the plot
% % %     for counter=1:numel(componentNames)
% % %         if isempty(find(ismember(simulation.variables,componentNames(counter))))
% % %             temp=find(ismember(simulation.states,componentNames(counter)));
% % %         else
% % %             %         temp=find(ismember(simulation.variables,componentNames(counter)))+k-2; %-1
% % %             temp=find(ismember(simulation.variables,componentNames(counter)))+k-numel(IQMvariables(model))+numel(plotvariables); %-1
% % %         end
% % %         subplot(3,4,temp), hold on, errorbar(time,values(:,counter),values(:,counter)-minvalues(:,counter),maxvalues(:,counter)-values(:,counter),'r.','MarkerSize',20), hold off
% % %         hold on
% % %         
% % %         sgtitle(figureTitles(indExpToOpt(counter2)))
% % %     end
    
end

size(plotDataAll)
save WORKSPACE

%% median + 68% confidence interval
% close all

fontSize=14
lineWidth=1
printAs='-djpeg'  %'-depsc' %'-djpeg'
titleOn=1

% figure
figure('color','white')
for counter=1:numel(indExpToOpt)
    counter
    plotData=plotDataAll(:,:,:,counter);
    mPlotData=median(plotData,3);
    
% % %     collectData1=[0 0 0 simulationTime];
% % %     collectData2=[0 0 0 simulationTime];
% % %     collectData3=[0 0 0 simulationTime];
% % %     collectData4=[0 0 0 simulationTime];
%     figure
    for counter6=[10:1:68] %conf levels
        %             fadingColor=[1 (counter6+20)/100 (counter6+20)/100];
        % % %         fadingColor=[1 (counter6-10)/58*0.8+0.1 (counter6-10)/58*0.8+0.1]; %red %for range [10, 68]
        if counter==1
            fadingColor=[(counter6-10)/58*0.8+0.1 (counter6-10)/58*0.8+0.1 (counter6-10)/58*0.8+0.1]*0.9; %red %for range [10, 68]
        else
            fadingColor=[(counter6-10)/58*0.8+0.1 (counter6-10)/58*0.8+0.1 1]; %red %for range [10, 68]
        end
        lbPlotData=prctile(plotData,50-counter6/2,3); %68% range of fits / predictions (similar standard deviation)
        ubPlotData=prctile(plotData,50+counter6/2,3);
% % %         for counter2=1:length(plotstates)
% % %             plotState=plotstates(counter2);
% % %             subplot(3,4,counter2), hold on, plot(simulationTime,lbPlotData(:,counter2),'-','LineWidth',lineWidth,'Color',fadingColor), hold off
% % %             subplot(3,4,counter2), hold on, plot(simulationTime,ubPlotData(:,counter2),'-','LineWidth',lineWidth,'Color',fadingColor), hold off
% % %         end
        for counter3=1:length(plotvariables)
            plotState=plotvariables(counter3);
            subplot(4,4,(counter-1)*4+counter3), hold on, plot(simulationTime,lbPlotData(:,counter3),'-','LineWidth',lineWidth,'Color',fadingColor), hold off
            subplot(4,4,(counter-1)*4+counter3), hold on, plot(simulationTime,ubPlotData(:,counter3),'-','LineWidth',lineWidth,'Color',fadingColor), hold off
% % %             if counter3==1
% % %                collectData1=[collectData1; fadingColor, lbPlotData(:,counter3)'];
% % %                collectData1=[collectData1; fadingColor, ubPlotData(:,counter3)'];
% % %             end
% % %             if counter3==2
% % %                collectData2=[collectData2; fadingColor, lbPlotData(:,counter3)'];
% % %                collectData2=[collectData2; fadingColor, ubPlotData(:,counter3)'];
% % %             end
% % %             if counter3==3
% % %                collectData3=[collectData3; fadingColor, lbPlotData(:,counter3)'];
% % %                collectData3=[collectData3; fadingColor, ubPlotData(:,counter3)'];
% % %             end
% % %             if counter3==4
% % %                collectData4=[collectData4; fadingColor, lbPlotData(:,counter3)'];
% % %                collectData4=[collectData4; fadingColor, ubPlotData(:,counter3)'];
% % %             end
        end
    end
% % %     for counter2=1:length(plotstates)
% % %         plotState=plotstates(counter2);
% % %         subplot(3,4,counter2), hold on, plot(simulationTime,mPlotData(:,plotState),'r--','LineWidth',lineWidth*3)
% % %         if titleOn, title(simulation.states(plotState)), end
% % %         set(gca,'FontSize',fontSize,'XTickLabel',{'-2','0','2','4','6','8'}), hold off
% % %     end
for counter3=1:length(plotvariables)
    plotState=plotvariables(counter3);
    % % %         subplot(4,4,(counter-1)*4+counter3), hold on, plot(simulationTime,mPlotData(:,counter3),'r--','LineWidth',lineWidth*3)
    if counter==1
        subplot(4,4,(counter-1)*4+counter3), hold on, plot(simulationTime,mPlotData(:,counter3),'k--','LineWidth',lineWidth*3)
    else
        subplot(4,4,(counter-1)*4+counter3), hold on, plot(simulationTime,mPlotData(:,counter3),'b--','LineWidth',lineWidth*3)
    end
    if titleOn, title(simulation.variables(plotState)), end
    %         set(gca,'FontSize',fontSize,'XTickLabel',{'0','2','4','6'}), hold off
    set(gca,'FontSize',fontSize), hold off
    xlim([0 5])
    switch counter3
        case 1
            ylim([0 2])
                switch counter
                    case 1
                        ylabel('Ctrl','FontSize',fontSize+6,'fontweight','bold')
                    case 2
                        ylabel('Kras','FontSize',fontSize+6,'fontweight','bold')
                    case 3
                        ylabel('Nras','FontSize',fontSize+6,'fontweight','bold')
                    case 4
                        ylabel('Hras','FontSize',fontSize+6,'fontweight','bold')
                end
            case 2
                ylim([75 100])
            case 3
                ylim([0 20])
            case 4
                ylim([0 90])
                temp=mPlotData(:,counter3);
                disp(' ')
                disp(['For experiment: ' num2str(indExpToOpt(counter))])
                disp('pCil on Day 3: ')
                disp(temp(find(simulationTime==3)))
        end
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
        
        subplot(4,4,(counter-1)*4+temp), hold on, errorbar(time,values(:,counter4),values(:,counter4)-minvalues(:,counter4),maxvalues(:,counter4)-values(:,counter4),'k.','MarkerSize',30), hold off
        hold on
    end
    sgtitle('Diff with serum switch (Fig.3): Experiments 3 - 4 - 5 - 6')
end

h=gcf
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% print(h,'Fig_Exp3456.pdf','-dpdf','-r0')

%% svg export ??
exportgraphics(gcf, 'Fig_Exp3456.pdf', 'ContentType', 'Vector')
% exportgraphics(gcf, 'XXX.eps', 'ContentType', 'Vector')
% exportgraphics(gcf, 'XXX.tiff', 'ContentType', 'Vector')

%% legend
figure
plot(0.7,1,'k.','MarkerSize',30)
text(1.1,1.01,'Measurement','FontSize',20)
hold on
plot([0.4 1],[0.7 0.7],'r--','LineWidth',lineWidth*3)
text(1.1,0.71,'Model (median +/- 68%)','FontSize',20)
axis([0 5 0 2])
axis off
hold off
