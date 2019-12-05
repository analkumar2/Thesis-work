function [] = plot_features(cell_id, protocol_name, rep_num)
% Plot analysis feature for a given cell, protocol and repetition.
% rCell#.mat has to be in ../Cells/
% aCell#.nwb has to be in ../Analysis/
%
% Args:
%     cell_id (int): the ID of the cell.
%     protocol_name (str): the name of the protocol. Should be one of
%         'Activation', 'Ramp', 'Deactivation', 'AP', 'Inactivation'.
%     rep_num (int): the repetition number, e.g. 1 for repetition #1.
%
% Example:
%     Plot analysis feature for repetition 1 of protocol 'Activation' for cell ID 2:
%
%     >> plotFeatures(6240, 'Activation', 1, 1);

dwnSample = 1;

protocols = {'Activation', 'Ramp', 'Deactivation', 'AP', 'Inactivation'};
if ~ismember(protocol_name, protocols)
    fprintf(1,'\nProtocol unknown or not supported.\n')
    return
end

file_name = ['../Analysis/aCell', num2str(cell_id), '.mat'];
x = load(file_name);
aCell = x.aCell;
if isempty(aCell)
    fprintf(1, '\nError aCell %d not found [%s]', cell_id, file_name);
    return
end
if ~isfield(aCell.(protocol_name),'repetition')
    fprintf(1,'\nWARNING: no protocol %s in aCell %d ',protocol_name, cell_id);
    return
end

[time, data] = nwbGetProtocolTraces(cell_id, protocol_name, rep_num);
if isempty(data)
    fprintf(1, '\nError timeseries %d not found [%s]', cell_id, file_name);
    return
end
% convert to milliseconds
if time(1) ~= 1
  time = time * 1000;
end

hfig = figure();    
[clr, fitClr] = getTempColor(aCell.Temp);
[nPtTotal, nTrace] = size(data);
if isfield(aCell.(protocol_name).repetition(rep_num), 'mean_base')
    meanBase = aCell.(protocol_name).repetition(rep_num).mean_base;
    for i=1:nTrace
        data(:,i) = data(:,i) - meanBase(i);
    end
    data = data ./ aCell.(protocol_name).repetition(rep_num).max_I;
end
time1   = downsample(time, dwnSample);
Data1   = downsample(data, dwnSample); 
plot(time1, Data1, 'color', clr);
hold on;
switch protocol_name
    case 'Activation'
        plotPhaseFeatures(aCell, time, protocol_name, rep_num, 4);
    case 'Inactivation'
        plotPhaseFeatures(aCell, time, protocol_name, rep_num, 4);
        plotPhaseFeatures(aCell, time, protocol_name, rep_num, 5);
    case 'Deactivation'
        plotPhaseFeatures(aCell, time, protocol_name, rep_num, 5);
    case 'Ramp'
        plotRampIFeatures(aCell, rep_num, time, data);
    case 'AP'
        plotAPFeatures(aCell, rep_num, time, data);
end
ha = gca();
set(ha, 'Box', 'off',  'Units','normalized', 'clipping' , 'off', 'XMinorTick', 'off', 'YMinorTick', 'off');
ylabel('I', 'FontSize', 14, 'FontName', 'Microsoft Sans Serif');
xlabel('Time(ms)', 'FontSize', 14, 'FontName', 'Microsoft Sans Serif');
strTitle = ['cell ID ', num2str(aCell.cellDBID), ', ', aCell.ionChannel, ', Temp. ', aCell.Temp];
title(strTitle,'interpreter', 'none', 'FontSize', 12, 'FontWeight', 'bold');
set(hfig, 'color', [1 1 1]);
set(hfig, 'Position', [150, 300, 800, 600]);
axis('tight')

if exist('ha')
    ha.XTickLabel = str2double(ha.XTickLabel)/1000;
end



function [] = plotAPFeatures(aCell, rep_num, time, Data)
    [clr, fitClr] = getTempColor(aCell.Temp);
    rep   = aCell.AP.repetition(rep_num);
    for i = 1:length(rep.NI.peakTime)
        plot(rep.NI.peakTime{i}, rep.NI.peakVal{i}, '*','color', fitClr);
        plot(rep.phaseTime{i}, rep.phaseVal{i}, '*','color', fitClr);
    end


function [] = plotRampIFeatures(aCell, rep_num, time, Data)
    [clr, fitClr] = getTempColor(aCell.Temp);
    phaseIDs = [4, 5, 7, 8, 10, 11, 13, 14];
%      for i = 1:length(phaseIDs)
%          p  =  aCell.Ramp.repetition(rep_num).phase(phaseIDs(i));
%          NI =  p.NI;
%          plot(time(p.idStart+NI.minIDs), Data(p.idStart+NI.minIDs, :), '.k');
%          plot(time(p.idStart+NI.maxIDs), Data(p.idStart+NI.maxIDs, :), '.k');
%      end
    [nPt, nTrace] = size(Data);
    clr1 = [0.75, 0, 0.75];
    for i = 1:length(phaseIDs)
        p  =  aCell.Ramp.repetition(rep_num).phase(phaseIDs(i));
        NI =  p.NI;
        for j = 1:nTrace
            x = time(p.idStart:p.idEnd); y = Data(p.idStart:p.idEnd, j);
            figure(2); hold on; 
            plot(x, y, 'color', clr);
            x1 = x(1); x2 = x(end);
            y1 = mean(y(1:50));  y2 = mean(y(end-50:end));
            m  = (y2-y1)/(x2-x1);
            line = (y1 + (m*(x'-x1)));
            newData = y - line;
            plot(x, line, 'color', clr1);
            figure(3); hold on;
            plot(x, newData, 'color', clr);
            plot(x, zeros(size(x)), 'color', clr1);
            plot(x(NI.minIDs(j)), newData(NI.minIDs(j)), '*', 'color', fitClr );
            if(~isempty(find(phaseIDs(i)==[4, 7, 10, 13], 1)))
                plot(x(NI.maxIDs(j)), newData(NI.maxIDs(j)), '*', 'color', fitClr );
            end
        end
        %plot(time(p.idStart+NI.maxIDs), Data(NI.maxIDs, :), '*', 'color', fitClr);
        %plot(time(p.idStart+NI.minIDs), Data(NI.minIDs, :), 'g*');
    end
    
    
function [] = plotPhaseFeatures(aCell, time, protocol, rep_num,  phaseID)
        [clr, fitClr] = getTempColor(aCell.Temp);
        NI = aCell.(protocol).repetition(rep_num).phase(phaseID).NI;
        tStart  = time(aCell.(protocol).repetition(rep_num).phase(phaseID).idStart);
        idStart = aCell.(protocol).repetition(rep_num).phase(phaseID).idStart;
        idEnd   = aCell.(protocol).repetition(rep_num).phase(phaseID).idEnd;
        for i = 1: length(NI.Volt)
            x = 0:0.1:NI.peak_time(i);
            if(isfield(NI, 'Act_Tau'))
                [x, y, strEq] = singleExp(x, NI.peak_value(i), NI.Act_Tau(i), 4);  %[x, y, strEq] = singleExp(x, NI.Act_A1(i), NI.Act_Tau1(i), 4);
            else
                [x, y, strEq] = singleExp(x, NI.peak_value(i), NI.Act_Tau1(i), 4);
            end
            ids = find(y>1.2); if(ids) y(ids)=1.2; end
            plot(tStart+NI.peak_time(i), NI.peak_value(i), '*', 'MarkerSize', 8, 'color', [1, 0, 1]);
            plot(tStart+x, y, 'g--');
            %plot(tStart+time(NI.id5pc), NI.peak_value_5pc(i),   '+', 'color', fitClr);
            %plot(tStart+time(NI.id25pc), NI.peak_value_25pc(i), '+', 'color', fitClr);
            %plot(tStart+time(NI.id50pc), NI.peak_value_50pc(i), '+', 'color', fitClr);
            %plot(tStart+time(NI.id75pc), NI.peak_value_75pc(i), '+', 'color', fitClr);
            %plot(tStart+time(NI.id95pc), NI.peak_value_95pc(i), '+', 'color', fitClr);
            x = NI.peak_time(i)+tStart:0.1:time(idEnd); 
            x = x -x(1); 
            %%[x, y, strEq] = doubleExp(x, NI.Inact_2A1(i), NI.Inact_2Tau1(i), NI.Inact_2A2(i), NI.Inact_2Tau2(i), 1);
            %[x, y, strEq] = doubleExp(x, NI.Inact_A1(i), NI.Inact_Tau1(i), NI.Inact_A2(i), NI.Inact_Tau2(i), 1);
            [x, y, strEq] = offsetSingleExp(x, NI.Inact_C(i), NI.Inact_A(i), NI.Inact_Tau(i), 1);
            ids = find(y<-0.1); if(ids) y(ids)=-0.1; end
            plot(tStart+NI.peak_time(i)+x, y, '--', 'color', fitClr);
        end
    
