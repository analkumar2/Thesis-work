function nwbPlotProtocolTraces(cell_id, protocol_name, rep_num)
% Plot a figure for a given cell, protocol and repetition.
% rCell#.mat has to be in ../Cells/
%
% Args:
%     cell_id (int): the ID of the cell.
%     protocol_name (str): the name of the protocol. Should be one of
%         'VRest', 'Activation', 'Ramp', 'Deactivation', 'AP',
%         'Inactivation' or 'Recovery'.
%     rep_num (int): the repetition number, e.g. 1 for repetition #1.
%
% Example:
%     Plot repetition 1 of protocol 'Activation' for cell ID 2:
%
%     >> nwbPlotProtocolTraces(6240, 'Activation', 1);

[timestamps, data] = nwbGetProtocolTraces(cell_id, protocol_name, rep_num);

if strcmp(protocol_name, 'Recovery')
    hold on;
    cellfun(@plot, timestamps, data),
else
    plot(timestamps, data);
end

xlabel('t [ms]');
ylabel('I [nA]');
file_name = ['../Cells/rCell', num2str(cell_id), '.nwb'];
id = num2str(h5read(file_name, '/general/cell_id'));
ion_channel = h5read(file_name, '/general/channel_info/ion_channel');
temp = h5read(file_name, '/general/experiment/temp');
figTitle = sprintf('cell ID %s, %s, Temp. %sC', id, ion_channel{1}, temp{1});
figTitle = strrep(figTitle,'_','\_'); % replace does not work in matlab r2015
title(figTitle);


