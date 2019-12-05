function [timestamps, data] = nwbGetProtocolTraces(cell_id, protocol_name, rep_num)
% Get the traces data for a given cell, protocol and repetition.
%
% Args:
%     cell_id (int): the ID of the cell.
%     protocol_name (str): the name of the protocol. Should be one of
%         'VRest', 'Activation', 'Ramp', 'Deactivation', 'AP',
%         'Inactivation' or 'Recovery'.
%     rep_num (int): the repetition number, e.g. 1 for repetition #1.
%
% Returns:
%     For protocols with a single dataset per repetition, returns two vectors
%     containing the timestamps and data.
%     For protocols with multiple dataset per repetition, returns two cell
%     arrays containing the timestamps and data for each dataset.
%
% Examples:
%     - Get and plot the data (for normal protocols):
%
%     >> [timestamps, data] = nwbGetProtocolTraces(1234, 'Activation', 1);
%     >> plot(timestamps, data)
%
%     - Get and plot the data (for Recovery, which has multiple datasets
%       per repetition):
%
%     >> [timestamps, data] = nwbGetProtocolTraces(1234, 'Recovery', 1);
%     >> hold on
%     >> cellfun(@plot, timestamps, data)

if rep_num < 1
    error('Repetition number must be postiive.')
end

file_name = ['../Cells/rCell', num2str(cell_id), '.nwb'];
if exist(file_name, 'file') ~= 2
    error('Cell %d not found. Please make sure you placed it in the Cells/ directory.', cell_id)
end

info = h5info(file_name, '/stimulus/presentation/');
group_names = cell(1, length(info.Groups));
for i = 1:length(info.Groups)
    path = split(string(info.Groups(i).Name), '/');
    group_names{i} = path{4};
end

PROTOCOLS = {'VRest', 'Activation', 'Ramp', 'Deactivation', 'AP', 'Inactivation', 'Recovery'};
if ~any(strcmp(group_names, protocol_name))
    err_msg = ['Protocol "', protocol_name, '" not found. Available protocols: '];
    for i = 1:length(PROTOCOLS)
        protocol = PROTOCOLS{i};
        if any(strcmp(group_names, protocol))
            err_msg = [err_msg, protocol, ', '];
        end
    end
    error(err_msg(1:end-2))
end

info = h5info(file_name,  ['/acquisition/timeseries/', protocol_name, '/repetitions/']);
if rep_num > length(info.Groups)
    error('Repetition %d not found. There are %d repetitions in %s.', rep_num, length(info.Groups), protocol_name)
end

protPath  = ['/acquisition/timeseries/', protocol_name, '/', 'repetitions/repetition', num2str(rep_num)];
x_interval = h5read(file_name, [protPath, '/x_interval']);  
x_interval = x_interval(1);
x_start = h5read(file_name, [protPath, '/x_start']);
x_start = x_start(1);
nPoints = h5read(file_name, [protPath, '/n_points']);
hdfPath = [protPath, '/data'];

try % whether it's a non-KvRecovery data format (alternatively use h5info)
    data = h5read(file_name, hdfPath)'; %transposed as in readData!
    maxnPoints = max(nPoints); % to avoid crash when first trace has 0 points
    x_end = x_start+(x_interval * double(maxnPoints-1));
    timestamps = linspace(x_start, x_end, maxnPoints);
    % testing new format, same for all protocols
    [nRow, nCol] = size(data);
    if nCol ~= length(nPoints)
        fprintf(1,'\nWARNING: this cell had a corrupted traces (which was not stored) in %s \n', hdfPath);
    end
    if nRow ~= nPoints(1)
        fprintf(1,'\nERROR: this cell has a corrupted traces which was not discarded!!! %s \n', hdfPath);
    end
catch % if it gives error it's a structured data format
    % cellfun needs a cell array
    timestamps = {};
    data = {};
    for i= 1:length(nPoints)
        dataPath = [hdfPath, '/data', num2str(i), '/data'];
        try % trick to avoid missing (corrupted) traces
            data{end+1} = h5read(file_name, dataPath);
            x_end = x_start+(x_interval * double(nPoints(i)-1));
            timestamps{end+1} = linspace(x_start, x_end, nPoints(i)); % linspace avoid rounding issues of : :
            if size(data{i}) ~= nPoints(i)
                fprintf(1,'\nERROR: this cell has a corrupted traces which was not discarded!!! %s \n', hdfPath);
            end
        catch
            fprintf(1,'\nWARNING: this cell had a corrupted traces (which was not stored) in %s \n', hdfPath);
        end
    end
end
