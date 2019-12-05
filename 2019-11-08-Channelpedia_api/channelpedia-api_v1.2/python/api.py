#!/usr/bin/env python3
'''Channelpedia Python API.

Allows users to easily plot or extract data from raw or analyzed cell data.

Example:
    >>> import api

    >>> # Plot a protocol/repetition.
    >>> api.plot_protocol(2, 'Activation', 1)
    >>> # Save a protocol/repetition plot to a file (any supported format can
    >>> # be specified in the extension).
    >>> api.save_protocol(2, 'Activation', 1)
    >>> # A custom output path and format can be specified if necessary.
    >>> api.save_protocol(2, 'Activation', 1, '/home/lnmc/custom_name.png')
    >>> # save_protocol(2, 'Activation', 1, 'custom_name.pdf')

    >>> # Get and plot the data (for normal protocols).
    >>> timestamps, data = api.get_protocol_traces(2, 'Activation', 10)
    >>> import matplotlib.pyplot as plt
    >>> plt.plot(timestamps, data)
    >>> plt.show()

    >>> # Get and plot the data (for Recovery, which has multiple datasets
    >>> # per repetition).
    >>> for timestamps, data in api.get_protocol_traces(2, 'Recovery', 1):
    >>>     plt.plot(timestamps, data)  # Plot each trace.
    >>> plt.show()  # Once done, show the resulting plot.
'''

import hashlib
import logging
import os

import h5py
import matplotlib.pyplot as plt
import numpy as np


_PROTOCOLS = ['VRest', 'Activation', 'Ramp', 'Deactivation', 'AP',
              'Inactivation', 'Recovery']
_NWB_FILES_CACHE = {}
logging.basicConfig(level=logging.DEBUG, format='%(message)s')
_LOGGER = logging.getLogger(__name__)


def _reset_ax_params():
    ax = plt.axes()
    ax.margins(0)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')

def _get_nwb_file(cell_id):
    if cell_id in _NWB_FILES_CACHE and _NWB_FILES_CACHE[cell_id]:
        return _NWB_FILES_CACHE[cell_id]
    else:
        file_name = '../Cells/rCell{}.nwb'.format(cell_id)
        if not os.path.isfile(file_name):
            logging.warning(('Cell %d not found. Please make sure you placed '
                             'it in the Cells/ directory.'), cell_id)
            return None
        _NWB_FILES_CACHE[cell_id] = h5py.File(file_name, 'r')
        return _NWB_FILES_CACHE[cell_id]

def _get_timestamps(rep, data, period):
    starting_time = rep['starting_time']
    start = starting_time.value
    stop = start + period * (len(data)-1)
    timestamps = np.linspace(start, stop, len(data), endpoint=False)
    timestamps = np.transpose(timestamps)
    return timestamps

def _get_protocol_traces_generator(data, period):
    ret_val = [];
    for trace in data.values():
        d = trace['data'].value
        timestamps = _get_timestamps(trace, d, period)
        #yield timestamps, d # this would prevent to use nwb_f.close()
        ret_val.append((timestamps, d))
    return ret_val

def _get_repetition(nwb_f, protocol_name, rep_num):
    timeseries = nwb_f['acquisition']['timeseries']
    if protocol_name not in timeseries:
        _LOGGER.error('Protocol "%s" not found. Available protocols:',
                      protocol_name)
        for protocol in _PROTOCOLS:
            if protocol in timeseries:
                num_rep = len(timeseries[protocol]['repetitions'])
                _LOGGER.error('- %s (%d repetition%s)', protocol, num_rep,
                              '' if num_rep == 1 else 's')
        return None
    protocol = timeseries[protocol_name]
    rep_name = 'repetition{}'.format(rep_num)
    if rep_name not in protocol['repetitions']:
        num_rep = len(protocol['repetitions'])
        _LOGGER.error(('Repetition %d not found. There %s %d repetition%s in '
                       '%s.'),
                      num_rep, 'is' if num_rep == 1 else 'are', num_rep,
                      '' if num_rep == 1 else 's', protocol_name)
        return None
    rep = protocol['repetitions'][rep_name]
    return rep

def _get_temp_color(temp):
    if temp == '15c':
        return (0, 0, 1), (1, 0, 0)
    elif temp == '25c':
        return (0, 0, 0), (1, 0, 0)
    elif temp == '35c':
        return (1, 0, 0), (0, 0, 1)
    elif temp == 'RT':
        return (0.1647, 0.3843, 0.2745), (0, 0, 1)
    else:
        return (0.5843, 0.3882, 0.3882), (0, 0, 1)

def _prepare_plot(cell_id, protocol_name, rep_num):
    nwb_f = _get_nwb_file(cell_id)
    if not nwb_f: return False
    rep = _get_repetition(nwb_f, protocol_name, rep_num)
    if not rep:
        nwb_f.close()
        return False
    data = rep['data']
    temp = str(nwb_f['general']['experiment']['temp'].value, 'utf-8')
    temp_color, _ = _get_temp_color(temp)
    period = rep['x_interval'][0]
    if isinstance(data, h5py._hl.group.Group):
        for trace in data.values():
            d = trace['data']
            starting_time = trace['starting_time']
            start = starting_time.value
            stop = start + period * (len(d)-1)
            timestamps = np.linspace(start, stop, len(d), endpoint=False)
            timestamps = np.transpose(timestamps)
            plt.plot(timestamps, d, color=temp_color)
    else:
        starting_time = rep['starting_time']
        start = starting_time.value
        stop = start + period * (len(data)-1)
        timestamps = np.linspace(start, stop, len(data), endpoint=False)
        timestamps = np.transpose(timestamps)
        plt.plot(timestamps, data, color=temp_color)
    plt.xlabel('Time (ms)', fontsize=16)
    plt.ylabel('I', fontsize=16)
    id_ = nwb_f['general']['cell_id'].value
    ion_channel = str(nwb_f['general']['channel_info']['ion_channel'].value,
                      'utf-8')
    data_name = str(nwb_f['identifier'].value, 'utf-8')
    title = 'ID : {}: {} ({}) Temp = {}'.format(
        id_, ion_channel, data_name, temp)
    plt.title(title, fontsize=14, fontweight='bold')
    nwb_f.close()
    return True

def get_protocol_traces(cell_id, protocol_name, rep_num):
    '''Get the traces data for a given cell, protocol and repetition.

    Args:
        cell_id (int): the ID of the cell.
        protocol_name (str): the name of the protocol. Should be one of
            'VRest', 'Activation', 'Ramp', 'Deactivation', 'AP',
            'Inactivation' or 'Recovery'.
        rep_num (int): the repetition number, e.g. 1 for repetition #1.

    Returns:
        For protocols with a single dataset per repetition, returns a tuple
        containing the timestamps and data.
        For protocols with multiple dataset per repetition, returns a generator
        yielding tuples containing the timestamps and data.

    Examples:
        - Get and plot the data (for normal protocols):

        >>> import matplotlib.pyplot as plt
        >>> timestamps, data = get_protocol_traces(2, 'Activation', 1)
        >>> plt.plot(timestamps, data)
        >>> plt.show()

        - Get and plot the data (for Recovery, which has multiple datasets
          per repetition):

        >>> import matplotlib.pyplot as plt
        >>> for timestamps, data in get_protocol_traces(2, 'Recovery', 1):
        >>>     plt.plot(timestamps, data)  # Plot each trace.
        >>> plt.show()  # Once done, show the resulting plot.
    '''
    nwb_f = _get_nwb_file(cell_id)
    if not nwb_f: return None
    rep = _get_repetition(nwb_f, protocol_name, rep_num)
    if not rep:
        return None
    data = rep['data']
    period = rep['x_interval'][0]
    if isinstance(data, h5py._hl.group.Group):
        ret_val = _get_protocol_traces_generator(data, period)
    else:
        data = data.value
        timestamps = _get_timestamps(rep, data, period)
        ret_val = (timestamps, data)
    # the following crashes with Recovery since _get_protocol_traces_generator is a generator
    nwb_f.close()
    return ret_val

def plot_protocol(cell_id, protocol_name, rep_num):
    '''Plot a figure for a given cell, protocol and repetition.

    Args:
        cell_id (int): the ID of the cell.
        protocol_name (str): the name of the protocol. Should be one of
            'VRest', 'Activation', 'Ramp', 'Deactivation', 'AP',
            'Inactivation' or 'Recovery'.
        rep_num (int): the repetition number, e.g. 1 for repetition #1.

    Example:
        Plot repetition 1 of protocol 'Activation' for cell ID 2:

        >>> plot_protocol(2, 'Activation', 1)
    '''
    if _prepare_plot(cell_id, protocol_name, rep_num):
        _reset_ax_params()
        plt.show()
        plt.gcf().clear()
        _reset_ax_params()

def save_protocol(cell_id, protocol_name, rep_num, output_path=None, dpi=600):
    '''Save a figure of a given cell, protocol and repetition to disk.

    Args:
        cell_id (int): the ID of the cell.
        protocol_name (str): the name of the protocol. Should be one of
            'VRest', 'Activation', 'Ramp', 'Deactivation', 'AP',
            'Inactivation' or 'Recovery'.
        rep_num (int): the repetition number, e.g. 1 for repetition #1.
        output_path (str, optional): the path the figure will be saved to. The
            output format is deduced from the extension of the filename. If not
            specified, the figure is saved as
            'rCell{cell_id}_{protocol_name}_rep{rep_num}.png' in the current
            directory (e.g. 'rCell2_Activation_rep2.png').
        dpi (int, optional): the resolution in dots per inch. Defaults to 600.

    Examples:
        - Save plot of repetition 1 of protocol 'Activation' for cell ID 2:

        >>> save_protocol(2, 'Activation', 1)

        - Save plot of repetition 1 of protocol 'Activation' for cell ID 2 as
          'custom_name.png' with a resolution of 1200:

        >>> save_protocol(2, 'Activation', 1, 'custom_name.png', 1200)

        - Save plot of repetition 1 of protocol 'Activation' for cell ID 2 as
          'custom_name.pdf' (using the PDF output format):

        >>> save_protocol(2, 'Activation', 1, 'custom_name.pdf')
    '''
    if _prepare_plot(cell_id, protocol_name, rep_num):
        if not output_path:
            output_path = 'rCell{}_{}_rep{}.png'.format(
                cell_id, protocol_name, rep_num)
        _reset_ax_params()
        plt.savefig(output_path, dpi=600)
        plt.gcf().clear()
        _reset_ax_params()


plt.rc('path', simplify=False)
plt.rc('figure', facecolor='white')
plt.rcParams['axes.spines.top'] = False
plt.rcParams['axes.spines.right'] = False
plt.rcParams['xtick.major.size'] = 10
plt.rcParams['ytick.major.size'] = 10
plt.rcParams['lines.linewidth'] = 0.5
_reset_ax_params()

if __name__ == '__main__':
    print(('This script cannot be run directly. Please import it from your '
           'Python interpreter.'))
