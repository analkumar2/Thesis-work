*********************
Channelpedia API Help
*********************

Version 1.1: compatible with Channelpedia data file nwb_version 1.1

Table of Contents
=================

-  `General Instructions`_
-  `MATLAB API`_
-  `Python API`_

-------------------------------------------------------------------------------

General Instructions
====================

-  Download and unzip the Channelpedia API.
-  Place raw data downloaded from Channelpedia (e.g. ``rCell1234.nwb``)
   in the ``Cells/`` directory,
   and any downloaded analysis files (e.g. ``aCell1234.nwb.gz``) in
   the ``Analysis/`` directory.
-  You can find the Python API scripts in the ``python/`` directory, and the
   MATLAB API scripts in the ``matlab/`` directory.

-------------------------------------------------------------------------------

MATLAB API
==========

Pre-requisites
--------------

-  You need to have a valid MATLAB installation. `Please refer to the Mathworks
   support page for help <https://www.mathworks.com/help/install/>`__.

Usage
-----

-  Start MATLAB.
-  In MATLAB, make sure you're in the directory where the Channelpedia MATLAB
   API script files are (``matlab/``). You can do so using the ``cd`` command,
   for example::

       >> cd 'C:\Users\Me\Download\channelpedia-api\matlab\'

-  You can now use the provided functions to easily access or plot trace data.

Available API functions
-----------------------

``nwbGetProtocolTraces(cell_id, protocol_name, rep_num)``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Get the traces data for a given cell, protocol and repetition.

Parameters
""""""""""

- **cell_id** (*int*): the ID of the cell.
- **protocol_name** (*str*): the name of the protocol. Should be one of
  ‘VRest’, ‘Activation’, ‘Ramp’, ‘Deactivation’, ‘AP’, ‘Inactivation’
  or ‘Recovery’.
- **rep_num** (*int*): the repetition number, e.g. 1 for repetition #1.                                                                                                                                                                     |

Returns
"""""""

For protocols with a single dataset per repetition, returns two vectors
containing the timestamps and data.
For protocols with multiple dataset per repetition, returns two cell arrays
containing the timestamps and data for each dataset.

Examples
""""""""

- Get and plot the data (for normal protocols)::

    >> [timestamps, data] = nwbGetProtocolTraces(1234, 'Activation', 1);
    >> plot(timestamps, data)

- Get and plot the data (for Recovery, which has multiple datasets per
  repetition)::

    >> [timestamps, data] = nwbGetProtocolTraces(1234, 'Recovery', 1);
    >> hold on
    >> cellfun(@plot, timestamps, data)

``nwbPlotProtocolTraces(cell_id, protocol_name, rep_num)``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Plot the traces data for a given cell, protocol and repetition.

Parameters
""""""""""

- **cell_id** (*int*): the ID of the cell.
- **protocol_name** (*str*): the name of the protocol. Should be one of
  ‘VRest’, ‘Activation’, ‘Ramp’, ‘Deactivation’, ‘AP’, ‘Inactivation’
  or ‘Recovery’.
- **rep_num** (*int*): the repetition number, e.g. 1 for repetition #1.                                                                                                                                                                     |

Example
"""""""

Plot repetition 1 of protocol ‘Activation’ for cell ID 1234::

    >> nwbPlotProtocolTraces(1234, 'Activation', 1)

``plot_features(cell_id, protocol_name, rep_num)``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Plot analysis feature for a given cell, protocol and repetition.

Parameters
""""""""""

- **cell_id** (*int*): the ID of the cell.
- **protocol_name** (*str*): the name of the protocol. Should be one of
  ‘VRest’, ‘KvActivation’, ‘KvRamp’, ‘KvDeactivation’, ‘AP’, ‘KvInactivation’
  or ‘KvRecovery’.
- **rep_num** (*int*): the repetition number, e.g. 1 for repetition #1.                                                                                                                                                                     |

Example
"""""""

Plot analysis feature for repetition 1 of protocol ‘KvActivation’ for cell ID
1234::

    >>> plot_features(1234, 'KvActivation', 1)

-------------------------------------------------------------------------------

Python API
==========

Pre-requisites
--------------

-  Python 3.5.x (or later) is required (tested with 3.6). If you have Python 2 installed
   (check with running ``python -V`` in a terminal) your operative system
   likely needs it, and you should not upgrade it: it is recommended instead
   to install Anaconda with Python 3, https://www.anaconda.com/distribution/,
   through which you can open a special terminal with Python 3, when needed.
-  The Channelpedia Python API requires a few Python modules to be installed;
   the easiest way to install them is using pip.
   If you've installed Python from the official website, you probably already
   have pip installed. To verify, run ``pip -V`` in a terminal. If not,
   to install pip run the following commands in a terminal:

       $ curl https://bootstrap.pypa.io/get-pip.py -o get-pip.py
       $ python get-pip.py --user

   **Note:** on certain Linux distributions, the pip version found in the
   distro's repositories may be outdated. In case of any errors when installing
   dependencies, it is recommended to install pip from the official source as
   shown above.
-  Once pip is installed, you can have it fetch and install the required
   dependencies for the Channelpedia Python API.

   -  Open a terminal inside the ``python/`` directory.
   -  Run the following command to have pip fetch and install the required
      dependencies::

          $ pip install -r requirements.txt --user

      **Note:** on some systems (e.g. Ubuntu), you may need to use ``pip3``
      instead of ``pip``, as the latter may default to the Python 2.x version.

-  (Optional) If you want to analyze raw data cells yourself, you'll also need
   to install the scipy Python module. For Linux/macOS users, simply run
   ``pip install scipy --user`` in a terminal to do so.
   For Windows users, download and install the following pre-built packages:

   -  **Note:** download the appropriate package versions depending on the
      Python version you have installed ("cp35" if you have Python 3.5, "cp36"
      if you have Python 3.6) and architecture ("win32" for 32-bit Python,
      "win\_amd64" for 64-bit Python).
      To find out your Python version, run the following command in a
      terminal::

          $ python -V

     To verify the architecture version, run the following command in a
     terminal::

          $ python -c "import platform; print(platform.architecture())"

   -  Numpy+MKL: `download the appropriate
      package <http://www.lfd.uci.edu/~gohlke/pythonlibs/#numpy>`__ then
      install it using pip. You can do so by opening a terminal in the
      directory where you downloaded the ``.whl`` file, then running the
      following command (modify the file name to point to the version you
      downloaded)::

          $ pip install "numpy‑1.13.1+mkl‑cp35‑cp35m‑win32.whl" --user

   -  SciPy: `download the appropriate
      package <http://www.lfd.uci.edu/~gohlke/pythonlibs/#scipy>`__ then
      install it using pip as above, for example::

          $ pip install "scipy‑0.19.1‑cp35‑cp35m‑win32.whl" --user

Usage
-----

-  Open a terminal in the ``python/`` directory and start the Python
   interpreter::

       $ python

-  Import the Channelpedia Python API::

       >>> import api

-  You can now use all the functions the API provides. For example, to show the
   plot for the first repetition of the Activation protocol for cell 1234,
   run the following command::

       >>> api.plot_protocol(1234, 'Activation', 1)

Available API functions
-----------------------

**Note:** You can view this documentation at any time from the Python
interpreter::

    >>> import api
    >>> help(api)


``api.get_protocol_traces(cell_id, protocol_name, rep_num)``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Get the traces data for a given cell, protocol and repetition.

Parameters
""""""""""

- **cell_id** (*int*): the ID of the cell.
- **protocol_name** (*str*): the name of the protocol. Should be one of
  ‘VRest’, ‘Activation’, ‘Ramp’, ‘Deactivation’, ‘AP’, ‘Inactivation’
  or ‘Recovery’.
- **rep_num** (*int*): the repetition number, e.g. 1 for repetition #1.                                                                                                                                                                     |

Returns
"""""""

For protocols with a single dataset per repetition, returns a tuple containing
the timestamps and data.
For protocols with multiple dataset per repetition, returns a generator
yielding tuples containing the timestamps and data.

Examples
""""""""

- Get and plot the data (for normal protocols)::

    >>> import matplotlib.pyplot as plt
    >>> timestamps, data = api.get_protocol_traces(1234, 'Activation', 1)
    >>> plt.plot(timestamps, data)
    >>> plt.show()

- Get and plot the data (for Recovery, which has multiple datasets per
  repetition)::

    >>> import matplotlib.pyplot as plt
    >>> for timestamps, data in api.get_protocol_traces(1234, 'Recovery', 1):
    >>>     plt.plot(timestamps, data)  # Plot each trace.
    >>> plt.show()  # Once done, show the resulting plot.

``api.plot_protocol(cell_id, protocol_name, rep_num)``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Plot a figure for a given cell, protocol and repetition.

Parameters
""""""""""

- **cell_id** (*int*): the ID of the cell.
- **protocol_name** (*str*): the name of the protocol. Should be one of
  ‘VRest’, ‘Activation’, ‘Ramp’, ‘Deactivation’, ‘AP’, ‘Inactivation’
  or ‘Recovery’.
- **rep_num** (*int*): the repetition number, e.g. 1 for repetition #1.                                                                                                                                                                     |

Example
"""""""

Plot repetition 1 of protocol ‘Activation’ for cell ID 1234::

    >>> api.plot_protocol(1234, 'Activation', 1)

``api.save_protocol(cell_id, protocol_name, rep_num, output_path=None, dpi=600)``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Save a figure of a given cell, protocol and repetition to disk.

Parameters
""""""""""

- **cell_id** (*int*): the ID of the cell.
- **protocol_name** (*str*): the name of the protocol. Should be one of
  ‘VRest’, ‘Activation’, ‘Ramp’, ‘Deactivation’, ‘AP’, ‘Inactivation’
  or ‘Recovery’.
- **rep_num** (*int*): the repetition number, e.g. 1 for repetition #1.                                                                                                                                                                     |
- **output_path** (*str*, optional): the path the figure will be saved to. The
  output format is deduced from the extension of the filename. If not
  specified, the figure is saved as
  ‘rCell{cell_id}_{protocol_name}_rep{rep_num}.png’ in the current directory
  (e.g. ‘rCell1234_Activation_rep1.png’).
- **dpi** (*int*, optional): the resolution in dots per inch. Defaults to 600.

Examples
""""""""

- Save plot of repetition 1 of protocol ‘Activation’ for cell ID 1234::

    >>> api.save_protocol(1234, 'Activation', 1)

- Save plot of repetition 1 of protocol ‘Activation’ for cell ID 1234 as
  ‘custom_name.png’ with a resolution of 1200::

    >>> api.save_protocol(1234, 'Activation', 1, 'custom_name.png', 1200)

- Save plot of repetition 1 of protocol ‘Activation’ for cell ID 1234 as
  ‘custom_name.pdf’ (using the PDF output format)::

    >>> api.save_protocol(1234, 'Activation', 1, 'custom_name.pdf')
