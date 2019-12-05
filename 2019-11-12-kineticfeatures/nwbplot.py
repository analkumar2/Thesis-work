import h5py
import matplotlib.pyplot as plt
reader = h5py.File('../../Raw_data/Channelpedia/Kv1.1/DataKv1.1RatCHO/rCell10070.nwb')
reader['acquisition']['timeseries'].keys() #Find what all protocols are available
data = reader['acquisition']['timeseries']['Activation']['repetitions']['repetition1']['data']
plt.plot(data[:,4:15])
plt.show()

reader.keys()
# <KeysViewHDF5 ['acquisition', 'analysis', 'checksums', 'epochs', 'file_create_date', 'general', 'identifier', 'nwb_version', 'processing', 'session_description', 'session_start_time', 'stimulus']>
reader['acquisition']['timeseries'].keys() # give activation or AP or inactivation etc
reader['general']['cell_id'][:] # gives cell_id
reader['general']['cell_info'] #More info on the cell
reader['general']['channel_info'] #More info on the ion channel expressed
# The general key has more metadata

reader['session_description'][:] #Description of the data
reader['stimulus']['presentation'] #Has the protocal data but my brain does not know the syntax # Now it knows.
# Example. reader['stimulus']['presentation']['Activation']['command'][:] will give array([b'-80:0:-80:40;-90:0:-90:10;-80:0:-80:50;-90:10:80:500;-80:0:-80:100;'],dtype='|S67'). There are 5 terms. Each term is of the form A:B:C:D. D is the duration to hold at various levels from A to C with B deveisions


#So, finally, The actual data is in 'acquisition'. 'general', 'session_description' and ['stimulus']['presentation'] provides some metadata. Other keys are useless.
# The actual data is in reader['acquisition']['timeseries']['Activation']['repetitions']['repetition1']['data']
reader['acquisition']['timeseries']['Activation']['repetitions']['repetition1']['amp'] #Provides some capacitance, v_offset etc
reader['acquisition']['timeseries']['Activation']['repetitions']['repetition1']['n_points'][:] #gives number of points

# y axis is in nA. Data acquisition is at 10kHz
