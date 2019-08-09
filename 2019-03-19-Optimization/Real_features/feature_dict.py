# exec(open('Somatic model/feature_dict.py').read())
import pandas as pd

feature_range = {}

feature_range['E_rest_150.0_min'] = -72.75347392
feature_range['E_rest_150.0_max'] = -60.64017436
# scale = 2/(feature_range['E_rest_150.0_max'] - feature_range['E_rest_150.0_min'])
# feature_range['E_rest_150.0_scale'] = -60.64017436

feature_range['AP1_amp_150.0_min'] = 91.05240275
feature_range['AP1_amp_150.0_max'] =123.3158946

feature_range['APp_amp_150.0_min'] =89.25186564
feature_range['APp_amp_150.0_max'] =114.1301036

feature_range['AP1_width_150.0_min'] =0.00080002
feature_range['AP1_width_150.0_max'] =0.001400035

feature_range['APp_width_150.0_min'] =0.000950024
feature_range['APp_width_150.0_max'] =0.001700043

feature_range['AP1_thresh_150.0_min'] =-54.9621582
feature_range['AP1_thresh_150.0_max'] =-44.43359375

feature_range['APp_thresh_150.0_min'] =-48.5534668
feature_range['APp_thresh_150.0_max'] =-35.55297852

feature_range['AP1_lat_150.0_min'] =0.001102063
feature_range['AP1_lat_150.0_max'] =0.043507304

feature_range['ISI1_150.0_min'] =0.005500138
feature_range['ISI1_150.0_max'] =0.035450886

feature_range['ISIl_150.0_min'] =0.025201008
feature_range['ISIl_150.0_max'] =0.099902498

feature_range['ISIavg_150.0_min'] =0.017354005
feature_range['ISIavg_150.0_max'] =0.057626441

feature_range['freq_150.0_min'] =16
feature_range['freq_150.0_max'] =58

feature_range['Adptn_id_150.0_min'] =0.373612824
feature_range['Adptn_id_150.0_max'] =0.907872697

feature_range['fAHP_AP1_amp_150.0_min'] =4.974965413
feature_range['fAHP_AP1_amp_150.0_max'] =26.44668986

feature_range['fAHP_APp_amp_150.0_min'] =11.32262166
feature_range['fAHP_APp_amp_150.0_max'] =32.06192424

feature_range['mAHP_AP1_amp_150.0_min'] =13.89707642
feature_range['mAHP_AP1_amp_150.0_max'] =19.96658427

feature_range['mAHP_APp_amp_150.0_min'] =10.19347127
feature_range['mAHP_APp_amp_150.0_max'] =26.87393595

feature_range['mAHP_AP1_dur_150.0_min'] =0.282738095
feature_range['mAHP_AP1_dur_150.0_max'] =0.488188976

feature_range['mAHP_APp_dur_150.0_min'] =0.24024024
feature_range['mAHP_APp_dur_150.0_max'] =0.821428571

feature_range['ADP_AP1_amp_150.0_min'] =0
feature_range['ADP_AP1_amp_150.0_max'] =0

feature_range['ADP_APp_amp_150.0_min'] =17.4032959
feature_range['ADP_APp_amp_150.0_max'] =25.42174581

feature_range['mAHP_stimend_amp_150.0_min'] =-2.074595133
feature_range['mAHP_stimend_amp_150.0_max'] =2.136272176

feature_range['sAHP_stimend_amp_150.0_min'] =-2.181726074
feature_range['sAHP_stimend_amp_150.0_max'] = 0.030559285


feature_range_df = pd.DataFrame(feature_range, index=[0])
feature_range_df = feature_range_df.transpose()
feature_range_df = feature_range_df.rename(index=str, columns={0:'Min'})
a = list(feature_range_df['Min'])[1:]
a.append(0)
feature_range_df.insert(1, 'Max', a)

maxidx = []
for i in feature_range_df.index:
    if 'max' in i:
        maxidx.append(i)

feature_range_df = feature_range_df.drop(index=maxidx)
feature_range_df['feature'] = feature_range_df.index
feature_range_df['feature'] = [i[:-4] for i in feature_range_df['feature']]
feature_range_df = feature_range_df.set_index('feature')
