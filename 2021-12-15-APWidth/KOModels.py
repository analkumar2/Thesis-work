Models = {}

Models['Model1'] = {
    "Parameters": {
        "notes": "",
        "Morphology": {"sm_len": 6.73077545020806e-05, "sm_diam": 6.73077545020806e-05},
        "Passive": {
            "Cm": 1.19e-10,
            "Rm": 609705941.9371204,
            "Em": 0.06022479405782892,
        },
        "Channels": {
            "Na_Chan": {
                "Gbar": 0.00012480475550619003,
                "Erev": 0.06,
                "Kinetics": "../../Compilations/Kinetics/Na_Chan_Custom4",
                "KineticVars": {
                    "m_vhalf_inf": -0.0316,
                    "h_vhalf_inf": -0.066,
                    "s_vhalf_inf": -0.033,
                },
            },
            "K_DR_Chan": {
                "Gbar": 1.0502259538910637e-7,
                "Erev": -0.09,
                "Kinetics": "../../Compilations/Kinetics/K_DR_Chan_Custom3",
                "KineticVars": {"n_F": 0.00306},
            },
            "K_A_Chan": {
                "Gbar": 1.008422244061249e-06,
                "Erev": -0.09,
                "Kinetics": "../../Compilations/Kinetics/K_A_Chan_Custom3",
                "KineticVars": {},
            },
            "K_M_Chan": {
                "Gbar": 4.016076778584836e-09,
                "Erev": -0.09,
                "Kinetics": "../../Compilations/Kinetics/K_M_Chan_Custom1",
                "KineticVars": {"factor": 3.3e-05},
            },
            "h_Chan": {
                "Gbar": 5.3739087243907273e-11,
                "Erev": -0.04,
                "Kinetics": "../../Compilations/Kinetics/h_Chan_Custom1",
                "KineticVars": {},
            },
            "Ca_L_Chan": {
                "Gbar": 1.0330355973445023e-10,
                "Erev": 0.14,
                "Kinetics": "../../Compilations/Kinetics/Ca_L_Chan_Custom1",
            },
            "K_SK_Chan": {
                "Gbar": 2.0605218247275846e-09,
                "Erev": -0.09,
                "Kinetics": "../../Compilations/Kinetics/K_SK_Chan_Custom4",
            },
        },
        "Ca_Conc": {
            "Ca_B": 75427936887.46373,
            "Ca_tau": 0.038,
            "Ca_base": 8e-05,
            "Kinetics": "../../Compilations/Kinetics/Ca_Conc_(Common)",
        },
    },
}