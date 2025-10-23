#%%

GD_occupancy = 61

undamaged_EDR = 0.22592592592592592
# AC4DC_LF_EDR = 0.2201203502568035
# AC4DC_HF_EDR = 0.14757419134927216

# AC4DC_LF_EDR = 0.22448807690021502 # get from charge_contrast.py
# AC4DC_HF_EDR = 0.14357845509910588

#full

AC4DC_LF_EDR = 0.22010737717332937
AC4DC_HF_EDR = 0.14100013805308534

# no Gd or salt
# AC4DC_LF_EDR = 0.21526242046484392
# AC4DC_HF_EDR = 0.12258268770743477

charge_contrast_expected = (AC4DC_LF_EDR - AC4DC_HF_EDR)/undamaged_EDR*GD_occupancy
print(charge_contrast_expected)

#%% Carbon-based
# OLD
# LF_C_charge = 0.4384777598337092
# HF_C_charge = 2.6847904609202007
# LF_Gd_charge = 6.7197708159447505
# HF_Gd_charge =  39.816266334290525

# I-avged charges
LF_Q = dict(
    C = 0.29485160740356786,
    N = 0.21761973087515285,
    O = 0.1670088532164954,
    Gd = 6.691835350303332
)
HF_Q = dict(
    C = 2.5068978063978844,
    N = 2.7146018013601503,
    O = 2.7767624248616536,
    Gd = 39.737147616938515
)

LF_naive_Q = dict(
    C = 0.130034681641397,
    N = 0.09005350642415147,
    O = 0.06987722724765072,
    Gd = 6.691835350303332
) 

HF_naive_Q = dict(
    C = 1.9235137032919358,
    N = 1.9886329207592346,
    O = 1.9642125169346185,
    Gd = 39.737147616938515
)


def ratio(Q_dict,ignore_light_charge= False):
    light_atoms = dict(
        C = [6,20],
        N = [7,10],
        O = [8,10],
    )  
    heavy_atoms = dict(Gd=[64,1])
    light_density = heavy_density = 0
    for k, v in light_atoms.items():
        Z, N = v
        if ignore_light_charge:
            light_density += N*(Z) 
        else:
            light_density += N*(Z - Q_dict[k]) 
    for k, v in heavy_atoms.items():
        Z, N = v
        heavy_density += N*(Z - Q_dict[k]) 
    return heavy_density/light_density


def Gd_charge_contrast(dictLF,dictHF,ignore_light_atom_charges=False):
    return (ratio(dictLF,ignore_light_atom_charges) - ratio(dictHF,ignore_light_atom_charges))/ratio(dictLF,ignore_light_atom_charges)*61  # 61 because in limit of HF_Q = 0 we have max charge contrast equal to initial occupancy

Gd_effective_charge = Gd_charge_contrast(LF_Q,HF_Q,True)
charge_contrast_expected = Gd_charge_contrast(LF_Q,HF_Q)
# If we compute without considering cascades instigated by Gd and salt (of course, we still need to model primary ionization of Gd)
charge_contrast_expected_naive = Gd_charge_contrast(LF_naive_Q,HF_naive_Q)
print("Gd",Gd_effective_charge)
print("naive correction",charge_contrast_expected_naive)
print("correction", charge_contrast_expected)




# OLD (no salt no H 50% solvent)
# AC4DC_LF_EDR = 0.22448807690021502
# AC4DC_HF_EDR = 0.14357845509910588
# %%
