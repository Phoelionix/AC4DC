#%%
#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy
import scipy.optimize




def read_validation_file(cif_file):
    data_dict = dict (h='_refln.index_h',
                    k='_refln.index_k',
                    l='_refln.index_l',
                    I_meas='_refln.intensity_meas',
                    I_sigma='_refln.intensity_sigma'
                )
    rows = []
    with open(cif_file, 'r') as f:
        lines = f.readlines()
        j = 0
        # Move to part of file containing reflections 
        for line in lines:
            if line.strip().startswith('_refln.wavelength_id'): # TODO not all start with this...
                break
            j+=1
            
        headers = []
        # Collect header lines
        for line in lines[j:]:   
            if not line.strip().startswith('_'):
                break
            headers.append(line.strip())
            j+=1
        
        for line in lines[j:]:
            #assert ('_refln.intensity_meas' in headers and '_refln.intensity_sigma' in headers)

            # Read data lines
            if line[0].strip().isnumeric():
                data = line.split()
                row = []
                for k,v in data_dict.items():
                    idx = headers.index(v)
                    row.append(float(data[idx]))
                rows.append(row)
            else:
                break
    return pd.DataFrame(rows, columns=list(data_dict.keys()))




def plot_linear(df):
    plt.figure()
    plt.scatter(df['I_meas'], df['I_sigma'], s=1)
    plt.xlabel('I_meas')
    plt.ylabel('$\sigma_I$')
    plot_curve_fit()
    plt.tight_layout()
    plt.xlim(0,100)
    plt.ylim(0,100)


def plot_sqrt(df):
    plt.figure()
    plt.scatter(np.sqrt(df['I_meas']), df['I_sigma'], s=1)
    plt.xlabel('sqrt I_meas')
    plt.ylabel('$\sigma_I$')
    plt.tight_layout()

def plot_photon_count_error_analysis(df):
    plt.figure()
    y = df['I_sigma']/np.sqrt(df['I_meas'])
    plt.scatter(df['I_meas'],y, s=10)
    plt.xlabel('I_meas')
    plt.ylabel('$\sigma_I$')
    plt.tight_layout()
    plt.xscale('log')
    plt.yscale('log')


def plot_log(df,show_labels=False):
    plt.figure()
    plt.scatter(df['I_meas'], df['I_sigma'], s=1)
    if show_labels:
        for i in range(len(df['k'])):
            plt.annotate(f"{df['h'][i].astype(str)},{df['k'][i].astype(str)},{df['l'][i].astype(str)}",(df['I_meas'][i],df['I_sigma'][i]))
    #plt.annotate(df['k'].astype(str),(df['I_meas'],df['I_sigma']))

    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('I_meas')
    plt.ylabel('$\sigma_I$')
    plot_curve_fit(log=True)
    #plt.grid(True, which="both", ls="--")
    plt.tight_layout()

# def noise_and_photon_count_curve(x,a,b,c):
#     # if np.any(a*np.sqrt(x)+b*x +c < 0):
#     #     return -9999999999999
#     # if np.any(c < 0):
#     #     return -9999999999999
#     return a*np.sqrt(x)+c+10

def noise_and_photon_count_curve(x,a,b,c,d):
    if a < 0:  # Enforce positive photon counting error. Probably not how you're meant to do it.
        return -9999999
    if b < 0 or d < 0:  # Enforce that other sources of error are positive
        return -9999999
    y = a * x**0.5 +b*x**c + d 
    return y
# def noise_and_photon_count_curve(x,a,b,c,d):

#     #y = a * x**0.5 +b*x**c + d 
#     y = a * x**0.5 +20
#     return y
# def noise_and_photon_count_curve(x,a,b,c,d,e):
#     y = a * x**b +d*x**c + e
#     return y
def plot_curve_fit(log=False):
    curve = noise_and_photon_count_curve
    popt, pcov = scipy.optimize.curve_fit(curve,df['I_meas'],df['I_sigma'],sigma=df['I_meas'], absolute_sigma=True)
    print(popt,pcov)
    #fit = np.polynomial.Polynomial.fit(df['I_meas'],df['I_sigma'],10)
    if log:
        x = np.logspace(0,5,100)
    else:
        x = np.linspace(0,5e4,100000)
    plt.scatter(x,curve(x,*popt),s=1)

if __name__ == "__main__":
    cif_file = "/home/speno/PhenixWorkspace/data/4et8-sf.cif"
    df = read_validation_file(cif_file)
    plot_log(df)
    plot_linear(df)
#plot_sqrt()
#plot_photon_count_error_analysis()




#plt.show()
# %%
