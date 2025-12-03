from matplotlib import pyplot
import numpy as np
def generate_charges(num_to_sample,stime,dQ_arrays,Q_array,lowest_charge,plot=True,plot_tag=""):
    # dQ : dQ[i][j] is an array of densities contributed to charge state j by charge state i for each time step 
    state_transitions = np.array(dQ_arrays)



    #N = max(num_to_sample,10000)
    N = max(num_to_sample,1)

    nstates_from,nsteps,nstates_to=dQ_arrays.shape
    assert nstates_from==nstates_to
    nstates=nstates_to
    assert nsteps==len(stime)

    for i in nstates:
        state_transitions[i,:,i] = Q_array # Include the probability of "no change" 


    charges = np.zeros((nsteps, N), dtype=np.int64) # matrix to store the charge of each the atom
    rng = np.random.default_rng()
    for step in range(1, nsteps):
        _, counts = np.unique(charges[step-1],return_counts=True)
        for from_charge_idx, dQ in enumerate(dQ_arrays):
            num_from=counts[from_charge_idx]
            transitions = rng.multinomial(num_from, state_transitions)  # how many atoms of from_charge_idx transitioned to each charge state
            
            new_vals = np.array([[i,]*count for i,count in enumerate(transitions)])
            new_vals.flatten()
            new_vals=np.random.shuffle(new_vals)


            atom_idxes = (charges[step-1]==from_charge_idx)
            charges[step][atom_idxes]=new_vals
    charges+=lowest_charge 


    if plot:
        
        fig, ax = pyplot.subplots()
        ax.plot(stime, Q_array)
        ax.set(xlabel="Time (s)", ylabel="Fractional population")
        ax.legend(["+%d" % cs for cs in range(nstates)], fontsize=10)
        pyplot.savefig(f"original_charge_dist{plot_tag}.png")
        pyplot.close()
        #####

        fig, ax = pyplot.subplots()
        ax.plot(stime, np.dot(Q_array, range(nstates)), label="cretin")
        ax.plot(stime, np.mean(charges, axis=1), label="discretized")
        ax.set(xlabel="Time (s)", ylabel="Average ionization")
        ax.legend(fontsize=10)
        pyplot.savefig(f"avg_ionization{plot_tag}.png")
        pyplot.close()
        #####

        plot_tag = "_"+plot_tag if plot_tag is not "" else ""
        fig, ax = pyplot.subplots()
        for cs in range(nstates):
            ax.plot(stime, np.sum(charges == cs, axis=1))
        ax.set(xlabel="Time (s)", ylabel="Charge state distribution among atoms")
        ax.legend(["+%d" % cs for cs in range(nstates)], fontsize=10)
        pyplot.savefig(f"modified_charge_dist{plot_tag}.png")
        pyplot.close()
        #####

    charges = charges[:,:num_to_sample]

    return np.swapaxes(charges,0,1)