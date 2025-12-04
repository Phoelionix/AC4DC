from matplotlib import pyplot
import numpy as np
def generate_charges(num_to_sample,stime,dQ_arrays,charge_data,Z,plot=True,plot_tag="",species_name_for_plot="atoms"):
    # dQ : dQ[i][t][j] is an array of densities contributed to states with j bound electrons by states with i bound electrons over time step t-1 to t
    # Z: Atomic number
    state_transitions = np.array(dQ_arrays) 
    # So charge in ascending order
    state_transitions = state_transitions[::-1] 
    state_transitions = state_transitions[...,::-1]

    #N = max(num_to_sample,10000)
    N = max(num_to_sample,1)

    nstates_from,nsteps,nstates_to=state_transitions.shape
    assert nstates_from==nstates_to
    nstates=nstates_to
    assert nsteps==len(stime)

    og_state_transitions=np.copy(state_transitions)
    for i in range(nstates):
        assert charge_data.shape==state_transitions.shape[1:]
        state_transitions[i,:,i] = (charge_data[:,i]-np.sum(state_transitions[i,:],axis=1)) # Include the probability of "no change" 
        force_min_zero=True
        if force_min_zero:
            state_transitions[i,:,i]=np.maximum(state_transitions[i,:,i],0)
        #state_transitions[i,:,i] = charge_data[:,i] # Include the probability of "no change" 
    del i

    charges = np.full((nsteps, N),99, dtype=np.int32) # matrix to store the charge of each the atom
    charges[0]*=0
    rng = np.random.default_rng()
    #for step in range(1, nsteps):
    min_charge=Z+1-nstates
    for step in range(1, nsteps):
        unique, counts = np.unique(charges[step-1],return_counts=True)
        for from_charge, num_from in zip(unique,counts):
            TEMPORARY_BUG_PATCH=True
            if TEMPORARY_BUG_PATCH and from_charge==Z:
                atom_idxes = (charges[step-1]==from_charge)
                charges[step][atom_idxes]=from_charge
                continue
            pvals = state_transitions[from_charge][step]/np.sum(state_transitions[from_charge][step])

            # handle roundoff error
            if np.sum(pvals)!=1:
                j = np.argmax(pvals)
                mask = np.full(pvals.shape,True,dtype=bool)
                mask[j]=False
                pvals[j]=1-np.sum(pvals[mask])
            try:
                transitions = rng.multinomial(num_from, pvals)  # how many atoms of from_charge_idx transitioned to each charge state
            except Exception as e:
                print("-----")
                print(from_charge,num_from,step,Z)
                print()
                print(charge_data[step-1])
                print("Charge data")
                print(charge_data[step])
                print("dQ arrays")
                print(dQ_arrays[-from_charge-1][step])
                print("state transitions")
                print(state_transitions[from_charge][step])
                print("state transitions pre change")
                print(og_state_transitions[from_charge][step-1])
                print(og_state_transitions[from_charge][step])
                print("no change")
                print((charge_data[:,from_charge]-np.sum(og_state_transitions[from_charge,:],axis=1))[step])
                print("from charge state")
                print(np.sum(og_state_transitions[from_charge,:],axis=1)[step])
                print("sum of changes")
                print(np.sum(og_state_transitions[from_charge,step]))
                print("pvals sum")
                print(np.sum(pvals))

                raise e
            #print(transitions)
            #new_vals = [[i,]*count for i,count in enumerate(transitions)]
            #new_vals = [v for sublist in new_vals for v in sublist]
            new_vals = [v for i, count in enumerate(transitions) for v in [i,]*count]
            new_vals= min_charge + np.array(new_vals)
            np.random.shuffle(new_vals)
            atom_idxes = (charges[step-1]==from_charge)
            charges[step][atom_idxes]= new_vals
    


    if plot:
        
        fig, ax = pyplot.subplots()
        ax.plot(stime, charge_data/np.sum(charge_data[0]))
        ax.set(xlabel="Time (s)", ylabel="Fractional population")
        ax.legend(["+%d" % cs for cs in range(nstates)], fontsize=10)
        pyplot.savefig(f"original_charge_dist{plot_tag}.png")
        pyplot.close()
        #####

        fig, ax = pyplot.subplots()
        ax.plot(stime, np.dot(charge_data/np.sum(charge_data[0]), range(nstates)), label="AC4DC")
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
        ax.set(xlabel="Time (s)", ylabel=f"Charge state distribution among {species_name_for_plot}")
        ax.legend(["+%d" % cs for cs in range(nstates)], fontsize=10)
        pyplot.savefig(f"modified_charge_dist{plot_tag}.png")
        pyplot.close()
        #####
        fig, ax = pyplot.subplots()
        ax.plot(stime, charges[:, :5])
        ax.set(xlabel="Time (s)", ylabel="Charge", title="the charge for a couple of atoms")
        pyplot.savefig(f"sample_states{plot_tag}.png")
        pyplot.close()

    charges = charges[:,:num_to_sample]

    return np.swapaxes(charges,0,1)