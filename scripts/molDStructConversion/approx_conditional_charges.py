from matplotlib import pyplot
import numpy
def generate_charges(num_to_sample,stime,P_charge,plot=True,plot_tag=""):
    N = max(num_to_sample,10000)

    nsteps,nstates=P_charge.shape
    assert nsteps==len(stime)
    charges = numpy.zeros((nsteps, N), dtype=numpy.int64)  # matrix to store the charge of each the atom
    rng = numpy.random.default_rng()
    for step in range(1, nsteps):
        prev = charges[step - 1].copy()
        target = rng.multinomial(N, P_charge[step])  # how many atoms of each charge state should exist
        current = numpy.bincount(prev, minlength=P_charge.shape[-1])  # how many atoms of each charge state exist now
        diff = target - current  # discrepancy

        for cs in range(1, nstates):
            if diff[cs] > 0:  # if we are missing atoms for this charge state
                sel = numpy.where(prev < cs)[0]
                if len(sel) > 0:
                    chosen = rng.choice(sel, size=diff[cs])
                    prev[chosen] += 1
        charges[step] = prev


    if plot:
        
        fig, ax = pyplot.subplots()
        ax.plot(stime, P_charge)
        ax.set(xlabel="Time (s)", ylabel="Fractional population")
        ax.legend(["+%d" % cs for cs in range(nstates)], fontsize=10)
        pyplot.savefig(f"original_charge_dist{plot_tag}.png")
        pyplot.close()
        #####

        fig, ax = pyplot.subplots()
        ax.plot(stime, numpy.dot(P_charge, range(nstates)), label="cretin")
        ax.plot(stime, numpy.mean(charges, axis=1), label="discretized")
        ax.set(xlabel="Time (s)", ylabel="Average ionization")
        ax.legend(fontsize=10)
        pyplot.savefig(f"avg_ionization{plot_tag}.png")
        pyplot.close()
        #####

        plot_tag = "_"+plot_tag if plot_tag is not "" else ""
        fig, ax = pyplot.subplots()
        for cs in range(nstates):
            ax.plot(stime, numpy.sum(charges == cs, axis=1))
        ax.set(xlabel="Time (s)", ylabel="Charge state distribution among atoms")
        ax.legend(["+%d" % cs for cs in range(nstates)], fontsize=10)
        pyplot.savefig(f"modified_charge_dist{plot_tag}.png")
        pyplot.close()
        #####

    charges = charges[:,:num_to_sample]

    return numpy.swapaxes(charges,0,1)