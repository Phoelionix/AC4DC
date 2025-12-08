from core_functions import get_sim_params
import sys
import io
import contextlib

sim_handle=sys.argv[1]
param_keys=sys.argv[2:]

# silence output
with contextlib.redirect_stdout(io.StringIO()) as f:
    sim_params=get_sim_params(sim_handle)[0]
sim_params["dt_ps"]=f'{sim_params["dt"]/1e3:.9f}'
print(
    ' '.join([ str(sim_params[k]) for k in param_keys ])
)