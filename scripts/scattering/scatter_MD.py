
from scatter import XFEL, Crystal,Results,RESULTS_LOCAL_PATH
from core_functions import get_sim_params
import os.path as path



    # def set_ff_calculator_snapshot_time(self,time):
    #     self.xfel.target.ff_calculator.initialise_form_factor_params(time,time,self.max_q,self.photon_energy,t_fineness=self.t_fineness)

class MD_Crystal:
    def __init__(self,num_times, md_struct_path, allowed_atoms, positional_stdv = 0, is_damaged=True, include_symmetries = False, rocking_angle = 0.3, cell_packing = "SC", CNO_to_N = False, supercell_scale = 1,num_supercells=1, supercell_simulations = 1, S_to_N=False,convert_excluded_elements_to_N=False,random_waters=None,use_bfactors=True,zero_bfactors=False):
        self.md_struct_path = md_struct_path
        assert not include_symmetries, "symmetries not implemented for MD input"
        self.crystal_kwargs =  {k: v for k, v in locals().items() if k not in ("self","num_times","md_struct_path")}

        assert path.exists(self.md_struct_path)
        T = self.read_times(num_times, self.md_struct_path)


        #t_fineness = num_times-1
        # in picoseconds
        self.times = [T[0] + n/(num_times-1)*(T[-1]-T[0]) for n in range(num_times)]
        print("set times:", [t_pico*1e3 for t_pico in self.times])
         
    @staticmethod
    def read_times(num_times, md_struct_path):
        times:list[float] = []
        with open(md_struct_path) as f:
            for line in f:
                if line.startswith("TITLE"):
                    times.append(float(line.split()[-1]))
        return times



    def set_crystal_snapshot(self,t:float):
        tmp_file_path = path.abspath(path.join(__file__ ,"../","snapshot.pdb"))

        assert t in self.times, f"{t} not found in times ({self.times})"
        snapshot_lines:list[str] = []

        with open(self.md_struct_path) as f:
            for line in f:
                if line.startswith("TITLE") and float(line.split()[-1]) == t:
                    snapshot_lines.append(line)
                    break

            for line in f:
                snapshot_lines.append(line)
                if line.startswith("TER"):
                    continue
                if line.startswith("ENDMDL"):
                    break
        with open(tmp_file_path,'w') as f_snap:
            f_snap.writelines([f"{l}\n" for l in snapshot_lines])

        self.crystal_snapshot = Crystal(tmp_file_path,**self.crystal_kwargs)
    # def set_base_crystal(self, struct_file_path, allowed_atoms, positional_stdv = 0, is_damaged=True, include_symmetries = None, rocking_angle = 0.3, cell_packing = "SC", CNO_to_N = False, supercell_scale = 1,num_supercells=1, supercell_simulations = 1, S_to_N=False,convert_excluded_elements_to_N=False,random_waters=None,use_bfactors=True,zero_bfactors=False):
    #     args = {k: v for k, v in locals().items() if k != "self"}
    #     self.base_crystal = Crystal(**args) 
    # def set_time(self):
    #     self.base_crystal.positi



# MolDStruct or any dynamical structure input
class MD_XFEL:
    def __init__(self, experiment_name, photon_energy, detector_distance_mm=100, q_minimum = None, q_cutoff = None, max_miller_idx = None, screen_type = "hemisphere", num_orients_crys=1, orientation_axis_crys = None, x_orientations = 1, y_orientations = 1, pixels_per_ring = 400, num_rings = 50,t_fineness=100,SPI_y_rotation = 0,SPI_x_rotation = 0,SPI_z_rotation = 0,all_miller_indices=False, custom_cell_dims_for_miller_indices=None,override_max_q = False,miller_indices_override=None,spot_fraction_per_orient=None):
        args = {k: v for k, v in locals().items() if k != "self"}
        self.xfel = XFEL(**args)
    
    #def spooky_laser(self, start_time, end_time, sim_data_handle, sim_parent_dir_path, target : Crystal, SPI_resolution = None, results_parent_dir = RESULTS_LOCAL_PATH, circle_grid = False, pixels_across = 10, clear_output = False, random_orientation = False, SPI=False,do_not_integrate_times=False):
    def spooky_laser(self, sim_data_handle, sim_parent_dir_path, md_target : MD_Crystal, SPI_resolution = None, results_parent_dir = RESULTS_LOCAL_PATH, circle_grid = False, pixels_across = 10, clear_output = False, random_orientation = False, SPI=False,do_not_integrate_times=False):
        laser_kwargs = {k: v for k, v in locals().items() if k not in (
            "self",
            "sim_data_handle",
            "sim_parent_dir_path",
            "md_target"
            )}

        param,_,_ = get_sim_params(sim_data_handle)

        I = None
        for t_pico in md_target.times:  # gromacs output is in picoseconds
            t = t_pico*1e3 + param["start_t"]  # AC4DC time
            print(f"Snapshot t = {t} fs")
            
            md_target.set_crystal_snapshot(t_pico)
            
            
            results:Results = self.xfel.spooky_laser(t,t,sim_data_handle,sim_parent_dir_path,
                                   md_target.crystal_snapshot,
                                   **laser_kwargs)
            I = I + results.I if I is not None else I

        
    def set_orientation_set(self,orientation_set):
        self.xfel.set_orientation_set(orientation_set)