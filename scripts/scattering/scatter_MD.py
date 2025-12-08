
from scatter import XFEL, Crystal,Results,RESULTS_LOCAL_PATH
from core_functions import get_sim_params
import os.path as path
import numpy as np
import pickle


    # def set_ff_calculator_snapshot_time(self,time):
    #     self.xfel.target.ff_calculator.initialise_form_factor_params(time,time,self.max_q,self.photon_energy,t_fineness=self.t_fineness)




class MD_Crystal:
    def __init__(self,num_times, md_struct_path, allowed_atoms, t_cutoff_frac=None, positional_stdv = 0, is_damaged=True, include_symmetries = False, rocking_angle = 0.3,
                 cell_packing = "SC", CNO_to_N = False, supercell_scale = 1,num_supercells=1, supercell_simulations = 1,
                 S_to_N=False,convert_excluded_elements_to_H=False,convert_excluded_elements_to_N=False,allow_skip_species=False,random_waters=None,
                 use_bfactors=True,zero_bfactors=False):
        self.current_snapshot_time=None
        self.md_struct_path = md_struct_path
        self.is_damaged = is_damaged
        assert not include_symmetries, "symmetries not implemented for MD input"
        self.crystal_kwargs =  {k: v for k, v in locals().items() if k not in ("self","num_times","md_struct_path","is_damaged","t_cutoff_frac")}

        assert path.exists(self.md_struct_path), f"{self.md_struct_path} not found!" 
        T = self.read_times(num_times, self.md_struct_path)

        if t_cutoff_frac is not None: 
            truncted_T = []
            assert 0 <= t_cutoff_frac <= 1
            for t in T:
                if t <= t_cutoff_frac*T[-1]:
                    truncted_T.append(t)
            T = truncted_T


        #t_fineness = num_times-1
        # in picoseconds

        if num_times is None:
            if not is_damaged:
                self.times=[T[0]]
            elif len(T)>0:
                self.times=T[1:]
        else:
            times_to_aim_for = [T[0] + n/(num_times-1)*(T[-1]-T[0]) for n in range(num_times)]  # test: self.times = T[0:2]
            self.times = self.get_nearest_time(times_to_aim_for,T,tol_fs=1)
            if not is_damaged:
                self.times = [T[0]]
        print("chose times:", [t_pico*1e3 for t_pico in self.times])
            
         
    @staticmethod
    def get_nearest_time(times,allowed_times,tol_fs=None):
        allowed_times=np.array(allowed_times)
        if tol_fs is None:
            tol_fs = min(5e-1,(times[-1]-times[0])/100) 
        
        tol = tol_fs/1000 # convert to ps
        #n = np.argmin(np.abs(self.timeData[None,:] - time[:,None]))
        #assert np.all(np.abs(self.timeData[n]-time)<tol) , f"would use time at {self.timeData[n]} fs not {t} fs" 

        n = []
        for t in times:
            n.append(np.argmin(np.abs(t - allowed_times)))
            assert np.abs(allowed_times[n[-1]]-t)<tol , f"would use time at {allowed_times[n[-1]]*1e3} fs not {t*1e3} fs" 

        return allowed_times[np.array(n)]

    @staticmethod
    def read_times(num_times, md_struct_path,t_cutoff=None):
        times:list[float] = []
        with open(md_struct_path) as f:
            for line in f:
                if line.startswith("TITLE"):
                    t = float(line.split()[-1])
                    if t_cutoff is None or t <= t_cutoff:
                        times.append(t)
        return times

    
    EXCLUDE_WATER=True
    def set_crystal_snapshot(self,t:float,exclude_water=False,exclude_SOL=True,exclude_water_H=True):
        # if not self.is_damaged:
        #     assert self.current_snapshot_time is None
        #     t = self.times[0]
        #     print(t)

        # if self.current_snapshot_time == t and skip_if_time_unchanged:
        #     if t!= self.times[0]:
        #         print("Warning: Reusing mid-dynamics snapshot")
        #     return self.crystal_snapshot
            

        tmp_file_path = path.abspath(path.join(__file__ ,"../","snapshot.pdb"))

        assert t in self.times, f"{t} not found in times ({self.times})"
        snapshot_lines:list[str] = []

        with open(self.md_struct_path) as f:
            reading_block=False
            true_resnum=0
            last_read_resnum=-1 
            last_resname = ""
            for line in f:
                if line.startswith("TITLE") and float(line.split()[-1]) == t:
                    reading_block = True
                if not reading_block:
                    continue

                # Make res numbers unique
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    read_resnum=int(line[22:26])
                    resname = line[17:20]
                    name = line[12:16]
                    if resname == "HOH" and exclude_water:
                        continue
                    if name in ("HW1","HW2") and exclude_water_H:
                        continue
                    if resname == "SOL" and exclude_SOL:
                        continue
                    if read_resnum != last_read_resnum or resname!=last_resname:
                        true_resnum+=1
                    
                    
                    line = line[:22]+str(true_resnum) + ' ' * (4 - len(str(true_resnum))) + line[26:]
                    
                    last_read_resnum = read_resnum
                    last_resname = resname

                snapshot_lines.append(line)
                if line.startswith("TER"):
                    continue
                if line.startswith("ENDMDL"):
                    break

        assert reading_block
        with open(tmp_file_path,'w') as f_snap:
            #f_snap.writelines([f"{l}\n" for l in snapshot_lines])
            f_snap.writelines(snapshot_lines)

        self.crystal_snapshot = Crystal(tmp_file_path,is_damaged=self.is_damaged,**self.crystal_kwargs)
        self.current_snapshot_time = t
        #self.crystal_snapshot.plot_me()
    # def set_base_crystal(self, struct_file_path, allowed_atoms, positional_stdv = 0, is_damaged=True, include_symmetries = None, rocking_angle = 0.3, cell_packing = "SC", CNO_to_N = False, supercell_scale = 1,num_supercells=1, supercell_simulations = 1, S_to_N=False,convert_excluded_elements_to_N=False,random_waters=None,use_bfactors=True,zero_bfactors=False):
    #     args = {k: v for k, v in locals().items() if k != "self"}
    #     self.base_crystal = Crystal(**args) 
    # def set_time(self):
    #     self.base_crystal.positi



# MolDStruct or any dynamical structure input
class MD_XFEL:
    def __init__(self, experiment_name, photon_energy, detector_distance_mm=100, q_minimum = None, q_cutoff = None, max_miller_idx = None, screen_type = "hemisphere", num_orients_crys=1, orientation_axis_crys = None, x_orientations = 1, y_orientations = 1, pixels_per_ring = 400, num_rings = 50,t_fineness=0,SPI_y_rotation = 0,SPI_x_rotation = 0,SPI_z_rotation = 0,all_miller_indices=False, custom_cell_dims_for_miller_indices=None,override_max_q = False,miller_indices_override=None,spot_fraction_per_orient=None):
        assert t_fineness==0
        args = {k: v for k, v in locals().items() if k != "self"}
        self.xfel = XFEL(**args)
    
    def laser_multi_target(self, sim_data_handle, sim_parent_dir_path, md_target_list : list[MD_Crystal], SPI_resolution = None, results_parent_dir = RESULTS_LOCAL_PATH, circle_grid = False, pixels_across = 10, clear_output = False, random_orientation = False, SPI=False,do_not_integrate_times=False):
        #TODO THIS IS TERRIBLE, WE NEED THE INTERFERENCE!! 
        
        kwargs = {k: v for k, v in locals().items() if k not in (
            "self",
            "sim_data_handle",
            "sim_parent_dir_path",
            "md_target_list",
            )}
        I_tot = None
        for md_target in md_target_list:
            results = self.spooky_laser(sim_data_handle,sim_parent_dir_path,md_target,**kwargs)
            I_tot = I_tot + results.I if I_tot is not None else results.I
        out_results = results
        out_results.I = I_tot/len(md_target_list)
        with open(out_results.save_path,"wb") as pickle_out:
            pickle.dump(out_results,pickle_out)
        return out_results
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
            I = I + results.I if I is not None else results.I

            if not md_target.crystal_snapshot.is_damaged:
                break
        out_results = results #XXX 
        out_results.I = I

        with open(out_results.save_path,"wb") as pickle_out:
            pickle.dump(out_results,pickle_out)
        
        return out_results
        
    def set_orientation_set(self,orientation_set):
        self.xfel.set_orientation_set(orientation_set)
    def used_orientations(self):
        return self.xfel.used_orientations