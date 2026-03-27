
from scatter import XFEL, Crystal,Results,RESULTS_LOCAL_PATH,USE_PHENIX
import scatter
from core_functions import get_sim_params
import os.path as path
import os
import numpy as np
import pickle
import sys, pathlib
import traceback
sys.path.append(str(pathlib.Path(__file__).parent.parent.parent))
from pipeline_scripting.UntanglerStuff.UntangleFunctions import prepare_pdb

    # def set_ff_calculator_snapshot_time(self,time):
    #     self.xfel.target.ff_calculator.initialise_form_factor_params(time,time,self.max_q,self.photon_energy,t_fineness=self.t_fineness)




class MD_Crystal:
    def __init__(self,name,num_times, md_struct_path, allowed_atoms, 
                 charges_path,debye_path,start_t,end_t, 
                 t_cutoff_frac=None, positional_stdv = 0, is_damaged=True, include_symmetries = False, rocking_angle = 0.3,
                 cell_packing = "SC", CNO_to_N = False, supercell_scale = 1,num_supercells=1, supercell_simulations = 1,
                 S_to_N=False,convert_excluded_elements_to_H=False,convert_excluded_elements_to_N=False,allow_skip_species=False,random_waters=None,
                 use_bfactors=True,zero_bfactors=False,ignore_water_H=True,electronic_damage=True,nuclear_damage=True):
        self.name=name
        self.current_snapshot_time=None
        self.md_struct_path = md_struct_path
        self.electronic_damage = is_damaged and electronic_damage
        self.nuclear_damage = is_damaged and nuclear_damage
        assert not include_symmetries, "symmetries not implemented for MD input"
        self.crystal_kwargs =  {k: v for k, v in locals().items() if k not in ("self","name","num_times","md_struct_path","charges_path","debye_path","start_t","end_t",
                                                                               "charges","is_damaged",
                                                                               "t_cutoff_frac","electronic_damage","nuclear_damage")}


        assert path.exists(self.md_struct_path), f"{self.md_struct_path} not found!" 
        T = self.read_times(self.md_struct_path)
        if self.electronic_damage:
            assert len(T)>0, T
            all_charges=read_charges_binary(charges_path,debye_path)
            all_times=np.linspace(start_t,end_t,all_charges.shape[0])
            self.charges=all_charges[np.searchsorted(all_times,T)]



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
                self.times=T

        else:
            if num_times == 2:
                self.times=[T[1],T[-1]]
            else:
                assert False, "bugged"
                times_to_aim_for = [T[0] + (n)/(num_times-1)*(T[-1]-T[0]) for n in range(num_times)]  # test: self.times = T[0:2]
                self.times = list(self.get_nearest_time(times_to_aim_for,T,tol_fs=1))
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
    def read_times(md_struct_path,t_cutoff=None):
        times:list[float] = []
        with open(md_struct_path) as f:
            for line in f:
                if line.startswith("TITLE"):
                    t = float(line.split()[-1])
                    if t_cutoff is None or t <= t_cutoff:
                        times.append(t)
        return times

    
    def set_crystal_snapshot(self,t_pico:float,exclude_water=False,exclude_SOL=True,exclude_water_H=True):
        # if not self.is_damaged:
        #     assert self.current_snapshot_time is None
        #     t = self.times[0]
        #     print(t)

        # if self.current_snapshot_time == t and skip_if_time_unchanged:
        #     if t!= self.times[0]:
        #         print("Warning: Reusing mid-dynamics snapshot")
        #     return self.crystal_snapshot
            
        if not self.nuclear_damage:
            t_pico=self.times[0]

        tmp_file_path = path.abspath(path.join(__file__ ,"../",f"{self.name}_snapshot-{t_pico*1e3}fs.pdb"))

        assert t_pico in self.times, f"{t_pico} not found in times ({self.times})"
        snapshot_lines:list[str] = []

        print(f"Loading snapshot t = {t_pico*1e3} fs from {self.md_struct_path}")
        with open(self.md_struct_path) as f:
            reading_block=False
            true_resnum=0
            last_read_resnum=-1 
            last_resname = ""
            for line in f:
                if line.startswith("TITLE") and float(line.split()[-1]) == t_pico:
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
                    
                    resnum_to_use=int(true_resnum%1e4)
                    line = line[:22]+str(resnum_to_use) + ' ' * (4 - len(str(resnum_to_use))) + line[26:]
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

        prepare_pdb(tmp_file_path,tmp_file_path,allow_no_altloc=True,repeated_names_altlocs_are_new_residues=True) # Otherwise Bio.PDB.PDBParser will silently ignore repeat residues!!!

        self.crystal_snapshot = Crystal(tmp_file_path,is_damaged=self.electronic_damage,
            use_intensity_for_time=t_pico*1e3 if (self.nuclear_damage and not self.electronic_damage) else None,
            charge_states=(None if not self.electronic_damage else self.charges[self.times.index(t_pico)]),
            **self.crystal_kwargs
        )
        os.remove(tmp_file_path)
        self.current_snapshot_time = t_pico
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
    
    def laser_multi_target(self, sim_data_handle, sim_parent_dir_path, md_target_list : list[MD_Crystal], SPI_resolution = None, results_parent_dir = RESULTS_LOCAL_PATH, circle_grid = False, pixels_across = 10, clear_output = False, random_orientation = False, SPI=False,do_not_integrate_times=False,
                           artificial_I_scale=1,reflections_file_symmetry_override=None,update_reflection_file_every_target=True,ground_truth_pdb=None):
        #TODO Add interference between targets.
        
        kwargs = {k: v for k, v in locals().items() if k not in (
            "self",
            "sim_data_handle",
            "sim_parent_dir_path",
            "md_target_list",
            "artificial_I_scale",
            "reflections_file_symmetry_override",
            "update_reflection_file_every_target",
            "ground_truth_pdb",
            )}
        def save_results():
            out_results = results # XXX nonlocal
            assert len(md_target_list)>0
            out_results.I = I_tot/len(md_target_list)
            out_results.I_snapshots = I_tot_snapshots/len(md_target_list)
            with open(out_results.save_path,"wb") as pickle_out:
                pickle.dump(out_results,pickle_out)

            if kwargs["SPI"]:
                assert self.xfel.num_x_orientations==self.xfel.num_y_orientations==1
                rotation_str = '_'.join([f"{vars(self.xfel)['input_SPI_'+s+'_rotation']}" for s in ('x','y','z') ])
                cell_intensity_log_path = path.abspath(path.join(__file__ ,"../","SPI_out",self.xfel.experiment_name+"_"+rotation_str+".csv"))
                os.makedirs(os.path.dirname(cell_intensity_log_path),exist_ok=True)
                print(f"Writing intensities to {cell_intensity_log_path}")
                with open(cell_intensity_log_path,'w') as f:
                    #f.write(f"q, I, pixel_idx_x, pixel_idx_y\n")
                    assert all(mdt.times==md_target_list[0].times for mdt in md_target_list)
                    f.write(f"row, col, resolution, intensity at t={', '.join([str(t) for t in md_target_list[0].times])} \n")
                    #for q_val, I_val in zip(np.array(results.q).flatten(), np.array(results.I).flatten()):
                    #assert results.q.shape == results.I.shape
                    for x_idx in range(len(results.q)):
                        for y_idx in range(len(results.q[x_idx])):
                            resolution=scatter.q_to_res(results.q[x_idx,y_idx])*scatter.ang_per_bohr
                            f.write(f"{x_idx}, {y_idx}, {resolution}, "+ ', '.join([f"{I[x_idx,y_idx]:.3e}" for I in results.I_snapshots])+"\n")

                #cutoff_log_intensity = -1
                #scatter.scatter_scatter_plot(SPI_result1=out_results,SPI_result2=None,radial_lim=radial_lim,plot_against_q = use_q,log_radial=log_radial,cmap=cmap,cmap2=cmap2,**kwargs)
                #for log_range in (10,20,30):
                for log_range in (10,20,30):
                    scatter.stylin(self.xfel.experiment_name,None,self.xfel.max_q,
                                SPI=True,SPI_result1=results,results_parent_dir=results_parent_dir,
                                spi_full_rings_only=False,
                                log_range=log_range,
                                show_grid=True)
                results.I = results.I_snapshots[0]
                for log_range in (10,15,20):
                    scatter.stylin(self.xfel.experiment_name,None,self.xfel.max_q,
                                SPI=True,SPI_result1=results,results_parent_dir=results_parent_dir,
                                spi_full_rings_only=False,
                                log_range=log_range,
                                show_grid=True)
                results.I = results.I_snapshots[-1]
                for log_range in (10,15,20):
                    scatter.stylin(self.xfel.experiment_name,None,self.xfel.max_q,
                                SPI=True,SPI_result1=results,results_parent_dir=results_parent_dir,
                                spi_full_rings_only=False,
                                log_range=log_range,
                                show_grid=True)
                                
            else: # Reflections at Miller indices
                scatter.create_reflection_file(self.xfel.experiment_name,results_parent_dir=results_parent_dir,
                                    artificial_I_scale=artificial_I_scale,symmetry_override=reflections_file_symmetry_override)
                _, mtz_file = scatter.rfl_to_sca(self.xfel.experiment_name,create_mtz=USE_PHENIX)
                if USE_PHENIX and ground_truth_pdb is not None:
                    gen_true_phases=False
                    if gen_true_phases: # for e. dens. map making.
                        high_res=scatter.q_to_res(self.xfel.max_q)
                        cplx_data = scatter.phenix_fcalc(ground_truth_pdb,high_res,real=False) 
                    fcalc = scatter.phenix_fcalc_from_file(ground_truth_pdb,mtz_file,real=True)
                    scatter.phenix_R(ground_truth_pdb,mtz_file)
                    scatter.phenix_R(ground_truth_pdb,fcalc) # Should be ~0

            return out_results
        
        I_tot = None
        I_tot_snapshots=None
        for i, md_target in enumerate(md_target_list):
            print(f"Capturing trajectory {i+1}/{len(md_target_list)}")
            results = self.fire_laser(sim_data_handle,sim_parent_dir_path,md_target,**kwargs)
            I_tot = I_tot + results.I if I_tot is not None else results.I
            I_tot_snapshots = I_tot_snapshots + results.I_snapshots if I_tot_snapshots is not None else results.I_snapshots
            if update_reflection_file_every_target and i < len(md_target_list)-1:
                save_results()
        return save_results()


    #def fire_laser(self, start_time, end_time, sim_data_handle, sim_parent_dir_path, target : Crystal, SPI_resolution = None, results_parent_dir = RESULTS_LOCAL_PATH, circle_grid = False, pixels_across = 10, clear_output = False, random_orientation = False, SPI=False,do_not_integrate_times=False):
    def fire_laser(self, sim_data_handle, sim_parent_dir_path, md_target : MD_Crystal, SPI_resolution = None, results_parent_dir = RESULTS_LOCAL_PATH, circle_grid = False, pixels_across = 10, clear_output = False, random_orientation = False, SPI=False,do_not_integrate_times=False):
        laser_kwargs = {k: v for k, v in locals().items() if k not in (
            "self",
            "sim_data_handle",
            "sim_parent_dir_path",
            "md_target"
            )}

        param,_,_ = get_sim_params(sim_data_handle)

        I = None
        I_snapshots=[]
        for k, t_pico in enumerate(md_target.times):  # gromacs output is in picoseconds
            t = t_pico*1e3 + param["start_t"]  # AC4DC time
            print(f"Snapshot t = {t} fs ({k+1}/{len(md_target.times)})")
            
            try:
                if md_target.nuclear_damage:
                    md_target.set_crystal_snapshot(t_pico)
                elif k == 0:
                    md_target.set_crystal_snapshot(t_pico)
            except Exception as e:
                print(traceback.format_exc())
                print(e)
                print(f"Unexpected error for {sim_parent_dir_path} at snapshot {t} fs - skipping")
                continue
            
            
            results:Results = self.xfel.fire_laser(t,t,sim_data_handle,sim_parent_dir_path,
                                   md_target.crystal_snapshot,
                                   **laser_kwargs)            
            assert not np.any(np.isnan(results.I)), results.I
            I_snapshots.append(results.I)
            I = I + results.I if I is not None else results.I


        out_results = results #XXX 
        out_results.I = I
        out_results.I_snapshots=np.array(I_snapshots)

        return out_results
        
    def set_orientation_set(self,orientation_set):
        self.xfel.set_orientation_set(orientation_set)
    def get_used_orientations(self):
        return self.xfel.get_used_orientations()
    


def read_charges_binary(charges_path,debye_path):
    # NOTE debye_path is purely used to get times #XXX
    data_1d=np.fromfile(debye_path,dtype=np.float32)
    num_timesteps = len(data_1d)

    charges = np.fromfile(charges_path, dtype=np.ushort)    
    charges=charges.reshape((num_timesteps,-1))
    print(f"Charges shape={charges.shape}")
    return charges