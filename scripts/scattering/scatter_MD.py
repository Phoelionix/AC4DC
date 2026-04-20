
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
#NUM_THREADS=1



class MD_Crystal:
    def __init__(self,name,num_times, md_struct_path, allowed_atoms, 
                 charges_path,debye_path,start_t,end_t, timespan_ps_MD,
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
        self.crystal_kwargs =  {k: v for k, v in locals().items() if k not in ("self","name","num_times","md_struct_path","charges_path","debye_path","start_t","end_t","timespan_ps_MD",
                                                                               "charges","is_damaged",
                                                                               "t_cutoff_frac","electronic_damage","nuclear_damage", "start_time_ps")}


        assert path.exists(self.md_struct_path), f"{self.md_struct_path} not found!" 
        T_ps = self.read_times(self.md_struct_path)



        if t_cutoff_frac is not None: 
            truncated_T = []
            assert 0 <= t_cutoff_frac <= 1
            for t in T_ps:
                if t <= t_cutoff_frac*T_ps[-1]:
                    truncated_T.append(t)
            T_ps = truncated_T


        #t_fineness = num_times-1
        # in picoseconds

        if num_times is None:
            if not is_damaged:
                self.times_ps=[T_ps[0]]
            elif len(T_ps)>0:
                self.times_ps=T_ps[1:]
            else:
                self.times_ps=T_ps

        else:
            if num_times == 2:
                self.times_ps=[T_ps[1],T_ps[-1]]
            else:
                assert False, "bugged"
                times_to_aim_for = [T_ps[0] + (n)/(num_times-1)*(T_ps[-1]-T_ps[0]) for n in range(num_times)]  # test: self.times_ps = T_ps[0:2]
                self.times_ps = list(self.get_nearest_time(times_to_aim_for,T_ps,tol_fs=1))
                if not is_damaged:
                    self.times_ps = [T_ps[0]]
        print("chose times:", [t_pico*1e3 for t_pico in self.times_ps])
            
        self.times_ac4dc_scale = np.array([t_pico*1e3 for t_pico in self.times_ps])
        self.times_ac4dc_scale = self.times_ac4dc_scale -timespan_ps_MD*1e3 + end_t
        if self.electronic_damage:
            assert len(T_ps)>0, T_ps
            all_charges=read_charges_binary(charges_path,debye_path)
            all_times=np.linspace(0,timespan_ps_MD,all_charges.shape[0])
            print(np.searchsorted(all_times,self.times_ps))
            print(all_times, self.times_ps)
            self.charges=all_charges[np.searchsorted(all_times,self.times_ps)]
         
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
        #     t = self.times_ps[0]
        #     print(t)

        # if self.current_snapshot_time == t and skip_if_time_unchanged:
        #     if t!= self.times_ps[0]:
        #         print("Warning: Reusing mid-dynamics snapshot")
        #     return self.crystal_snapshot
            
        if not self.nuclear_damage:
            t_pico=self.times_ps[0]


        assert t_pico in self.times_ps, f"{t_pico} not found in times ({self.times_ps})"
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

        assert reading_block, (t_pico, self.md_struct_path)
        tmp_file_path = path.abspath(path.join(__file__ ,"../",f"{self.name}_snapshot-{t_pico*1e3}fs.pdb"))
        with open(tmp_file_path,'w') as f_snap:
            #f_snap.writelines([f"{l}\n" for l in snapshot_lines])
            f_snap.writelines(snapshot_lines)

        prepare_pdb(tmp_file_path,tmp_file_path,allow_no_altloc=True,repeated_names_altlocs_are_new_residues=True) # Otherwise Bio.PDB.PDBParser will silently ignore repeat residues!!!

        self.crystal_snapshot = Crystal(tmp_file_path,is_damaged=self.electronic_damage,
            use_intensity_for_time=t_pico*1e3 if (self.nuclear_damage and not self.electronic_damage) else None,
            charge_states=(None if not self.electronic_damage else self.charges[self.times_ps.index(t_pico)]),
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
                upper_resolution=scatter.q_to_res(self.xfel.max_q)*scatter.ang_per_bohr
                out_handle="SPI_out",self.xfel.experiment_name+"_"+rotation_str+"_"+f"{upper_resolution}_A"
                cell_intensity_log_path = path.abspath(path.join(__file__ ,"../",out_handle+".csv"))
                os.makedirs(os.path.dirname(cell_intensity_log_path),exist_ok=True)
                print(f"Writing intensities to {cell_intensity_log_path}")
                with open(cell_intensity_log_path,'w') as f:
                    #f.write(f"q, I, pixel_idx_x, pixel_idx_y\n")
                    assert all(mdt.times_ps==md_target_list[0].times_ps for mdt in md_target_list)
                    f.write(f"row, col, resolution, intensity at t={', '.join([str(t) for t in md_target_list[0].times_ps])} \n")
                    #for q_val, I_val in zip(np.array(results.q).flatten(), np.array(results.I).flatten()):
                    #assert results.q.shape == results.I.shape
                    for x_idx in range(len(results.q)):
                        for y_idx in range(len(results.q[x_idx])):
                            resolution=scatter.q_to_res(results.q[x_idx,y_idx])*scatter.ang_per_bohr
                            f.write(f"{x_idx}, {y_idx}, {resolution}, "+ ', '.join([f"{I[x_idx,y_idx]:.3e}" for I in results.I_snapshots])+"\n")

                #cutoff_log_intensity = -1
                #scatter.scatter_scatter_plot(SPI_result1=out_results,SPI_result2=None,radial_lim=radial_lim,plot_against_q = use_q,log_radial=log_radial,cmap=cmap,cmap2=cmap2,**kwargs)
                #for log_range in (10,20,30):
                #for log_range in (10,20,30):
                plot_kwargs=dict(
                    SPI=True,
                    SPI_result1=results,
                    results_parent_dir=results_parent_dir,
                    spi_full_rings_only=False,
                    show_grid=True,
                    show_plot=False,
                )
                for log_range in (15,):
                    scatter.stylin(self.xfel.experiment_name,None,self.xfel.max_q,
                                **plot_kwargs,
                                log_range=log_range,
                                plot_handle="integrated_"+out_handle,
                                )
                for log_range in (20,):
                    scatter.stylin(self.xfel.experiment_name,None,self.xfel.max_q,
                                **plot_kwargs,
                                log_range=log_range,
                                plot_handle="integratedBrighter_"+out_handle,
                                )
                results.I = results.I_snapshots[0]
                #for log_range in (10,15,20):
                for log_range in (15,):
                    scatter.stylin(self.xfel.experiment_name,None,self.xfel.max_q,
                                log_range=log_range,
                                plot_handle="undamaged_"+out_handle,
                                **plot_kwargs,
                                )
                results.I = results.I_snapshots[-1]
                #for log_range in (10,15,20):
                for log_range in (15,):
                    scatter.stylin(self.xfel.experiment_name,None,self.xfel.max_q,
                                log_range=log_range,
                                plot_handle="damaged_"+out_handle,
                                **plot_kwargs,
                                )
                                
            else: # Reflections at Miller indices
                scatter.create_reflection_file(self.xfel.experiment_name,results_parent_dir=results_parent_dir,
                                    artificial_I_scale=artificial_I_scale,symmetry_override=reflections_file_symmetry_override)
                _, mtz_file = scatter.rfl_to_sca(self.xfel.experiment_name,create_mtz=USE_PHENIX)
                if USE_PHENIX and ground_truth_pdb is not None:
                    gen_true_phases=False
                    if gen_true_phases: # for e. dens. map making.
                        high_res=scatter.q_to_res(self.xfel.max_q)
                        cplx_data = scatter.phenix_fcalc(ground_truth_pdb,high_res,real=False)
                    FCALC_COMPARISON=False
                    if FCALC_COMPARISON: 
                        fcalc = scatter.phenix_fcalc_from_file(ground_truth_pdb,mtz_file,real=True)
                        scatter.phenix_R(ground_truth_pdb,fcalc) # Should be ~0
                    scatter.phenix_R(ground_truth_pdb,mtz_file)

            return out_results
        
        I_tot = None
        I_tot_snapshots=None
        
        # def process(i):
        #     md_target = md_target_list[i]
        #     print(f"Capturing trajectory {i+1}/{len(md_target_list)}")
        #     return self.fire_laser(sim_data_handle,sim_parent_dir_path,md_target,**kwargs)
        #     #I_tot = I_tot + results.I if I_tot is not None else results.I
        #     #I_tot_snapshots = I_tot_snapshots + results.I_snapshots if I_tot_snapshots is not None else results.I_snapshots

        # with Pool(12) as p:
        #     results_list = p.map(process,range(len(md_target_list)))
        # I_tot = np.zeros(results_list[0].I.shape)
        # I_tot_snapshots = np.zeros(results_list[0].I_snapshots.shape)
        # for result in results_list:
        #     I_tot+=results.I
        #     I_tot_sn
                
        num_snapshots = max(len(md_target.times_ps) for md_target in md_target_list)
        for i, md_target in enumerate(md_target_list):
            if len(md_target.times_ps)!=num_snapshots:
                print(f"skipping target {i}, has only {len(md_target.times_ps)} snapshots (expected {num_snapshots})")
                continue
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

        #param,_,_ = get_sim_params(sim_data_handle)

        global process
        def process(k):
            t_pico = md_target.times_ps[k] # gromacs output is in picoseconds
            t = md_target.times_ac4dc_scale[k]
            print(f"Snapshot t = {t} fs ({k+1}/{len(md_target.times_ps)})")
            
            try:
                if md_target.nuclear_damage:
                    md_target.set_crystal_snapshot(t_pico)
                elif k == 0:
                    md_target.set_crystal_snapshot(t_pico)
            except Exception as e:
                print(traceback.format_exc())
                print(e)
                print(f"Unexpected error for {sim_parent_dir_path} at snapshot {t} fs - skipping")
                return
            
            
            results:Results = self.xfel.fire_laser(t,t,sim_data_handle,sim_parent_dir_path,
                                   md_target.crystal_snapshot,
                                   **laser_kwargs)            
            print(f"Finished snapshot({k+1})")
            assert not np.any(np.isnan(results.I)), (results.I, k)
            #I_snapshots.append(results.I)
            #I = I + results.I if I is not None else results.I
            return results
        # with Pool(NUM_THREADS) as p:
        #      results_list = p.map(process,range(len(md_target.times_ps)))
        results_list = [process(_i) for _i in range(len(md_target.times_ps))]
        results_list = [r for r in results_list if r is not None]
        out_results = results_list[0]
        out_results.I_snapshots = np.array([result.I for result in results_list])
        out_results.I = np.sum(out_results.I_snapshots,axis=0)
        if len(results_list)>1:
            assert out_results.I.shape == results_list[1].I.shape
        else:
            print(f"Warning, number of snapshots is {len(results_list)}")

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