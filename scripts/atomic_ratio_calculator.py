#%%
import numpy as np
from diffpy.structure.spacegroups import GetSpaceGroup

# no hydrogen
#space_group_pdb_name = "C 1 2 1" 

class Solvent:
    molecule_idx_iterator=0
    def __init__(self,VS,VM=None,solution_density=None):
        self.molecules = []
        self.VS = VS # Solvent content %(v/v)
        self.VM = VM # Matthews Coefficient, Ang^3/Da  (note 1 Da -> 1 g/mol)
        self.solution_density=solution_density
    def add_molecule(self, molecule,M = None,solvent_v_on_v = None,solvent_w_on_v=None): # v_on_v and w_on_v given as percentages
        measure_check=solvent_v_on_v if (solvent_v_on_v is not None) else (solvent_w_on_v if (solvent_w_on_v is not None) else None)
        if measure_check is not None and measure_check < 1: 
            print(f"Warning: v_on_v or w_on_v is {measure_check}% - make sure it is given as a percentage")
        if solvent_w_on_v is not None:
            assert solvent_v_on_v is None
            assert molecule.density is not None
            assert False, "Not coded"
            #solvent_v_on_v = solvent_w_on_v*
        molecule.set_concentration(M,solvent_v_on_v,self.solution_density)
        self.molecules.append(molecule)
    def clear_molecules(self):
        self.molecules = []
    

class Molecule:
    def __init__(self,molar_mass=None,density_UNUSED=None,undiluted_molarity=None,**atom_dict):
        if undiluted_molarity is not None:
            print("WARNING this ignores density differences") # TODO
        self.idx = Solvent.molecule_idx_iterator
        Solvent.molecule_idx_iterator+=1
        self.undiluted_molarity = undiluted_molarity
        self.molar_mass=molar_mass
        #self.density=density_in_water
        self.atom_dict = atom_dict

        # Properties of solvent
        self.fraction_of_solvent_volume=None
        self.M = None
    def set_concentration(self,M = None,solvent_v_on_v = None,solution_density=None):
        have_measure = False
        for elem in solvent_v_on_v,M:
            if elem is not None:
                assert not have_measure, "Given multiple variables for concentration, should only pass one!"
                have_measure = True
        assert have_measure, "Missing concentration information"

        if M is not None:
            assert(self.molar_mass is not None and solution_density is not None)
            self.M = M
            self.fraction_of_solvent_volume = M*self.molar_mass/(solution_density*1e3)
        elif solvent_v_on_v is not None:
            self.fraction_of_solvent_volume = solvent_v_on_v/100
            if solution_density is not None and self.molar_mass is not None:
                self.M = self.fraction_of_solvent_volume*solution_density*1e3/self.molar_mass
            elif self.undiluted_molarity is not None:
                self.M = self.undiluted_molarity*self.fraction_of_solvent_volume
            else:
                print(f"Can't determine molarity of molecule in solution") 

        #print("M",self.M)
        #print("%(v/v)",self.fraction_of_solvent_volume*100)
        assert self.fraction_of_solvent_volume is not None



def calculate(solvent,protein_light_atoms,protein_heavy_atoms,lengths,angles,num_asymm_units):

    ## Protein ##

    num_protein_CNO = 0
    for val in protein_light_atoms.values():
        num_protein_CNO += val
    protein_hydrogens = {"H":1.01*num_protein_CNO}

    protein_atoms = protein_light_atoms | protein_heavy_atoms | protein_hydrogens

    ## Solvent ##

    # sg = GetSpaceGroup(space_group_pdb_name)
    # num_asymm_units = sg.num_primitive_sym_equiv #sg.symop_list[:sg.num_primitive_sym_equiv]

    non_water_solvent_volume_fraction = 0
    for molecule in solvent.molecules:
        non_water_solvent_volume_fraction += molecule.fraction_of_solvent_volume 
    assert non_water_solvent_volume_fraction <= 1, non_water_solvent_volume_fraction
    #
    #solvent_molecules = non_water_solvent_molecules.append([dict(H = 2, O = 1),1 - non_water_solvent_volume_fraction])
    solvent.add_molecule(water,solvent_v_on_v=(1-non_water_solvent_volume_fraction)*100)



    def cos(x):
        return np.cos(np.radians(x))

    V_cell = np.prod(lengths)*np.sqrt(1+2*np.prod([cos(x) for x in angles]) + np.sum([-np.square(cos(x)) for x in angles]))

    # Combined
    asym_atoms = dict(protein_atoms)
    #print(asym_atoms)

    # Num molecules in asymmetric unit
    def get_num_asymm_molecules(volume_fraction,M):
        return  (V_cell/num_asymm_units * M * volume_fraction
            * 1e-27*6.02214076e23)

    for molecule in solvent.molecules:
        N = get_num_asymm_molecules(solvent.VS/100,molecule.M)
        for element, num in molecule.atom_dict.items():
            if element in asym_atoms:
                asym_atoms[element] += N*num
            else:
                asym_atoms[element] = N*num



    print("-----------")
    for k,v in asym_atoms.items():
        print(k,f"{v:.3f}")
    print(f"\nV: {V_cell/num_asymm_units:.1f} Ang^3")

# MOLECULE DEFINITIONS

def PEG(peg_molar_mass):
    peg_n = (peg_molar_mass -18.02)/44.05
    return Molecule(peg_molar_mass,1.125,
        C=2*peg_n,
        H=4*peg_n+2,
        O=peg_n+1,
    ); 

PEG_8000 = PEG(8000)


sodium_cacodylate = Molecule(137.9977,1.1,
    C=2,
    H=7,
    As=1,
    O=2
)
glycerol = Molecule(92.09382,1.26,
    C=3,
    H=8,
    O=3
)

water = Molecule(18.02, #undiluted_molarity=55.56,
    H = 2,
    O = 1
)


NaCl = Molecule(58.44,# 1.02,
    Na = 1,
    Cl = 1                
)

CsCl = Molecule(168.36,# 1.02,
    Cs = 1,
    Cl = 1                
)

Gadoteridol = Molecule(558.69,#1.3,
    Gd = 1,
    C = 17,
    H = 29,
    N = 4,
    O = 7
)

sodium_acetate = Molecule(82.0343,#1.02,
C = 2,
H = 3,
Na = 1,
O = 2,
)

KI=Molecule(166.0028, #1.31,  # https://advancedthermo.com/electrolytes/density_KI.html
K=1,
I=1
)

PEG_6000=PEG(6000)

########################################

#%%
if __name__ == "__main__":

    # unit cell dimensions angstrom (paralleliped)
    lengths =  117.896, 64.199, 74.585
    angles = 90, 125.81, 90
    num_asymm_units = 8
    protein_light_atoms = dict(
        C = 946,
        N = 242,
        O = 282
    )
    protein_heavy_atoms = dict(
        S = 8
    )
    solvent = Solvent(55.1,2.74)
    solvent.add_molecule(PEG_8000,solvent_v_on_v=12)
    solvent.add_molecule(sodium_cacodylate,M=0.1)
    solvent.add_molecule(glycerol,solvent_v_on_v=6.25)


    non_water_solvent_molecules = [PEG_8000,sodium_cacodylate,glycerol]

    # non-water molecules and fractions in solution


    calculate(solvent,protein_light_atoms,protein_heavy_atoms,lengths,angles,num_asymm_units)





# %%
if __name__ == "__main__":
    lengths = 79.470,  79.470,   38.320  
    angles = 90, 90, 90
    num_asymm_units=8

    protein_light_atoms = dict(
        C = 632,
        N = 197,
        O = 193,
    )
    protein_heavy_atoms = dict(
        S=10
    )

    solvent = Solvent(41.73)
    #solvent.add_molecule(NaOAc,M=0.05)
    solvent.add_molecule(PEG_6000,solvent_v_on_v=16.7)
    solvent.add_molecule(NaCl,M=1.7)
    
    calculate(solvent,protein_light_atoms,protein_heavy_atoms,lengths,angles,num_asymm_units)

# %% 
if __name__ == "__main__":
    lengths = 79.470,  79.470,   38.320  
    angles = 90, 90, 90
    num_asymm_units=8

    protein_light_atoms = dict(
        C = 632*35.1/75,
        N = 197*35.1/75,
        O = 193*35.1/75,
    )
    protein_heavy_atoms = dict(
        S=10*35.1/75,
        #Gd=2
    )
    #solvent = Solvent(35.1,solution_density=1.1)
    solvent = Solvent(75,solution_density=1.1)
    #solvent.add_molecule(NaCl,M=1.71)
    #solvent.add_molecule(sodium_acetate,M=0.1)
    #solvent.add_molecule(Gadoteridol,M=0.1)
    
    calculate(solvent,protein_light_atoms,protein_heavy_atoms,lengths,angles,num_asymm_units)
# %%


# %% more solvent
if __name__ == "__main__":
    lengths = 79.470,  79.470,   38.320  
    angles = 90, 90, 90
    num_asymm_units=3

    protein_light_atoms = dict(
        C = 632,
        N = 197,
        O = 193,
    )
    protein_heavy_atoms = dict(
        S=10,
        Gd=2
    )

    solvent = Solvent(78.14)
    #solvent.add_molecule(NaOAc,M=0.05)
    #solvent.add_molecule(PEG_6000,solvent_v_on_v=16.7)
    #solvent.add_molecule(NaCl,solvent_v_on_v=10)
    #solvent.add_molecule(NaCl,M=1.71)
    solvent.add_molecule(NaCl,M=1.71)
    #solvent.add_molecule(NaCl,solvent_w_on_v=10)
    #solvent.add_molecule(NaCl,solvent_v_on_v=9.7)
    
    calculate(solvent,protein_light_atoms,protein_heavy_atoms,lengths,angles,num_asymm_units)
# %% 
# 7W6B
if __name__ == "__main__":
    lengths = 52.271,  422.928,   48.391
    angles = 90, 90, 90
    num_asymm_units=2

    protein_light_atoms = dict(
        C = 3635,
        N = 1022,
        O = 1216,
    )
    protein_heavy_atoms = dict(
        S = 4,
        Mg=1,
        Ca=1,
    )

    # https://advancedthermo.com/electrolytes/density_KI.html
    # But PEG lighter

    
    solvent = Solvent(59.48,solution_density=1.1) 
    solvent.add_molecule(KI,M=1)
    solvent.add_molecule(PEG(3350),solvent_v_on_v=25)
    #TODO 100 mM HEPES
    
    calculate(solvent,protein_light_atoms,protein_heavy_atoms,lengths,angles,num_asymm_units)
# %% 
# 7W6BnoKI
if __name__ == "__main__":
    lengths = 52.271,  422.928,   48.391
    angles = 90, 90, 90
    num_asymm_units=2

    protein_light_atoms = dict(
        C = 3635,
        N = 1022,
        O = 1216,
    )
    protein_heavy_atoms = dict(
        S = 4,
        Mg=1,
        Ca=1,
    )

    # https://advancedthermo.com/electrolytes/density_KI.html
    # But PEG lighter

    
    solvent = Solvent(59.48,solution_density=1.1) 
    #solvent.add_molecule(KI,M=1)
    solvent.add_molecule(PEG(3350),solvent_v_on_v=25)
    #TODO 100 mM HEPES
    
    calculate(solvent,protein_light_atoms,protein_heavy_atoms,lengths,angles,num_asymm_units)

# %%
if __name__ == "__main__":
    lengths = 29.53,  30,   30
    angles = 90, 90, 90
    num_asymm_units=1

    protein_light_atoms = dict(
    )
    protein_heavy_atoms = dict(
    )

    # https://advancedthermo.com/electrolytes/density_KI.html
    # But PEG lighter

    
    solvent = Solvent(100,solution_density=1.732) 
    #solvent.add_molecule(KI,M=1)
    solvent.add_molecule(CsCl,M=6)
    #TODO 100 mM HEPES
    
    calculate(solvent,protein_light_atoms,protein_heavy_atoms,lengths,angles,num_asymm_units)# %%
# %% Hemoglobin 3PEL
# REMARK 280 CRYSTALLIZATION CONDITIONS: 1.7 M AMSO4, 100 MM GLYCINE PH 9.0,      
# REMARK 280  0.6 M 3-(1-PYRIDINO)-1-PROPANE SULFONATE [NDSB-201] AND 21%         
# REMARK 280  GLYCEROL, VAPOR DIFFUSION, SITTING DROP, TEMPERATURE 298.0K, PH     
# REMARK 280  8.5  
if __name__ == "__main__":
    lengths = 87.973,   88.045,   53.073  
    angles = 90.00, 103.37,  90.00
    num_asymm_units=4

    protein_light_atoms = dict(
        C=1418,
        N=380,
        O=405,
        H=2230
    )
    protein_heavy_atoms = dict(
        S=5
    )

    # https://advancedthermo.com/electrolytes/density_KI.html
    # But PEG lighter

    
    solvent = Solvent(61.54,solution_density=1.1) 
    #solvent.add_molecule(KI,M=1)
    solvent.add_molecule()
    #TODO 100 mM HEPES
    
    calculate(solvent,protein_light_atoms,protein_heavy_atoms,lengths,angles,num_asymm_units)# %%



# %% Hemoglobin 2QSP (using same atoms as 3PEL)
# REMARK 280 CRYSTALLIZATION CONDITIONS: 0.4 M NA CACODYLATE 12-15% PEG 3350,     
# REMARK 280  PH 5.7, VAPOR DIFFUSION, HANGING DROP, TEMPERATURE 298K  
if __name__ == "__main__":
    lengths = 65.033,   78.273,   109.085  
    angles = 90.00, 90.00,  90.00
    num_asymm_units=4

    protein_light_atoms = dict(
        C=1418,
        N=380,
        O=405,
        H=2230
    )
    protein_heavy_atoms = dict(
        S=5
    )

    # https://advancedthermo.com/electrolytes/density_KI.html
    # But PEG lighter

    
    solvent = Solvent(44.97,solution_density=1.1) 
    #solvent.add_molecule(KI,M=1)
    solvent.add_molecule(PEG(3350),solvent_v_on_v=13.5)
    solvent.add_molecule(sodium_cacodylate,M=0.4)
    #TODO 100 mM HEPES
    
    calculate(solvent,protein_light_atoms,protein_heavy_atoms,lengths,angles,num_asymm_units)# %%


# %%
