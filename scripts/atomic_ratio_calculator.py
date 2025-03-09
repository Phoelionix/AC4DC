#%%
import numpy as np
from diffpy.structure.spacegroups import GetSpaceGroup

# no hydrogen
#space_group_pdb_name = "C 1 2 1" 

class Solvent:
    molecule_idx_iterator=0
    def __init__(self,VS,VM=None):
        self.molecules = []
        self.VS = VS # Solvent content %(v/v)
        self.VM = VM # Matthews Coefficient, Ang^3/Da  (note 1 Da -> 1 g/mol)
    def add_molecule(self, molecule,M = None,solvent_v_on_v = None):
        molecule.set_concentration(M,solvent_v_on_v)
        self.molecules.append(molecule)
    def clear_molecules(self):
        self.molecules = []
    

class Molecule:
    def __init__(self,molar_mass=None,density=None,undiluted_molarity=None,**atom_dict):
        self.idx = Solvent.molecule_idx_iterator
        Solvent.molecule_idx_iterator+=1
        self.undiluted_molarity = undiluted_molarity
        self.molar_mass=molar_mass
        self.density=density
        self.atom_dict = atom_dict

        # Properties of solvent
        self.fraction_of_solvent_volume=None
        self.M = None
    def set_concentration(self,M = None,solvent_v_on_v = None):
        have_measure = False
        for elem in solvent_v_on_v,M:
            if elem is not None:
                assert not have_measure, "Given multiple variables for concentration, should only pass one!"
                have_measure = True
        assert have_measure, "Missing concentration information"

        if M is not None:
            assert(self.molar_mass is not None and self.density is not None)
            self.M = M
            self.fraction_of_solvent_volume = M*self.molar_mass/(self.density*1e3)
        elif solvent_v_on_v is not None:
            self.fraction_of_solvent_volume = solvent_v_on_v/100
            if self.density is not None and self.molar_mass is not None:
                self.M = self.fraction_of_solvent_volume*self.density*1e3/self.molar_mass
                print(self.density*1e3/self.molar_mass)
            elif self.undiluted_molarity is not None:
                self.M = self.undiluted_molarity*self.fraction_of_solvent_volume
            else:
                print(f"Can't determine molarity of molecule in solution") 

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
    print(asym_atoms)

    # Num molecules in asymmetric unit
    def get_num_asymm_molecules(volume_fraction,M):
        return  (V_cell/num_asymm_units * M * volume_fraction
            * 1e-27*6.02214076e23)

    for molecule in solvent.molecules:
        N = get_num_asymm_molecules(solvent.VS/100,molecule.M)
        print(N)
        for element, num in molecule.atom_dict.items():
            if element in asym_atoms:
                asym_atoms[element] += N*num
            else:
                asym_atoms[element] = N*num



    print("-----------")
    for k,v in asym_atoms.items():
        print(k,f"{v:.3f}")
    print(f"\nV: {V_cell/num_asymm_units:.1f}")

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

water = Molecule(undiluted_molarity=55.56,
    H = 2,
    O = 1
)

NaCl = Molecule(58.44,2.16,
    Na = 1,
    Cl = 1                
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
