# Fix hydrogen names for gromacs (pdb2gmx)
#%%
# file = "scripts/scattering/targets/4et8H_full_struct.pdb"
# out_file = "scripts/scattering/targets/4et8H_full_struct_Hfix.pdb"
file = "../scripts/scattering/targets/2qspH_unit.pdb"
out_file = "../scripts/scattering/targets/2qsp_unit_Hfix.pdb"

OPLS_FORMAT = False

# For CHARMM36 (2020)

with open(file) as f_in, open(out_file, 'w') as f_out:
    last_name = None
    fix_next=False
    lines = f_in.readlines()
    last_res_num=None
    res_start_idx=0
    new_lines=[]
    for i, line in enumerate(lines):
        if not line.startswith("ATOM") and not line[17:20].strip() in ["HEM","HEME"]:
            new_lines.append(line)
            continue
        # terminal res name fix
        # res_num=int(line[22:26].strip())
        # res_name = line[17:20].strip()
        # name = line[12:16].strip()
        # if last_res_num!=res_num:
        #     res_start_idx=i
        # if name == "OT2" and len(res_name)!=4:
        #     for j,line_to_change in enumerate(lines[res_start_idx:i+1]):
        #         lines[res_start_idx+j]=line_to_change[:17]+"C"+res_name+line_to_change[21:]

        # last_res_num=res_num
        res_num=int(line[22:26].strip())
        res_name = line[17:20].strip()
        name = line[12:16].strip()
        new_lines.append(line)
        if name == "OT2" and not lines[i+1].startswith("TER"):
                new_lines.append("TER\n")
    lines=new_lines
    for i, line in enumerate(lines):
        if not line.startswith("ATOM") and not line[17:20].strip() in ["HEM","HEME"]:
            f_out.write(line)
            continue




        name = line[12:16].strip()
        res_name = line[17:21].strip()


        #new_convert_dict=dict(HT1="H ",HT2="HA",HT3="HB")
        new_convert_dict=dict(HT1="H ",HT2="HA",HT3="HB",OT1="O1",OT2="O2",CT1="C1",CT2="C2",CT3="C3")
        if name in new_convert_dict:
            #line = line[:12]+ ' ' + name[0]+name[-1] + ' '  + line[16:]  
            line = line[:12]+ ' ' + new_convert_dict[name] + ' '  + line[16:]  

        

        if name == "HG" and res_name in ["SER","CYS"]:
            line = line[:12]+ ' ' + "HG1"  + line[16:]  
        his_to_hisd_convert_dict=dict(HD2="HE1",HE1="HD1",HE2="HD2")
        if name == "HE2" and res_name == "HISD":
            line = line[:12]+ ' ' + "HD2"  + line[16:]  



        HEME_convert_dict=dict(HHA="HA",HHB="HB",HHC="HC",HHD="HD")
        if res_name in ["HEM","HEME"] and name in HEME_convert_dict:  # ????? why do they do this?
            line = line[:12]+ ' ' + HEME_convert_dict[name] + ' '  + line[16:]  


        ####
        #ignore = not OPLS_FORMAT

        ignore = False
        if res_name in "TRP" and name[:2]=="HZ" or name == "HH2" or res_name in ["HEM","HEME"]:  # ????? why do they do this?
            ignore = True


        modified_line = line
        if not ignore:
            if fix_next and name[0] == "H" and name[-1] == "3":
                fix_next=False
                modified_line = line[:12]+ ' '*(4-len(name)) + name[:-1]+"2" + line[16:]    
            elif name[0] == "H" and name[-1] == "2" and last_name != name[:-1]+"1":
                modified_line = line[:12]+ ' '*(4-len(name)) + name[:-1]+"1"  + line[16:]    
                fix_next = True

        f_out.write(modified_line)
        last_name = name

# %%
