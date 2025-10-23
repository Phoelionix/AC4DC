# Fix hydrogen names for gromacs (pdb2gmx)

file = "4et8H_full_struct.pdb"
out_file = "4et8H_full_struct_Hfix.pdb"

with open(file) as f_in, open(out_file, 'w') as f_out:
    last_name = None
    fix_next=False
    for line in f_in:
        if not line.startswith("ATOM"):
            f_out.write(line)
            continue
        name = line[12:16].strip()
        res_name = line[17:20]

        ignore = False
        if res_name == "TRP" and name[:2]=="HZ" or name == "HH2":  # ????? why do they do this?
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
