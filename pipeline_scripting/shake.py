import random
import sys


if len(sys.argv)!=4:
    print("Usage: python shake.py [input] [output] [shake_amount]")
input_file = sys.argv[1]
output_file = sys.argv[2]
shake = float(sys.argv[3])  

with open(input_file) as f:
    lines = f.readlines()

header=lines[:2]
coords = lines[2:-1]
box = lines[-1]

modified_lines = []
for line in coords:
    x = float(line[20:28]) + random.uniform(-shake, shake)
    y = float(line[28:36]) + random.uniform(-shake, shake)
    z = float(line[36:44]) + random.uniform(-shake, shake)

    modified_lines.append(f"{line[:20]}{x:8.3f}{y:8.3f}{z:8.3f}{line[44:]}")

with open(output_file, "w") as f:
    f.writelines(header)
    f.writelines(modified_lines)
    f.write(box)