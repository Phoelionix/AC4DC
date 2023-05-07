set -x
#!/usr/bin/env

for f in tetrapeptide/tetra-{1..16}.mol; do cp tetrapeptide/tetra_template.mol "$f"; done
# cp tetrapeptide/tetra_template.mol tetrapeptide/test.mol
