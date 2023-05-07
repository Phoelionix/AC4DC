set -x
#!/usr/bin/env

for f in input/tetra/tetra-{15..27}.mol; do ./ac4dc "$f"; done
