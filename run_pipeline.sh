set -u 


id_start=1
id_end=9
NUM_MD=15

for ((i=0; i<NUM_MD; i++)) {
    for ((runid=id_start; runid<=id_end; runid+=1)) {
        bash hemoglobin_MD_pipeline.sh $runid 1
    }
}