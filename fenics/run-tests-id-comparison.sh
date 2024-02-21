source functions.sh

AAM=3 # Anderson acceleration parameter

#### SLAB
NX=24
NY=8
NZ=8
BS="10"
LOADS="1 10 25 50 75 100"
for LOAD in $LOADS; do
    # inverse
    OUT="${OUTDIR}/comparison/slab/id_${LOAD}"
    runSlab $NX $NY $NZ $LOAD inverse $BS $OUT
    # forward after inverse
    OUT="${OUTDIR}/comparison/slab/id-forw_${LOAD}"
    runSlab $NX $NY $NZ $LOAD inverse $BS $OUT "--forward-steps ${BS}"
    # forward
    OUT="${OUTDIR}/comparison/slab/forw_${LOAD}"
    runSlab $NX $NY $NZ $LOAD inverse $BS $OUT "--backward-steps 0 --forward-steps ${BS}"
done

#### LV
MESH="h4_v2_ASCII"
BS="10"
ASS="0  1000 5000 10000 20000" # Active stresses, only those that converge
ENDOS="0 100 1000 2500 5000" # Endo pressures, only those that converge
for AS in $ASS; do
for ENDO in $ENDOS; do
    if [ "${AS}" = "0" -o "${ENDO}" = "0" ]; then # test separately
        # inverse
        OUT="${OUTDIR}/comparison/LV/id_${ENDO}_${AS}_${BS}"
        runLV $MESH $ENDO $AS inverse $BS $OUT

        # inverse then forward
        OUT="${OUTDIR}/comparison/LV/id-forw_${ENDO}_${AS}_${BS}"
        runLV $MESH $ENDO $AS inverse $BS $OUT "--forward-steps ${BS}"

        # only forward
        OUT="${OUTDIR}/comparison/LV/forw_${ENDO}_${AS}_${BS}"
        runLV $MESH $ENDO $AS inverse $BS $OUT "--backward-steps 0 --forward-steps ${BS}"
    fi
done
done
