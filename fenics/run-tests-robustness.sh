source functions.sh

AAM=3 # Anderson acceleration parameter

#### SLAB
NX=24
NY=8
NZ=8
BSS="1 10 100"
LOADS="1 10 25 50 75 100"
for BS in $BSS; do
    for LOAD in $LOADS; do
        # inverse
        OUT="${OUTDIR}/robustness/slab/inverse_${LOAD}_${BS}"
        runSlab $NX $NY $NZ $LOAD inverse $BS $OUT

        # sellier
        OUT="${OUTDIR}/robustness/slab/sellier_${LOAD}_${BS}"
        runSlab $NX $NY $NZ $LOAD sellier $BS $OUT

        # sellier armijo
        OUT="${OUTDIR}/robustness/slab/sellier-armijo_${LOAD}_${BS}"
        runSlab $NX $NY $NZ $LOAD sellier $BS $OUT "--sellier-armijo"

        # sellier anderson 5
        OUT="${OUTDIR}/robustness/slab/sellier-aa_${LOAD}_${BS}"
        runSlab $NX $NY $NZ $LOAD sellier $BS $OUT "--sellier-aa-m ${AAM}"
    done
done

#### LV
MESH="h4_v2_ASCII"
BSS="1 10 100"
ASS="0  1000 5000 10000 20000 40000" # Active stresses
ENDOS="0 100 1000 2500 5000 10000" # Endo pressures
for BS in $BSS; do
    for AS in   $ASS; do
    for ENDO in $ENDOS; do
        if [ "${AS}" = "0" -o "${ENDO}" = "0" ]; then # test separately
            # inverse
            OUT="${OUTDIR}/robustness/LV/inverse_${ENDO}_${AS}_${BS}"
            runLV $MESH $ENDO $AS inverse $BS $OUT

            # sellier
            OUT="${OUTDIR}/robustness/LV/sellier_${ENDO}_${AS}_${BS}"
            runLV $MESH $ENDO $AS sellier $BS $OUT

            # sellier armijo
            OUT="${OUTDIR}/robustness/LV/sellier-armijo_${ENDO}_${AS}_${BS}"
            runLV $MESH $ENDO $AS sellier $BS $OUT "--sellier-armijo"

            # sellier anderson 5
            OUT="${OUTDIR}/robustness/LV/sellier-aa_${ENDO}_${AS}_${BS}"
            runLV $MESH $ENDO $AS sellier $BS $OUT "--sellier-aa-m ${AAM}"
        fi
    done
    done
done
