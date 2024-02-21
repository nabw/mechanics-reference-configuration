source functions.sh 

AAM=3

#### Slab
BS="1"
LOADS="10 20"
NS="2 4 8 16 24 32 40"
for LOAD in $LOADS; do
    for N in $NS; do
        NY=$N
        NZ=$N
        NX=$((3 * ${N}))
        # inverse
        OUT="${OUTDIR}/optimality/slab/inverse_${LOAD}_${BS}_${N}"
        runSlab $NX $NY $NZ $LOAD inverse $BS $OUT

        # sellier
        OUT="${OUTDIR}/optimality/slab/sellier_${LOAD}_${BS}_${N}"
        runSlab $NX $NY $NZ $LOAD sellier $BS $OUT

        # sellier armijo
        OUT="${OUTDIR}/optimality/slab/sellier-armijo_${LOAD}_${BS}_${N}"
        runSlab $NX $NY $NZ $LOAD sellier $BS $OUT "--sellier-armijo"

        # sellier anderson 5
        OUT="${OUTDIR}/optimality/slab/sellier-aa_${LOAD}_${BS}_${N}"
        runSlab $NX $NY $NZ $LOAD sellier $BS $OUT "--sellier-aa-m ${AAM}"
    done
done

#### LV
MESHES="4mm h4_v2_ASCII 4mm_ref h4_v2_ASCII_ref"
BS=10
ASMAX="5000"
ENDOMAX="200"
for MESH in ${MESHES}; do
    for AS in 0 ${ASMAX}; do
    for ENDO in 0 ${ENDOMAX}; do
        if [ "$AS" = "0" -a "$ENDO" = "$ENDOMAX" -o "$AS" = "$ASMAX" -a "$ENDO" = "0" ]; then
            # inverse
            OUT="${OUTDIR}/optimality/LV/inverse_${ENDO}_${AS}_${MESH}"
            runLV $MESH $ENDO $AS inverse $BS $OUT

            # sellier
            OUT="${OUTDIR}/optimality/LV/sellier_${ENDO}_${AS}_${MESH}"
            runLV $MESH $ENDO $AS sellier $BS $OUT

            # sellier armijo
            OUT="${OUTDIR}/optimality/LV/sellier-armijo_${ENDO}_${AS}_${MESH}"
            runLV $MESH $ENDO $AS sellier $BS $OUT "--sellier-armijo"

            # sellier anderson 5
            OUT="${OUTDIR}/optimality/LV/sellier-aa_${ENDO}_${AS}_${MESH}"
            runLV $MESH $ENDO $AS sellier $BS $OUT "--sellier-aa-m ${AAM}"
        fi
    done
    done
done
