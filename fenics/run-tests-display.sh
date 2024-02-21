source functions.sh

# Params
BS="10" # Backward steps
METHOD="inverse"

# Generate solutions for slab
NX=24
NY=8
NZ=8
LOADS="5 10 15 20 25"
for LOAD in ${LOADS}; do
    OUT="${OUTDIR}/display/slab_${LOAD}"
    ARGS="--output-file ${OUT}"
    runSlab $NX $NY $NZ $LOAD $METHOD $BS $OUT "$ARGS"
done


# Generate solutions for LV
MESH="h4_v2_ASCII"
ASS="1000 5000 10000 20000" # Active stresses
ENDOS="100 1000 2500 5000 10000" # Endocardial pressures

# First for active stress
for AS in ${ASS}; do
    ENDO="0.0"
    OUT="${OUTDIR}/display/LV_as_${AS}"
    ARGS="--output-file ${OUT}"
    runLV $MESH $ENDO $AS $METHOD $BS $OUT "$ARGS"
done

# Then for endocardial pressures
for ENDO in ${ENDOS}; do
    AS="0.0"
    OUT="${OUTDIR}/display/LV_endo_${ENDO}"
    ARGS="--output-file ${OUT}"
    if [ "${ENDO}" = "10000" ]; then
        BS=50
    fi
    runLV $MESH $ENDO $AS $METHOD $BS $OUT "$ARGS"
done

