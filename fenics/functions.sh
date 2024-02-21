SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

# Create output folders
OUTDIR="output" #Library uses relative dir and prepends output
mkdir -p $OUTDIR
mkdir -p ${OUTDIR}/display
mkdir -p ${OUTDIR}/optimality/slab
mkdir -p ${OUTDIR}/optimality/LV
mkdir -p ${OUTDIR}/robustness/slab
mkdir -p ${OUTDIR}/robustness/LV
mkdir -p ${OUTDIR}/comparison/slab
mkdir -p ${OUTDIR}/comparison/LV

# Solver methods create a .done file if the test has finished successfully.
# Remove them to rerun a test

# Run Slab solver
runSlab (){
    NX=$1
    NY=$2
    NZ=$3
    LOAD=$4
    TYPE=$5
    BS=$6
    OUT=$7
    ARGS="$8"
    if [ ! -e "${OUT}.done" ]; then 
        printf "" > $OUT
        echo "Running main_prestress_slab.py --nx $NX --ny $NY --nz $NZ --z-load-scale $LOAD --solver-type $TYPE --backward-steps $BS $ARGS | tee $OUT"
        tsp sh -c "python main_prestress_slab.py --nx $NX --ny $NY --nz $NZ --z-load-scale $LOAD --solver-type $TYPE --backward-steps $BS $ARGS | tee $OUT; touch ${OUT}.done"
    fi
}

# Run left ventricle solver
runLV (){
    MESH=$1
    ENDO=$2
    AS=$3
    TYPE=$4
    BS=$5
    OUT=$6
    ARGS="$7"
    if [ ! -e "${OUT}.done" ]; then 
        printf "" > $OUT
        echo "Running main_prestress_LV.py --heart-geometry-name $MESH --z-load-scale 0 --solver-type $TYPE --backward-steps $BS --endo-pressure $ENDO --active-peak $AS $ARGS | tee $OUT"
        tsp sh -c "python main_prestress_LV.py --heart-geometry-name $MESH --z-load-scale 0 --solver-type $TYPE --backward-steps $BS --endo-pressure $ENDO --active-peak $AS $ARGS | tee $OUT; touch ${OUT}.done"
    fi
}

