#!/bin/bash
BooleanFunction=("Mux" "Maj");
ARRAY_AGING=(50 300 800 1500)
ARRAY_dr=(0.01 0.02 0.04 0.06 0.08 0.10 0.30 0.50 0.70)
ARRAY_mrr=(0.01)
ARRAY_reps=(001 002 003 004 005 006 007 008 009 010)

path=$(pwd)

for((b=0;b<2;b++));
do
    for((r=0;r<10;r++));
    do
        for ((a=0;a<4;a++));
        do
            for ((d=0;d<9;d++));
            do
                for ((m=0;m<1;m++));
                do
                    EXPE_PATH=$path/${ARRAY_reps[$r]}_${BooleanFunction[$b]}_dr${ARRAY_dr[$d]}_mrr${ARRAY_mrr[$m]}_a${ARRAY_AGING[$a]}
                    if [ ! -d ${EXPE_PATH} ]; then
                        mkdir ${EXPE_PATH}
                        mkdir ${EXPE_PATH}/NetworksData/
                        mkdir ${EXPE_PATH}/DOT/
                        mkdir ${EXPE_PATH}/ParetoNonPareto/
                    else
                        echo "${EXPE_PATH} already exists."
                    fi
                done;
            done;
        done;
    done;
done;
