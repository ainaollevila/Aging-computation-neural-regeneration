#!/bin/bash
BooleanFunction=(0 1); # 0=Mux; 1=Maj;
ARRAY_AGING=(50 300 800 1500)
ARRAY_dr=(0.01 0.02 0.04 0.06 0.08 0.10 0.30 0.50 0.70)
ARRAY_mrr=(0.01)
ARRAY_reps=(1 2 3 4 5 6 7 8 9 10)

for((b=0;b<2;b++));
do
    for((r=0;r<10;r++));
    do
        for ((a=0;a<4;a++));
        do
            for ((d=0;d<9;d++));
            do
                dr="$(printf '%.2f\n' "${ARRAY_dr[$d]}")"
                for ((m=0;m<1;m++));
                do
                    mrr="$(printf '%.2f\n' "${ARRAY_mrr[$m]}")"
                    ./exe ${mrr} ${dr} ${ARRAY_AGING[$a]} ${BooleanFunction[$b]} ${ARRAY_reps[$r]}
                done;
            done;
        done;
    done;
done;
