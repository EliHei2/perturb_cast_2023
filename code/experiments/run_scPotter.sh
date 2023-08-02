# bsub -q production-rh74 -P gpu -gpu  - python3 scPotter/script_py.py 

for seed in 0 1 2 3 4 5 6 7 8 9; do
    for tag_tarin in embryo1_2; do
        for tag_test in embryo1_2 embryo1_5 embryo2_2 embryo2_5 embryo3_2 embryo3_5; do
            bsub -q production-rh74 -M 8000 -R rusage[mem=8000] \
            -P gpu -gpu - "python3 scPotter/run_scPotter.py \
            -i ../../data_tidy -tag seqfish_mouse_embryo \
            -tag_train $tag_tarin \
            -tag_test  $tag_test \
            --oo ../../output \
            --seed $seed" 
        done
    done
done