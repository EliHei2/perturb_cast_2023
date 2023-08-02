for seed in 0 1 2 3 4 5 6 7 8 9; do
    for tag in emgryo1_5 emrbryo2_2 embryo2_5 embryo3_2 embryo3_5; do
        bsub -W 1:00 -R "rusage[mem=6000]"
        "python run_novosparc.py \
        --data_dir data_tidy/seqfish_mouse_embryo
        --out_dir output/R
        --tag tag \
        --atlas_locs_f cell_meta_embryo1_2.txt \
        --atlas_mtx_f data_embryo1_2.txt \
        --expr_mtx_f data_$tag.txt \  
        --nns 5 \
        --nnt 5 \
        --alpha 0.5 \
        --seed $seed"
    done
done