import os.path as path
from glob import glob


input_files = glob(path.join(config["root_dir"], config["samples_csv"], '*.csv'))


localrules: all


rule all:
    input:
        [path.join(config["root_dir"], config["samples_adjusted_swaths"],
                          path.split(path.splitext(fname)[0] + '_adjusted.csv')[1])
         for fname in input_files]


rule adjust_file:
    input:
        path.join(config["root_dir"], config["samples_csv"], "{sample}.csv")
    output:
        path.join(config["root_dir"], config["samples_adjusted_swaths"], "{sample}_adjusted.csv")
    run:
        for in_csv, out_csv in zip(input, output):
            shell("Rscript scripts/util/adjust_swaths.R -i {in_csv} -o {out_csv}")

