from grim.RunGrim import run_original_grim

conf_file = "conf/minimal-configuration.json"
run_original_grim(conf_file, True, True, True)


# TODO: Check if don.problem gets only short gls or if it's really 2 and 3 loci inputs
# TODO: multithread for 9loci. There's already a skeleton code in scripts/parallel-imputation.py
