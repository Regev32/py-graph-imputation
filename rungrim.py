from grim.RunGrim import run_original_grim

conf_file = "conf/minimal-configuration.json"
# processes=0 imputes one subject per core. The workers share the graph this
# process builds, so the core count is the limit, not the memory.
run_original_grim(conf_file, True, True, True, processes=0)


# TODO: Check if don.problem gets only short gls or if it's really 2 and 3 loci inputs
