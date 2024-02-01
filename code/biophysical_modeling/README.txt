Thomas C. Day
2024 January

There are 4 scripts in this folder that are meant to run and analyze simulations where snowflake yeast of varying cell size, shape, and bond strength, are grown and measured.

kai_sim_sweep_MaxStress.m sweeps a range of overlap thresholds (a proxy for cell bond strength) and simulates N groups at each overlap threshold.

Results_overlap_sweep.m compares the simulated cluster sizes from kai_sim_sweep_MaxStress to measured cluster size values to help choose a relevant threshold.

kai_phashemap_snowflake_gen.m is the main simulation script, and it contains many subfunctions. This script simulates N clusters of a certain size and shape, given a particular overlap thresh defined from the scripts above.

kai_phasemap_gen.m collects the simulations generated from kai_phasemap_snowflake_gen.m and measures/plots cluster sizes.