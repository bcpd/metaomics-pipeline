# Debugging error in atlas_genome rule

Ocassionally, an error in the atlas_genome rule can appear and will stop the process. The error (see below) appears when atlas create a rooted tree and use it to classify the genomes using the GTDB scheme. It is unpredictable when it happens so it is hard to debug. We hope the developers of Atlas will fix it in the future. Nevertheless, there is a way to bypass this error (see below).


*Error:*
```
Error in rule atlas_genomes:
    jobid: 5
    output: ./finished_genomes
    log: ./logs/atlas_genomes.log (check log file(s) for error message)
    conda-env: atlas
    shell:
        sed -i "s/^- dram/# - dram/" ./config.yaml;sed -i "s/^- kegg_modules/# - kegg_modules/" ./config.yaml;(atlas run genomes -j 10 -w ./ --max-mem 128 --omit-from run_dram --config-file ./config.yaml) 2> ./logs/atlas_genomes.log;
        (one of the commands exited with non-zero exit code; note that snakemake uses bash strict mode!)

Shutting down, this might take some time.
Exiting because a job execution failed. Look above for error message
```

The error has been previously reported [here](https://github.com/metagenome-atlas/atlas/issues/530) and it is bypassed by copying the unrooted tree and renaming the copy. The file can be found in the working directory in the *genomes* folder inside a folder called *tree*

*Solution*:
```
cd genomes/tree
cp gtdbtk.bac120.unrooted.nwk  gtdbtk.bac120.nwk
cd ../../

# Resume the pipeline
M5P -c M5p_config.yaml
```
