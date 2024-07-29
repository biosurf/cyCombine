# cyCombine 0.2.20

* Batch correction can now be performed iteratively by submitting xdim and ydim as vectors. E.g., `xdim = c(4,8)` and `ydim = c(1,8)`. This will first correct using a 4x1 grid and then an 8x8 grid.
* You can now adjust the bin size for EMD calculations.
* Added \ to removed characters in clean_colnames, but - and + are no longer removed.
* Added markers as input in detect_batch_effects_express
* Added compensation in prepare_data
* Experimental Seurat support
* Added mode selection in SOM clustering
* Added binSize argument to detect_batch_effect* functions
* batch_correct can now be parallelized with `mc.cores` parameter using `pbmcapply`
* Data normalization prior to clustering is now also parallelizable
* Added option for alternative clustering methods, such as `FlowSOM` and `kmeans`
* Vectorize and parallelize EMD computation

# cyCombine 0.2.19

* Revised data export functions
* Experimental support of ComBat_seq

# cyCombine 0.2.16

* Added description for installing with conda. Closes #37
* Added plot_density functionality to create X plots per file. Closes #34
* Fixed wrong function names in README 

# cyCombine < 0.2.16

* Check commit messages for changelog
