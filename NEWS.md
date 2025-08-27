# genepicoloc 0.9.3

## Major changes
* Restructuring of pipeline architecture for improved scalability
  - Introduced three-level hierarchy: `genepicoloc_wrapper()` → `genepicoloc_job()` → `genepicoloc_run()`
  - Implemented two-level parallelization strategy (regions × datasets)
  - Added job-based processing with configurable `max_regions_per_job` parameter
  - Improved memory efficiency: each job needs <1Gb when setting <100 regions (flexible)

## New features
* Added tar archive storage for efficient file management
  - Summary statistics saved to `sumstats.tar`
  - Colocalization results saved to `coloc.tar`
* Added `save_sumstats` parameter (default: FALSE) to save processed summary statistcs
* Added `p_min_save` parameter for saving only significant summary statistics
* Added `p_filt` parameter for additional p-value filtering
* Added `batch_size` parameter for controlling Level 2 parallelization
* New `filter_significant_regions()` function to reduce storage requirements
* New `consolidate_coloc_results()` function for study-wise result aggregation
* Enhanced gene annotation with configurable `nearest` parameter
* Updated documentation for pipeline functions
* Removed dependency on `writexl` package for core functionality

# genepicoloc 0.9.2

## Documentation improvements
* Removed progress bar with processx
* Added progress report using base R
* Added name_by_position function

## Bug fixes
* Added proper call of data.table::copy

# genepicoloc 0.9.1

## Documentation improvements
* Enhanced documentation for all exported functions
* Added comprehensive parameter descriptions
* Added examples where appropriate
* Marked internal functions with `@keywords internal`

## Bug fixes
* Removed debug print statement in `get_coloc_regions()`
* Fixed duplicate function definitions in query_sumstats_custom.R

## Minor improvements
* Cleaned up code comments
* Improved function descriptions
* Added proper S3 method documentation

# genepicoloc 0.9.0

* Initial release

