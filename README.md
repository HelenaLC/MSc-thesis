# MSc thesis on "Differential Analysis of scRNA-seq data with complex experimental designs"

### contents

- **scripts:** R code to reproduce all analyses & figures.
- **results:** Figures & data produced by scripts.

### packages

- **pkg:** R package containing method wrappers and utilities for plotting & evaluation.
- **fda:** `FDA` package fork w/ modified code of `tperm.fd()` to decrease runtime.
- **scDD:** `scDD` package fork w/ modified version of `simulateSet()` to prevent repeated running of `findIndex()` & make simulated counts non-continuous.

### scripts

- **dd_patterns.R**: <br/>
    Generates a schematic of differential distribution patterns <br/>
    (reproduces: *dd_patterns*)
    
- **ECDFs.R**: <br/>
    Generates an exemplary set of ECDFs for a 3 vs. 3 sample comparison. <br/>
    (reproduces: *ECDFs*)
    
- **scDD-sim_ex.R**: <br/>
    Visualises an exemplary `scDD` simulation <br/>
    (reproduces: *scDD_sim_ex-med_exprs*, *scDD_sim_ex-expr_profiles*)
  
- **simDD-seurat.R** <br/>
    Evaluates `Seurat` clustering performance on 
    10 simulation replicates w/ randomised parameters <br/>
    (reproduces: *simDD-seurat_scores*)
 
- **simDD-sim_qc.R**: <br/>
    Generates basic quality control plots for `simDD` simulated data
    and data from Koh et al. <br/>
    (reproduces: *qc_var_explained*, *qc_lib_sizes*, 
        *qc_top_expr*, *qc_expr_freq_vs_mean*, *qc_disp_vs_mean*)
        
- **scDD-null_sim.R**: <br/>
    Evaluates method performances on 3 replicates of a null simulation <br/>
    (reproduces: *scDD_null_sim*)
  
- **diffcyt_runmodes.R**: <br/>
    Evaluates the performance of `diffcyt` for varying data inputs & summary statistics <br/>
    (reproduces: *diffcyt_runmodes*)

- **kang-data_prep.R**: <br/>
    Performed `Seurat` preprocessing & constructs a `daFrame` from
    the Kang et al. raw data available at accession # GSE96583
  
- **kang-data_overview.R**: <br/>
    Generates general data overview plot for the Kang et al. data set <br/>
    (reproduces: *kang_cluster_props*, *kang_tsne*, *kang_nb_cells*)

- **kang-DS_analysis.R**: <br/>
    Performs differential analysis using `diffcyt` & `edgeR` methods
    & compares obtained results with thouse published
    (reproduces: *kang_nb_de_gs*, *kang_overlap*, *kang_pvals*, 
        *kang_top_undetected*, *kang_highest_pvals_xxx*)

- **runtimes-nb_gs.R**: <br/>
    Measures method runtimes for increasing numbers of genes
    (reproduces: *runtimes*)
    
- **runtimes-FDA_reso.R**, **runtimes-FDA_nperm.R**: <br/>
    Measure `FDA` runtimes for increasing `reso` and `n_perm` parameters
    (reproduces: *runtimes_FDA_reso/nperm*)
