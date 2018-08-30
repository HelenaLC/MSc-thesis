# MSc thesis on "Differential Analysis of scRNA-seq data with complex experimental designs"

### contents

- **scripts:** R code to reproduce all analyses & figures.
- **results:** Figures produced by scripts.

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
- **scDD_sim_ex.R**: 
<br/>Visualises an exemplary `scDD` simulation <br/>
  (reproduces: *scDD_sim_ex-med_exprs*, *scDD_sim_ex-expr_profiles*)
- **scDD_null_sim.R**: <br/>
  Evaluates method performances on 3 replicates of a null simulation <br/>
  (reproduces: *scDD_null_sim*)
- **diffcyt_runmodes.R**: <br/>
  Evaluates the performance of `diffcyt` for varying data inputs & summary statistics <br/>
  (reproduces: *diffcyt_runmodes*)

