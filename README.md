# MSc thesis on "Differential Analysis of scRNA-seq data with complex experimental designs"

### contents

- **fda:** forked repository of the `FDA` package w/ modified code of `tperm.fd()` to decrease runtime.
- **scDD:** forked repository of the `scDD` package w/ modified version of `simulateSet()` to prevent repeated running of `findIndex()`.
- **scripts:** R code to reproduce all analyses & figures.
- **results:** figures produced by scripts.

### scripts

- **dd_patterns.R**: <br/>
  Generates a schematic of differential distribution patterns <br/>
  (reproduces: *dd_patterns*)
- **scDD_sim_ex**: 
<br/>Visualises an exemplary `scDD` simulation <br/>
  reproduces: *scDD_sim_ex-med_exprs*, *scDD_sim_ex-expr_profiles*)
- **scDD_null_sim**: <br/>
  Evaluates method performances on 3 replicates of a null simulation <br/>
    *scDD_null_sim*
- **diffcyt_runmodes**: <br/>
  Evaluates the performance of `diffcyt` for varying data inputs & summary statistics <br/>
  *diffcyt_runmodes*

