# MSc thesis on "Differential Analysis of scRNA-seq data with complex experimental designs"

### contents

- **fda:** forked repository of the `FDA` package w/ modified code of `tperm.fd()` to decrease runtime.
- **scDD:** forked repository of the `scDD` package w/ modified version of `simulateSet()` to prevent repeated running of `findIndex()`.
- **scripts:** R code to reproduce all analyses & figures.
- **results:** figures produced by scripts.

### scripts

- **dd_patterns.R**: <br/>
  Generates a schematic of differential distribution patterns \linebreak
  (reproduces: *dd_paterns*)
- **scDD_sim_ex**: Visualises an exemplary `scDD` simulation (**scDD_sim_ex-med_exprs**, **scDD_sim_ex-expr_profiles**)
- **scDD_null_sim**: Evaluates method performances on 3 replicates of a null simulation (**scDD_null_sim**)
- **diffcyt_runmodes**: Evaluates the performance of `diffcyt` for varying data inputs & summary statistics (**diffcyt_runmodes**)

