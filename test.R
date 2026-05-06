source("load_reservoirs.R")
source("network_generators.R")
source("method.R")
source("evaluation.R")
source("run_stress_test.R")

res_obj <- load_reservoirs("coloc-simulation-results.csv.gz")

res <- run_stress_test(
  n_reps = 10,
  scenarios = c("triangle", "bowtie", "ld_trap", "local_bridge_trap"),
  reservoir_shared = res_obj$reservoir_shared,
  reservoir_shared_weak = res_obj$reservoir_shared_weak,
  reservoir_distinct_bg = res_obj$reservoir_distinct_bg,
  reservoir_distinct_trap = res_obj$reservoir_distinct_trap,
  reservoir_ambiguous = res_obj$reservoir_ambiguous,
  cpm_tau_grid = seq(0.1, 0.7, by = 0.05)
)

summary_tab <- summarise_results(res)

print(summary_tab)

if (!dir.exists("output")) dir.create("output", recursive = TRUE)

readr::write_csv(res, "output/stress_test_results.csv")
readr::write_csv(summary_tab, "output/stress_test_summary.csv")