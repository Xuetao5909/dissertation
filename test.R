source("load_reservoirs.R")
source("network_generators.R")
source("method.R")
source("evaluation.R")
source("run_stress_test.R")

res_obj <- load_reservoirs("coloc-simulation-results.csv.gz")

res <- run_stress_test(
  n_reps = 20,
  scenarios = c("triangle", "bowtie", "ld_trap", "local_bridge_trap"),
  reservoir_shared = res_obj$reservoir_shared,
  reservoir_shared_weak = res_obj$reservoir_shared_weak,
  reservoir_shared_moderate = res_obj$reservoir_shared_moderate,
  reservoir_shared_strong = res_obj$reservoir_shared_strong,
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
names(res)
library(dplyr)
library(tidyr)
library(ggplot2)

ari_long <- res %>%
  filter(scenario %in% c("ld_trap", "local_bridge_trap")) %>%
  select(
    scenario,
    replicate,
    louvain_ari,
    cpm_ari,
    spinglass_ari
  ) %>%
  pivot_longer(
    cols = c(louvain_ari, cpm_ari, spinglass_ari),
    names_to = "method",
    values_to = "ARI"
  ) %>%
  mutate(
    method = recode(
      method,
      louvain_ari = "Louvain",
      cpm_ari = "CPM",
      spinglass_ari = "Signed Spinglass"
    ),
    scenario = recode(
      scenario,
      ld_trap = "LD trap",
      local_bridge_trap = "Local bridge trap"
    )
  )

p_ari_box <- ggplot(ari_long, aes(x = method, y = ARI)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +
  geom_jitter(width = 0.12, alpha = 0.6, size = 1.8) +
  facet_wrap(~ scenario) +
  labs(
    title = "Adjusted Rand Index across stress-test replicates",
    x = "Method",
    y = "Adjusted Rand Index"
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 25, hjust = 1),
    plot.title = element_text(face = "bold")
  )

print(p_ari_box)
