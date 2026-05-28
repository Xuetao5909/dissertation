source("load_reservoirs.R")
source("network_generators.R")
source("method.R")
source("evaluation.R")
source("run_stress_test.R")
source("run_scaling_test.R")

res_obj <- load_reservoirs("coloc-simulation-results.csv.gz")

res <- run_stress_test(
  n_reps = 100,
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
# Runtime summary table

# -------------------------
runtime_total_by_scenario <- res %>%
  group_by(scenario) %>%
  summarise(
    louvain_total_time_sec = sum(louvain_time_sec, na.rm = TRUE),
    cpm_total_time_sec = sum(cpm_time_sec, na.rm = TRUE),
    spinglass_total_time_sec = sum(spinglass_time_sec, na.rm = TRUE),
    .groups = "drop"
  )

print(runtime_total_by_scenario)

if (!dir.exists("output")) dir.create("output", recursive = TRUE)

readr::write_csv(res, "output/stress_test_results.csv")
readr::write_csv(summary_tab, "output/stress_test_summary.csv")
readr::write_csv(runtime_tab, "output/runtime_total_by_scenario.csv")
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
runtime_total_long <- res %>%
  group_by(scenario) %>%
  summarise(
    Louvain = sum(louvain_time_sec, na.rm = TRUE),
    CPM = sum(cpm_time_sec, na.rm = TRUE),
    `Signed Spinglass` = sum(spinglass_time_sec, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  tidyr::pivot_longer(
    cols = c(Louvain, CPM, `Signed Spinglass`),
    names_to = "method",
    values_to = "total_time_sec"
  )

print(runtime_total_long)
p_runtime_total <- ggplot(runtime_total_long, aes(x = method, y = total_time_sec)) +
  geom_col(width = 0.65) +
  facet_wrap(~ scenario) +
  labs(
    title = "Total runtime by method and scenario",
    x = "Method",
    y = "Total runtime seconds"
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 25, hjust = 1),
    plot.title = element_text(face = "bold")
  )

print(p_runtime_total)





scaling_res <- run_scaling_test(
  n_reps = 5,
  node_sizes = c(10, 20, 30, 50, 100),
  reservoir_shared_strong = res_obj$reservoir_shared_strong,
  reservoir_shared_moderate = res_obj$reservoir_shared_moderate,
  reservoir_distinct_bg = res_obj$reservoir_distinct_bg,
  reservoir_distinct_trap = res_obj$reservoir_distinct_trap,
  cpm_tau_grid = seq(0.1, 0.7, by = 0.1)
)
scaling_summary <- summarise_scaling_runtime(scaling_res)

print(scaling_summary)
library(tidyr)
library(ggplot2)

scaling_long <- scaling_summary %>%
  select(
    scaling_type,
    scenario,
    n_nodes,
    true_cluster_size,
    louvain_mean_time_sec,
    cpm_mean_time_sec,
    spinglass_mean_time_sec
  ) %>%
  pivot_longer(
    cols = c(louvain_mean_time_sec, cpm_mean_time_sec, spinglass_mean_time_sec),
    names_to = "method",
    values_to = "mean_time_sec"
  ) %>%
  mutate(
    method = recode(
      method,
      louvain_mean_time_sec = "Louvain",
      cpm_mean_time_sec = "CPM",
      spinglass_mean_time_sec = "Signed Spinglass"
    ),
    scaling_type = recode(
      scaling_type,
      background_expansion = "Background-node expansion",
      expanded_ldtrap = "Expanded LD-trap"
    ),
    scenario = recode(
      scenario,
      triangle = "Triangle",
      ld_trap = "LD trap",
      local_bridge_trap = "Local bridge trap"
    )
  )

p_scaling_runtime <- ggplot(
  scaling_long,
  aes(x = n_nodes, y = mean_time_sec, group = method, colour = method)
) +
  geom_line(aes(linetype = method), linewidth = 0.9) +
  geom_point(size = 2.4) +
  facet_grid(scaling_type ~ scenario, scales = "free_y") +
  scale_colour_manual(
    values = c(
      "Louvain" = "#4E79A7",
      "CPM" = "#F28E2B",
      "Signed Spinglass" = "#59A14F"
    )
  ) +
  labs(
    title = "Runtime scaling with number of nodes",
    x = "Number of nodes",
    y = "Mean runtime per replicate (seconds)",
    colour = "Method",
    linetype = "Method"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 25, hjust = 1),
    legend.position = "right"
  )

print(p_scaling_runtime)
p_scaling_runtime_log <- ggplot(
  scaling_long,
  aes(x = n_nodes, y = mean_time_sec, group = method, colour = method)
) +
  geom_line(aes(linetype = method), linewidth = 0.9) +
  geom_point(size = 2.4) +
  scale_y_log10() +
  facet_grid(scaling_type ~ scenario, scales = "free_y") +
  scale_colour_manual(
    values = c(
      "Louvain" = "#4E79A7",
      "CPM" = "#F28E2B",
      "Signed Spinglass" = "#59A14F"
    )
  ) +
  labs(
    title = "Runtime scaling with number of nodes, log scale",
    x = "Number of nodes",
    y = "Mean runtime per replicate (seconds, log scale)",
    colour = "Method",
    linetype = "Method"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 25, hjust = 1),
    legend.position = "right"
  )

print(p_scaling_runtime_log)
