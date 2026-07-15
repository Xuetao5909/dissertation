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

library(readr)


write_csv(summary_tab, "simulation_summary_results.csv")
library(dplyr)
library(tidyr)
library(ggplot2)

success_long <- summary_tab %>%
  select(
    scenario,
    louvain_success_rate,
    infomap_success_rate,
    cpm_success_rate,
    spinglass_success_rate
  ) %>%
  pivot_longer(
    cols = -scenario,
    names_to = "method",
    values_to = "success_rate"
  ) %>%
  mutate(
    method = recode(
      method,
      louvain_success_rate = "Louvain",
      infomap_success_rate = "Infomap",
      cpm_success_rate = "CPM",
      spinglass_success_rate = "Signed Spinglass"
    ),
    method = factor(
      method,
      levels = c("Louvain", "Infomap", "CPM", "Signed Spinglass")
    ),
    scenario = recode(
      scenario,
      triangle = "Triangle",
      bowtie = "Bowtie overlap",
      ld_trap = "LD trap",
      local_bridge_trap = "Local bridge trap"
    ),
    scenario = factor(
      scenario,
      levels = c("Triangle", "Bowtie overlap", "LD trap", "Local bridge trap")
    )
  )

p_success <- ggplot(
  success_long,
  aes(x = scenario, y = success_rate * 100, fill = method)
) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  scale_y_continuous(limits = c(0, 100)) +
  labs(
    title = "Figure 3.1. Success rate across stress-test scenarios",
    x = "Scenario",
    y = "Success rate (%)",
    fill = "Method"
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 25, hjust = 1),
    plot.title = element_text(face = "bold"),
    legend.position = "bottom"
  )

print(p_success)

ggsave(
  "output/Figure_3_1_success_rate.png",
  p_success,
  width = 10,
  height = 5.5,
  dpi = 300
)
# Runtime summary table

# -------------------------
runtime_total_by_scenario <- res %>%
  group_by(scenario) %>%
  summarise(
    louvain_total_time_sec = sum(louvain_time_sec, na.rm = TRUE),
    infomap_total_time_sec = sum(infomap_time_sec, na.rm = TRUE),
    cpm_total_time_sec = sum(cpm_time_sec, na.rm = TRUE),
    spinglass_total_time_sec = sum(spinglass_time_sec, na.rm = TRUE),
    .groups = "drop"
  )

print(runtime_total_by_scenario)

if (!dir.exists("output")) dir.create("output", recursive = TRUE)

readr::write_csv(res, "output/stress_test_results.csv")
readr::write_csv(summary_tab, "output/stress_test_summary.csv")
readr::write_csv(runtime_total_by_scenario, "output/runtime_total_by_scenario.csv")
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
    infomap_ari,
    cpm_ari,
    spinglass_ari
  ) %>%
  pivot_longer(
    cols = c(louvain_ari, infomap_ari, cpm_ari, spinglass_ari),
    names_to = "method",
    values_to = "ARI"
  ) %>%
  mutate(
    method = recode(
      method,
      louvain_ari = "Louvain",
      infomap_ari = "Infomap",
      cpm_ari = "CPM",
      spinglass_ari = "Signed Spinglass"
    ),
    method = factor(
      method,
      levels = c("Louvain", "Infomap", "CPM", "Signed Spinglass")
    ),
    scenario = recode(
      scenario,
      ld_trap = "LD trap",
      local_bridge_trap = "Local bridge trap"
    )
  )

p_ari_box <- ggplot(ari_long, aes(x = method, y = ARI)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +
  geom_jitter(width = 0.12, height = 0, alpha = 0.6, size = 1.8) +
  facet_wrap(~ scenario, nrow = 1) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    title = "Figure 3.3. Adjusted Rand Index across stress-test replicates",
    x = "Method",
    y = "Adjusted Rand Index"
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 25, hjust = 1),
    plot.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold")
  )

print(p_ari_box)

ggsave(
  filename = "output/Figure_3_2_ARI_boxplot.png",
  plot = p_ari_box,
  width = 9,
  height = 5,
  dpi = 300
)

print(p_ari_box)
runtime_total_long <- res %>%
  group_by(scenario) %>%
  summarise(
    Louvain = sum(louvain_time_sec, na.rm = TRUE),
    Infomap = sum(infomap_time_sec, na.rm = TRUE),
    CPM = sum(cpm_time_sec, na.rm = TRUE),
    `Signed Spinglass` = sum(spinglass_time_sec, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  tidyr::pivot_longer(
    cols = c(Louvain, Infomap, CPM, `Signed Spinglass`),
    names_to = "method",
    values_to = "total_time_sec"
  ) %>%
  mutate(
    method = factor(
      method,
      levels = c("Louvain", "Infomap", "CPM", "Signed Spinglass")
    )
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
    infomap_mean_time_sec,
    cpm_mean_time_sec,
    spinglass_mean_time_sec
  ) %>%
  pivot_longer(
    cols = c(
      louvain_mean_time_sec,
      infomap_mean_time_sec,
      cpm_mean_time_sec,
      spinglass_mean_time_sec
    ),
    names_to = "method",
    values_to = "mean_time_sec"
  ) %>%
  mutate(
    method = recode(
      method,
      louvain_mean_time_sec = "Louvain",
      infomap_mean_time_sec = "Infomap",
      cpm_mean_time_sec = "CPM",
      spinglass_mean_time_sec = "Signed Spinglass"
    ),
    method = factor(
      method,
      levels = c("Louvain", "Infomap", "CPM", "Signed Spinglass")
    ),
    scaling_type_label = recode(
      scaling_type,
      background_expansion = "Background-node expansion",
      expanded_ldtrap = "Expanded LD-trap"
    ),
    scenario_label = recode(
      scenario,
      triangle = "Triangle",
      ld_trap = "LD trap",
      local_bridge_trap = "Local bridge trap"
    ),
    panel_label = paste0(scaling_type_label, "\n", scenario_label)
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
      "Infomap" = "#B07AA1",
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
scaling_plot_data <- scaling_long %>%
  filter(
    scaling_type == "expanded_ldtrap",
    scenario == "ld_trap"
  )

method_cols <- c(
  "Louvain" = "#4E79A7",
  "Infomap" = "#B07AA1",
  "CPM" = "#F28E2B",
  "Signed Spinglass" = "#59A14F"
)

p_scaling_runtime_log <- ggplot(
  scaling_plot_data,
  aes(
    x = n_nodes,
    y = mean_time_sec,
    group = method,
    colour = method,
    linetype = method
  )
) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 3) +
  scale_y_log10() +
  scale_colour_manual(values = method_cols, drop = FALSE) +
  labs(
    title = "Figure 3.4. Runtime scaling in the expanded LD-trap scenario",
    x = "Number of nodes",
    y = "Mean runtime per replicate (seconds, log scale)",
    colour = "Method",
    linetype = "Method"
  ) +
  theme_bw(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 25, hjust = 1),
    legend.position = "bottom"
  )

print(p_scaling_runtime_log)

ggsave(
  filename = "output/Figure_3_3_runtime_scaling_expanded_ldtrap_log.png",
  plot = p_scaling_runtime_log,
  width = 8.5,
  height = 5,
  dpi = 300
)
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
      "Infomap" = "#B07AA1",
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

library(CliquePercolation)
# -------------------------
# -------------------------
# Runtime plot: one representative scenario for presentation
# -------------------------

runtime_plot_data <- runtime_total_long %>%
  filter(scenario == "ld_trap") %>%
  mutate(
    method = factor(
      method,
      levels = c("Signed Spinglass", "CPM", "Infomap", "Louvain")
    )
  )

p_runtime_total_one <- ggplot(
  runtime_plot_data,
  aes(x = total_time_sec, y = method)
) +
  geom_col(width = 0.55, fill = "grey45") +
  labs(
    title = "Total runtime by method: LD trap scenario",
    x = "Total runtime across 20 replicates (seconds)",
    y = "Method"
  ) +
  theme_bw(base_size = 16) +
  theme(
    plot.title = element_text(face = "bold", size = 18),
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 14),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )

print(p_runtime_total_one)

ggsave(
  filename = "output/runtime_ld_trap_one_scenario.png",
  plot = p_runtime_total_one,
  width = 8.5,
  height = 4.8,
  dpi = 300
)
# -------------------------
# Runtime scaling plot for presentation
# -------------------------

library(dplyr)
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
    scaling_type_label = recode(
      scaling_type,
      background_expansion = "Background-node expansion",
      expanded_ldtrap = "Expanded LD-trap"
    ),
    scenario_label = recode(
      scenario,
      triangle = "Triangle",
      ld_trap = "LD trap",
      local_bridge_trap = "Local bridge trap"
    )
  )

# Choose one representative result for the presentation
scaling_plot_data <- scaling_long %>%
  filter(
    scaling_type_label == "Expanded LD-trap",
    scenario_label == "LD trap"
  )

p_scaling_runtime_slide <- ggplot(
  scaling_plot_data,
  aes(
    x = n_nodes,
    y = mean_time_sec,
    group = method,
    colour = method,
    linetype = method
  )
) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 3) +
  labs(
    title = "Runtime scaling with number of nodes",
    subtitle = "Representative result: expanded LD-trap scenario",
    x = "Number of nodes",
    y = "Mean runtime per replicate (seconds)",
    colour = "Method",
    linetype = "Method"
  ) +
  theme_bw(base_size = 16) +
  theme(
    plot.title = element_text(face = "bold", size = 19),
    plot.subtitle = element_text(size = 14),
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 13),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12),
    legend.position = "right",
    panel.grid.minor = element_blank()
  )

print(p_scaling_runtime_slide)

if (!dir.exists("output")) dir.create("output", recursive = TRUE)

ggsave(
  filename = "output/runtime_scaling_time_vs_nodes_slide.png",
  plot = p_scaling_runtime_slide,
  width = 9,
  height = 5.2,
  dpi = 300
)

scaling_res_test <- run_scaling_test(
  n_reps = 2,
  node_sizes = c(10, 20),
  reservoir_shared_strong = res_obj$reservoir_shared_strong,
  reservoir_shared_moderate = res_obj$reservoir_shared_moderate,
  reservoir_distinct_bg = res_obj$reservoir_distinct_bg,
  reservoir_distinct_trap = res_obj$reservoir_distinct_trap,
  cpm_tau_grid = seq(0.1, 0.7, by = 0.1)
)

scaling_summary_test <- summarise_scaling_runtime(scaling_res_test)

print(scaling_summary_test)
names(scaling_summary_test)
unique(scaling_long$method)
names(scaling_summary)
library(ggplot2)
library(dplyr)
library(tibble)
library(patchwork)

# Fixed bowtie node layout: two triangles share node 3
nodes <- tibble(
  node = c("1", "2", "3", "4", "5"),
  x = c(-1.2, -1.2, 0, 1.2, 1.2),
  y = c(0.8, -0.8, 0, 0.8, -0.8)
)

edges <- tibble(
  from = c("1", "1", "2", "3", "3", "4"),
  to   = c("2", "3", "3", "4", "5", "5")
)

make_edges <- function(edges, nodes) {
  edges %>%
    left_join(nodes %>% rename(from = node, x_from = x, y_from = y), by = "from") %>%
    left_join(nodes %>% rename(to = node, x_to = x, y_to = y), by = "to")
}

edges_plot <- make_edges(edges, nodes)

# Method-specific representative solutions
solution_nodes <- bind_rows(
  nodes %>%
    mutate(
      method = "Louvain",
      community = case_when(
        node %in% c("1", "2", "3") ~ "Community 1",
        node %in% c("4", "5") ~ "Community 2"
      )
    ),
  nodes %>%
    mutate(
      method = "Infomap",
      community = "Merged community"
    ),
  nodes %>%
    mutate(
      method = "CPM",
      community = case_when(
        node %in% c("1", "2") ~ "Community 1",
        node == "3" ~ "Overlap hub",
        node %in% c("4", "5") ~ "Community 2"
      )
    ),
  nodes %>%
    mutate(
      method = "Signed Spinglass",
      community = case_when(
        node %in% c("1", "2", "3") ~ "Community 1",
        node %in% c("4", "5") ~ "Community 2"
      )
    )
)

solution_edges <- bind_rows(
  edges_plot %>% mutate(method = "Louvain"),
  edges_plot %>% mutate(method = "Infomap"),
  edges_plot %>% mutate(method = "CPM"),
  edges_plot %>% mutate(method = "Signed Spinglass")
)

solution_nodes$method <- factor(
  solution_nodes$method,
  levels = c("Louvain", "Infomap", "CPM", "Signed Spinglass")
)

solution_edges$method <- factor(
  solution_edges$method,
  levels = c("Louvain", "Infomap", "CPM", "Signed Spinglass")
)

community_cols <- c(
  "Community 1" = "#4E79A7",
  "Community 2" = "#F28E2B",
  "Overlap hub" = "#EAC54F",
  "Merged community" = "#59A14F"
)

p_bowtie_solutions <- ggplot() +
  geom_segment(
    data = solution_edges,
    aes(x = x_from, y = y_from, xend = x_to, yend = y_to),
    linewidth = 0.8,
    colour = "grey35"
  ) +
  geom_point(
    data = solution_nodes,
    aes(x = x, y = y, fill = community),
    shape = 21,
    size = 10,
    colour = "grey25",
    stroke = 0.8
  ) +
  geom_text(
    data = solution_nodes,
    aes(x = x, y = y, label = node),
    size = 4.5,
    fontface = "bold"
  ) +
  facet_wrap(~ method, nrow = 1) +
  scale_fill_manual(values = community_cols) +
  coord_equal(xlim = c(-1.8, 1.8), ylim = c(-1.2, 1.2), clip = "off") +
  labs(
    title = "Figure 3.2. Representative clustering outputs for one bowtie example",
    fill = "Detected community"
  ) +
  theme_void(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    strip.text = element_text(face = "bold", size = 11),
    legend.position = "bottom"
  )

print(p_bowtie_solutions)

ggsave(
  filename = "output/Figure_3_2_bowtie_example_solutions.png",
  plot = p_bowtie_solutions,
  width = 11,
  height = 4.5,
  dpi = 300
)