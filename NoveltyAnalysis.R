library(readr)
library(dplyr)
library(stringr)
library(tidyr)
library(purrr)
library(igraph)
library(tidygraph)
library(ggraph)
library(ggrepel)
library(graphlayouts)

path_string <- "stringhi_edges.csv"
path_hidg   <- "hidgidb_edges.csv"

# --------- Load ----------
string_edges <- read_csv(path_string, show_col_types = FALSE)
hidg_edges   <- read_csv(path_hidg,   show_col_types = FALSE)

# Expect columns: color, from, to, width, title
stopifnot(all(c("color","from","to","width") %in% names(string_edges)))
stopifnot(all(c("color","from","to","width") %in% names(hidg_edges)))

# --------- Helpers ----------
node_type <- function(x) {
  case_when(
    str_starts(x, "@GENE_")     ~ "gene",
    str_starts(x, "@CHEMICAL_") ~ "chemical",
    str_starts(x, "@DISEASE_")  ~ "disease",
    TRUE                        ~ "other"
  )
}

node_label <- function(x) {
  # removes @GENE_ / @CHEMICAL_ / etc
  str_remove(x, "^@[^_]+_")
}

make_vertices <- function(edges_df) {
  tibble(name = unique(c(edges_df$from, edges_df$to))) %>%
    mutate(
      type  = node_type(name),
      label = node_label(name)
    )
}


make_graph <- function(edges_df) {
  verts <- make_vertices(edges_df)
  
  # make sure optional columns exist
  if (!"title" %in% names(edges_df)) edges_df$title <- NA_character_
  if (!"width" %in% names(edges_df)) edges_df$width <- NA_real_
  
  g <- graph_from_data_frame(
    d = edges_df %>%
      transmute(
        from,
        to,
        edge_class = color,
        title = title,
        width = width
      ),
    directed = FALSE,
    vertices = verts
  )
  
  # start with explicit per-edge evidence bookkeeping
  E(g)$dup_n       <- 1L
  E(g)$kg_red_n    <- as.integer(E(g)$edge_class == "red")
  E(g)$kg_green_n  <- as.integer(E(g)$edge_class == "green")
  
  # simplest weight = one unit per occurrence
  #E(g)$weight <- 1
  
  g
}


combine_edge_class <- function(x) {
  ux <- unique(na.omit(x))
  if (length(ux) == 0) return(NA_character_)
  if (length(ux) == 1) return(ux)
  "mixed"
}

combine_title <- function(x) {
  x <- unique(na.omit(x))
  if (length(x) == 0) return(NA_character_)
  paste(x, collapse = " | ")
}

combine_width <- function(x) {
  mean(x, na.rm = TRUE)
}


zscore <- function(x) {
  (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
}


# ============================================================
# 2) Reproduce basic summaries
# ============================================================

cat("\n--- Fig 2c edge list summary (stringhi_edges.csv) ---\n")
print(string_edges %>% count(color) %>% arrange(desc(n)))
cat("Unique nodes:", length(unique(c(string_edges$from, string_edges$to))), "\n")
cat("Unique genes:", length(unique(c(string_edges$from, string_edges$to))[str_starts(unique(c(string_edges$from, string_edges$to)), "@GENE_")]), "\n")

cat("\n--- Fig 3c edge list summary (hidgidb_edges.csv) ---\n")
print(hidg_edges %>% count(color) %>% arrange(desc(n)))

hidg_etypes <- hidg_edges %>%
  mutate(
    from_type = node_type(from),
    to_type   = node_type(to),
    etype     = paste(pmin(from_type, to_type), pmax(from_type, to_type), sep = "-")
  ) %>%
  count(etype) %>%
  arrange(desc(n))

print(hidg_etypes)

all_nodes <- unique(c(string_edges$from, string_edges$to, hidg_edges$from, hidg_edges$to))
cat("\n--- Union node counts ---\n")
cat("Total nodes:", length(all_nodes), "\n")
cat("Genes:", sum(str_starts(all_nodes, "@GENE_")), "\n")
cat("Chemicals:", sum(str_starts(all_nodes, "@CHEMICAL_")), "\n")

# ============================================================
# 3) Build KG-derived combined network (red + green edges only)
#    - red   = KG-only
#    - green = overlap (KG + external)
# ============================================================

kg_edges2 <- bind_rows(
  string_edges %>% filter(color %in% c("red", "green")),
  hidg_edges   %>% filter(color %in% c("red", "green"))
)

if (!"title" %in% names(kg_edges2)) kg_edges2$title <- NA_character_
if (!"width" %in% names(kg_edges2)) kg_edges2$width <- NA_real_

kg_edges_agg <- kg_edges2 %>%
  mutate(
    a = pmin(from, to),
    b = pmax(from, to),
    dup_n      = 1L,
    kg_red_n   = as.integer(color == "red"),
    kg_green_n = as.integer(color == "green"),
    weight     = 1
  ) %>%
  group_by(a, b) %>%
  summarise(
    dup_n      = sum(dup_n),
    kg_red_n   = sum(kg_red_n),
    kg_green_n = sum(kg_green_n),
    weight     = sum(weight),
    edge_class = case_when(
      kg_green_n > 0 & kg_red_n == 0 ~ "green",
      kg_red_n > 0 & kg_green_n == 0 ~ "red",
      TRUE ~ "mixed"
    ),
    title = {
      x <- unique(na.omit(title))
      if (length(x) == 0) NA_character_ else paste(x, collapse = " | ")
    },
    width = mean(width, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  rename(from = a, to = b)

kg_edges_agg <- kg_edges_agg %>%
  rename(color = edge_class)

g_kg <- make_graph(kg_edges_agg)

# Giant component
comp <- components(g_kg)
giant_id <- which.max(comp$csize)
g <- induced_subgraph(g_kg, vids = V(g_kg)[comp$membership == giant_id])

cat("\n--- KG combined network (giant component) ---\n")
cat("Nodes:", vcount(g), "Edges:", ecount(g), "\n")
cat("Genes:", sum(V(g)$type == "gene"), "Chemicals:", sum(V(g)$type == "chemical"), "\n")

# Plot
set.seed(78)

V(g)$degree <- degree(g)

g_tbl <- as_tbl_graph(g) %>%
  activate(nodes) %>%
  mutate(
    type = factor(type),
    deg  = degree,
    # tweak these two to taste
    node_size = scales::rescale(deg, to = c(3, 14))
  )

p <- ggraph(
  g_tbl,
  layout = "centrality",
  centrality = centrality_betweenness(directed = FALSE),  
  scale = TRUE,
  niter = 500,
  tseq = seq(0, 1, 0.2)
) +
  geom_edge_link0(aes(color = edge_class), alpha = 0.5, width = 0.3) +
  scale_edge_colour_manual(values = c("green" = "darkgreen", "red" = "darkred")) +
  geom_node_point(aes(shape = type, color = type,size = node_size),alpha = 0.5, stroke = 1) +
  geom_node_text(
    aes(label = label),
    size = 3,
    repel = TRUE,
    max.overlaps = Inf,
    box.padding = 0.5,nudge_y = 0.5,
    point.padding = 0.2
  ) +
  scale_size_identity() +
  theme_void()

print(p)

ggsave("kg_giant_component.png", p, width = 14, height = 10, dpi = 300)


# ============================================================
# 4) Centrality + composite KG-score
# ============================================================

deg <- degree(g, mode = "all")
bet <- betweenness(g, directed = FALSE, normalized = TRUE)
clo <- closeness(g, normalized = TRUE)

V(g)$degree      <- deg
V(g)$betweenness <- bet
V(g)$closeness   <- clo
V(g)$kg_score    <- zscore(deg) + zscore(bet) + zscore(clo)

node_tbl <- tibble(
  name       = V(g)$name,
  label      = V(g)$label,
  type       = V(g)$type,
  degree     = V(g)$degree,
  betweenness= V(g)$betweenness,
  closeness  = V(g)$closeness,
  kg_score   = V(g)$kg_score
) %>% arrange(desc(kg_score))


# ============================================================
# 5) Community detection (fast_greedy)
# ============================================================

library(ideanet)
library(dplyr)

greedy <- cluster_fast_greedy(
  g,
  merges = TRUE,
  modularity = TRUE,
  membership = TRUE
  
)

# store on graph
V(g)$community <- membership(greedy)

# summary
K <- length(unique(membership(greedy)))
Q <- modularity(g, membership(greedy), weights = E(g)$weight)

cat("Chosen partition: leiden(res=)",
    "K=", K,
    "Q=", Q, "\n")

# continue exactly like before
comm_summary <- node_tbl %>%
  mutate(community = V(g)$community[match(name, V(g)$name)]) %>%
  group_by(community, type) %>%
  summarise(n = n(), .groups = "drop") %>%
  arrange(community, desc(n))

cat("\n--- Community sizes (by type) ---\n")
print(comm_summary)

top_gene_by_comm <- node_tbl %>%
  filter(type == "gene") %>%
  mutate(community = V(g)$community[match(name, V(g)$name)]) %>%
  group_by(community) %>%
  slice_max(order_by = kg_score, n = 100, 
            with_ties = FALSE) %>%
  ungroup() %>%
  arrange(community, desc(kg_score))

all_by_comm <- node_tbl %>%
  #filter(type == "gene") %>%
  mutate(community = V(g)$community[match(name, V(g)$name)]) %>%
  group_by(community) %>%
  #slice_max(order_by = kg_score, n = 100, 
  #          with_ties = FALSE) %>%
  ungroup() %>%
  arrange(community, desc(kg_score))

write.csv(all_by_comm,file = "all_by_comm.csv")

# ============================================================
# 6) Gene enrichment: gProfile
# ============================================================

library(gprofiler2)
for (com in c(1:max(unique(membership(greedy))))) {
  genes <- all_by_comm %>% filter(type == "gene" & community == com)
  
  gostres <- gost(query = as.vector(genes[,"label"])[[1]],
                  organism = "hsapiens")
  df <- gostres$result
  
  df[] <- lapply(df, function(col) {
    if (is.list(col)) {
      sapply(col, function(x) paste(x, collapse = "; "))
    } else {
      col
    }
  })
  
  write.csv(df, paste0("GoEnrichment_forCom_", com, ".csv"), row.names = FALSE)
  saveRDS(gostres,file = paste("GoEnrichment_forCom_",com,".RDS"))
  out_plot <- gostplot(gostres, capped = FALSE, interactive = F)
  
  highlight_terms <- NA
  for (db in c("CORUM","GO:BP","KEGG","REAC")) {
    highlight_terms <- c(highlight_terms,gostres$result[gostres$result$source == db,]$term_id[1:6])
  }
   
  out_plot <- publish_gostplot(out_plot, highlight_terms = highlight_terms)
  
  pdf(paste("GoEnrichment_forCom_",com,".pdf"),width = 10, height = 20)
  plot(out_plot)
  dev.off()
}

g_tbl2 <- g_tbl %>%
  activate(nodes) %>%
  left_join(
    all_by_comm %>%
      select(name, betweenness, closeness, kg_score, community),
    by = "name"
  )

my_labels <- c("FAH","PAH", "SOST", "DKK1","CTNNB1", "TNFSF11", "TNFSF11A","SYK","PLCG2",
               "VAV1","NFATC1","TNF","IL1B","NFE2L2","KEAP1","SOD2","CAT","FOXO3",
               "IL6", "IL1A","Metformin","Sirolimus","Ascorbic_Acid",
               "romosozumab","Denosumab","Teriparatide","alpha_Tocopherol","Thioctic_Acid")

g_tbl3 <- g_tbl2 %>%
  activate(nodes) %>%
  mutate(label_plot = ifelse(label %in% my_labels, label, NA))

p_com <- ggraph(
  g_tbl3,
  #layout = "stress",
  layout = "centrality",
  centrality = centrality_betweenness(directed = FALSE),  
  scale = TRUE,
  #niter = 500,
  tseq = seq(0, 1, 0.2)
) +
    geom_edge_link0(aes(color = "lightgray"),#edge_class), 
                  alpha = 1, width = 0.3) +
  scale_edge_colour_manual(values = c("green" = "darkgreen", "red" = "darkred")) +
  geom_node_point(aes(shape = type, color = as.character(community),size = node_size/1.5),alpha = 0.8, stroke = 1) +
  geom_node_text(
    aes(label = label_plot),
    size = 3,
    repel = F,
    #max.overlaps = 1,
    #box.padding = 0.5,
    nudge_y = 1,
    #point.padding = 0.2
  ) +
  scale_size_identity() + 
  theme_void()+
  guides(color=guide_legend("Modules")) 

print(p_com)

ggsave("kg_giant_component_modules.png", p_com, width = 7, height = 5, dpi = 400)





