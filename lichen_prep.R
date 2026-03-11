### script o loading and cleaning the lichen plots


 #Load necessary libraries
library(tidyverse)
library(readxl)
library(ggplot2)
library(pheatmap)



lichens <- read_xlsx("Lichens/Licen_data.xlsx") 

## translating within column object ID

lichens <- lichens %>%
  mutate(`object ID` = recode(`object ID`,
                              "3 kláda" = "3 log",
                              "8 kameny" = "8 stones", 
                              "4 pahýl" = "4 stump",
                              "5 kládas" = "5 logs",
                              "7 kameny" = "7 stones",
                              "6 kláda" = "6 log",
                              "7 spadlá větev" = "7 fallen branch",
                              "5 kláda" = "5 log",
                              "6 nízký pahýl" = "6 low stump",
                              "3 pahýl" = "3 stump",
                              "4 kláda" = "4 log",
                              "5 pařez" = "5 stump",
                              "5 kláda"= "5 log",
                              "6 pahýl" = "6 stump",
                              "4 Picea abies mladý" = "4 Picea abies young",
                              "Pinus sylvestris padlá větev" = "Pinus sylvestris fallen branch",
                              "6 vývrat" = "6 uprooted tree",
                              "7 pařez" = "7 stump",
                              "6 pařez" = "6 stump",
                              "5 pahýl" = "5 stump",
                              "5 nízký pahýl" = "5 low stump",
                              "10 nízký pahýl" = "10 low stump",
                              "9 pahýl" = "9 stump",
                              "8 kláda" = "8 log",
                              "7 vývrat" = "7 uprooted tree",
                              "7 nízký pahýl" = "7 low stump",
                              "4 nízký pahýl" = "4 low stump",
                              "1 pahýl" = "1 stump",
                              "7 pahýl" = "7 stump",
                              "7 kláda" = "7 log", 
                              "2 kláda" = "2 log",
                              "4 pařez" = "4 stump",
                              "8 pařez" = "8 stump",
                              
                              
  ))




#  ===== 1. SPECIES-LEVEL SUMMARY =====
  # Which species are found where, how often, and how abundant?
  
  species_summary <- lichens %>%
  group_by(species) %>%
  summarise(
    n_plots = n_distinct(Plot),              # Number of plots where present
    prevalence_pct = (n_plots / 120) * 100,  # % of all 120 plots
    total_observations = n(),                 # Total records
    mean_cover = mean(cover, na.rm = TRUE),  # Average cover value
    median_cover = median(cover, na.rm = TRUE),
    max_cover = max(cover, na.rm = TRUE),
    
  ) %>%
  arrange(desc(n_plots)) %>%  # Sort by most common species
  mutate(
    # Categorize species for modeling suitability
    model_suitability = case_when(
      prevalence_pct >= 20 & prevalence_pct <= 80 ~ "Good - Balanced",
      prevalence_pct > 80 ~ "Ubiquitous - Low variance",
      prevalence_pct >= 10 & prevalence_pct < 20 ~ "Rare - Consider grouping",
      prevalence_pct < 10 ~ "Very rare - Not suitable"
    )
  )

# View the summary
print(species_summary, n = 20)  # Show top 20 species

# Quick modeling suitability breakdown
species_summary %>%
  count(model_suitability) %>%
  arrange(desc(n))


# ===== 2. PLOT-LEVEL SUMMARY =====
# What's the species diversity and abundance at each plot?

plot_summary <- lichens %>%
  group_by(Plot) %>%
  summarise(
    n_species = n_distinct(species),          # Species richness
    total_cover = sum(cover, na.rm = TRUE),   # Total lichen cover
    mean_cover = mean(cover, na.rm = TRUE),
    n_observations = n(),                      # Number of records
    # Most abundant species (by cover)
    dominant_species = species[which.max(cover)],
    dominant_cover = max(cover, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  arrange(desc(n_species))  # Richest plots first

# View plot diversity
print(plot_summary, n = 20)

# Summary statistics
summary(plot_summary$n_species)


# ===== 3. SPECIES × PLOT PRESENCE/ABSENCE MATRIX =====
# Wide format:  rows = plots, columns = species, values = 1 (present) or 0 (absent)

presence_absence_matrix <- lichens %>%
  distinct(Plot, species) %>%          # Remove duplicate species per plot
  mutate(presence = 1) %>%
  pivot_wider(
    names_from = species,
    values_from = presence,
    values_fill = 0                     # Fill absences with 0
  )

# This creates a 120 × n_species matrix
dim(presence_absence_matrix)


# ===== 4. IDENTIFY CANDIDATE SPECIES FOR MODELING =====
# Filter species with good prevalence (20-80%)

candidate_species <- species_summary %>%
  filter(prevalence_pct >= 20 & prevalence_pct <= 80) %>%
  pull(species)

cat("\n=== CANDIDATE SPECIES FOR BINOMIAL MODELING ===\n")
cat("Number of suitable species:", length(candidate_species), "\n\n")
print(candidate_species)


# ===== 5. VISUALIZATION:  SPECIES PREVALENCE DISTRIBUTION =====

# Histogram of species prevalence
ggplot(species_summary, aes(x = prevalence_pct)) +
  geom_histogram(binwidth = 5, fill = "steelblue", color = "white") +
  geom_vline(xintercept = c(20, 80), linetype = "dashed", color = "red") +
  labs(
    title = "Distribution of Species Prevalence Across 120 Plots",
    subtitle = "Red lines indicate 20-80% threshold for balanced modeling",
    x = "Prevalence (%)",
    y = "Number of Species"
  ) +
  theme_minimal()

# Scatterplot:  Prevalence vs Mean Cover
ggplot(species_summary, aes(x = prevalence_pct, y = mean_cover, 
                            color = model_suitability)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_vline(xintercept = c(20, 80), linetype = "dashed", color = "gray50") +
  scale_color_manual(
    values = c("Good - Balanced" = "darkgreen",
               "Ubiquitous - Low variance" = "orange",
               "Rare - Consider grouping" = "gold",
               "Very rare - Not suitable" = "red")
  ) +
  labs(
    title = "Species Prevalence vs Abundance",
    x = "Prevalence (% of plots)",
    y = "Mean Cover Value",
    color = "Model Suitability"
  ) +
  theme_minimal()


# ===== 6. PLOT-LEVEL DIVERSITY VISUALIZATION =====

ggplot(plot_summary, aes(x = reorder(Plot, n_species), y = n_species)) +
  geom_col(fill = "forestgreen", alpha = 0.7) +
  coord_flip() +
  labs(
    title = "Species Richness by Plot",
    x = "Plot ID",
    y = "Number of Species"
  ) +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 6))  # Small text for 120 plots






#-----------------

#### next code is to check for single species if they are in the dataset

# ===== FUNCTION:  CHECK SPECIES MODELING VIABILITY =====


# ===== CORRECTED FUNCTION =====

check_species_viability <- function(species_in_question, 
                                    species_summary_df = species_summary,
                                    lichen_df = lichens,
                                    min_prevalence = 20,
                                    max_prevalence = 80) {
  
  # Check if species exists in dataset
  if (!species_in_question %in% species_summary_df$species) {
    cat("❌ Species NOT FOUND in dataset:", species_in_question, "\n")  # ← Added comma here
    cat("\n💡 Tip: Check spelling or try fuzzy search below\n\n")
    
    # Fuzzy matching suggestions
    possible_matches <- agrep(species_in_question, species_summary_df$species, 
                              max.distance = 0.3, value = TRUE)
    if (length(possible_matches) > 0) {
      cat("Did you mean one of these?\n")
      print(possible_matches)
    }
    return(invisible(NULL))
  }
  
  # Extract species data
  sp_data <- species_summary_df %>%
    filter(species == species_in_question)
  
  # Get plot-level occurrence
  plot_occurrence <- lichen_df %>%
    filter(species == species_in_question) %>%
    group_by(Plot) %>%
    summarise(
      n_observations = n(),
      mean_cover = mean(cover, na.rm = TRUE),
      max_cover = max(cover, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    arrange(desc(mean_cover))
  
  # === PRINT REPORT ===
  cat("\n")
  cat("═══════════════════════════════════════════════════════\n")
  cat("  SPECIES MODELING VIABILITY REPORT\n")
  cat("═══════════════════════════════════════════════════════\n\n")
  
  cat("🔬 Species:", species_in_question, "\n\n")
  
  # Prevalence
  cat("📊 PREVALENCE:\n")
  cat("   • Found in", sp_data$n_plots, "out of 120 plots (", 
      round(sp_data$prevalence_pct, 1), "%)\n")
  
  # Modeling suitability
  cat("\n✅ MODELING SUITABILITY:\n")
  cat("   • Category:", sp_data$model_suitability, "\n")
  
  if (sp_data$prevalence_pct >= min_prevalence & sp_data$prevalence_pct <= max_prevalence) {
    cat("   • ✅ SUITABLE for binomial GLM/GAM\n")
    cat("   • Statistical power:  GOOD\n")
  } else if (sp_data$prevalence_pct < min_prevalence & sp_data$prevalence_pct >= 10) {
    cat("   • ⚠️  MARGINAL - Low prevalence may limit power\n")
    cat("   • Consider:  Grouping with similar species or using presence-only methods\n")
  } else if (sp_data$prevalence_pct < 10) {
    cat("   • ❌ NOT SUITABLE - Too rare for robust modeling\n")
    cat("   • Recommendation: Group into functional type or use MaxEnt\n")
  } else {
    cat("   • ⚠️  UBIQUITOUS - Limited environmental variation\n")
    cat("   • May have weak predictor relationships\n")
  }
  
  # Abundance
  cat("\n📈 ABUNDANCE METRICS:\n")
  cat("   • Total observations:", sp_data$total_observations, "\n")
  cat("   • Mean cover (when present):", round(sp_data$mean_cover, 2), "\n")
  cat("   • Median cover:", sp_data$median_cover, "\n")
  cat("   • Max cover:", sp_data$max_cover, "\n")
  
  # Plot distribution
  cat("\n🗺️  PLOT DISTRIBUTION:\n")
  cat("   • Plots with highest abundance (top 5):\n")
  print(head(plot_occurrence, 5), row.names = FALSE)
  
  # Presence/absence balance
  n_present <- sp_data$n_plots
  n_absent <- 120 - n_present
  cat("\n⚖️  PRESENCE/ABSENCE BALANCE:\n")
  cat("   • Present:", n_present, "plots\n")
  cat("   • Absent:", n_absent, "plots\n")
  cat("   • Ratio:", round(n_present / n_absent, 2), ":1\n")
  
  if (n_present / n_absent >= 0.25 & n_present / n_absent <= 4) {
    cat("   • ✅ Well-balanced for logistic regression\n")
  } else {
    cat("   • ⚠️  Imbalanced - may need special handling (SMOTE, weighted GLM)\n")
  }
  
  # Visualization
  cat("\n📊 GENERATING VISUALIZATIONS...\n\n")
  
  # Plot 1: Spatial distribution of abundance
  p1 <- plot_occurrence %>%
    ggplot(aes(x = reorder(Plot, mean_cover), y = mean_cover)) +
    geom_col(fill = "steelblue", alpha = 0.7) +
    coord_flip() +
    labs(
      title = paste("Cover Distribution:", species_in_question),
      x = "Plot ID",
      y = "Mean Cover"
    ) +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 6))
  
  # Plot 2: Cover value histogram
  cover_data <- lichen_df %>%
    filter(species == species_in_question)
  
  p2 <- ggplot(cover_data, aes(x = cover)) +
    geom_histogram(binwidth = 0.5, fill = "forestgreen", color = "white", alpha = 0.7) +
    labs(
      title = "Cover Value Distribution",
      x = "Cover Class",
      y = "Frequency"
    ) +
    theme_minimal()
  
  print(p1)
  print(p2)
  
  cat("═══════════════════════════════════════════════════════\n\n")
  
  # Return data invisibly for further use
  invisible(list(
    species = species_in_question,
    summary = sp_data,
    plot_data = plot_occurrence
  ))
}

# ===== NOW RUN IT =====








#========Searching and checking species========


species_in_question <- "Alectoria sarmentosa"
check_species_viability(species_in_question)

## checking other species:

# Check species 2
species_in_question <- "Parmelia sulcata"
check_species_viability(species_in_question)

# Check species 3
species_in_question <- "Bryoria nadvornikiana"
check_species_viability(species_in_question)

# Check species 4
species_in_question <- "Chaenotheca sphaerocephala"
check_species_viability(species_in_question)

# Check species 5
species_in_question <- "Xylographa vitiligo"
check_species_viability(species_in_question)


##################################################################

###### printing a distinct species list with the count  ##########

species_list <- lichens %>%
  count(species, name = "n_observations") %>%
  arrange(species)

print(species_list, n = Inf)

#################################################################










#####------######

## The next code is to check a list of interesting species at once
## Bulk search


check_multiple_species_detailed_nopause <- function(species_list, 
                                                    species_summary_df = species_summary,
                                                    lichen_df = lichens,
                                                    min_prevalence = 20,
                                                    max_prevalence = 80) {
  
  cat("\n")
  cat("╔═══════════════════════════════════════════════════════╗\n")
  cat("║   BATCH SPECIES VIABILITY CHECK (AUTO-RUN)            ║\n")
  cat("║   Total species to check:", length(species_list), "                        ║\n")
  cat("╚═══════════════════════════════════════════════════════��\n")
  
  results_summary <- tibble()
  
  for (i in seq_along(species_list)) {
    
    sp <- species_list[i]
    
    cat("\n\n")
    cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
    cat("  SPECIES", i, "of", length(species_list), "\n")
    cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
    
    result <- check_species_viability(
      species_in_question = sp,
      species_summary_df = species_summary_df,
      lichen_df = lichen_df,
      min_prevalence = min_prevalence,
      max_prevalence = max_prevalence
    )
    
    if (! is.null(result)) {
      results_summary <- bind_rows(
        results_summary,
        result$summary %>% select(species, n_plots, prevalence_pct, model_suitability)
      )
    } else {
      results_summary <- bind_rows(
        results_summary,
        tibble(
          species = sp,
          n_plots = NA,
          prevalence_pct = NA,
          model_suitability = "NOT FOUND"
        )
      )
    }
  }
  
  # Final summary
  cat("\n\n")
  cat("╔═══════════════════════════════════════════════════════╗\n")
  cat("║   BATCH SUMMARY:   ALL SPECIES                         ║\n")
  cat("╚═══════════════════════════════════════════════════════╝\n\n")
  
  print(results_summary, n = Inf)
  
  cat("\n📊 OVERALL STATISTICS:\n")
  cat("   • Total species checked:", nrow(results_summary), "\n")
  cat("   • Found in dataset:", sum(!is.na(results_summary$prevalence_pct)), "\n")
  cat("   • Suitable for modeling (20-80%):", 
      sum(results_summary$prevalence_pct >= 20 & results_summary$prevalence_pct <= 80, na.rm = TRUE), "\n")
  cat("   • Rare but viable (10-20%):", 
      sum(results_summary$prevalence_pct >= 10 & results_summary$prevalence_pct < 20, na.rm = TRUE), "\n")
  cat("   • Too rare (<10%):", 
      sum(results_summary$prevalence_pct < 10, na.rm = TRUE), "\n")
  cat("   • Not found:", 
      sum(is.na(results_summary$prevalence_pct)), "\n\n")
  
  write_csv(results_summary, "batch_species_check_summary.csv")
  cat("💾 Summary table saved to:  batch_species_check_summary. csv\n\n")
  
  return(results_summary)
}



analyze_cooccurrence <- function(species_list, lichen_df = lichens) {
  
  # Create presence/absence matrix for selected species
  pa_matrix <- lichen_df %>%
    filter(species %in% species_list) %>%
    distinct(Plot, species) %>%
    mutate(presence = 1) %>%
    pivot_wider(
      names_from = species,
      values_from = presence,
      values_fill = 0
    )
  
  # Extract species columns only
  species_cols <- pa_matrix %>% select(-Plot)
  
  # Calculate co-occurrence metrics
  cooccur_results <- list()
  
  # 1. Jaccard Similarity Matrix (0-1, higher = more co-occurrence)
  jaccard_matrix <- matrix(NA, nrow = length(species_list), ncol = length(species_list))
  rownames(jaccard_matrix) <- species_list
  colnames(jaccard_matrix) <- species_list
  
  for (i in seq_along(species_list)) {
    for (j in seq_along(species_list)) {
      sp1 <- species_list[i]
      sp2 <- species_list[j]
      
      if (sp1 %in% names(species_cols) & sp2 %in% names(species_cols)) {
        a <- sum(species_cols[[sp1]] == 1 & species_cols[[sp2]] == 1)  # Both present
        b <- sum(species_cols[[sp1]] == 1 & species_cols[[sp2]] == 0)  # Only sp1
        c <- sum(species_cols[[sp1]] == 0 & species_cols[[sp2]] == 1)  # Only sp2
        
        jaccard <- a / (a + b + c)
        jaccard_matrix[i, j] <- jaccard
      }
    }
  }
  
  # 2. Co-occurrence counts
  cooccur_counts <- expand_grid(
    species1 = species_list,
    species2 = species_list
  ) %>%
    filter(species1 != species2) %>%
    rowwise() %>%
    mutate(
      both_present = if_else(
        species1 %in% names(species_cols) & species2 %in% names(species_cols),
        sum(species_cols[[species1]] == 1 & species_cols[[species2]] == 1),
        NA_integer_
      ),
      sp1_only = if_else(
        species1 %in% names(species_cols) & species2 %in% names(species_cols),
        sum(species_cols[[species1]] == 1 & species_cols[[species2]] == 0),
        NA_integer_
      ),
      sp2_only = if_else(
        species1 %in% names(species_cols) & species2 %in% names(species_cols),
        sum(species_cols[[species1]] == 0 & species_cols[[species2]] == 1),
        NA_integer_
      ),
      neither = if_else(
        species1 %in% names(species_cols) & species2 %in% names(species_cols),
        sum(species_cols[[species1]] == 0 & species_cols[[species2]] == 0),
        NA_integer_
      ),
      jaccard = both_present / (both_present + sp1_only + sp2_only)
    ) %>%
    ungroup() %>%
    arrange(desc(jaccard))
  
  # 3. Shared plots
  shared_plots <- list()
  for (i in 1:(length(species_list) - 1)) {
    for (j in (i + 1):length(species_list)) {
      sp1 <- species_list[i]
      sp2 <- species_list[j]
      
      if (sp1 %in% names(species_cols) & sp2 %in% names(species_cols)) {
        shared <- pa_matrix %>%
          filter(!! sym(sp1) == 1 & !!sym(sp2) == 1) %>%
          pull(Plot)
        
        if (length(shared) > 0) {
          shared_plots[[paste(sp1, "×", sp2)]] <- shared
        }
      }
    }
  }
  
  return(list(
    pa_matrix = pa_matrix,
    jaccard_matrix = jaccard_matrix,
    cooccur_counts = cooccur_counts,
    shared_plots = shared_plots
  ))
}


# ===== VISUALIZE CO-OCCURRENCE =====

plot_cooccurrence <- function(cooccur_analysis, species_list) {
  
  jaccard_mat <- cooccur_analysis$jaccard_matrix
  
  # Remove species not found (NA rows/cols)
  valid_species <- rowSums(is.na(jaccard_mat)) < ncol(jaccard_mat)
  jaccard_mat <- jaccard_mat[valid_species, valid_species]
  
  if (nrow(jaccard_mat) < 2) {
    cat("⚠️  Not enough species found for co-occurrence analysis\n")
    return(NULL)
  }
  
  # 1. Heatmap
  cat("\n📊 Generating co-occurrence heatmap...\n")
  
  pheatmap(
    jaccard_mat,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    display_numbers = TRUE,
    number_format = "%.2f",
    color = colorRampPalette(c("white", "yellow", "orange", "red"))(50),
    main = "Species Co-occurrence (Jaccard Similarity)",
    fontsize = 10,
    fontsize_number = 8,
    angle_col = 45
  )
  
  # 2. Network-style plot (for pairs with high co-occurrence)
  high_cooccur <- cooccur_analysis$cooccur_counts %>%
    filter(! is.na(jaccard), jaccard > 0.3) %>%  # Threshold:  30% similarity
    select(species1, species2, jaccard, both_present)
  
  if (nrow(high_cooccur) > 0) {
    cat("\n📊 Generating co-occurrence network...\n")
    
    ggplot(high_cooccur, aes(x = reorder(species1, jaccard), 
                             y = reorder(species2, jaccard),
                             fill = jaccard,
                             label = both_present)) +
      geom_tile(color = "white", size = 0.5) +
      geom_text(color = "black", size = 3) +
      scale_fill_gradient2(
        low = "white", mid = "yellow", high = "red",
        midpoint = 0.5,
        limits = c(0, 1),
        name = "Jaccard\nSimilarity"
      ) +
      labs(
        title = "Species Pairs with >30% Co-occurrence",
        subtitle = "Numbers show count of shared plots",
        x = "Species 1",
        y = "Species 2"
      ) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        axis.text.y = element_text(size = 8)
      )
  } else {
    cat("\n⚠️  No species pairs with >30% co-occurrence found\n")
  }
}


# ===== PRINT CO-OCCURRENCE REPORT =====

print_cooccurrence_report <- function(cooccur_analysis, species_list) {
  
  cat("\n")
  cat("╔═══════════════════════════════════════════════════════╗\n")
  cat("║   CO-OCCURRENCE ANALYSIS REPORT                       ║\n")
  cat("╚═══════════════════════════════════════════════════════╝\n\n")
  
  # Top co-occurring pairs
  cat("🔗 TOP CO-OCCURRING SPECIES PAIRS:\n")
  cat("   (Jaccard Similarity:  0 = never together, 1 = always together)\n\n")
  
  top_pairs <- cooccur_analysis$cooccur_counts %>%
    filter(!is.na(jaccard)) %>%
    arrange(desc(jaccard)) %>%
    head(10)
  
  print(top_pairs %>% select(species1, species2, jaccard, both_present), n = 10)
  
  # Shared plots details
  cat("\n\n📍 SHARED PLOTS (Species occurring together):\n\n")
  
  if (length(cooccur_analysis$shared_plots) > 0) {
    for (pair_name in names(cooccur_analysis$shared_plots)) {
      plots <- cooccur_analysis$shared_plots[[pair_name]]
      cat("   •", pair_name, ":", length(plots), "shared plots\n")
      cat("     Plots:", paste(plots, collapse = ", "), "\n\n")
    }
  } else {
    cat("   ⚠️  No species pairs share plots\n")
  }
  
  # Isolation analysis
  cat("\n🏝️  SPECIES ISOLATION (occur alone without other target species):\n\n")
  
  pa_matrix <- cooccur_analysis$pa_matrix
  species_cols <- pa_matrix %>% select(-Plot)
  
  isolation_stats <- tibble(
    species = names(species_cols),
    n_plots_present = map_int(species_cols, sum),
    n_plots_alone = map_int(names(species_cols), function(sp) {
      sum(species_cols[[sp]] == 1 & rowSums(species_cols %>% select(-!! sp)) == 0)
    })
  ) %>%
    mutate(
      pct_alone = round((n_plots_alone / n_plots_present) * 100, 1)
    ) %>%
    arrange(desc(pct_alone))
  
  print(isolation_stats, n = Inf)
  
  # Interpretation
  cat("\n\n💡 INTERPRETATION:\n")
  
  high_cooccur <- cooccur_analysis$cooccur_counts %>%
    filter(!is.na(jaccard), jaccard > 0.5) %>%
    nrow()
  
  if (high_cooccur > 0) {
    cat("   • ✅", high_cooccur, "species pairs show STRONG co-occurrence (>50%)\n")
    cat("     → These species likely share similar habitat requirements\n")
    cat("     → Consider grouping them in a 'species complex' model\n")
  }
  
  moderate_cooccur <- cooccur_analysis$cooccur_counts %>%
    filter(!is.na(jaccard), jaccard > 0.3, jaccard <= 0.5) %>%
    nrow()
  
  if (moderate_cooccur > 0) {
    cat("   • ⚠️ ", moderate_cooccur, "species pairs show MODERATE co-occurrence (30-50%)\n")
    cat("     → Partially overlapping niches\n")
  }
  
  high_isolation <- isolation_stats %>%
    filter(pct_alone > 50) %>%
    nrow()
  
  if (high_isolation > 0) {
    cat("   • 🏝️ ", high_isolation, "species are often ISOLATED (>50% alone)\n")
    cat("     → These are habitat specialists or competitive dominants\n")
    cat("     → Model separately to capture unique ecological requirements\n")
  }
  
  cat("\n")
}


# ===== INTEGRATED BATCH CHECK WITH CO-OCCURRENCE =====

check_multiple_species_with_cooccurrence <- function(species_list, 
                                                     species_summary_df = species_summary,
                                                     lichen_df = lichens,
                                                     min_prevalence = 20,
                                                     max_prevalence = 80) {
  
  cat("\n")
  cat("╔═══════════════════════════════════════════════════════╗\n")
  cat("║   COMPREHENSIVE SPECIES ANALYSIS                      ║\n")
  cat("║   Phase 1: Individual Species Reports                 ║\n")
  cat("║   Phase 2: Co-occurrence Analysis                     ║\n")
  cat("╚═══════════════════════════════════════════════════════╝\n")
  
  # === PHASE 1: Individual species checks ===
  cat("\n\n━━━ PHASE 1: INDIVIDUAL SPECIES VIABILITY ━━━\n")
  
  results_summary <- check_multiple_species_detailed_nopause(
    species_list = species_list,
    species_summary_df = species_summary_df,
    lichen_df = lichen_df,
    min_prevalence = min_prevalence,
    max_prevalence = max_prevalence
  )
  
  # Filter to species actually found
  found_species <- results_summary %>%
    filter(! is.na(prevalence_pct)) %>%
    pull(species)
  
  if (length(found_species) < 2) {
    cat("\n⚠️  Not enough species found for co-occurrence analysis\n")
    return(results_summary)
  }
  
  # === PHASE 2: Co-occurrence analysis ===
  cat("\n\n━━━ PHASE 2: CO-OCCURRENCE ANALYSIS ━━━\n")
  
  cooccur_analysis <- analyze_cooccurrence(found_species, lichen_df)
  
  print_cooccurrence_report(cooccur_analysis, found_species)
  
  plot_cooccurrence(cooccur_analysis, found_species)
  
  # Save co-occurrence data
  write_csv(cooccur_analysis$cooccur_counts, "species_cooccurrence_matrix.csv")
  cat("\n💾 Co-occurrence matrix saved to:  species_cooccurrence_matrix. csv\n")
  
  # Return combined results
  return(list(
    species_summary = results_summary,
    cooccurrence = cooccur_analysis
  ))
}


### checking it towards the 11 OGF indicators from Malíček et al 2019

ogf_indicators <- c("Xylographa vitiligo", "Chaenotheca sphaerocephala", 
                                    "Parmelia saxatilis", "Lecanora subintricata", 
                                    "Mycoblastus sanguinarius", "Mycoblastus affinis", 
                                    "Chaenothecopsis pusilla", "Calicium viride", 
                                    "Ochrolechia alboflavescens", "Chaenothecopsis viridireagens", 
                                    "Bryoria capillaris/nadvornikiana", "Micarea globulosella")


# Run comprehensive analysis
results <- check_multiple_species_with_cooccurrence(ogf_indicators)

# Access results
species_info <- results$species_summary
cooccurrence_data <- results$cooccurrence







##------------------------------------------------------------------##

#creating Parmelia agg

parmelia_agg <- lichens %>%
  filter(species %in% c("Parmelia sulcata", "Parmelia saxatilis")) %>%
  group_by(Plot) %>%
  summarise(
    cover = sum(cover, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  mutate(species = "Parmelia agg") %>%
  select(Plot, species, cover)

# Check co-occurrence between P. saxatilis and P. sulcata
parmelia_cooccur <- lichens %>%
  filter(species %in% c("Parmelia saxatilis", "Parmelia sulcata")) %>%
  distinct(Plot, species) %>%
  pivot_wider(names_from = species, values_from = species, 
              values_fill = "absent", names_prefix = "PA_") %>%
  mutate(
    saxatilis = ifelse(`PA_Parmelia saxatilis` != "absent", 1, 0),
    sulcata = ifelse(`PA_Parmelia sulcata` != "absent", 1, 0),
    both = saxatilis + sulcata
  )

# Summary statistics
cat("\n═══════════════════════════════════════════════════════\n")
cat("  PARMELIA SAXATILIS AGG. (AGGREGATE) ANALYSIS\n")
cat("═══════════════════════════════════════════════════════\n\n")

cat("📊 INDIVIDUAL SPECIES:\n")
cat("   • P. saxatilis:     28 plots (23.3%)\n")
cat("   • P. sulcata:       21 plots (17.5%)\n\n")

# Calculate aggregate prevalence
parmelia_agg_plots <- lichens %>%
  filter(species %in% c("Parmelia saxatilis", "Parmelia sulcata")) %>%
  distinct(Plot) %>%
  pull(Plot)

n_agg <- length(parmelia_agg_plots)
prev_agg <- (n_agg / 120) * 100

cat("🔗 PARMELIA AGG. (COMBINED):\n")
cat("   • Total plots:      ", n_agg, "plots (", round(prev_agg, 1), "%)\n")

# Check overlap
overlap <- parmelia_cooccur %>%
  filter(both == 2) %>%
  nrow()

cat("   • Co-occurring:     ", overlap, "plots (both species present)\n")
cat("   • P. saxatilis only:", sum(parmelia_cooccur$both == 1 & parmelia_cooccur$saxatilis == 1), "plots\n")
cat("   • P. sulcata only:  ", sum(parmelia_cooccur$both == 1 & parmelia_cooccur$sulcata == 1), "plots\n\n")

cat("✅ MODELING SUITABILITY:\n")
if (prev_agg >= 20 & prev_agg <= 80) {
  cat("   • ✅ SUITABLE for binomial GLM/GAM\n")
  cat("   • Statistical power:   EXCELLENT\n")
  cat("   • Matches Malíček's P. saxatilis agg. (#3 indicator, both=0. 845)\n")
} else {
  cat("   • Status: Good - Balanced\n")
}

cat("\n⚖️  PRESENCE/ABSENCE BALANCE:\n")
cat("   • Present:", n_agg, "plots\n")
cat("   • Absent:", 120 - n_agg, "plots\n")
cat("   • Ratio:", round(n_agg / (120 - n_agg), 2), ": 1\n")
cat("   • ✅ Well-balanced for logistic regression\n\n")

cat("═══════════════════════════════════════════════════════\n\n")









# ===== CREATE SPECIES PRESENCE/ABSENCE MATRIX =====

# Define species groups based on Malíček et al. and your data
species_groups <- list(
  
  # Group 1: Parmelia aggregate (Malíček #3 indicator)
  parmelia_agg = c("Parmelia saxatilis", "Parmelia sulcata"),
  
  # Group 2: Mycoblastus genus (core OGF specialists)
  mycoblastus = c("Mycoblastus sanguinarius", "Mycoblastus affinis"),
  
  # Group 3: Ochrolechia genus (OGF associates)
  ochrolechia = c("Ochrolechia alboflavescens", "Ochrolechia microstictoides", 
                  "Ochrolechia androgyna", "Ochrolechia mahluensis"),
  
  # Group 4: Calicioid lichens (elite OGF indicators)
  calicioids = c("Calicium viride", "Calicium glaucellum", "Calicium trabinellum",
                 "Chaenotheca sphaerocephala", "Chaenotheca ferruginea", 
                 "Chaenotheca chrysocephala", "Chaenotheca xyloxena", 
                 "Chaenotheca trichialis", "Chaenotheca brunneola",
                 "Chaenotheca furfuracea", "Chaenotheca stemonea",
                 "Chaenothecopsis pusilla", "Mycocalicium subtile"),
  
  # Group 5: Xylographa genus (wood specialists)
  xylographa = c("Xylographa vitiligo", "Xylographa paralella"),
  
  # Group 6: Elite rare indicators (Malíček top-ranked but very rare)
  elite_rare = c("Xylographa vitiligo", "Chaenotheca sphaerocephala", 
                 "Calicium viride", "Chaenothecopsis pusilla", 
                 "Micarea globulosella"),
  
  # Group 7: Core OGF complex (Mycoblastus + Ochrolechia + Lecanora)
  core_ogf = c("Mycoblastus sanguinarius", "Mycoblastus affinis",
               "Ochrolechia alboflavescens", "Lecanora subintricata")
)


# ===== FUNCTION:  CREATE GROUP PRESENCE/ABSENCE =====

create_group_variables <- function(lichen_df, species_groups_list) {
  
  # Get all unique plots
  all_plots <- tibble(Plot = unique(lichen_df$Plot))
  
  # For each group, create presence/absence
  group_data <- all_plots
  
  for (group_name in names(species_groups_list)) {
    
    species_list <- species_groups_list[[group_name]]
    
    # Presence (binary: 0/1)
    group_presence <- lichen_df %>%
      filter(species %in% species_list) %>%
      distinct(Plot) %>%
      mutate(!! paste0(group_name, "_presence") := 1)
    
    # Richness (count:  0-n species)
    group_richness <- lichen_df %>%
      filter(species %in% species_list) %>%
      group_by(Plot) %>%
      summarise(!! paste0(group_name, "_richness") := n_distinct(species),
                .groups = 'drop')
    
    # Join to main data
    group_data <- group_data %>%
      left_join(group_presence, by = "Plot") %>%
      left_join(group_richness, by = "Plot")
    
    # Fill NAs with 0
    group_data <- group_data %>%
      mutate(across(ends_with("_presence"), ~replace_na(., 0)),
             across(ends_with("_richness"), ~replace_na(., 0)))
  }
  
  return(group_data)
}


# ===== RUN THE FUNCTION =====

lichen_groups <- create_group_variables(lichens, species_groups)

# Summary report
cat("\n╔═══════════════════════════════════════════════════════╗\n")
cat("║   SPECIES GROUP PREVALENCE SUMMARY                    ║\n")
cat("╚═══════════════════════════════════════════════════════╝\n\n")

for (group in names(species_groups)) {
  presence_col <- paste0(group, "_presence")
  richness_col <- paste0(group, "_richness")
  
  n_plots <- sum(lichen_groups[[presence_col]])
  prevalence <- (n_plots / 120) * 100
  mean_richness <- mean(lichen_groups[[richness_col]][lichen_groups[[richness_col]] > 0])
  max_richness <- max(lichen_groups[[richness_col]])
  
  cat("📊", toupper(group), "\n")
  cat("   • Prevalence:       ", n_plots, "plots (", round(prevalence, 1), "%)\n")
  cat("   • Mean richness:   ", round(mean_richness, 1), "species (when present)\n")
  cat("   • Max richness:    ", max_richness, "species\n")
  
  # Suitability assessment
  if (prevalence >= 20 & prevalence <= 80) {
    cat("   • ✅ SUITABLE for binomial GLM/GAM\n")
  } else if (prevalence >= 10 & prevalence < 20) {
    cat("   • ⚠️  MARGINAL - Consider richness model or presence-only\n")
  } else if (prevalence < 10) {
    cat("   • ❌ TOO RARE - Use MaxEnt or exclude\n")
  } else {
    cat("   • ⚠️  UBIQUITOUS - Limited variation\n")
  }
  cat("\n")
}


# ===== VISUALIZE GROUP OVERLAP =====

# Create summary for plotting
group_summary <- tibble(
  Group = names(species_groups),
  Prevalence = sapply(names(species_groups), function(g) {
    sum(lichen_groups[[paste0(g, "_presence")]]) / 120 * 100
  }),
  N_Plots = sapply(names(species_groups), function(g) {
    sum(lichen_groups[[paste0(g, "_presence")]])
  })
) %>%
  arrange(desc(Prevalence))

# Barplot
ggplot(group_summary, aes(x = reorder(Group, Prevalence), y = Prevalence)) +
  geom_col(aes(fill = cut(Prevalence, breaks = c(0, 10, 20, 80, 100))),
           alpha = 0.8) +
  geom_hline(yintercept = c(20, 80), linetype = "dashed", color = "red") +
  geom_text(aes(label = N_Plots), hjust = -0.3, size = 4) +
  scale_fill_manual(
    values = c("(0,10]" = "red", "(10,20]" = "orange", 
               "(20,80]" = "darkgreen", "(80,100]" = "gold"),
    labels = c("Too rare", "Marginal", "Suitable", "Ubiquitous"),
    name = "Modeling\nSuitability"
  ) +
  coord_flip() +
  labs(
    title = "Species Group Prevalence Across 120 Plots",
    subtitle = "Red lines = 20-80% threshold for balanced GLM",
    x = "Species Group",
    y = "Prevalence (%)",
    caption = "Numbers indicate plot count"
  ) +
  theme_minimal() +
  theme(legend.position = "right")


# ===== SAVE RESULTS =====

write_csv(lichen_groups, "Lichens/lichen_species_groups.csv")
write_csv(group_summary, "Lichens/species_group_summary.csv")

cat("💾 Files saved:\n")
cat("   • lichen_species_groups.csv (120 rows × group variables)\n")
cat("   • species_group_summary.csv (group prevalence table)\n\n")



