#Monitoringa atskaitei datu sagatavošana (sugu frequency utt un analīze)


library(readxl)
library(dplyr)
library(tidyr)
library(purrr)



setwd("D:/Users/solvita/Documents/ABC_Solvita/R_ZP/ZP") 
# excel contains all 3 year data for only experts (some data omited due to nepilnīga laika rinda)
#read.csv2 see semicolons as separators

#prepare file for calculating cumulative species number by selecting only appropriate rows ----
exp <- read.csv2 ("exp3g_atsk.csv")
str(exp)
nrow(exp)


#unique_values <- unique(farm_exp$Rel_type)
#print(unique_values)
#"basic_ZP"     "basic_ZP_Pia" "basic_ZP_PIa" "detailed_ZP"  "farmer_ZP"   
#use only basic_ZP and farmer_ZP, others omit from file

# keep only those rows
#farm_exp_b <- subset(farm_exp, Rel_type %in% c("farmer_ZP", "basic_ZP"))


# export to a new csv
#write.csv2(farm_exp_c, "farm_exp_c.csv", row.names = FALSE)
#str(farm_exp_c)
#nrow(farm_exp_c)


# calculating  species number per plot (priekš vid.sugu skaita ----
library(readr)
species_cols <- names(exp)[57:106]

exp1 <- exp %>%
  # 1) Count species with value >= 1 in each subplot (row)
  mutate(
    spp_richness_subplot = rowSums(
      across(all_of(species_cols), ~ suppressWarnings(as.numeric(.x)) >= 1),
      na.rm = TRUE
    )
  ) %>%
  # 2) Compute average richness across the 10 subplots within each plot (Uniq_Rprog)
  group_by(Uniq_Rprog) %>%
  mutate(
    avg_spp_richness_10subplots = mean(spp_richness_subplot, na.rm = TRUE)
  ) %>%
  ungroup()

write.csv2(exp1, "exp3g_atsk_with_avg_richness.csv")


# calculating frequency of each species  per parcel ----
exp1 <- read.csv2 ("exp3g_atsk_with_avg_richness.csv")
str(exp1)
nrow(exp1)


species_cols <- names(exp1)[58:107]
species_names <- names(exp1)[species_cols]

exp2 <- exp1 %>%
  # make sure species columns are numeric (safe if already numeric)
  mutate(across(all_of(species_cols), ~ suppressWarnings(as.numeric(.x)))) %>%
  group_by(Uniq_Rprog) %>%
  # Frequency across the 10 subplots within each plot:
  # - proportion of subplots where species is present (value >= 1)
  mutate(
    across(
      all_of(species_cols),
      ~ mean(.x >= 1, na.rm = TRUE),
      .names = "{.col}_freq10"
    ),
    # (optional) also add counts out of 10 subplots
    across(
      all_of(species_cols),
      ~ sum(.x >= 1, na.rm = TRUE),
      .names = "{.col}_n10"
    )
  ) %>%
  ungroup()

write.csv2(exp2, "exp3g_atsk_richness_freq.csv")

#calculating cover-weighted frequency for each species per parcel ----

exp2 <- read.csv2 ("exp3g_atsk_richness_freq.csv")
str(exp2)

species_cols <- names(exp2)[59:108]

exp3 <- exp2 %>%
  mutate(across(all_of(species_cols), ~ suppressWarnings(as.numeric(.x)))) %>%
  group_by(Uniq_Rprog) %>%
  mutate(
    # Frequency F (proportion of subplots with cover > 0)
    across(
      all_of(species_cols),
      ~ mean(.x > 0, na.rm = TRUE),
      .names = "{.col}_F10"
    ),
    
    # Mean cover when present C_pos (0 if never present)
    across(
      all_of(species_cols),
      ~ {
        x <- .x
        if (sum(x > 0, na.rm = TRUE) == 0) 0 else mean(x[x > 0], na.rm = TRUE)
      },
      .names = "{.col}_Cpos"
    ),
    
    # Weighted frequency WF = F * C_pos (0 if never present)
    across(
      all_of(species_cols),
      ~ {
        x <- .x
        F <- mean(x > 0, na.rm = TRUE)
        Cpos <- if (sum(x > 0, na.rm = TRUE) == 0) 0 else mean(x[x > 0], na.rm = TRUE)
        F * Cpos
      },
      .names = "{.col}_WF10"
    )
  ) %>%
  ungroup()
write.csv2(exp3, "exp3g_atsk_richness_freq_weightfr.csv")

#Notes:
  
#  *_F10 is 0–1.

#*_Cpos is mean cover among subplots where the species is present; it becomes NA if the species is absent in all 10 subplots.

#*_WF10 is F10 * Cpos (also NA when absent everywhere). If you prefer 0 instead of NA for “absent everywhere”, tell me and I’ll adjust one line.

# calculating cumulative species number per parcel ----

#calculate cumulative species number per parcel and add value into the new column 'cumul_spno'

exp3 <- read.csv2 ("exp3g_atsk_richness_freq_weightfr.csv")
str(exp3)

species_cols <- names(exp3)[60:109]

# 1) Calculate cumulative species number per plot (Uniq_Rprog)
cumul_tbl <- exp3 %>%
  group_by(Uniq_Rprog) %>%
  summarise(
    Cumul_spno = sum(sapply(pick(all_of(species_cols)), function(x) {
      x <- suppressWarnings(as.numeric(x))   # coerce safely if needed
      any(x > 0, na.rm = TRUE)               # present in at least one of 10 subplots?
    })),
    .groups = "drop"
  )

# 2) Add it back to every subplot row (original cover columns remain unchanged)
exp4 <- exp3 %>%
  left_join(cumul_tbl, by = "Uniq_Rprog")


write.csv2(exp4, "exp3g_atsk_richness_freq_weightfr_cumul.csv")

# calulate index where a plot scores “better” when it has many species and those species collectively have high cover (i.e., not just many tiny traces).“how many equally-covered species” times “how much cover there is”. This tends to reward plots where cover is spread across multiple species (not just one dominant).

exp4 <- read.csv2 ("exp3g_atsk_richness_freq_weightfr_cumul.csv")


species_cols <- names(exp4)[61:110]

# 1) Make a numeric + NA->0 copy for calculations ONLY (original df stays unchanged)
df_calc <- exp4 %>%
  mutate(across(all_of(species_cols), ~ suppressWarnings(as.numeric(.x)))) %>%
  mutate(across(all_of(species_cols), ~ replace_na(.x, 0)))

# 2) Aggregate species cover across the 10 subplots per plot (Uniq_Rprog)
plot_species_sum <- df_calc %>%
  group_by(Uniq_Rprog) %>%
  summarise(across(all_of(species_cols), ~ sum(.x, na.rm = TRUE)), .groups = "drop")

# 3) Compute Hill numbers (q=0,1,2), total cover, and a score (Hill_q1 * total cover)
plot_indices <- plot_species_sum %>%
  rowwise() %>%
  mutate(
    Total_cover = sum(c_across(all_of(species_cols))),
    Hill_q0 = sum(c_across(all_of(species_cols)) > 0),                 # richness
    Hill_q1 = if (Total_cover == 0) 0 else {                           # Shannon effective species
      p <- c_across(all_of(species_cols)) / Total_cover
      exp(-sum(p[p > 0] * log(p[p > 0])))
    },
    Hill_q2 = if (Total_cover == 0) 0 else {                           # Simpson effective species
      p <- c_across(all_of(species_cols)) / Total_cover
      1 / sum(p[p > 0]^2)
    },
    Score = Hill_q1 * Total_cover                                      # “many well-covered species” * cover
  ) %>%
  ungroup() %>%
  select(Uniq_Rprog, Hill_q0, Hill_q1, Hill_q2, Total_cover, Score)

# 4) Join indices back to original data (species cover columns remain unchanged)
exp5 <- exp4 %>%
  left_join(plot_indices, by = "Uniq_Rprog")

write.csv2(exp5, "exp3g_atsk_richness_freq_weightfr_cumul_weightcumul.csv")

#result is 
#Total_cover What it is: the sum of cover scores (1/2/3) across all species and all 10 subplots in that plot. Interpretation: a proxy for “how much vegetation cover / abundance” (in your ordinal units). Range: depends on how many species are recorded and how often; with 50 species × 10 subplots × max cover 3, the theoretical max would be 1500, but in practice much lower.

#Hill_q0 (q = 0). What it is: cumulative species richness across the 10 subplots. How computed: counts how many species have summed cover > 0 in the plot (present in at least one subplot). Interpretation: “how many species occur in this plot at least somewhere.” Range: 0 to number of species columns (here up to 50).

#Hill_q1 (q = 1, Shannon effective species). What it is: the effective number of species based on how evenly cover is distributed among species.How computed: convert each species’ plot-level cover sum to a proportion. Interpretation: “how many equally-covered species would give the same diversity.” Behavior:  high when many species share cover fairly evenly lower when cover is dominated by a few species, even if richness is high Range: 0 to Hill_q0 (it cannot exceed richness).

#Hill_q2 (q = 2, Simpson effective species). What it is: another effective species number, but more sensitive to dominance.Interpretation: “effective number of dominant species.”Behavior: decreases strongly if a few species have most of the cover. Range: 0 to Hill_q1 to Hill_q0 (typically Hill_q2 ≤ Hill_q1 ≤ Hill_q0).

#Score What it is: your combined “better = many well-covered species”. Interpretation: high when a plot has high total cover and that cover is spread across multiple species (since Hill_q1 drops when one/few species dominate). Note: it’s in “cover-units × effective species”; best used for relative comparisons among your plots, not as an absolute ecological quantity.

