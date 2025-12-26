#farmers expert scoring difference
#plus changes between years 2023 and 2024
#Wilcoxon and Hedges-Lehman
#Heatmap plot
library(readxl)
library(dplyr)
library(tidyr)
library(purrr)



setwd("D:/Users/solvita/Documents/ABC_Solvita/R_ZP/ZP") 
# excel contains all data from 2023 and 2024 and expert data for scores and species basic monit for 2025 (without structure data). Expert group for species counts grouped a bit differently than for Thesis (more experts in exp1 group)
#read.csv2 see semicolons as separators

#prepare file for calculating cumulative specie snumber by selecting only appropriate rows ----
farm_exp <- read.csv2 ("SR_farmerexp_preparation.csv")
str(farm_exp)
nrow(farm_exp)


unique_values <- unique(farm_exp$Rel_type)
print(unique_values)
#"basic_ZP"     "basic_ZP_Pia" "basic_ZP_PIa" "detailed_ZP"  "farmer_ZP"   
#use only basic_ZP and farmer_ZP, others omit from file

# keep only those rows
farm_exp_b <- subset(farm_exp, Rel_type %in% c("farmer_ZP", "basic_ZP"))

# keep only  rows where there is a value in scoring field 'Gr_Qualit_exp' (omit rows without monitoring results for species - excluded or stopped to participate) - value 'na' in the field 'Gr_Qualit_exp'
unique_values <- unique(farm_exp_b$Gr_qualit_exp)
print(unique_values)

# keep only those rows
farm_exp_c <- subset(farm_exp_b, Gr_qualit_exp %in% c("average",   "bad",       "good", "very_bad",  "excellent"))

# export to a new csv
write.csv2(farm_exp_c, "farm_exp_c.csv", row.names = FALSE)
str(farm_exp_c)
nrow(farm_exp_c)

# calculating cumulative species number per parcel ----

#calculate cumulative species number per parcel and add value into the new column 'cumul_spno'

# Define the species columns (from column 33 to 82)
species_cols <- names(farm_exp_c)[33:82]

# Convert cover values to presence/absence (1 if cover value is 1, 2 or 3; else 0)
farm_exp_c[species_cols] <- lapply(farm_exp_c[species_cols], function(x) {
  as.integer(as.character(x) %in% c("1", "2", "3"))
})

# Calculate cumulative species number per Uniq_no group
cumul_species <- farm_exp_c %>%
  group_by(Uniq_no) %>%
  summarise(cumul_spno = sum(colSums(across(all_of(species_cols))) > 0), .groups = "drop")

# Merge cumulative species count back to the main dataframe
farm_exp_c_result <- left_join(farm_exp_c, cumul_species, by = "Uniq_no")

# Save the result to a new CSV file
write.csv2(farm_exp_c_result, "farm_exp_c_result.csv", row.names = FALSE)

#end of species number calculation

# omit parcels not usable in farmer-expert comparison -  all 2025 data

# keep only needed rows
farm_exp_c_no25 <- subset(farm_exp_c_result, Year %in% c("2023", "2024"))
str(farm_exp_c_no25)
nrow(farm_exp_c_no25)

# omit parcels not usable in farmer-expert comparison -  controls, stopped and volunteers

# keep only needed rows
farm_exp_c_nocontr <- subset(farm_exp_c_no25, Status %in% c("Participant"))
str(farm_exp_c_nocontr)
nrow(farm_exp_c_nocontr)

# omit parcels not usable in farmer-expert comparison -  exp3 (exclude intern and 1st year experieced expert data)

# keep only needed rows
farm_exp_c_noexp3 <- subset(farm_exp_c_nocontr, Exp_gr %in% c("exp1", "farmer"))
str(farm_exp_c_noexp3)
nrow(farm_exp_c_noexp3)



# export to a new csv
write.csv(farm_exp_sub, "farm_exp_or_basic_ZP.csv", row.names = FALSE)










#1.2.1.	Are farmers consistently resulting in less/more points per polygon than experts do?
#Use real species scores – excel column Spp_score and not Spp_scor50
#pairs excpert1, farmer, each year seperatly
#shapiro
#wilcoxon

unique_values <- unique(scoring_data$Rel_type)
print(unique_values)
#unse only basic_ZP

scoring_data <- scoring_data %>%
  filter(Rel_type %in% c("basic_ZP", "farmer_ZP"))


#use only vegetation condition normal, check difference in dataset size afzerwards

unique_values <- unique(scoring_data$Veg_cond_sampling)
print(unique_values)




#use only basic_ZP
scoring_data <- scoring_data %>%
  filter(Veg_cond_sampling == "normal")
#8750 rows to 7160 rows

unique_values <- unique(scoring_data$Eksp_group)
print(unique_values)
#use only those who are farmer and Rusina as well as Medne and set last ones (Rusina and Medne) to exp1

scoring_data <- scoring_data %>%
  filter(Eksp_group %in% c("farmer", "Rusina", "Medne"))

scoring_data <- scoring_data %>%
  mutate(Eksp_group = ifelse(Eksp_group %in% c("Rusina", "Medne"), "exp1", Eksp_group))

# Kontrolle:
table(scoring_data$Eksp_group)

#use only 1st mont pnt
scoring_data <- scoring_data %>%
  filter(Mont_pnt == "1")

str(scoring_data)
#nur noch Paare 
filtered_paired <- scoring_data %>%
  filter(Eksp_group %in% c("farmer", "exp1")) %>%
  group_by(Uniq_Rprog) %>%
  filter(all(c("farmer", "exp1") %in% Eksp_group)) %>%
  ungroup()

# Kontrolle:
table(filtered_paired$Eksp_group)
length(unique(filtered_paired$Uniq_Rprog)) # Gibt die Anzahl der Paare aus
#185 Paare

#SR. export_Tabele_to_excel
library(openxlsx)
write.xlsx(filtered_paired, file = "filtered_paired.xlsx", rowNames = FALSE)

#shapiro
scoring_data$Spp_score <- as.numeric(as.character(scoring_data$Spp_score))
unique_values <- unique(scoring_data$Spp_score)
print(unique_values)
table(is.na(scoring_data$Spp_score))


normality_check <- scoring_data %>%
  filter(Eksp_group %in% c("exp1", "farmer")) %>%
  group_by(Eksp_group) %>%
  summarise(
    n = sum(!is.na(Spp_score)),
    W = if (n() >= 3) shapiro.test(Spp_score)$statistic else NA,
    p_value = if (n() >= 3) shapiro.test(Spp_score)$p.value else NA,
    .groups = "drop"
  )

print(normality_check)
#no normal distribution


#wilcoxon for pairs
# Preparation: Wide table with one row per pair
paired_data <- scoring_data %>%
  filter(Eksp_group %in% c("exp1", "farmer")) %>%
  select(Uniq_Rprog, Eksp_group, Spp_score) %>%
  tidyr::pivot_wider(names_from = Eksp_group, values_from = Spp_score)

# Oder (robuster, bei nicht-Normalverteilung):
wilcox.test(paired_data$exp1, paired_data$farmer, paired = TRUE)

#There is no statistically significant difference between the two groups with regard to the median values of


# Assuming your paired_data has columns 'farmer' and 'exp1' with scores per polygon
total_polygons <- nrow(paired_data)

# Count how many polygons have farmer > exp1
farmer_greater <- sum(paired_data$farmer > paired_data$exp1, na.rm = TRUE)

# Count how many polygons have farmer < exp1
farmer_lesser <- sum(paired_data$farmer < paired_data$exp1, na.rm = TRUE)

percent_farmer_greater <- (farmer_greater / total_polygons) * 100
percent_farmer_lesser <- (farmer_lesser / total_polygons) * 100

cat(sprintf("Percentage polygons where Farmer > Expert: %.2f%%\n", percent_farmer_greater))
cat(sprintf("Percentage polygons where Farmer < Expert: %.2f%%\n", percent_farmer_lesser))

#use ewith confodence inetrval
result <- wilcox.test(paired_data$farmer, paired_data$exp1, 
                      paired = TRUE, conf.int = TRUE, conf.level = 0.90)

print(result$conf.int)    # 90% confidence interval for the median difference (Hodges-Lehmann)
print(result$estimate)    # Hodges-Lehmann median difference (the point estimate)

# Absolute difference vector
abs_diff <- abs(paired_data$farmer - paired_data$exp1)

# Median absolute difference
median_abs_diff <- median(abs_diff, na.rm = TRUE)

# Percentage within ±1 point
percent_within_1 <- mean(abs_diff <= 1, na.rm = TRUE) * 100

cat(sprintf("Median |difference| = %.2f; %.1f%% of polygons within ±1 point.\n", 
            median_abs_diff, percent_within_1))







#seperate for years (CI, wilcoxon etc.)
# If the year is always the first two characters, use substr:
paired_data <- paired_data %>%
  mutate(Year = substr(Uniq_Rprog, 1, 2)) %>%
  mutate(Year = as.integer(Year) + 2000)  # '23' → 2023, '24' → 2024

# Kontrolle:
table(paired_data$Year)


# Jahr extrahieren und ergänzen
paired_data <- paired_data %>%
  mutate(Year = as.integer(substr(Uniq_Rprog, 1, 2)) + 2000)

# Jahre extrahieren für Schleife oder Einzel-Analyse
for (y in c(2023, 2024)) {
  cat("\n==== Ergebnisse für Jahr", y, "====\n")
  
  dat <- paired_data %>% filter(Year == y)
  
  # Wilcoxon-Test & Hodges-Lehmann
  result <- wilcox.test(dat$farmer, dat$exp1, paired = TRUE, conf.int = TRUE, conf.level = 0.90)
  
  cat("Wilcoxon-Test:\n")
  print(result)
  cat("\nHodges-Lehmann 90%-CI:", 
      sprintf("[%.2f, %.2f]", result$conf.int[1], result$conf.int[2]), "\n")
  cat("Point estimate (median diff):", round(as.numeric(result$estimate), 3), "\n")
  
  # Median |difference|
  abs_diff <- abs(dat$farmer - dat$exp1)
  med_abs <- median(abs_diff, na.rm = TRUE)
  pct_within1 <- mean(abs_diff <= 1, na.rm = TRUE) * 100
  
  cat("Median |difference|:", round(med_abs, 2), "\n")
  cat(sprintf("%% polygons within ±1 point: %.1f%%\n", pct_within1))
  
  # Prozent Farmer > und < Expert
  pct_farmer_greater <- mean(dat$farmer > dat$exp1, na.rm = TRUE) * 100
  pct_farmer_lesser <- mean(dat$farmer < dat$exp1, na.rm = TRUE) * 100
  cat(sprintf("%% Farmer > Expert: %.1f%%\n", pct_farmer_greater))
  cat(sprintf("%% Farmer < Expert: %.1f%%\n", pct_farmer_lesser))
}








#Wilcoxon for the difference between the differences between the two years
# Abbreviation: Polygon ID without year prefix (e.g. "23_987_1" → "987_1")
paired_data <- paired_data %>%
  mutate(Poly_ID = substring(Uniq_Rprog, 3))

# Unterschiedsvektor Farmer-Expert für jedes Jahr berechnen
paired_data <- paired_data %>%
  mutate(diff = farmer - exp1)

# Finde nur die Poly_IDs, die in beiden Jahren vorhanden sind
ids_both_years <- paired_data %>%
  group_by(Poly_ID) %>%
  summarise(n_years = n_distinct(Year)) %>%
  filter(n_years == 2) %>%
  pull(Poly_ID)

# Filtere den Datensatz auf diese gemeinsamen Poly_IDs
paired_data_both <- paired_data %>%
  filter(Poly_ID %in% ids_both_years)

# Erzeuge breite Tabelle: jede Zeile eine Poly_ID, je eine Diff-Spalte pro Jahr
diff_wide <- paired_data_both %>%
  select(Poly_ID, Year, diff) %>%
  pivot_wider(names_from = Year, values_from = diff, names_prefix = "diff_")

# Jetzt gibt es diff_2023 und diff_2024 pro Poly_ID

result <- wilcox.test(diff_wide$diff_2023, diff_wide$diff_2024, paired = TRUE, conf.int = TRUE, conf.level = 0.90)
print(result)

paired_data_both %>%
  group_by(Poly_ID, Year) %>%
  summarise(count = n())

summary(diff_wide)

#visualisierung
valid <- !is.na(diff_wide$diff_2023) & !is.na(diff_wide$diff_2024)
result <- wilcox.test(diff_wide$diff_2023[valid], diff_wide$diff_2024[valid], paired = TRUE, conf.int = TRUE, conf.level = 0.90)

paired_data_both %>%
  group_by(Poly_ID, Year) %>%
  summarise(count = n())




#	Do  Farmers detect individual species the same as experts – individual species occurence per polygon?

curves_monitoring_head_experts <- read_excel ("Species_R_aprekiniem_1dala_15032025_v3.xlsx")
str(curves_monitoring_head_experts)

#Filter only those Polygons recorded before mowing, farmer and exp1 in basic monitoring only
curves_experts_cond_normal <- subset(curves_monitoring_head_experts, Veg_cond_sampling == "normal")
print(curves_experts_cond_normal)
str(curves_experts_cond_normal)

unique_values <- unique(curves_experts_cond_normal$Rel_type)
print(unique_values)
df_filtered <- curves_experts_cond_normal %>%
  filter(Rel_type %in% c("basic_ZP", "farmer_ZP"), Mont_pnt == 1)

unique_values <- unique(df_filtered$Exp_fact)
print(unique_values)
table(df_filtered$Exp_fact)

df_filtered <- df_filtered %>%
  filter(Exp_fact %in% c("farmer", "exp1"))
unique_values <- unique(df_filtered$Exp_fact)
print(unique_values)

#only keep pairs

#nur noch Paare 
df_filtered <- df_filtered %>%
  filter(Exp_fact %in% c("farmer", "exp1")) %>%
  group_by(Uniq_rprog) %>%
  filter(all(c("farmer", "exp1") %in% Exp_fact)) %>%
  ungroup()

# Kontrolle:
table(df_filtered$Exp_fact)
length(unique(df_filtered$Uniq_rprog)) # Gibt die Anzahl der Paare aus
#185 Paare

str(df_filtered)
which(colnames(df_filtered) == "d19_Agrimoniaeupatoria")
which(colnames(df_filtered) == "d19_Viscariavulgaris")    

df_filtered[, 36:85] <-df_filtered[, 36:85] / 10




#shapiro
spalten <- names(df_filtered)[36:85]

shapiro_grouped <- df_filtered %>%
  select(Exp_fact, all_of(spalten)) %>%
  pivot_longer(cols = -Exp_fact, names_to = "Variable", values_to = "Wert") %>%
  group_by(Exp_fact, Variable) %>%
  filter(n_distinct(Wert) > 1) %>%                # nur Variablen mit mehr als einem eindeutigen Wert
  summarise(
    n = sum(!is.na(Wert)),
    W = if (n >= 3) shapiro.test(Wert)$statistic else NA,
    p_value = if (n >= 3) shapiro.test(Wert)$p.value else NA,
    .groups = "drop"
  )

print(shapiro_grouped)
#nichts int normalverteilt.
str(df_filtered)





#nun wilcoxon paarweise nach uniq_rprog
#wide <- paired_data %>%
#  select(Uniq_Rprog, Exp_fact, starts_with("d19_")) %>%
#  pivot_wider(
#    id_cols = Uniq_rprog,
#    names_from = Exp_fact,
#    values_from = starts_with("d19_"),
#    names_glue = "{.value}_{Exp_fact}"
#  )
#spalten_basis <- unique(gsub("_(exp1|farmer)$", "", 
#                             grep("_(exp1|farmer)$", names(wide), value = TRUE)
#))

# Annahme: df_filtered enthält die Spalte "Uniq_rprog", "Exp_fact" und die Art-Spalten (spalten)
spalten <- names(df_filtered)[36:85]

wilcoxon_results <- map_dfr(
  spalten,
  function(var) {
    # Datensatz für relevanten Artspalte aus beiden Gruppen:
    tab <- df_filtered %>%
      select(Uniq_rprog, Exp_fact, !!var) %>%
      filter(Exp_fact %in% c("farmer", "exp1"))
    
    # Pivotieren: Jeweils pro Uniq_rprog die Werte für beide Gruppen in eigene Spalten
    wide_tab <- tab %>%
      pivot_wider(names_from = Exp_fact, values_from = !!var)
    
    v_exp1   <- wide_tab$exp1
    v_farmer <- wide_tab$farmer
    
    n_exp1   <- sum(!is.na(v_exp1) & v_exp1 != 0)
    n_farmer <- sum(!is.na(v_farmer) & v_farmer != 0)
    valid <- !is.na(v_exp1) & !is.na(v_farmer)
    exp1_vec   <- v_exp1[valid]
    farmer_vec <- v_farmer[valid]
    diff_vec   <- farmer_vec - exp1_vec
    abs_diff   <- abs(diff_vec)
    n_pairs    <- sum(valid)
    
    if (length(exp1_vec) < 3 | all(diff_vec == 0 | is.na(diff_vec))) {
      return(tibble(
        Variable = var, n_exp1 = n_exp1, n_farmer = n_farmer,
        n_pairs = n_pairs, V = NA, p_value = NA, HL_CI_low = NA, HL_CI_high = NA,
        HL_est = NA, median_abs_diff = NA, pct_within_1 = NA,
        pct_farmer_gt_expert = NA, pct_farmer_lt_expert = NA
      ))
    }
    wt <- wilcox.test(farmer_vec, exp1_vec, paired = TRUE, conf.int = TRUE, conf.level = 0.90, exact = FALSE)
    tibble(
      Variable = var,
      n_exp1 = n_exp1,
      n_farmer = n_farmer,
      n_pairs = n_pairs,
      V = wt$statistic,
      p_value = wt$p.value,
      HL_CI_low = wt$conf.int[1],
      HL_CI_high = wt$conf.int[2],
      HL_est = as.numeric(wt$estimate),
      median_abs_diff = median(abs_diff, na.rm = TRUE),
      pct_within_1 = mean(abs_diff <= 1, na.rm = TRUE) * 100,
      pct_farmer_gt_expert = mean(diff_vec > 0, na.rm = TRUE) * 100,
      pct_farmer_lt_expert = mean(diff_vec < 0, na.rm = TRUE) * 100
    )
  }
)
wilcoxon_results <- wilcoxon_results %>%
  mutate(
    CI_in_1 = HL_CI_low >= -1 & HL_CI_high <= 1,
    CI_in_2 = HL_CI_low >= -2 & HL_CI_high <= 2
  )

print(wilcoxon_results)

#write.csv(wilcoxon_results, "wilcoxon_CI_IDspecies_farmers_experts.csv", row.names = FALSE)
# Alle Variablen, bei denen das 90%-CI innerhalb [-1, 1] liegt:
wilcoxon_results %>% 
  filter(HL_CI_low >= -5, HL_CI_high <= 5)

# Alle Variablen, bei denen das 90%-CI innerhalb [-2, 2] liegt:
wilcoxon_results %>%
  filter(HL_CI_low >= -2, HL_CI_high <= 2)

# Alle Variablen, bei denen das 90%-CI innerhalb [-2, 2] liegt:
wilcoxon_results %>%
  filter(HL_CI_low >= -1, HL_CI_high <= 1)




#visualisation
#afterwards species changes in years
#Manhattan Plot in viridis
library(ggplot2)
library(viridis)

wilcoxon_results <- wilcoxon_results %>%
  mutate(
    CI_band = case_when(
      HL_CI_low >= -1 & HL_CI_high <= 1 ~ "+/-1",
      HL_CI_low >= -2 & HL_CI_high <= 2 ~ "+/-2",
      HL_CI_low >= -3 & HL_CI_high <= 3 ~ "+/-3",
      TRUE ~ "mehr"
    ),
    CI_band = factor(CI_band, levels = c("+/-1", "+/-2", "+/-3", "mehr"))
  )

ggplot(wilcoxon_results, aes(x = reorder(Variable, HL_est), y = HL_est, 
                             color = CI_band, size = -log10(p_value))) +
  geom_point(alpha = 0.8) +
  scale_color_viridis_d(option = "D") +
  labs(x = "Art", y = "Median-Differenz", 
       color = "90%-CI Bandbreite", 
       title = "Manhattan Plot: Median-Differenz & CI-Band") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle=90, vjust=0.5))

changes <- left_join(summary_stats_2023, summary_stats_2024, by="Variable", suffix = c("_2023", "_2024")) %>%
  mutate(delta_median = HL_est_2024 - HL_est_2023)

ggplot(changes, aes(x=reorder(Variable, delta_median), y=HL_est_2024)) +
  geom_point(aes(color=delta_median), size=4) +
  scale_color_viridis_c(option = "D") +
  labs(x="Art", y="Median Differenz (2024)", color="Mediane Change\n(2024-2023)",
       title="Shift der Median-Differenz: Farmer–Expert pro Art") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle=90))





#Did the differences between the individual species differ between Farmer and wxp1 in the two years?
  
spalten <- names(df_filtered)[36:85]

# ID-Vorbereitung
df_filtered <- df_filtered %>%
  mutate(Year = as.integer(substr(Uniq_rprog, 1, 2)) + 2000,
         Poly_ID = substring(Uniq_rprog, 3))

# Nur Polygone, die in beiden Jahren vorkommen:
ids_both_years <- df_filtered %>%
  group_by(Poly_ID) %>%
  summarise(n_years = n_distinct(Year)) %>%
  filter(n_years == 2) %>%
  pull(Poly_ID)

df_both <- df_filtered %>%
  filter(Poly_ID %in% ids_both_years)

result_table <- map_dfr(
  spalten,
  function(var) {
    dat <- df_both %>%
      select(Poly_ID, Year, Exp_fact, !!var) %>%
      filter(Exp_fact %in% c("farmer", "exp1")) %>%
      pivot_wider(names_from = Exp_fact, values_from = !!var)
    dat <- dat %>%
      mutate(Diff = farmer - exp1)
    diff_wide <- dat %>%
      select(Poly_ID, Year, Diff) %>%
      pivot_wider(names_from = Year, values_from = Diff, names_prefix = "Diff_")
    v_23 <- diff_wide$Diff_2023
    v_24 <- diff_wide$Diff_2024
    valid <- !is.na(v_23) & !is.na(v_24)
    n_pairs <- sum(valid)
    # Die Veränderung berechnen
    change_vec <- v_24[valid] - v_23[valid]
    if(n_pairs < 3 | all(v_23[valid] == v_24[valid])) {
      return(tibble(
        Variable = var, n_pairs = n_pairs,
        V = NA, p_value = NA, HL_est = NA, HL_CI_low = NA, HL_CI_high = NA,
        median_change = NA
      ))
    }
    wt <- wilcox.test(v_24[valid], v_23[valid], paired = TRUE, conf.int = TRUE)
    tibble(
      Variable = var,
      n_pairs = n_pairs,
      V = wt$statistic,
      p_value = wt$p.value,
      HL_est = as.numeric(wt$estimate),
      HL_CI_low = wt$conf.int[1],
      HL_CI_high = wt$conf.int[2],
      median_change = median(change_vec, na.rm = TRUE)
    )
  }
)

print(result_table, n = 50)

write.csv(result_table, "wilcoxon_CI_IDspecies_farmers_experts_difference_change.csv", row.names = FALSE)

#HL_est and median_change are positive if the Farmer-Expert difference 
#in 2024 was TYPICALLY greater than in 2023
#(Farmer is further behind Expert than in the previous year).


library(ggplot2)
library(viridis)

result_table <- result_table %>%
  mutate(
    CI_band = case_when(
      HL_CI_low >= -1 & HL_CI_high <= 1 ~ "+/-1",
      HL_CI_low >= -2 & HL_CI_high <= 2 ~ "+/-2",
      HL_CI_low >= -3 & HL_CI_high <= 3 ~ "+/-3",
      TRUE ~ "mehr"
    ),
    CI_band = factor(CI_band, levels = c("+/-1", "+/-2", "+/-3", "mehr"))
  )

ggplot(result_table %>% arrange(median_change), 
       aes(x = reorder(Variable, median_change), y = median_change, color = CI_band)) +
  geom_pointrange(aes(ymin = HL_CI_low, ymax = HL_CI_high), size = 1) +
  scale_color_viridis_d(option = "D") +
  geom_hline(yintercept=0, linetype="dashed") +
  labs(x = "Art", y = "Median der Änderung (2024-2023)", color = "90%-CI Bandbreite",
       title = "Änderung der Bewertung (Farmer-Expert) von 2023 zu 2024") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle=90, vjust=0.5))



#d19_Carlinavulgaris, d19_Listeraovata, d19_Lythrumsalicaria, d19_Sedumacre


#manhattan plot mit differenzen und Jahrestrends
#install.packages("leftjoin")
#library(leftjoin)
# Benenne p_value in result_table um:
result_table <- result_table %>%
  rename(p_value_change = p_value)

combined_results <- left_join(
  wilcoxon_results, 
  result_table %>% select(Variable, median_change, p_value_change),
  by = "Variable"
) %>%
  mutate(sig_change = p_value_change < 0.05)

library(ggplot2)
library(viridis)

library(dplyr)
library(tidyr)

# Für alle Arten und Jahre Median Farmer–Expert-Differenz bestimmen
medians_per_year <- df_both %>%
  select(Poly_ID, Year, Exp_fact, all_of(spalten)) %>%
  filter(Exp_fact %in% c("farmer", "exp1")) %>%
  pivot_longer(cols = spalten, names_to = "Variable", values_to = "Value") %>%
  pivot_wider(names_from = Exp_fact, values_from = Value) %>%
  mutate(Diff = farmer - exp1) %>%
  group_by(Variable, Year) %>%
  summarise(median_diff = median(Diff, na.rm = TRUE), .groups = "drop")
#als tabelle pivot 
medians_wide <- medians_per_year %>%
  pivot_wider(names_from = Year, values_from = median_diff, names_prefix = "median_diff_")
combined <- left_join(
  result_table %>% select(Variable, median_change, p_value), 
  medians_wide,
  by = "Variable"
)

names(result_table)

library(ggplot2)
library(viridis)

ggplot(combined %>% filter(p_value < 0.05), 
       aes(x = reorder(Variable, median_diff_2023), y = median_diff_2023)) +
  geom_segment(aes(
    xend = reorder(Variable, median_diff_2023),
    yend = median_diff_2024,
    color = Richtung
  ),
  size = 1.5,
  arrow = arrow(type = "closed", length = unit(0.13, "inches"))) +
  geom_point(aes(y = median_diff_2024), shape = 21, fill = "white", size = 2) +
  geom_point(shape = 21, fill = "black", size = 2) +
  scale_color_manual(values = c("Annäherung" = "forestgreen", "Entfernung" = "orangered", "Gleich" = "gray")) +
  labs(x = "Art", y = "Median Farmer–Expert Differenz", color = "Trend Richtung",
       title = "Signifikante Änderung & Richtung pro Art") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle=90, vjust=0.5))



line_scale <- 0.45
ggplot(combined_results, aes(x = reorder(Variable, HL_est), y = HL_est)) +
  # Pfeile wie bisher ... (s. vorherigen Plot)
  # Punktplot
  geom_point(aes(color = CI_band, size = -log10(p_value)), alpha = 0.9) +
  # Füge Sterne bei signifikanter Änderung und median_change ~ 0 hinzu
  geom_text(
    data = subset(combined_results, sig_change & abs(median_change) < 1e-5),
    aes(label = "*"), 
    nudge_y = 0.4, color = "black", size = 5
  ) +
  scale_color_viridis_d(option = "D") +
  scale_size_continuous(range = c(2, 6)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle=90, vjust=0.5))



#next try manhattan with annual changes and directions
# Determine median farmer–expert difference for all types and years
medians_per_year <- df_both %>%
  select(Poly_ID, Year, Exp_fact, all_of(spalten)) %>%
  filter(Exp_fact %in% c("farmer", "exp1")) %>%
  pivot_longer(cols = spalten, names_to = "Variable", values_to = "Value") %>%
  pivot_wider(names_from = Exp_fact, values_from = Value) %>%
  mutate(Diff = farmer - exp1) %>%
  group_by(Variable, Year) %>%
  summarise(median_diff = median(Diff, na.rm = TRUE), .groups = "drop")
#als tabelle pivot 
medians_wide <- medians_per_year %>%
  pivot_wider(names_from = Year, values_from = median_diff, names_prefix = "median_diff_")

combined <- left_join(
  result_table %>% select(Variable, median_change, p_value_change), 
  medians_wide,
  by = "Variable"
)
combined <- combined %>%
  mutate(
    Richtung = case_when(
      abs(median_diff_2024) < abs(median_diff_2023) ~ "Annäherung",
      abs(median_diff_2024) > abs(median_diff_2023) ~ "Entfernung",
      TRUE ~ "Gleich"
    )
  )

combined <- left_join(
  result_table %>% select(Variable, median_change, p_value_change, median_diff_2023, median_diff_2024, CI_band), 
  medians_wide,
  by = "Variable"
)

library(ggplot2)
library(viridis)

# line_scale: wie weit der Pfeil gehen soll (optional anpassen)
line_scale <- 0.45

# Daten-Subset NUR für signifikante Änderungen
sig_arrow <- combined %>% filter(p_value_change < 0.05)

# Deine Hauptpunkte (manhattan plot) für alle Arten
ggplot(combined, aes(x = reorder(Variable, HL_est), y = HL_est)) +
  geom_point(aes(color = CI_band, size = -log10(p_value)), alpha = 0.9) +
  # Nur für signifikante Arten einen Richtungspfeil:
  geom_segment(
    data = sig_arrow,
    aes(
      xend = reorder(Variable, HL_est),
      yend = median_diff_2024,  # oder, wenn du den shift sichtbarer haben willst, yend = HL_est + (median_diff_2024 - median_diff_2023)
      color = ifelse(abs(median_diff_2024) < abs(median_diff_2023), "Annäherung", "Entfernung")
    ),
    size = 1.5,
    arrow = arrow(type = "closed", length = unit(0.13, "inches")),
    show.legend = TRUE
  ) +
  scale_color_manual(
    values = c("Annäherung" = "forestgreen", "Entfernung" = "orangered", "+/-1" = "#440154FF", "+/-2" = "#21908CFF", "+/-3" = "#FDE725FF", "mehr" = "gray"),
    breaks = c("+/-1", "+/-2", "+/-3", "mehr", "Annäherung", "Entfernung"),
    name = "90%-CI / Trend"
  ) +
  scale_size_continuous(range = c(2, 6)) +
  labs(
    x = "Art", y = "Median-Differenz (Farmer-Expert)",
    size = "-log10(p)",
    title = "Manhattan Plot: Farmer-Expert-Differenz\nmit Richtungspfeil für signifikante Jahresänderungen",
    subtitle = "Punkte: Wilcoxon/Mittelwert-CI, Pfeilrichtung/-farbe: Annäherung oder Entfernung"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle=90, vjust=0.5))





