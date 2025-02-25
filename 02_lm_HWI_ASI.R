# SOCIAL BOND STRENGTH AND ACOUSTIC SIMILARITY 
# 
# Author:             Maria Silvia Labriola
# Created:            12/06/2024
# 
# R version:          R version 4.3.3
#
#
# PURPOSE
# 
# Assessing correlation between a) social bond strength measured using photo-ID data
# and SW data of matched individuals b) social bond strength and acoustic similarity


# Load required libraries 

library(Rraven)
library(dplyr)      
library(scales)     
library(tidyr)
library(ggplot2)
library(reshape2)
library(lme4)
library(MuMIn)
library(car)
library(Rraven)
library(mgcv)
library(voxel)


# Load and Clean Data

matchings <-  readRDS("matchings.RData")
SWedges <- read.csv("Data/SWedges.csv") # SOCPROG/Gephi outputs
PHOTOedges <- read.csv("Data/PHOTOedges.csv")

SWedges <- SWedges[,c(1,2,6)]
PHOTOedges <- PHOTOedges[,c(1,2,6)]


# Assign HWI to matched dyads ----

# Generate all possible pairs of SW_IDs
dyads_sw <- matchings %>% 
  dplyr::select(SW_ID) %>% 
  expand_grid(SW_A = SW_ID, SW_B = SW_ID) %>% 
  filter(SW_A != SW_B)  

# Find the matching PHOTO_IDs 
dyads_m  <- dyads_sw %>% 
  left_join(matchings, by = c("SW_A" = "SW_ID")) %>% 
  rename(PHOTO_A = PHOTO_ID) %>%
  left_join(matchings, by = c("SW_B" = "SW_ID")) %>%
  rename(PHOTO_B = PHOTO_ID)

# Sort A and B in each row to get unique dyad 
dyads_m$sw_dyad <- apply(dyads_m, 1, function(row) {
  paste(sort(c(row["SW_A"], row["SW_B"])), collapse = "_")
})

dyads_m <- dyads_m[!duplicated(dyads_m$sw_dyad), ]

dyads_m$ph_dyad <- apply(dyads_m, 1, function(row) {
  paste(sort(c(row["PHOTO_A"], row["PHOTO_B"])), collapse = "_")
})

dyads_m <- dyads_m[,c("sw_dyad", "ph_dyad")]


# Extract HWI for matched dyads

SWedges$sw_dyad <- apply(SWedges, 1, function(row) {
  paste(sort(c(row["Source"], row["Target"])), collapse = "_")
})


edges_m <- dyads_m %>%
  left_join(SWedges, by = "sw_dyad")


edges_m <- edges_m[,c("sw_dyad", "ph_dyad", "Weight")]
colnames(edges_m)[ncol(edges_m)] <- "sw_HWI"


PHOTOedges$ph_dyad <- apply(PHOTOedges, 1, function(row) {
  paste(sort(c(row["Source"], row["Target"])), collapse = "_")
})



edges_m <- edges_m %>%
  left_join(PHOTOedges, by = "ph_dyad")


colnames(edges_m)[ncol(edges_m)] <- "ph_HWI"
edges_m <- edges_m[,c("sw_dyad", "ph_dyad", "sw_HWI", "ph_HWI")]



# Acoustic Similarity Index ----

# Load RavenPro data (selecion tables with acoustic parameters)

rvn.dat <- imp_raven(path = "Data/RavenPro/", all.data = TRUE, name.from.file = TRUE, ext.case = "lower")
rvn.dat <- rvn.dat %>%
  filter(Category == "SW", Quality >= 1.5)

# select acoustic parameters to use
rvn.dat <- rvn.dat %>%
  rename(Delta_Time = `Delta Time (s)`, N_max = `N max`, N_min = `N min`, Inflection_Points = `Inflection Points`)


## for matchings ----

SW.rvn_m <- rvn.dat  %>%
  dplyr::select(ID_Category, Type, Contour, N_max, N_min, Inflection_Points, Delta_Time) %>%
  filter(ID_Category %in% matchings$SW_ID) %>%
  na.omit() 

# remove rows with empty strings
rows_with_empty <- apply(SW.rvn_m, 1, function(row) any(row == ""))
SW.rvn_m <- SW.rvn_m[!rows_with_empty, ]


colnames(SW.rvn_m) <- c("SW_ID", "whistle_type" ,"contour_type", "local_maxima", 
                        "local_minima", "inflection_points", "duration")


# For each SW:
# 1) Calculate median values for numerical variables for each SW and round
# to nearest integer
# 2) Assign types based on most frequent (generally the whistle and contour type
# should be the same for all replicates of the same SW

get_mode <- function(v) {
  uniq_v <- unique(v)
  uniq_v[which.max(tabulate(match(v, uniq_v)))]
}

SW.rvn_m <- SW.rvn_m %>%
  group_by(SW_ID) %>%
  summarise(
    whistle_type = get_mode(whistle_type),
    contour_type = get_mode(contour_type),
    local_maxima = median(local_maxima),
    local_minima = median(local_minima),
    inflection_points = median(inflection_points),
    duration = median(duration)  
  )


# Normalize numerical parameters using min-max scaling
SW.rvn_m <- SW.rvn_m %>%
  mutate(
    local_maxima_scaled = rescale(local_maxima),
    local_minima_scaled = rescale(local_minima),
    inflection_points_scaled = rescale(inflection_points),
    duration_scaled = rescale(duration)
  )

# Create a matrix of the scaled numerical parameters
num_features <- SW.rvn_m %>%
  dplyr::select(local_maxima_scaled, local_minima_scaled, inflection_points_scaled, duration_scaled)

# Compute the Euclidean distance matrix for numerical parameters
num_dist_matrix <- as.matrix(dist(num_features, method = "euclidean"))

# Function to compute binary similarity (1 if same, 0 if different)
binary_similarity <- function(v1, v2) {
  outer(v1, v2, Vectorize(function(x, y) ifelse(x == y, 1, 0)))
}

# Calculate the similarity matrix for whistle_type and contour_type
whistle_type_sim <- binary_similarity(SW.rvn_m$whistle_type, SW.rvn_m$whistle_type)
contour_type_sim <- binary_similarity(SW.rvn_m$contour_type, SW.rvn_m$contour_type)

# Combine the two categorical similarity matrices 
categorical_sim <- 0.5 * whistle_type_sim + 0.5 * contour_type_sim

# Convert numerical distance into similarity
num_sim_matrix <- 1 / (1 + num_dist_matrix)

# Combine numerical similarity and categorical similarity 
final_asi_matrix <- 0.5 * num_sim_matrix + 0.5 * categorical_sim
colnames(final_asi_matrix) = rownames(final_asi_matrix) = SW.rvn_m$SW_ID

# Visualize ASI matrix as a heatmap

asi_melted <- melt(final_asi_matrix)

ggplot(asi_melted, aes(Var1, Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "seashell", high = "darkblue") +
  labs(title = "Acoustic Similarity Index (ASI)", x = "Whistle ID", y = "Whistle ID") +
  theme_minimal()


# Add ASI to dataset
asi_melted$sw_dyad <- apply(asi_melted, 1, function(row) {
  paste(sort(c(row["Var1"], row["Var2"])), collapse = "_")
})

edges_m <- edges_m %>%
  left_join(asi_melted, by = "sw_dyad") 

edges_m <- edges_m[!duplicated(edges_m$sw_dyad), ]

colnames(edges_m)[ncol(edges_m)] <- "ASI"
edges_m <- edges_m[,-c((ncol(edges_m)-1), (ncol(edges_m)-2))]


saveRDS(edges_m, "edges_m.RData")


## for all SWs----

# Build dataset with all info on associations between SWs

# Generate all possible pairs of SW_IDs
uniqueSWs <- unique(c(SWedges$Source, SWedges$Target))
SWedges_all <- expand_grid(SW_A = uniqueSWs, SW_B = uniqueSWs) %>% 
  filter(SW_A != SW_B)  

# Sort A and B in each row to get unique dyad 
SWedges_all$sw_dyad <- apply(SWedges_all, 1, function(row) {
  paste(sort(c(row["SW_A"], row["SW_B"])), collapse = "_")
})

SWedges_all <- SWedges_all[!duplicated(SWedges_all$sw_dyad), ]


SWedges_all <- SWedges_all %>%
  left_join(SWedges, by = "sw_dyad")


SWedges_all <- SWedges_all[,c("SW_A", "SW_B", "sw_dyad", "Weight")]
colnames(SWedges_all)[ncol(SWedges_all)] <- "sw_HWI"

SWedges_all$sw_HWI[is.na(SWedges_all$sw_HWI)] <- 0


# Get ASI for all analysed SWs

SW.rvn <- rvn.dat  %>%
  dplyr::select(ID_Category, Type, Contour, N_max, N_min, Inflection_Points, Delta_Time) %>%
  filter(ID_Category %in% unique(c(SWedges_all$SW_A, SWedges_all$SW_B))) %>%
  na.omit() 

rows_with_empty <- apply(SW.rvn, 1, function(row) any(row == ""))
SW.rvn <- SW.rvn[!rows_with_empty, ]


colnames(SW.rvn) <- c("SW_ID", "whistle_type" ,"contour_type", "local_maxima", 
                      "local_minima", "inflection_points", "duration")

SW.rvn <- SW.rvn %>%
  group_by(SW_ID) %>%
  summarise(
    whistle_type = get_mode(whistle_type),
    contour_type = get_mode(contour_type),
    local_maxima = median(local_maxima),
    local_minima = median(local_minima),
    inflection_points = median(inflection_points),
    duration = median(duration)  
  )

SW.rvn <- SW.rvn %>%
  mutate(
    local_maxima_scaled = rescale(local_maxima),
    local_minima_scaled = rescale(local_minima),
    inflection_points_scaled = rescale(inflection_points),
    duration_scaled = rescale(duration)
  )

num_features <- SW.rvn %>%
  dplyr::select(local_maxima_scaled, local_minima_scaled, inflection_points_scaled, duration_scaled)

num_dist_matrix <- as.matrix(dist(num_features, method = "euclidean"))

whistle_type_sim <- binary_similarity(SW.rvn$whistle_type, SW.rvn$whistle_type)
contour_type_sim <- binary_similarity(SW.rvn$contour_type, SW.rvn$contour_type)

categorical_sim <- 0.5 * whistle_type_sim + 0.5 * contour_type_sim

num_sim_matrix <- 1 / (1 + num_dist_matrix)

final_asi_matrix <- 0.5 * num_sim_matrix + 0.5 * categorical_sim
colnames(final_asi_matrix) = rownames(final_asi_matrix) = SW.rvn$SW_ID

asi_melted <- melt(final_asi_matrix)

asi_melted$sw_dyad <- apply(asi_melted, 1, function(row) {
  paste(sort(c(row["Var1"], row["Var2"])), collapse = "_")
})

SWedges_all <- SWedges_all %>%
  left_join(asi_melted, by = "sw_dyad") 

SWedges_all <- SWedges_all[!duplicated(SWedges_all$sw_dyad), ]

colnames(SWedges_all)[ncol(SWedges_all)] <- "ASI"
SWedges_all <- SWedges_all[,-c(1,2,(ncol(SWedges_all)-1), (ncol(SWedges_all)-2))]
SWedges_all <- na.omit(SWedges_all) #remove dyads without ASI (no replicates with measured parameters)

colnames(SWedges_all)[which(colnames(SWedges_all)=="Weight")] <- "sw_HWI"


# Investigate the relationship between SW and PHOTO HWI ----

edges_withzeros <- edges_m
edges_withzeros$sw_HWI[is.na(edges_withzeros$sw_HWI)] <- 0

# Fit linear models
lm1 <- lm(ph_HWI ~ sw_HWI, data = edges_withzeros)
summary(lm1)
coef1 = summary(lm1)$coefficients[2]
pv1 = summary(lm1)$coefficients[8]
# plot(lm1)

ggplot(edges_withzeros, aes(x = sw_HWI, y = ph_HWI)) +
  geom_jitter(alpha = 0.7, size = 3, colour="darkorange") +
  geom_smooth(method = "lm", colour="black", linewidth=1.5, linetype="dashed") + 
  theme_test() +
  labs(x = "SW HWI (Acoustic)", y = "PHOTO HWI (Visual)",
       title = "Relationship Between Acoustic and Visual Social Bond Strengths, with zeros") +
  annotate("text", 
           x = min(edges_withzeros$sw_HWI) + 0.055, 
           y = max(edges_withzeros$ph_HWI) * 1.1,
           label = paste0("β = ", round(coef1, 3), " (p = ", round(pv1, 3), ")"))



# Investigate the relationship between social bond strength and ASI ----
# all SWs

lmASI <-  lm(sw_HWI ~ ASI, data = SWedges_all)
summary(lmASI)
coef3 = summary(lmASI)$coefficients[2]
pv3 = summary(lmASI)$coefficients[8]


ggplot(SWedges_all, aes(x = sw_HWI, y = ASI)) +
  geom_jitter(alpha = 0.5, size = 3, colour="darkgreen") +
  geom_smooth(method = "lm", colour="black", linewidth=1.5, linetype="dashed") + 
  theme_test() +
  labs(x = "SW HWI (Acoustic)", y =  "ASI",
       title = "Relationship Between Acoustic Social Bond Strengths and ASI") +
  annotate("text", 
           x = max(SWedges_all$sw_HWI) * 0.8, 
           y = min(SWedges_all$ASI) * 1.1,
           label = paste0("β = ", round(coef3, 3), " (p = ", round(pv3, 3), ")"))
