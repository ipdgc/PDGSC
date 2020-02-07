library(tibble)
library(tidyr)
library(readr)
library(purrr)
library(dplyr)
library(stringr)
library(ggplot2)
library(RColorBrewer)
library(CGEN)

setwd("/home/dkia/clustering_results_keepmono_50depth_50anc_2to1")
# setwd("/home/dkia/clustering_results")

# dataset_to_load <- Sys.getenv("DATASET_TO_LOAD")
#  load(dataset_to_load)



# Functions -------------------------------------------------------------------------------------------

##### Check eigenvec is proportional to variance explained #####

define_plot_theme <- function(){
  
  theme_DZ <- 
    theme_minimal() +
    theme(text = element_text(color = "gray20", family = "sans"),
          plot.margin = margin(1, 1, 1, 1, "cm"),
          plot.title = element_text(size = rel(1.1), hjust = 0.5, vjust = 3, face = "bold" ), 
          axis.text.x = element_text(face = "italic", size = rel(0.75), angle = 30, vjust = 0.5),
          axis.text.y = element_text(face = "italic", size = rel(1), hjust = 0.5),
          axis.title.x = element_text(size = rel(0.9)), 
          axis.title.y = element_text(size = rel(0.9)),
          axis.line = element_line(color = "gray40", size = 0.5), 
          legend.text = element_text(size = rel(0.6)),
          legend.title = element_text(size = rel(0.75), face = "bold"))
  
  return(theme_DZ)
  
}

format_PCA_eigenvec <- function(ped){
  
  PCA_eigenvec_rightcols <- ped[,c("iid","CASE", "mean_depth", "n_sites", str_c("PC", 1:20))]
  PCA_eigenvec_rightcols$CASE <- as.numeric(PCA_eigenvec_rightcols$CASE)
  PCA_eigenvec_rightcols$mean_depth <- as.numeric(PCA_eigenvec_rightcols$mean_depth)
  PCA_eigenvec_rightcols$n_sites <- as.numeric(PCA_eigenvec_rightcols$n_sites)
  PCA_eigenvec_rightcols$PC1 <- as.numeric(PCA_eigenvec_rightcols$PC1)
  PCA_eigenvec_rightcols$PC2 <- as.numeric(PCA_eigenvec_rightcols$PC2)
  PCA_eigenvec_rightcols$PC3 <- as.numeric(PCA_eigenvec_rightcols$PC3)
  PCA_eigenvec_rightcols$PC4 <- as.numeric(PCA_eigenvec_rightcols$PC4)
  PCA_eigenvec_rightcols$PC5 <- as.numeric(PCA_eigenvec_rightcols$PC5)
  PCA_eigenvec_rightcols$PC6 <- as.numeric(PCA_eigenvec_rightcols$PC6)
  PCA_eigenvec_rightcols$PC7 <- as.numeric(PCA_eigenvec_rightcols$PC7)
  PCA_eigenvec_rightcols$PC8 <- as.numeric(PCA_eigenvec_rightcols$PC8)
  PCA_eigenvec_rightcols$PC9 <- as.numeric(PCA_eigenvec_rightcols$PC9)
  PCA_eigenvec_rightcols$PC10 <- as.numeric(PCA_eigenvec_rightcols$PC10)
  PCA_eigenvec_rightcols$PC11 <- as.numeric(PCA_eigenvec_rightcols$PC11)
  PCA_eigenvec_rightcols$PC12 <- as.numeric(PCA_eigenvec_rightcols$PC12)
  PCA_eigenvec_rightcols$PC13 <- as.numeric(PCA_eigenvec_rightcols$PC13)
  PCA_eigenvec_rightcols$PC14 <- as.numeric(PCA_eigenvec_rightcols$PC14)
  PCA_eigenvec_rightcols$PC15 <- as.numeric(PCA_eigenvec_rightcols$PC15)
  PCA_eigenvec_rightcols$PC16 <- as.numeric(PCA_eigenvec_rightcols$PC16)
  PCA_eigenvec_rightcols$PC17 <- as.numeric(PCA_eigenvec_rightcols$PC17)
  PCA_eigenvec_rightcols$PC18 <- as.numeric(PCA_eigenvec_rightcols$PC18)
  PCA_eigenvec_rightcols$PC19 <- as.numeric(PCA_eigenvec_rightcols$PC19)
  PCA_eigenvec_rightcols$PC20 <- as.numeric(PCA_eigenvec_rightcols$PC20)
  
  PCA_eigenvec_formatted <- 
    PCA_eigenvec_rightcols %>% 
    arrange(desc(CASE))
  
  
  return(PCA_eigenvec_formatted)
  
}

compare_mean_range_of_eigenvec_to_eigenval <- function(PCA_eigenvec_formatted, PCA_eigenval){
  
  calc_range <- function(x){
    
    stopifnot(is.integer(x) || is.numeric(x))
    
    range_diff <- range(x)[2] - range(x)[1]
    
    return(range_diff)
    
  }
  
  
  PCA_eigenval_with_mean_range <- 
    PCA_eigenval %>% 
    mutate(Component = 1:20,
           C_ranges = map_dbl(.x = PCA_eigenvec_formatted[, str_c("PC", 1:20)], .f = calc_range),
           C_mean = map_dbl(.x = PCA_eigenvec_formatted[, str_c("PC", 1:20)], .f = mean)) %>% 
    rename(eigenval = X1)
  
  return(PCA_eigenval_with_mean_range)
  
}

plot_lm_eigenvec_mean_range_vs_eigenval <- function(PCA_eigenval_with_mean_range){
  
  PCA_eigenval_vs_eigenvec_range <- 
    ggplot(PCA_eigenval_with_mean_range, aes(x = eigenval, y = C_ranges, colour = Component)) +
    geom_point() + 
    ggtitle("PDGSC Group 2 PCA eigenval vs eigenvec range")
  theme_DZ
  
  ggsave("group2_PCA_eigenval_vs_eigenvec_range.pdf", PCA_eigenval_vs_eigenvec_range, 
         width = 8, height = 6)
  
  PCA_eigenval_vs_eigenvec_mean <- 
    ggplot(PCA_eigenval_with_mean_range, aes(x = eigenval, y = C_mean, colour = Component)) +
    geom_point() + 
    geom_smooth(method = "lm") +
    scale_x_continuous(name = "Eigenval") +
    scale_y_continuous(name = "Eigenvec Mean") +
    ggtitle("PDGSC Group 2 PCA eigenval vs eigenvec mean") + 
    theme_DZ
  
  ggsave("group2_PCA_eigenval_vs_eigenvec_mean.pdf", PCA_eigenval_vs_eigenvec_mean, 
         width = 8, height = 6)  
  
  lm_eigenval_vs_eigenvec_mean <- lm(eigenval ~ C_mean, data = PCA_eigenval_with_mean_range)
  capture.output(summary(lm_eigenval_vs_eigenvec_mean), file = "group2_lm_eigenval_vs_eigenvec_mean.txt")
  
}

plot_distribution_of_eigenval <- function(PCA_eigenval_with_mean_range){
  
  eigenval_vs_components_plot <- 
    ggplot(PCA_eigenval_with_mean_range, aes(x = Component, y = eigenval)) +
    geom_point() +
    scale_y_continuous(name = "Eigenval") +
    ggtitle("Distribution of Eigenvalues over Components 1-20") +
    theme_DZ
  
  ggsave("Distribution_of_Eigenvalues_over_Components_1-20.pdf", eigenval_vs_components_plot, 
         width = 8, height = 6)
  
  
}


normalise_PCA_eigenvec <- function(PCA_eigenvec_formatted){
  
  PCA_eigenvec_formatted_normalised <- PCA_eigenvec_formatted[,c("iid", "CASE", "mean_depth", "n_sites", "PC1", "PC2")]
  PCA_eigenvec_formatted_normalised$mean_depth_range <- max(PCA_eigenvec_formatted_normalised$mean_depth)-min(PCA_eigenvec_formatted_normalised$mean_depth)
  PCA_eigenvec_formatted_normalised$n_sites_range <- max(PCA_eigenvec_formatted_normalised$n_sites)-min(PCA_eigenvec_formatted_normalised$n_sites)
  PCA_eigenvec_formatted_normalised$PC1_range <- max(PCA_eigenvec_formatted_normalised$PC1)-min(PCA_eigenvec_formatted_normalised$PC1)
  PCA_eigenvec_formatted_normalised$PC2_range <- max(PCA_eigenvec_formatted_normalised$PC2)-min(PCA_eigenvec_formatted_normalised$PC2)
  PCA_eigenvec_formatted_normalised$mean_depth_numerator <- PCA_eigenvec_formatted_normalised$mean_depth-min(PCA_eigenvec_formatted_normalised$mean_depth)
  PCA_eigenvec_formatted_normalised$n_sites_numerator <- PCA_eigenvec_formatted_normalised$n_sites-min(PCA_eigenvec_formatted_normalised$n_sites)
  PCA_eigenvec_formatted_normalised$PC1_numerator <- PCA_eigenvec_formatted_normalised$PC1-min(PCA_eigenvec_formatted_normalised$PC1)
  PCA_eigenvec_formatted_normalised$PC2_numerator <- PCA_eigenvec_formatted_normalised$PC2-min(PCA_eigenvec_formatted_normalised$PC2)
  PCA_eigenvec_formatted_normalised$mean_depth_normalised_unweighted <- PCA_eigenvec_formatted_normalised$mean_depth_numerator/PCA_eigenvec_formatted_normalised$mean_depth_range
  PCA_eigenvec_formatted_normalised$n_sites_normalised_unweighted <- PCA_eigenvec_formatted_normalised$n_sites_numerator/PCA_eigenvec_formatted_normalised$n_sites_range
  PCA_eigenvec_formatted_normalised$PC1_normalised_unweighted <- PCA_eigenvec_formatted_normalised$PC1_numerator/PCA_eigenvec_formatted_normalised$PC1_range
  PCA_eigenvec_formatted_normalised$PC2_normalised_unweighted <- PCA_eigenvec_formatted_normalised$PC2_numerator/PCA_eigenvec_formatted_normalised$PC2_range

  PCA_eigenvec_formatted_normalised$mean_depth_normalised <- 50*PCA_eigenvec_formatted_normalised$mean_depth_normalised_unweighted
#  PCA_eigenvec_formatted_normalised$n_sites_normalised <- 10*PCA_eigenvec_formatted_normalised$n_sites_normalised_unweighted
  PCA_eigenvec_formatted_normalised$PC1_normalised <- 25*PCA_eigenvec_formatted_normalised$PC1_normalised_unweighted
  PCA_eigenvec_formatted_normalised$PC2_normalised <- 25*PCA_eigenvec_formatted_normalised$PC2_normalised_unweighted
  
  PCA_eigenvec_formatted_normalised$dummy <- 0
  
  PCA_eigenvec_formatted_normalised$mean_depth_range <- NULL
  PCA_eigenvec_formatted_normalised$n_sites_range <- NULL
  PCA_eigenvec_formatted_normalised$PC1_range <- NULL
  PCA_eigenvec_formatted_normalised$PC2_range <- NULL
  PCA_eigenvec_formatted_normalised$mean_depth_numerator <- NULL
  PCA_eigenvec_formatted_normalised$n_sites_numerator <- NULL
  PCA_eigenvec_formatted_normalised$PC1_numerator <- NULL
  PCA_eigenvec_formatted_normalised$PC2_numerator <- NULL  
  PCA_eigenvec_formatted_normalised$mean_depth_normalised_unweighted <- NULL
  PCA_eigenvec_formatted_normalised$n_sites_normalised_unweighted <- NULL
  PCA_eigenvec_formatted_normalised$PC1_normalised_unweighted <- NULL
  PCA_eigenvec_formatted_normalised$PC2_normalised_unweighted <- NULL
  
  return(PCA_eigenvec_formatted_normalised)
  
}

match_case_controls <- function(PCA_eigenvec_formatted, size, fixed, dist.vars){
  
  matched_sets <- 
    getMatchedSets(as.data.frame(PCA_eigenvec_formatted), CC=TRUE, NN=TRUE, ccs.var = "CASE",
                   dist.vars = dist.vars, 
                   size = size, fixed = fixed)
  
  return(matched_sets)
  
}

get_cluster_groups <- function(PCA, matched_sets){
  
  case_indexes <- 
    PCA %>% 
    mutate(index = 1:length(PCA$iid)) %>% 
    filter(CASE == 1) %>% 
    .[["index"]]
  
  cluster_groups <- list()
  
  for(i in seq_along(case_indexes)){
    
    regex_group <- str_c("^", case_indexes[i], "$")
    group <- which(str_detect(matched_sets[["CC"]], regex_group))   
    
    cluster_groups[[i]] <- group
    
  }
  
  names(cluster_groups) <- str_c("group_", seq_along(case_indexes))
  
  return(cluster_groups)
  
}

check_cluster_groups <- function(cluster_groups){
  
  cluster_sizes <- map_int(cluster_groups, length)
  all_group_sizes <- levels(factor(cluster_sizes))
  
  if(length(all_group_sizes) == 1){
    
    print(str_c("All groups include ", all_group_sizes, " samples."))
    
  }else{
    
    print(str_c("Group sizes include: ", str_c(all_group_sizes, collapse = ", "), ""))
    
  }
  
}

add_groups_to_PCA_data <- function(PCA, cluster_groups){
  
  PCA_grouped <- 
    PCA %>% 
    mutate(group = 0)
  
  for(i in seq_along(cluster_groups)){
    
    PCA_grouped$group[cluster_groups[[i]]] <- i
  }
  
  return(PCA_grouped)
}


update_matched_cohorts <- function(PCA_grouped){
  
  PCA_grouped$GROUP_matched <- NA
  
  for(i in seq_along(PCA_grouped$group)){
    matched_group <- PCA_grouped$group[i]
    matched_cohort <- PCA_grouped$STUDY_corrected[i]
    
    if(matched_group == 0){
      next
    }
    else if(PCA_grouped$CASE[i] == 1){
      PCA_grouped$GROUP_matched[i] <- PCA_grouped$GROUP[i]
    }else if(PCA_grouped$CASE[i] == 0){
      corresponding_cases <- subset(PCA_grouped, CASE == 1 & group == matched_group)     
      corresponding_studies <- corresponding_cases$STUDY_corrected
      corresponding_random_study <- corresponding_studies[sample(1:length(corresponding_studies),1)]
      if(length(corresponding_random_study) != 0){
        PCA_grouped$GROUP_matched[i] <- str_c(corresponding_random_study, "_control")
      }
    }
  }
  
  return(PCA_grouped)
}

plot_mean_depth_n_sites <- function(PCA_grouped, plot_filename, directory){
  
  kept_data <- subset(PCA_grouped, !is.na(GROUP_matched_ucl))
  PCA_grouped$GROUP_matched_ucl_unmatched <- ifelse(is.na(PCA_grouped$GROUP_matched_ucl), "unmatched_ADSP", PCA_grouped$GROUP_matched_ucl)
  
  PCA_grouped$GROUP_matched_with_n <- NA
  
  for(i in seq_along(PCA_grouped$GROUP_matched)){
    PCA_grouped$GROUP_matched_with_n[i] <- str_c(PCA_grouped$GROUP_matched_ucl_unmatched[i], "_n=", sum(PCA_grouped$GROUP_matched_ucl_unmatched == PCA_grouped$GROUP_matched_ucl_unmatched[i]))
  }
  
  ncase <- sum(kept_data$CASE==1)
  ncontrol <- sum(kept_data$CASE==0)
  
  ad_plot <- ggplot(PCA_grouped, aes(x = as.factor(GROUP_matched_with_n), y = as.numeric(mean_depth))) +
    geom_boxplot() + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggtitle(str_c("Mean depth (", ncase, " cases & ", ncontrol, " controls)"))
  ggsave(str_c("mean_depth_boxplot_", plot_filename, ".pdf"), plot = ad_plot, device="pdf", dpi = 300, height = 7, width = 7, units = "in", path = directory)
  #ggsave(str_c("plots/mean_depth_boxplot_NN_", NN_PARAMETER, "_SIZE_", SIZE_PARAMETER, "_", paste(DIST_VARIABLES_PARAMETER,collapse="_"), ".pdf"), plot = ad_plot, device="pdf", dpi = 300, height = 7, width = 7, units = "in")
  
  nsites_plot <- ggplot(PCA_grouped, aes(x = as.factor(GROUP_matched_with_n), y = as.numeric(n_sites))) +
    geom_boxplot() + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggtitle(str_c("N sites (", ncase, " cases & ", ncontrol, " controls)"))
  ggsave(str_c("nsites_plot_boxplot_", plot_filename, ".pdf"), plot = nsites_plot, device="pdf", dpi = 300, height = 7, width = 7, units = "in", path = directory)
  #ggsave(str_c("plots/nsites_plot_boxplot_NN_", NN_PARAMETER, "_SIZE_", SIZE_PARAMETER, "_", paste(DIST_VARIABLES_PARAMETER,collapse="_"), ".pdf"), plot = ad_plot, device="pdf", dpi = 300, height = 7, width = 7, units = "in")
  
}

plot_PCA_data_grouped <- function(PCA_grouped, file_name_prefix, directory){
  
  PCA_grouped <- PCA_grouped %>% separate(GROUP_matched, c("matched_cohort", "status"), sep="_", remove=FALSE)
  
  studies_to_plot <- unique(PCA_grouped$matched_cohort)
  studies_to_plot <- studies_to_plot[!is.na(studies_to_plot)]
  
  for(i in seq_along(studies_to_plot)){
    study <- studies_to_plot[i]
    data_to_plot <- subset(PCA_grouped, matched_cohort == study)
    
    ncase <- sum(data_to_plot$CASE==1)
    ncontrol <- sum(data_to_plot$CASE==0)
    
    PCA_data_10_groups_plot <- 
      ggplot(data_to_plot, 
             aes(x = as.numeric(PC1), y = as.numeric(PC2), colour = factor(status))) +
      geom_point(size = rel(1)) +
      scale_x_continuous(name = "Principal Component 1") +
      scale_y_continuous(name = "Principal Component 2") +
      scale_color_manual(name = "Status", values = c("black", brewer.pal(10, "Set3"))) +
      ggtitle(str_c("PCA for ", study, " (", ncase, " cases & ", ncontrol, " controls)")) +
      theme_bw() 
    
    ggsave(filename = str_c("pca_plot_", study, "_", file_name_prefix, ".pdf"), 
           plot = PCA_data_10_groups_plot, 
           path = directory,
           width = 8, 
           height = 6)
  }
}


# Main ------------------------------------------------------------------------------------------------
PCA_eigenval_input_file <- Sys.getenv("EIGENVAL_INPUT_FILE")
ped_input_file <- Sys.getenv("PED_INPUT_FILE")
SIZE_PARAMETER <- as.numeric(Sys.getenv("SIZE_PARAMETER"))
DIST_VARIABLES_PARAMETER_FILE <- Sys.getenv("DIST_VARIABLES_PARAMETERS_FILE")
DIST_VARIABLES_PARAMETER_FILE_df <- read_delim(DIST_VARIABLES_PARAMETER_FILE, delim="\t", col_names=FALSE, col_types=cols(.default = "c"))
DIST_VARIABLES_PARAMETER <- DIST_VARIABLES_PARAMETER_FILE_df$X1

PCA_eigenval <- read_table2(str_c("../", PCA_eigenval_input_file), col_names = F)
ped <- read_delim(str_c("../", ped_input_file), delim="\t", col_names = T, col_types=cols(.default="c"))

theme_DZ <- define_plot_theme()
PCA_eigenvec_formatted <- format_PCA_eigenvec(ped)
PCA_eigenval_with_mean_range <- compare_mean_range_of_eigenvec_to_eigenval(PCA_eigenvec_formatted, PCA_eigenval)
plot_lm_eigenvec_mean_range_vs_eigenval(PCA_eigenval_with_mean_range)
plot_distribution_of_eigenval(PCA_eigenval_with_mean_range)
PCA_eigenvec_formatted_normalised <- normalise_PCA_eigenvec(PCA_eigenvec_formatted)

print(summary(as.numeric(PCA_eigenvec_formatted_normalised$n_sites)))

print(str_c("Ready to cluster! Size parameter is ", SIZE_PARAMETER, ", matching variables are ", DIST_VARIABLES_PARAMETER, "!"))

system.time(matched_sets <- match_case_controls(PCA_eigenvec_formatted_normalised, size = SIZE_PARAMETER, fixed = F, dist.vars = c(DIST_VARIABLES_PARAMETER)))

save.image(file = str_c("matched_pairs_", "_SIZE_", SIZE_PARAMETER, "_", paste(DIST_VARIABLES_PARAMETER,collapse="_"), ".RData"))

print("Uploaded RData file!")

# plot_filename <- dataset_to_load %>% str_replace("matched_pairs_", "") %>%
#    str_replace(".RData", "")

plot_filename <- str_c("SIZE_", SIZE_PARAMETER, "_", paste(DIST_VARIABLES_PARAMETER,collapse="_"))

plot_directory <- str_c("matching_plots/", plot_filename, "/")

dir.create(plot_directory, recursive = TRUE)

cohort_info <- ped[,c("iid", "GROUP", "STUDY_corrected")]

PCA_eigenvec_formatted_normalised_cohortinfo <- left_join(PCA_eigenvec_formatted_normalised, cohort_info, by = "iid")

cluster_groups <- get_cluster_groups(PCA_eigenvec_formatted_normalised_cohortinfo, matched_sets)

check_cluster_groups(cluster_groups) 

PCA_eigenvec_formatted_normalised_grouped <- add_groups_to_PCA_data(PCA_eigenvec_formatted_normalised_cohortinfo, cluster_groups)

PCA_updated_cohorts <- update_matched_cohorts(PCA_eigenvec_formatted_normalised_grouped)

PCA_updated_cohorts$GROUP_matched_ucl <- ifelse(PCA_updated_cohorts$GROUP_matched == "ucl_case", "UCL2_case",
                                                ifelse(PCA_updated_cohorts$GROUP_matched == "ucl_control", "UCL2_control",
                                                       ifelse(PCA_updated_cohorts$GROUP_matched == "UCL_case", "UCL1_case",
                                                              ifelse(PCA_updated_cohorts$GROUP_matched == "UCL_control", "UCL1_control",
                                                                     PCA_updated_cohorts$GROUP_matched))))

plot_mean_depth_n_sites(PCA_updated_cohorts, plot_filename, directory = plot_directory)

plot_PCA_data_grouped(PCA_updated_cohorts, file_name_prefix = plot_filename, directory = plot_directory)

mathced_samples <- PCA_updated_cohorts[,c("iid", "GROUP_matched_ucl")]
ped_matched <- left_join(ped, mathced_samples, by ="iid")
ped_matched <- subset(ped_matched, !is.na(GROUP_matched_ucl))
write.table(ped_matched, "caseControlAnalysisData_matched_group2_2to1_depth_ancestry_210818.ped", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

print("Mission complete!")
