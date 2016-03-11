# This script is used to read in all data files


# From IGCSimulationSummary.Rmd

rm(list=ls())  # clean up workspace
IGC.path <- "/Users/xji3/GitFolders/IGCSimulation/"
paml.path <- "/Users/xji3/GitFolders/IGCSimulation/"
IGC.geo.list <- c(3.0, 10.0, 50.0, 100.0, 200.0, 300.0, 400.0, 500.0)
#IGC.geo.list <- c(1.0)
num.sim = 100
paralog = "YDR418W_YEL054C"


# Now read in PAML analysis of these simulated data set
data.path <- paste(paralog, "",sep = "/")
for (IGC.geo in IGC.geo.list){
  summary_mat <- NULL
  IGC.geo.path <- paste("IGCgeo_", toString(IGC.geo), ".0/", sep = "")
  file.name <- paste("geo", paste(toString(IGC.geo), ".0", sep = ""), "paml", "unrooted", "summary.txt", sep = "_")
  for (sim.num in 0:(num.sim - 1)){
    summary_file <- paste(paml.path, file.name, sep = "")
    if (file.exists(summary_file)){
      all <- readLines(summary_file, n = -1)
      col.names <- strsplit(all[1], ' ')[[1]][-1]
      row.names <- strsplit(all[length(all)], ' ')[[1]][-1]
      summary_mat <- as.matrix(read.table(summary_file, 
                                          row.names = row.names, 
                                          col.names = col.names))
      
    }
  }
  assign(paste("PAML", paste(toString(IGC.geo), ".0", sep = ""), "summary", sep = "_"), summary_mat)
}

#Read.summary <- function(IGC.path, IGC.x3.path, IGC.geo.list, num.sim = 100, paralog = "YDR418W_YEL054C"){
## Now read simulation result

sim.path <- paste("SimulationSummary", paralog, "",sep = "/")
for (IGC.geo in IGC.geo.list){
  summary_mat <- NULL
  IGC.geo.path <- paste("IGCgeo_", toString(IGC.geo), ".0/", sep = "")
  file.prefix <- paste(paralog, "MG94_geo", paste(toString(IGC.geo), ".0", sep = ""), "Sim", sep = "_")
  file.suffix <- "summary.txt"
  for (sim.num in 0:(num.sim - 1)){
    file.name <- paste(file.prefix, toString(sim.num), file.suffix, sep = "_")
    summary_file <- paste(IGC.path, sim.path, IGC.geo.path, file.name, sep = "")
    if (file.exists(summary_file)){
      all <- readLines(summary_file, n = -1)
      col.names <- paste("geo", paste(toString(IGC.geo), ".0", sep = ""), "sim", toString(sim.num), sep = "_")
      row.names <- strsplit(all[length(all)], ' ')[[1]][-1]
      summary_mat <- cbind(summary_mat, as.matrix(read.table(summary_file, 
                                                             row.names = row.names, 
                                                             col.names = col.names)))        
    }
  }
  assign(paste("geo", paste(toString(IGC.geo), ".0", sep = ""), "summary", sep = "_"), summary_mat)
}


# Now read in 'real'  % changes due to IGC in the simulation
data.path <- paste(paralog, "",sep = "/")
for (IGC.geo in IGC.geo.list){
  summary_mat <- NULL
  IGC.geo.path <- paste("IGCgeo_", toString(IGC.geo), ".0/", sep = "")
  file.prefix <- paste(paralog, "MG94_geo", paste(toString(IGC.geo), ".0", sep = ""), "Sim", sep = "_")
  file.suffix <- "short.log"
  for (sim.num in 0:(num.sim - 1)){
    file.name <- paste(file.prefix, toString(sim.num), file.suffix, sep = "_")
    file.name <- paste("sim", paste(toString(sim.num), "/", file.prefix, sep = ""),
                       toString(sim.num), file.suffix, sep = "_")
    summary_file <- paste(IGC.path, data.path, IGC.geo.path, file.name, sep = "")
    if (file.exists(summary_file)){
      all <- readLines(summary_file, n = -1)
      col.names <- paste("geo", paste(toString(IGC.geo), ".0", sep = ""), "sim", toString(sim.num), sep = "_")
      row.names <- strsplit(all[length(all)], ' ')[[1]][-1]
      summary_mat <- cbind(summary_mat, as.matrix(read.table(summary_file, 
                                                             row.names = row.names, 
                                                             col.names = col.names)))        
    }
  }
  assign(paste("pIGC", paste(toString(IGC.geo), ".0", sep = ""), "summary", sep = "_"), summary_mat)
}

###############################################################
###############################################################
###############################################################
###############################################################

# From CheckBlenAsymmetry.Rmd
for (IGC.geo in IGC.geo.list){
  for (localtree in 1:3){
    summary_mat <- NULL
    IGC.geo.path <- paste("IGCgeo_", toString(IGC.geo), ".0/", sep = "")
    file.name <- paste("geo", paste(toString(IGC.geo), ".0", sep = ""), "estimatedTau", "paml", "unrooted", "LocalTree", toString(localtree), "summary.txt", sep = "_")
    for (sim.num in 0:(num.sim - 1)){
      summary_file <- paste(paml.path, file.name, sep = "")
      if (file.exists(summary_file)){
        all <- readLines(summary_file, n = -1)
        col.names <- strsplit(all[1], ' ')[[1]][-1]
        row.names <- strsplit(all[length(all)], ' ')[[1]][-1]
        summary_mat <- as.matrix(read.table(summary_file, 
                                            row.names = row.names, 
                                            col.names = col.names))
        
      }
    }
    assign(paste("PAML", "estimatedTau", paste(toString(IGC.geo), ".0", sep = ""), "LocalTree", toString(localtree), "summary", sep = "_"), summary_mat)}
}


# Read in new PAML results
data.path <- paste(paralog, "",sep = "/")
for (IGC.geo in IGC.geo.list){
  summary_mat <- NULL
  IGC.geo.path <- paste("IGCgeo_", toString(IGC.geo), ".0/", sep = "")
  file.name <- paste("geo", paste(toString(IGC.geo), ".0", sep = ""), "estimatedTau", "paml", "unrooted", "1stTree", "summary.txt", sep = "_")
  for (sim.num in 0:(num.sim - 1)){
    summary_file <- paste(paml.path, file.name, sep = "")
    if (file.exists(summary_file)){
      all <- readLines(summary_file, n = -1)
      col.names <- strsplit(all[1], ' ')[[1]][-1]
      row.names <- strsplit(all[length(all)], ' ')[[1]][-1]
      summary_mat <- as.matrix(read.table(summary_file, 
                                          row.names = row.names, 
                                          col.names = col.names))
      
    }
  }
  assign(paste("PAML", "estimatedTau", paste(toString(IGC.geo), ".0", sep = ""), "1stTree", "summary", sep = "_"), summary_mat)
}

for (IGC.geo in IGC.geo.list){
  summary_mat <- NULL
  IGC.geo.path <- paste("IGCgeo_", toString(IGC.geo), ".0/", sep = "")
  file.name <- paste("geo", paste(toString(IGC.geo), ".0", sep = ""), "estimatedTau", "paml", "unrooted", "2ndTree", "summary.txt", sep = "_")
  for (sim.num in 0:(num.sim - 1)){
    summary_file <- paste(paml.path, file.name, sep = "")
    if (file.exists(summary_file)){
      all <- readLines(summary_file, n = -1)
      col.names <- strsplit(all[1], ' ')[[1]][-1]
      row.names <- strsplit(all[length(all)], ' ')[[1]][-1]
      summary_mat <- as.matrix(read.table(summary_file, 
                                          row.names = row.names, 
                                          col.names = col.names))
      
    }
  }
  assign(paste("PAML", "estimatedTau", paste(toString(IGC.geo), ".0", sep = ""), "2ndTree", "summary", sep = "_"), summary_mat)
}

# No asymmetry in the 0 branch length estimates
convert.summary <- function(localTree, summary, vector.names){
  new.summary <- vector("numeric", length(vector.names))
  names(new.summary) <- vector.names
  if (localTree == 1){
    new.summary[1:3] <- summary[1:3]
    new.summary["N0_N6"] <- summary["N0_N5"]
    new.summary["N0_kluyveriYDR418W"] <- summary["N0_kluyveriYDR418W"]
    new.summary["N0_N1"] <- 0.0
    new.summary["N1_N2"] <- summary["N0_N1"]
    new.summary["N1_castelliiYDR418W"] <- summary["N0_castelliiYDR418W"]
    new.summary["N2_bayanusYDR418W"] <- summary["N1_bayanusYDR418W"]
    new.summary["N2_N3"] <- summary["N1_N2"]
    new.summary["N3_N4"] <- summary["N2_N3"]
    new.summary["N3_kudriavzeviiYDR418W"] <- summary["N2_kudriavzeviiYDR418W"]
    new.summary["N4_N5"] <- summary["N3_N4"]
    new.summary["N4_mikataeYDR418W"] <- summary["N3_mikataeYDR418W"]
    new.summary["N5_paradoxusYDR418W"] <- summary["N4_paradoxusYDR418W"]
    new.summary["N5_cerevisiaeYDR418W"] <- summary["N4_cerevisiaeYDR418W"]
    new.summary["N6_N7"] <- summary["N5_N6"]
    new.summary["N6_castelliiYEL054C"] <- summary["N5_castelliiYEL054C"]
    new.summary["N7_N8"] <- summary["N6_N7"]
    new.summary["N7_bayanusYEL054C"] <- summary["N6_bayanusYEL054C"]
    new.summary["N8_N9"] <- summary["N7_N8"]
    new.summary["N8_kudriavzeviiYEL054C"] <- summary["N7_kudriavzeviiYEL054C"]
    new.summary["N9_mikataeYEL054C"] <- summary["N8_mikataeYEL054C"]
    new.summary["N9_N10"] <- summary["N8_N9"]
    new.summary["N10_paradoxusYEL054C"] <- summary["N9_paradoxusYEL054C"]
    new.summary["N10_cerevisiaeYEL054C"] <- summary["N9_cerevisiaeYEL054C"]
    
  }
  else if (localTree == 2){
    new.summary[1:3] <- summary[1:3]
    new.summary["N0_N6"] <- 0.0
    new.summary["N0_kluyveriYDR418W"] <- summary["N0_kluyveriYDR418W"]
    new.summary["N0_N1"] <- summary["N0_N1"]
    new.summary["N1_N2"] <- summary["N1_N2"]
    new.summary["N1_castelliiYDR418W"] <- summary["N1_castelliiYDR418W"]
    new.summary["N2_bayanusYDR418W"] <- summary["N2_bayanusYDR418W"]
    new.summary["N2_N3"] <- summary["N2_N3"]
    new.summary["N3_N4"] <- summary["N3_N4"]
    new.summary["N3_kudriavzeviiYDR418W"] <- summary["N3_kudriavzeviiYDR418W"]
    new.summary["N4_N5"] <- summary["N4_N5"]
    new.summary["N4_mikataeYDR418W"] <- summary["N4_mikataeYDR418W"]
    new.summary["N5_paradoxusYDR418W"] <- summary["N5_paradoxusYDR418W"]
    new.summary["N5_cerevisiaeYDR418W"] <- summary["N5_cerevisiaeYDR418W"]
    new.summary["N6_N7"] <- summary["N0_N6"]
    new.summary["N6_castelliiYEL054C"] <- summary["N0_castelliiYEL054C"]
    new.summary["N7_N8"] <- summary["N6_N7"]
    new.summary["N7_bayanusYEL054C"] <- summary["N6_bayanusYEL054C"]
    new.summary["N8_N9"] <- summary["N7_N8"]
    new.summary["N8_kudriavzeviiYEL054C"] <- summary["N7_kudriavzeviiYEL054C"]
    new.summary["N9_mikataeYEL054C"] <- summary["N8_mikataeYEL054C"]
    new.summary["N9_N10"] <- summary["N8_N9"]
    new.summary["N10_paradoxusYEL054C"] <- summary["N9_paradoxusYEL054C"]
    new.summary["N10_cerevisiaeYEL054C"] <- summary["N9_cerevisiaeYEL054C"]
  }
  else if (localTree == 3){
    new.summary[1:3] <- summary[1:3]
    new.summary["N0_N6"] <- 0.0
    new.summary["N0_kluyveriYDR418W"] <- summary["N0_kluyveriYDR418W"]
    new.summary["N0_N1"] <- 0.0
    new.summary["N1_N2"] <- summary["N0_N1"]
    new.summary["N1_castelliiYDR418W"] <- summary["N0_castelliiYDR418W"]
    new.summary["N2_bayanusYDR418W"] <- summary["N1_bayanusYDR418W"]
    new.summary["N2_N3"] <- summary["N1_N2"]
    new.summary["N3_N4"] <- summary["N2_N3"]
    new.summary["N3_kudriavzeviiYDR418W"] <- summary["N2_kudriavzeviiYDR418W"]
    new.summary["N4_N5"] <- summary["N3_N4"]
    new.summary["N4_mikataeYDR418W"] <- summary["N3_mikataeYDR418W"]
    new.summary["N5_paradoxusYDR418W"] <- summary["N4_paradoxusYDR418W"]
    new.summary["N5_cerevisiaeYDR418W"] <- summary["N4_cerevisiaeYDR418W"]
    new.summary["N6_N7"] <- summary["N0_N5"]
    new.summary["N6_castelliiYEL054C"] <- summary["N0_castelliiYEL054C"]
    new.summary["N7_N8"] <- summary["N5_N6"]
    new.summary["N7_bayanusYEL054C"] <- summary["N5_bayanusYEL054C"]
    new.summary["N8_N9"] <- summary["N6_N7"]
    new.summary["N8_kudriavzeviiYEL054C"] <- summary["N6_kudriavzeviiYEL054C"]
    new.summary["N9_mikataeYEL054C"] <- summary["N7_mikataeYEL054C"]
    new.summary["N9_N10"] <- summary["N7_N8"]
    new.summary["N10_paradoxusYEL054C"] <- summary["N8_paradoxusYEL054C"]
    new.summary["N10_cerevisiaeYEL054C"] <- summary["N8_cerevisiaeYEL054C"]
  }
  else if (localTree == 4){
    new.summary[vector.names] <- summary[vector.names]
  }
  return (new.summary)
}

for (IGC.geo in IGC.geo.list){
  max.lnL.summary <- NULL
  vector.names = rownames(get(paste("PAML_estimatedTau_", paste(toString(IGC.geo), ".0", sep = ""), "_1stTree_summary", sep = "")))
  for (sim.num in 1:100){
    max.pos <- which.max(c(get(paste("PAML_estimatedTau_", paste(toString(IGC.geo), ".0", sep = ""), "_LocalTree_1_summary", sep = ""))["ll", sim.num],
                           get(paste("PAML_estimatedTau_",paste(toString(IGC.geo), ".0", sep = ""), "_LocalTree_2_summary", sep = ""))["ll", sim.num],
                           get(paste("PAML_estimatedTau_",paste(toString(IGC.geo), ".0", sep = ""), "_LocalTree_3_summary", sep = ""))["ll", sim.num],
                           get(paste("PAML_estimatedTau_", paste(toString(IGC.geo), ".0", sep = ""), "_1stTree_summary", sep = ""))["ll", sim.num]))
    if (max.pos < 4){
      new.summary <- convert.summary(max.pos, get(paste("PAML_estimatedTau", paste(toString(IGC.geo), ".0", sep = ""), "LocalTree", toString(max.pos), "summary", sep = "_"))[, sim.num], vector.names)
    }
    else{
      new.summary <- get(paste("PAML_estimatedTau_", paste(toString(IGC.geo), ".0", sep = ""), "_1stTree_summary", sep = ""))[, sim.num]
    }
    max.lnL.summary <- cbind(max.lnL.summary, new.summary)
  }
  
  rownames(max.lnL.summary) <- rownames(get(paste("PAML_estimatedTau_", paste(toString(IGC.geo), ".0", sep = ""), "_1stTree_summary", sep = "")))
  colnames(max.lnL.summary) <- paste("Sim", 1:dim(max.lnL.summary)[2], sep = "_")
  assign(paste("geo", paste(toString(IGC.geo), ".0", sep = ""), "PAML", "maxlnL", "summary", sep = "_"), max.lnL.summary)
}

###############################################################
###############################################################
###############################################################
###############################################################

# From PAMLsummary.Rmd
library("ape")
for (IGC.geo in IGC.geo.list){
  PAML.local.tree.1 <- NULL
  PAML.local.tree.2 <- NULL
  PAML.local.tree.3 <- NULL
  tree.local.1 <- read.tree(file = paste( "/Users/xji3/GitFolders/IGCSimulation/YDR418W_YEL054C_estimatedTau/IGCgeo_", toString(IGC.geo), ".0/sim_0/unrooted_MG94_geo_", toString(IGC.geo), ".0_Sim_0_localTree_1_codeml_est.newick", sep = ""))
  edge.local.1 <- tree.local.1$edge
  check.1 <- TRUE
  tree.local.2 <- read.tree(file = paste( "/Users/xji3/GitFolders/IGCSimulation/YDR418W_YEL054C_estimatedTau/IGCgeo_", toString(IGC.geo), ".0/sim_0/unrooted_MG94_geo_", toString(IGC.geo), ".0_Sim_0_localTree_2_codeml_est.newick", sep = ""))
  edge.local.2 <- tree.local.2$edge
  check.2 <- TRUE
  tree.local.3 <- read.tree(file = paste( "/Users/xji3/GitFolders/IGCSimulation/YDR418W_YEL054C_estimatedTau/IGCgeo_", toString(IGC.geo), ".0/sim_0/unrooted_MG94_geo_", toString(IGC.geo), ".0_Sim_0_localTree_3_codeml_est.newick", sep = ""))
  edge.local.3 <- tree.local.3$edge
  check.3 <- TRUE
  
  for (sim.num in 0:99){
    newick.tree.file <- paste("/Users/xji3/GitFolders/IGCSimulation/YDR418W_YEL054C_estimatedTau/IGCgeo_", toString(IGC.geo), ".0/sim_", toString(sim.num),"/unrooted_MG94_geo_", toString(IGC.geo), ".0_Sim_", toString(sim.num), "_localTree_1_codeml_est.newick", sep = "") 
    tree.local.1 <- read.tree(file = newick.tree.file)
    check.1 <- check.1 & all(edge.local.1 == tree.local.1$edge)
    PAML.local.tree.1 <- cbind(PAML.local.tree.1, tree.local.1$edge.length)
    
    newick.tree.file <- paste("/Users/xji3/GitFolders/IGCSimulation/YDR418W_YEL054C_estimatedTau/IGCgeo_", toString(IGC.geo), ".0/sim_", toString(sim.num),"/unrooted_MG94_geo_", toString(IGC.geo), ".0_Sim_", toString(sim.num), "_localTree_2_codeml_est.newick", sep = "") 
    tree.local.2 <- read.tree(file = newick.tree.file)
    check.2 <- check.2 & all(edge.local.2 == tree.local.2$edge)
    PAML.local.tree.2 <- cbind(PAML.local.tree.2, tree.local.2$edge.length)
    
    newick.tree.file <- paste("/Users/xji3/GitFolders/IGCSimulation/YDR418W_YEL054C_estimatedTau/IGCgeo_", toString(IGC.geo), ".0/sim_", toString(sim.num),"/unrooted_MG94_geo_", toString(IGC.geo), ".0_Sim_", toString(sim.num), "_localTree_3_codeml_est.newick", sep = "") 
    tree.local.3 <- read.tree(file = newick.tree.file)
    check.3 <- check.3 & all(edge.local.3 == tree.local.3$edge)
    PAML.local.tree.3 <- cbind(PAML.local.tree.3, tree.local.3$edge.length)
  }
  
  # Now assign rownames
  node.names.1 <- c("cerevisiaeYDR418W", "paradoxusYDR418W", "mikataeYDR418W", "kudriavzeviiYDR418W",
                    "bayanusYDR418W", "castelliiYDR418W", "cerevisiaeYEL054C", "paradoxusYEL054C", 
                    "mikataeYEL054C", "kudriavzeviiYEL054C", "bayanusYEL054C", "castelliiYEL054C",
                    "kluyveriYDR418W",
                    "N0", "N1", "N2", "N3", "N4", "N5",
                    "N6", "N7", "N8", "N9")
  node.names.2 <- c("cerevisiaeYDR418W", "paradoxusYDR418W", "mikataeYDR418W", "kudriavzeviiYDR418W",
                    "bayanusYDR418W", "castelliiYDR418W", "cerevisiaeYEL054C", "paradoxusYEL054C", 
                    "mikataeYEL054C", "kudriavzeviiYEL054C", "bayanusYEL054C", "castelliiYEL054C",
                    "kluyveriYDR418W",
                    "N0", "N1", "N2", "N3", "N4", "N5",
                    "N6", "N7", "N8", "N9")
  node.names.3 <- c("cerevisiaeYDR418W", "paradoxusYDR418W", "mikataeYDR418W", "kudriavzeviiYDR418W",
                    "bayanusYDR418W", "castelliiYDR418W", "cerevisiaeYEL054C", "paradoxusYEL054C", 
                    "mikataeYEL054C", "kudriavzeviiYEL054C", "bayanusYEL054C", "castelliiYEL054C",
                    "kluyveriYDR418W",
                    "N0", "N1", "N2", "N3", "N4", "N5", "N6",
                    "N7", "N8")
  local.row.names.1 <- NULL
  for (i in 1:dim(edge.local.1)[1]){
    local.row.names.1 <- c(local.row.names.1, paste(node.names.1[edge.local.1[i, 1]], node.names.1[edge.local.1[i, 2]], sep = "_"))
  }
  local.row.names.2 <- NULL
  for (i in 1:dim(edge.local.2)[1]){
    local.row.names.2 <- c(local.row.names.2, paste(node.names.2[edge.local.2[i, 1]], node.names.2[edge.local.2[i, 2]], sep = "_"))
  }
  local.row.names.3 <- NULL
  for (i in 1:dim(edge.local.3)[1]){
    local.row.names.3 <- c(local.row.names.3, paste(node.names.3[edge.local.3[i, 1]], node.names.3[edge.local.3[i, 2]], sep = "_"))
  }
  rownames(PAML.local.tree.1) <- local.row.names.1
  rownames(PAML.local.tree.2) <- local.row.names.2
  rownames(PAML.local.tree.3) <- local.row.names.3
  assign(paste("ape_paml_localTree_1_geo_", toString(IGC.geo), ".0_blen_summary", sep = ""), PAML.local.tree.1)
  assign(paste("ape_paml_localTree_2_geo_", toString(IGC.geo), ".0_blen_summary", sep = ""), PAML.local.tree.2)
  assign(paste("ape_paml_localTree_3_geo_", toString(IGC.geo), ".0_blen_summary", sep = ""), PAML.local.tree.3)
}

for (IGC.geo in IGC.geo.list){
  PAML.tree.1 <- NULL
  PAML.tree.2 <- NULL
  
  tree.1 <- read.tree(file = paste("/Users/xji3/GitFolders/IGCSimulation/YDR418W_YEL054C_estimatedTau/IGCgeo_", toString(IGC.geo), ".0/sim_0/unrooted_MG94_geo_", toString(IGC.geo), ".0_Sim_0_codeml_tree1_est.newick", sep = ""))
  edge.1 <- tree.1$edge
  check.1 <- TRUE
  
  tree.2 <- read.tree(file = paste("/Users/xji3/GitFolders/IGCSimulation/YDR418W_YEL054C_estimatedTau/IGCgeo_", toString(IGC.geo), ".0/sim_0/unrooted_MG94_geo_", toString(IGC.geo), ".0_Sim_0_codeml_tree2_est.newick", sep = ""))
  edge.2 <- tree.2$edge
  check.2 <- TRUE
  
  for (sim.num in 0:99){
    newick.tree.file <- paste("/Users/xji3/GitFolders/IGCSimulation/YDR418W_YEL054C_estimatedTau/IGCgeo_", toString(IGC.geo), ".0/sim_", toString(sim.num),"/unrooted_MG94_geo_", toString(IGC.geo), ".0_Sim_", toString(sim.num), "_codeml_tree1_est.newick", sep = "") 
    tree.1 <- read.tree(file = newick.tree.file)
    check.1 <- check.1 & all(edge.1 == tree.1$edge)
    PAML.tree.1 <- cbind(PAML.tree.1, tree.1$edge.length)
    
    newick.tree.file <- paste("/Users/xji3/GitFolders/IGCSimulation/YDR418W_YEL054C_estimatedTau/IGCgeo_", toString(IGC.geo), ".0/sim_", toString(sim.num),"/unrooted_MG94_geo_", toString(IGC.geo), ".0_Sim_", toString(sim.num), "_codeml_tree2_est.newick", sep = "") 
    tree.2 <- read.tree(file = newick.tree.file)
    check.2 <- check.2 & all(edge.2 == tree.2$edge)
    PAML.tree.2 <- cbind(PAML.tree.2, tree.2$edge.length)
  }
  # Now assign rownames
  node.names.1 <- c("cerevisiaeYDR418W", "paradoxusYDR418W", "mikataeYDR418W", "kudriavzeviiYDR418W",
                    "bayanusYDR418W", "castelliiYDR418W", "cerevisiaeYEL054C", "paradoxusYEL054C", 
                    "mikataeYEL054C", "kudriavzeviiYEL054C", "bayanusYEL054C", "castelliiYEL054C",
                    "kluyveriYDR418W",
                    "N0", "N1", "N2", "N3", "N4", "N5",
                    "N6", "N7", "N8", "N9", "N10")
  node.names.2 <- c("kluyveriYDR418W", "castelliiYEL054C", "bayanusYEL054C", "kudriavzeviiYEL054C",
                    "mikataeYEL054C","cerevisiaeYEL054C", "paradoxusYEL054C", "castelliiYDR418W",
                    "bayanusYDR418W", "kudriavzeviiYDR418W", "mikataeYDR418W", "cerevisiaeYDR418W", 
                    "paradoxusYDR418W",
                    "N0", 
                    "N1", "N2", "N3", "N4", "N5",
                    "N6", "N7", "N8", "N9", "N10")
  
  row.names.1 <- NULL
  for (i in 1:dim(edge.1)[1]){
    row.names.1 <- c(row.names.1, paste(node.names.1[edge.1[i, 1]], node.names.1[edge.1[i, 2]], sep = "_"))
  }
  row.names.2 <- NULL
  for (i in 1:dim(edge.2)[1]){
    row.names.2 <- c(row.names.2, paste(node.names.2[edge.2[i, 1]], node.names.2[edge.2[i, 2]], sep = "_"))
  }
  rownames(PAML.tree.1) <- row.names.1
  rownames(PAML.tree.2) <- row.names.2
  assign(paste("ape_paml_Tree_1_geo_", toString(IGC.geo), ".0_blen_summary", sep = ""), PAML.tree.1)
  assign(paste("ape_paml_Tree_2_geo_", toString(IGC.geo), ".0_blen_summary", sep = ""), PAML.tree.2)
}


################################################################################
# New
true.blen <- setNames(c(0.0197240946542, 0.215682181791, 0.20925129872, 0.171684721483,
                        0.0257112589202, 0.0266075664688, 0.0321083243449, 0.0853588718458,
                         0.024947887926, 0.0566627496729, 0.0581451177847, 0.0218788166581), 
                      c("N0_N1", "N0_kluyveri", "N1_N2", "N1_castellii", "N2_N3",
                        "N2_bayanus", "N3_N4", "N3_kudriavzevii", "N4_N5",
                        "N4_mikatae", "N5_cerevisiae", "N5_paradoxus") )

IGC.branch.list <- c("(N0,N1)",
                 "(N1,N2)",
                 "(N1,castellii)",
                 "(N2,N3)",
                 "(N2,bayanus)",
                 "(N3,N4)",
                 "(N3,kudriavzevii)",
                 "(N4,N5)",
                 "(N4,mikatae)",
                 "(N5,cerevisiae)",
                 "(N5,paradoxus)")

IGC.assign.branch.list <- c("N0.N1",
                     "N1.N2",
                     "N1.castellii",
                     "N2.N3",
                     "N2.bayanus",
                     "N3.N4",
                     "N3.kudriavzevii",
                     "N4.N5",
                     "N4.mikatae",
                     "N5.cerevisiae",
                     "N5.paradoxus")

true.branch.list <- c("N0_N1",
                            "N1_N2",
                            "N1_castellii",
                            "N2_N3",
                            "N2_bayanus",
                            "N3_N4",
                            "N3_kudriavzevii",
                            "N4_N5",
                            "N4_mikatae",
                            "N5_cerevisiae",
                            "N5_paradoxus")

PAML.paralog1.branch.list <- c("N0_N1",
                     "N1_N2",
                     "N1_castelliiYDR418W",
                     "N2_N3",
                     "N2_bayanusYDR418W",
                     "N3_N4",
                     "N3_kudriavzeviiYDR418W",
                     "N4_N5",
                     "N4_mikataeYDR418W",
                     "N5_cerevisiaeYDR418W",
                     "N5_paradoxusYDR418W")

PAML.paralog2.branch.list <- c("N0_N6",
                               "N6_N7",
                               "N6_castelliiYEL054C",
                               "N7_N8",
                               "N7_bayanusYEL054C",
                               "N8_N9",
                               "N8_kudriavzeviiYEL054C",
                               "N9_N10",
                               "N9_mikataeYEL054C",
                               "N10_cerevisiaeYEL054C",
                               "N10_paradoxusYEL054C")

## Mar 11 2016 change
## from IGCSimulationSummary.Rmd
#Now summarized % changes due to IGC in simulation estimates
#
#It was estimated 0.2131 (0.02884)
percent.IGC.geo.3.0 <- colSums(geo_3.0_summary[34:45, ] + geo_3.0_summary[46:57, ]) / (colSums(geo_3.0_summary[34:45, ] + geo_3.0_summary[46:57, ] + geo_3.0_summary[58:69, ]))
percent.IGC.geo.10.0 <- colSums(geo_10.0_summary[34:45, ] + geo_10.0_summary[46:57, ]) / (colSums(geo_10.0_summary[34:45, ] + geo_10.0_summary[46:57, ] + geo_10.0_summary[58:69, ]))
percent.IGC.geo.50.0 <- colSums(geo_50.0_summary[34:45, ] + geo_50.0_summary[46:57, ]) / (colSums(geo_50.0_summary[34:45, ] + geo_50.0_summary[46:57, ] + geo_50.0_summary[58:69, ]))
percent.IGC.geo.100.0 <- colSums(geo_100.0_summary[34:45, ] + geo_100.0_summary[46:57, ]) / (colSums(geo_100.0_summary[34:45, ] + geo_100.0_summary[46:57, ] + geo_100.0_summary[58:69, ]))
percent.IGC.geo.200.0 <- colSums(geo_200.0_summary[34:45, ] + geo_200.0_summary[46:57, ]) / (colSums(geo_200.0_summary[34:45, ] + geo_200.0_summary[46:57, ] + geo_200.0_summary[58:69, ]))
percent.IGC.geo.300.0 <- colSums(geo_300.0_summary[34:45, ] + geo_300.0_summary[46:57, ]) / (colSums(geo_300.0_summary[34:45, ] + geo_300.0_summary[46:57, ] + geo_300.0_summary[58:69, ]))
percent.IGC.geo.400.0 <- colSums(geo_400.0_summary[34:45, ] + geo_400.0_summary[46:57, ]) / (colSums(geo_400.0_summary[34:45, ] + geo_400.0_summary[46:57, ] + geo_400.0_summary[58:69, ]))
percent.IGC.geo.500.0 <- colSums(geo_500.0_summary[34:45, ] + geo_500.0_summary[46:57, ]) / (colSums(geo_500.0_summary[34:45, ] + geo_500.0_summary[46:57, ] + geo_500.0_summary[58:69, ]))