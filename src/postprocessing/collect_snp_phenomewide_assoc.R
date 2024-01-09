library(tidyverse)
source("/home/rodrigo/01_repos/CardiacGWAS/analysis/loci_names.R")

g <- glue::glue

BASEDIR = "/home/rodrigo/tmp/snp_assocs/"

PARTITIONS = c("10n90", "15n85", "20n80", "25n75", "30n70")

RUNS_SETS = c(
  "original_mlruns", 
  "finetuned_mlruns", 
  "finetuned_2_mlruns"
)

# RUNS_SETS = c("finetuned_2_mlruns")

SEEDS = 1:10

process_snp_assocs <- function(path) {
  
  assocs = lapply(
    list.files(path, full.names = TRUE),
    function(file) {
      p_vals = read.table(gzfile(file), header=TRUE, sep=" ") %>% select(contains(".log10p")) 
      p_vals
    }
  )
  
  assocs_df <- bind_rows(lapply(assocs, function(x) x[1,]))
  
  rownames_df = sapply(strsplit(list.files(path), split = "-"), function(x) x[2])
  rownames_df = sapply(rownames_df, str_replace, ".gz", "")
  # rownames_df = sapply(rownames_df, function(x) ifelse(x %in% names(all_loci), all_loci[x], x))
  
  rownames(assocs_df) <- rownames_df
  colnames(assocs_df) <- sapply(
    colnames(assocs_df), 
    function(x) str_replace(x, ".log10p", "")
  )
  
  assocs_df
  
}

assocs_discovery_lst <- assocs_replication_lst <- list()
names_ <- character()

for (SEED in SEEDS) { for (RUNS in RUNS_SETS) { for (PARTITION in PARTITIONS) {
      #for (DATASET in c("replication", "discovery")) {
        
        logging::loginfo("seed: {SEED}, run {RUNS}, partition {PARTITION}" %>% g)
        # THRESHOLD = ifelse(DATASET=="discovery", 6, 1.3)
        d_results_dir = glue::glue("{BASEDIR}/seed_{SEED}/{RUNS}_{PARTITION}/discovery")
        r_results_dir = glue::glue("{BASEDIR}/seed_{SEED}/{RUNS}_{PARTITION}/replication")
        #assocs_discovery_df <- process_snp_assocs(path=d_results_dir)
        #saveRDS(assocs_discovery_df, glue::glue("{BASEDIR}/seed_{SEED}/{RUNS}_{PARTITION}/discovery.rds"))
        
        #assocs_replication_df <- process_snp_assocs(path=r_results_dir)
        #saveRDS(assocs_replication_df, glue::glue("{BASEDIR}/seed_{SEED}/{RUNS}_{PARTITION}/replication.rds"))
        
        assocs_replication_df <- readRDS("{r_results_dir}.rds" %>% g)
        assocs_discovery_df <- readRDS("{d_results_dir}.rds" %>% g)
        
        name = g("{RUNS}_{PARTITION}_seed_{SEED}")
        names_ <- c(names_, name)
        
        assocs_discovery_lst <- c(assocs_discovery_lst, list(assocs_discovery_df))        
        assocs_replication_lst <- c(assocs_replication_lst, list(assocs_replication_df))        
        
      #}
    }
  }
}

# names(assocs_discovery_lst) <- names_
# names(assocs_replication_lst) <- names_

# assocs_lst <- lapply(assocs_lst, function(x) { colnames(x) <- c("counts", "best_p"); x})

names_ <- names_ %>% gsub(pattern = "original_mlruns", replacement = "ENS1")
names_ <- names_ %>% gsub(pattern = "finetuned_2_mlruns", replacement = "ENS2")
names_ <- names_ %>% gsub(pattern = "finetuned_mlruns", replacement = "ENS3")

names(assocs_replication_lst) <- names_
names(assocs_discovery_lst) <- names_

# apply(X = assocs_df, MARGIN = 1, function(x) sum(x>6)) %>% .[order(-.)] %>% data.frame %>% View
get_counts_over_threshold <- function(df, threshold) {
  df %>% 
    apply(1, function(x) x > threshold) %>% 
    apply(2, sum) %>% 
    .[order(-.)] # order decreasingly
}

assoc_counts <- lapply(
  assocs_discovery_lst, 
  function(df) { 
    get_counts_over_threshold(df, 7.3) # %>% 
    # .[. > 5] %>% 
    # names %>% .[names(.) %in% sws_regions]
  }
)

assoc_best <- lapply(
  assocs_discovery_lst, 
  function(df) {
    apply(df, 1, max) %>% .[order(-.)]
  }    
) #%>% .[sapply(., length) > 3] %>% sapply(sum)

# THRESHOLD = 7.3
# 
# SEED = 2
# cbind(
#   assocs_discovery_lst[[g("{RUNS_SETS[1]}_15n85_seed_{SEED}")]],
#   assocs_discovery_lst[[g("{RUNS_SETS[2]}_15n85_seed_{SEED}")]],
#   assocs_discovery_lst[[g("{RUNS_SETS[3]}_15n85_seed_{SEED}")]]
# ) %>% apply(1, 
#   function(row) { sum(row > THRESHOLD) }
# ) %>% .[order(-.)] %>% .[sapply(names(.), function(x) {x %in% sws_loci})]

replace_region_name <- function(region_name) {
  ifelse(region_name %in% names(all_loci), all_loci[region_name], region_name)
}

get_results_df<- function(ensemble="original_mlruns", partition="30n70" ,seed="1", phase="discovery") {
  filename <- g("{BASEDIR}/seed_{seed}/{ensemble}_{partition}/{phase}.rds")
  # print(filename)sws
  readRDS(filename)
}

# ensemble="original_mlruns", partition="30n70" ,seed="1",
is_replicating <- function(discovery_df, replication_df, locus, n_best_runs=5) {
  
  
  phenotypes = sort(as.numeric(discovery_df["chr1_124",]), index.return=TRUE, decreasing=TRUE)$ix
  phenotypes = phenotypes[1:n_best_runs]
  # phenotype = names(which.max()); phenotype
  # print(phenotypes)  
  replication_df[locus, phenotypes]
  
}


which_loci_replicate <- function(ensemble="all", partition="30n70" ,seed="1", n_best_runs=5, threshold=1.3) {

  if (ensemble == "all") {
    discovery_df <- RUNS_SETS %>% lapply( function(x) get_results_df(ensemble = x, partition=partition, seed=seed, phase = "discovery")) %>% bind_cols()
    replication_df <- RUNS_SETS %>% lapply( function(x) get_results_df(ensemble = x, partition=partition, seed=seed, phase = "replication")) %>% bind_cols()
  } else {
    discovery_df <- get_results_df(ensemble, partition, seed, "discovery")
    replication_df <- get_results_df(ensemble, partition, seed, "replication")
  }
  
  sapply(
    rownames(discovery_df), 
    function(locus) {
      is_replicating(discovery_df, replication_df, locus, n_best_runs) %>% { . > threshold } %>% apply(1, any) 
    }
  )
}

# assoc_counts

assocs_discovery_lst %>% 
  lapply( function(x) get_counts_over_threshold(x, 7.3)) %>%  # get counts for each locus
  { nm = names(.); dd = bind_rows(.); rownames(dd) <- nm; dd } %>% # bind the list's elements into a dataframe
  { colnames(.) <- replace_region_name(colnames(.)); .} -> kk # replace region names for loci names when known
  # .[sapply(., function(x) dim(x)[1] != 0)] %>% # discard empty elements
  
  

# THRES_REPL = 1.3

replication_counts_per_seed <- function(ensemble, partition) {
  
  THRES_REPL = 1
  
  kk <- lapply(
   1:10, 
   function(seed) { 
     which_loci_replicate(
       ensemble=ensemble,
       partition=partition, 
       seed=seed, n_best_runs = 5, threshold = THRES_REPL
     ) %>% bind_rows
   }
  )
  
  kk <- bind_rows(kk) %>% t %>% { rownames(.) <- strsplit(rownames(.), "\\.")  %>% sapply(`[`, 1) %>% replace_region_name(); .} 
  kk <- kk[apply(kk, 1, sum) %>% .[order(-.)] %>% names, ] %>% .[all_loci,]
  kk
    
}
# kk[apply(kk, 1, sum) %>% .[order(-.)] %>% names, ]

replication_counts_per_seed("finetuned_mlruns", "10n90") %>% 
  .[setdiff(rownames(.), c("WAC", "LMO7", "CNOT7")),]

# Names of runs with MSE < THRESHOLD
MSE_THRESHOLD = 4
lapply(RUNS_SETS, function(x) g("~/01_repos/CardiacGWAS/results/mse-{x}-old_meshes.csv") %>% read.csv)[[1]] %>% 
  apply(2, mean) %>% 
  .[. < MSE_THRESHOLD] %>% 
  names %>% sapply(function(x) strsplit(x, split = "_")[[1]][1])


# lapply(RUNS_SETS, function(x) g("mse-{x}-old_meshes.csv") %>% read.csv)

# assocs_lst <- lapply(assocs_lst, function(x) { colnames(x) <- c("counts", "best_p"); x})