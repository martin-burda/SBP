
options(repos = c(CRAN="http://cran.rstudio.com"))
R_packages <- c("librarian", "bookdown", "copula", "evd", "ggplot2", "gridExtra", "ks",
                "plotly", "dplyr", "ggpubr", "grid", "cowplot", "parallel", "rmarkdown")
librarian::shelf(R_packages)


## ----------------------------------------------------------------------------------------------------------------------------
theta <- c(5, 15, 10) 
w <- c(0.4, 0.4, 0.2)
gsz <- 100 # grid size
dim <- 3 # dimensions of copula
rname <- "mixture-hd31"


## ----------------------------------------------------------------------------------------------------------------------------
cop1 <- frankCopula(param = theta[1], dim = dim) 
cop2 <- claytonCopula(param = theta[2], dim = dim)
cop3 <- gumbelCopula(param = theta[3], dim = dim)


## ----------------------------------------------------------------------------------------------------------------------------
# sequence of points in each dimension
seqs <- replicate(dim, seq(0.01, 0.99, length.out = gsz), simplify = FALSE)

# create names for the dimensions
names(seqs) <- paste0("u", 1:dim)

# construct the grid
grid <- expand.grid(seqs)
grid_mat <- as.matrix(grid)
rm(grid)


## ----------------------------------------------------------------------------------------------------------------------------
compute_chunk <- function(idx) {
  subgrid <- grid_mat[idx, , drop = FALSE]
  d1 <- dCopula(subgrid, cop1)
  d2 <- dCopula(subgrid, cop2)
  d3 <- dCopula(subgrid, cop3)
  w[1] * d1 + w[2] * d2 + w[3] * d3
}


## ----------------------------------------------------------------------------------------------------------------------------
n <- nrow(grid_mat)
core_max <- 8 # maximum number of cores used
nchunks <- min(detectCores(), core_max)  # cap at core_max
chunks <- split(seq_len(n), cut(seq_len(n), nchunks, labels = FALSE))


## ----------------------------------------------------------------------------------------------------------------------------
start <- proc.time()
if (.Platform$OS.type == "unix") {
  mix_list <- mclapply(chunks, compute_chunk, mc.cores = nchunks)
} else {
  cl <- makeCluster(nchunks)
  clusterEvalQ(cl, library(copula))  # load copula package on workers
  clusterExport(cl, c("grid_mat", "cop1", "cop2", "cop3", "w", "compute_chunk"),
                envir = environment())
  mix_list <- parLapply(cl, chunks, compute_chunk)
  stopCluster(cl)
}
end <- proc.time()
elapsed <- end - start
elapsed

#user – CPU time spent in R
#system – CPU time in system calls
#elapsed – wall-clock time


## ----------------------------------------------------------------------------------------------------------------------------
mix_dens <- unlist(mix_list)
rm(mix_list)


## ----------------------------------------------------------------------------------------------------------------------------
n  <- 1000 # Sample size
n1 <- round(n*w[1])
n2 <- round(n*w[2])
n3 <- n - (n1 + n2)

set.seed(123456)
u1 <- rCopula(n1, cop1) # Sample from each copula
u2 <- rCopula(n2, cop2)
u3 <- rCopula(n3, cop3)

cs <- rbind(u1, u2, u3) # Combine the samples to form the mixture


## ----------------------------------------------------------------------------------------------------------------------------
plot_lower_triangular <- function(cs_df, title = NULL,
                                  point_alpha = 0.4, point_size = 0.2) {
  d <- ncol(cs_df)
  vars <- colnames(cs_df)
  if (d < 2) stop("Need at least 2 columns for pairwise plots.")
  
  my_scatter <- function(data, xvar, yvar) {
    ggplot(data, aes(x = {{xvar}}, y = {{yvar}})) +
      geom_point(alpha = point_alpha, size = point_size) +
      coord_fixed(ratio = 1) +
      labs(x = "", y = "") +
      theme_minimal() +
      theme(
        plot.margin  = unit(c(0,0,0,0), "in"),
        panel.grid   = element_blank(),
        axis.text    = element_blank(),
        axis.ticks   = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5)
      )
  }
  
  empty_panel <- ggplot() + theme_void()
  
  # build full (d-1) x (d-1) matrix
  plot_list <- vector("list", (d-1)*(d-1))
  
  for (i in 2:d) {
    for (j in 1:(d-1)) {
      idx <- (i-2)*(d-1) + j
      if (j < i) {
        plot_list[[idx]] <- my_scatter(cs_df, .data[[vars[j]]], .data[[vars[i]]])
      } else {
        plot_list[[idx]] <- empty_panel
      }
    }
  }
  
  combined_plot <- ggarrange(
    plotlist = plot_list,
    ncol = d-1,
    nrow = d-1,
    align = "hv"
  )
  
  if (!is.null(title)) {
    combined_plot <- annotate_figure(
      combined_plot,
      top = text_grob(title, size = 8)
    )
  }
  
  return(combined_plot)
}


## ----fig.width=4, fig.height=4, fig.cap = "&nbsp;", warning=FALSE------------------------------------------------------------
cs_df <- as.data.frame(cs)
colnames(cs_df) <- paste0("U", 1:dim)
csplot <- plot_lower_triangular(cs_df, title="Mixture 1") 
csplot <- ggdraw(csplot) +
          theme(plot.background = element_rect(color = "gray", fill = NA, size = 0.5),
                plot.margin = unit(c(0., 0.1, -0.1, -0.1), "in"))

plot_file = paste0("../graphs_sim/",rname,"_true.png")
ggsave(plot_file, plot = csplot, width = 4, height = 4, units = "in", dpi = 300)

#csplot


## ----------------------------------------------------------------------------------------------------------------------------
true_dens <- mix_dens
sample_sizes <- c(100, 200, 300, 400, 500, 1000, 10000)
kmax <- 10 # number of batches for each sample size

## ----------------------------------------------------------------------------------------------------------------------------
set.seed(123456)

for (n in sample_sizes) {

  n1 <- round(n*w[1])
  n2 <- round(n*w[2])
  n3 <- n - (n1 + n2)
  
  for (k in 1:kmax) {
  
    u1 <- rCopula(n1, cop1) # Sample from each copula
    u2 <- rCopula(n2, cop2)
    u3 <- rCopula(n3, cop3)
    
    u_mix <- rbind(u1, u2, u3) # Combine the samples to form the mixture
    
    # Write out for Fortran input
    write.table(u_mix, file = paste0("../data_sim/",rname,"_dat_",n,"_",k,".csv"),
                sep = ",", row.names = FALSE, col.names = FALSE)
  }
}

# Run Fortran code

## ----------------------------------------------------------------------------------------------------------------------------
n  <- 10000 # for sample size
k <- 1 # for data from batch k
csim_SBP <- read.table(file = paste0("../output_sim/",rname,"_csim_",n,"_",k,".out"), header = FALSE)
csim_SBP_df <- as.data.frame(csim_SBP)
colnames(csim_SBP_df) <- paste0("U", seq_len(ncol(csim_SBP_df)))

csplot_sim <- plot_lower_triangular(csim_SBP_df, title="SBP copula") 
csplot_sim <- ggdraw(csplot_sim) +
                theme(plot.background = element_rect(color = "gray", fill = NA, size = 0.5),
                plot.margin = unit(c(0., 0.1, -0.1, -0.1), "in"))

plot_file_sim = paste0("../graphs_sim/",rname,"_SBP.png")
ggsave(plot_file_sim, plot = csplot_sim, width = 4, height = 4, units = "in", dpi = 300)

#csplot


## ----------------------------------------------------------------------------------------------------------------------------
start <- proc.time()
KL <- numeric() 
KLse <- numeric()
epsilon <- 1e-10 # Add small epsilon to avoid division by 0

# Make a list of tasks (one per n,k)
task_list <- expand.grid(n = sample_sizes, k = 1:kmax, KEEP.OUT.ATTRS = FALSE)

compute_KL <- function(task) {
  n <- task$n
  k <- task$k
  
  # Read density estimate from file
  fname <- paste0("../output_sim/", rname, "_cdens_", n, "_", k, ".out")
  cdens_SBP <- as.numeric(unlist(read.table(file = fname, header = FALSE)))
  
  # Compute KL contribution
  ratio <- (true_dens + epsilon) / (cdens_SBP + epsilon)
  KL_k_val <- log(ratio) * true_dens
  mean(KL_k_val)
}

# Split tasks across cores
ncores <- min(detectCores(), 8)
tasks_split <- split(seq_len(nrow(task_list)), cut(seq_len(nrow(task_list)), ncores, labels = FALSE))

# Parallel execution
if (.Platform$OS.type == "unix") {
  KL_vals <- unlist(mclapply(tasks_split, function(idx) {
    sapply(idx, function(i) compute_KL(task_list[i,]))
  }, mc.cores = ncores))
} else {
  cl <- makeCluster(ncores)
  clusterExport(cl, c("task_list", "true_dens", "rname", "epsilon", "compute_KL"),
                envir = environment())
  KL_vals <- unlist(parLapply(cl, tasks_split, function(idx) {
    sapply(idx, function(i) compute_KL(task_list[i,]))
  }))
  stopCluster(cl)
}

# Combine results back by sample size
for (n in sample_sizes) {
  KL_k <- KL_vals[task_list$n == n]
  KL    <- c(KL, mean(KL_k))
  KLse  <- c(KLse, sd(KL_k) / sqrt(length(KL_k)))
}

end <- proc.time()
elapsed <- end - start
elapsed

print(round(KL, 2))
print(round(KLse, 2))

# Knit 'myfile.Rmd' into HTML
render(paste0(rname,".R"), output_format = "html_document")
