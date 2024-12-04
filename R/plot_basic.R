library(data.table)

args <- commandArgs(trailingOnly = TRUE)
file <- args[1]
libraryid <- args[2]
output_dir <- args[3]

if (!dir.exists(output_dir)) dir.create(output_dir)

df <- fread(file, header=T, sep="\t", fill=T, nThread=6)

# get the indices for different values we want to plot 

# k and n 
fwd_k_ix <- grep("^fwK", names(df))
fwd_n_ix <- grep("^fwN", names(df))
bwd_k_ix <- grep("^bwK", names(df))
bwd_n_ix <- grep("^bwN", names(df))

# x and f 
fwd_x_ix <- grep("^fwdx", names(df))
fwd_x_ix <- fwd_x_ix[!grepl("Conf", names(df)[fwd_x_ix])]
bwd_x_ix <- grep("^bwdx", names(df))
bwd_x_ix <- bwd_x_ix[!grepl("Conf", names(df)[bwd_x_ix])]
fwd_f_ix <- grep("^fwf", names(df))
bwd_f_ix <- grep("^bwf", names(df))

# x axis (how many positions)
positions <- seq_along(bwd_n_ix)

# damage params 
a_ix <- grep("^A_b", names(df))
q_ix <- grep("^q_b", names(df))
c_ix <- grep("^c_b", names(df))
phi_ix <- grep("^phi_b", names(df))

# lca info 
name_ix <- grep("name", names(df))
rank_ix <- grep("rank", names(df))
taxid_ix <- grep("taxid", names(df))

# reads and length 
nreads_ix <- grep("nreads", names(df))
rlen_ix <- grep("mean_rlen", names(df))


# do row by row 
for (i in 1:nrow(df)) {
  row_data <- df[i, ]
  
  # fancy title 
  title <- paste0(
    row_data[, ..rank_ix], ": ", row_data[, ..name_ix], 
    " (", row_data[, ..taxid_ix], ")\n",
    "A_b: ", row_data[, ..a_ix], ", q_b: ", row_data[, ..q_ix], 
    ", c_b: ", row_data[, ..c_ix], ", phi_b: ", row_data[, ..phi_ix], 
    "\nnreads: ", row_data[, ..nreads_ix], ", readlen: ", row_data[, ..rlen_ix]
  )

  # fancy filename 
  pdf_filename <- paste0(
    output_dir, "/", 
    libraryid, "_", 
    gsub("/| |\\(|\\)|\\.|:", "_",  row_data[, ..name_ix]), "_", 
    row_data[, ..taxid_ix], "_", 
    gsub("/| ", "_", row_data[, ..rank_ix]), ".pdf"
  )


  if (!is.na(row_data[,..a_ix]) # data must exist 
    & (row_data[,..rank_ix] %in% c("species","genus") | row_data[,..name_ix == "root"]) # only leaf or leafish nodes 
    & row_data[,..nreads_ix] >= 10) { # at least some reads 

    pdf(pdf_filename, width = 8, height = 6)

    par(mfrow = c(2, 1))

    # plot n and k 
    # n = total number of (eg) C 
    # k = total number of (eg) C-T
    k_fwd <- row_data[,..fwd_k_ix]
    n_fwd <- row_data[,..fwd_n_ix]
    k_bwd <- row_data[,..bwd_k_ix]
    n_bwd <- row_data[,..bwd_n_ix]
    
    plot(positions, k_fwd, type = "l", col = "red", xlab = "Position", ylab = "Number of nucleotides", ylim = c(0,max(n_fwd, n_bwd)), main=title)
    lines(positions, k_bwd, type = "l", col = "blue")
    lines(positions, n_fwd, type = "o", col = "red")
    lines(positions, n_bwd, type = "o", col = "blue")
    legend("right",              
           legend = c("k_fwd", "k_bwd", "n_fwd", "n_bwd"),
           col = c("red", "blue", "red", "blue"),        
           lty = c(1,1,1,1),                                      
           pch = c(NA,NA,1,1),
           cex = 0.8)

    # plot damage 
    # x is smoothed damage est 
    # f is point estimates  
    x_fwd <- row_data[,..fwd_x_ix]
    x_bwd <- row_data[,..bwd_x_ix]
    f_fwd <- row_data[,..fwd_f_ix]
    f_bwd <- row_data[,..bwd_f_ix]

    plot(positions, f_fwd, type = "o", col = "red", xlab = "Position", ylab = "Damage estimate", ylim = c(0,max(x_fwd, x_bwd)))
    lines(positions, f_bwd, type = "o", col = "blue")  
    lines(positions, x_bwd, type = "l", col = "red")  
    lines(positions, x_bwd, type = "l", col = "blue")  
    legend("topright",              
           legend = c("f_fwd", "f_bwd", "x_fwd", "x_bwd"),
           col = c("red", "blue", "red", "blue"),        
           lty = c(1,1,1,1),                                      
           pch = c(1,1,NA,NA),
           cex = 0.8)  

    dev.off()  
  }

}
