
library(vegan)
library(ggplot2)
library(plotly)
library(reshape2)
library(stringr)

##rarefaction curve
rarecurve_mine <- function(otu, step = 10, sample = 400) {
  otu <- t(otu)
  out <- rarecurve(otu,
                   step = step,
                   sample = sample,
                   label = FALSE)
  names(out) <- paste("species", 1:dim(otu)[1], sep = "")
  
  # Coerce data into "long" form.
  protox <- mapply(
    FUN = function(x, y) {
      mydf <- as.data.frame(x)
      colnames(mydf) <- "value"
      mydf$species <- y
      mydf$subsample <- attr(x, "Subsample")
      mydf
    },
    x = out,
    y = as.list(names(out)),
    SIMPLIFY = FALSE
  )
  
  xy <- do.call(rbind, protox)
  rownames(xy) <- NULL  # pretty
  
  # Plot.
  #ggplot(xy, aes(x = subsample, y = value, color = species)) +theme_bw() +scale_color_discrete(guide = FALSE) +   geom_line()
  
  
  ggplotly(ggplot(xy, aes(
    x = subsample, y = value, color = species
  )) +
    theme_bw() +
    theme(legend.position = "none") +  # ggplotly doesn't respect scales?
    geom_line())
}





### niche breadth and niche overlap
Bcom<-function(otu, filter=2e-5){
  ##as described in From surviving to thriving, the assembly processes of microbial communities in stone biodeterioration: A case study of the West Lake UNESCO World Heritage area in China
  ##doi.org/10.1016/j.scitotenv.2021.150395
  ##otu_table, row as species, col as community
  t_otu <- t(otu)
  t_otu <- t_otu[, colMeans(t_otu / (rowSums(t_otu))) > filter]
  otu_niche_breadth <- spaa::niche.width(t_otu, "levins")
  otu_binary <- t(t_otu)
  otu_binary[otu_binary > 0] <- 1
  otu_binary <- otu_binary * as.numeric(otu_niche_breadth)
  return(colMeans(otu_binary, na.rm = T))
}
niche.overlap.noparal<-function (t_otu) {
  mat<-t_otu
  match.arg(method)
  mat <- na.omit(mat)
  result <- matrix(0, as.numeric(ncol(mat)) , as.numeric(ncol(mat)))

  for (i in 1:(ncol(mat) - 1)) {
    for (j in (i + 1):ncol(mat)) {
                  result[j, i] <- spaa::niche.overlap.pair(mat[,i], mat[, j], method = "morisita")
    }
  }
  rownames(result) <- colnames(mat)
  colnames(result) <- colnames(mat)
  return(result)
}

niche.overlap.paral<-function(t_otu_table, cores=12){ 
  t_otu_table <- na.omit(t_otu_table)
  otu_number<-dim(t_otu_table)[2]
  
  library("foreach")
  library("doParallel")
  cl <- makeCluster(cores)
  registerDoParallel(cl)    #进行进程注册
  niche_overlap_matrix <- foreach(i = 1:(otu_number - 1),
                                  #输入等待请求的参数
                                  .combine = rbind,
                                  .verbose = F) %dopar% {
                                    mat <- (matrix(ncol = 3, nrow = otu_number))  # construct first, add element then 1.9s
                                    colnames(mat) <- c("i", "j", "value")
                                    print(i)
                                    for (j in (i + 1):otu_number) {
                                      mat[j, ] <-
                                        c(i,
                                          j,
                                          spaa::niche.overlap.pair(t_otu_table[, i], t_otu_table[, j], method = "morisita"))
                                    }
                                    
                                    mat <- mat[!is.na(mat[, 3]), ]
                                    # for (j in 1:i ){
                                    #   mat[j,]<-c(i,j,  0) ## provide
                                    # }
                                    if (i != (otu_number - 1)) {
                                      mat <- as.data.frame(mat)
                                    } else {
                                      mat <- as.data.frame(t(mat))
                                    }
                                    
                                    mat[, 1:2] <- apply(mat[, 1:2], 2, as.integer)
                                    gc()
                                    return(mat)
                                  }
  gc()
  stopCluster(cl)
  niche_overlap_matrix<-Matrix::sparseMatrix(i = niche_overlap_matrix[, 2], j = niche_overlap_matrix[, 1], x =
                                                 niche_overlap_matrix[, 3] , triangular = T)
  rownames(niche_overlap_matrix) <- colnames(t_otu_table)
  colnames(niche_overlap_matrix) <- colnames(t_otu_table)
  gc()
  return(niche_overlap_matrix)
}

Ocom<-function (otu_table, method = "morisita",filter=2e-5, paral=T, cores=12) { #dist_matrix=NULL, get_matrix=F, 
  ##as described in From surviving to thriving, the assembly processes of microbial communities in stone biodeterioration: A case study of the West Lake UNESCO World Heritage area in China
  ##doi.org/10.1016/j.scitotenv.2021.150395
  ##otu_table, row as species, col as community  
  t_otu <- t(otu_table)
  t_otu <- t_otu[, colMeans(t_otu / (rowSums(t_otu))) > filter]
  
  if (paral) {
    niche_overlap_matrix <- niche.overlap.paral(t_otu, cores = cores)
    gc()
  } else {
    niche_overlap_matrix <- niche.overlap.noparal(t_otu)
  }
  print("niche_overlap_matrix calculation complete")

  overlap_list <- list()
  niche_over <- vector()
  for (i in 1:dim(t_otu)[1]) {
    # print(i) #######the number of calculated samples
    tmp <-
      niche_overlap_matrix[as.vector(t_otu[i, ] != 0), as.vector(t_otu[i, ] != 0)]
    overlap_list[[i]] <- tmp[lower.tri(tmp)]
    niche_over[i] <- mean(overlap_list[[i]])
  }
  names(niche_over) <- rownames(t_otu)
  return(niche_over)
}



###DOC calculation using the DOC r package
library(DOC)
results <- DOC(otu_bac[, target_sample], R = 1000)
(sum(results$LME$Slope >= 0) + 1) / (length(results$LME$Slope) + 1)
plot(results) + theme_bw()



###### betaNTI, RCI and diversity response to betaNTI, function codes in references
  Beta_NTI_bac <- Beta_NTI(phy_bac, (otu_bac))###get a dist about the betaNTI of pairwise community
  rcbray_bac <- raup_crick((otu_bac))###get a dist about the RCI of pairwise community
  
  dist_melt <- function(dist_a) {
    a <- dist_a
    a <- as.matrix(a)
    a[upper.tri(a)] <- 10000
    diag(a) <- 10000
    betamat <- melt(a)
    betamat <- betamat[!(betamat$value %in% c(10000)), ]
    return(betamat)
  }
  
  t_otu <- t(otu_bac)
  null_mat <- cbind(dist_melt(Beta_NTI_bac), dist_melt(rcbray_bac)[, c(3, 2)])[, 1:4]
  
  null_mat$assembly[null_mat$value <= -2] <- "Homo Selection"
  null_mat$assembly[null_mat$value >= 2] <- "Hetro Selection"
  null_mat$assembly[(abs(null_mat$value) < 2) &
                      null_mat$value.1 >= 0.95] <-
    "Dispersal Limitation"
  null_mat$assembly[(abs(null_mat$value) < 2) &
                      null_mat$value.1 <= -0.95] <-
    "Homogenizing Dispersal"
  null_mat$assembly[(abs(null_mat$value) < 2) &
                      null_mat$value.1 > -0.95 &
                      null_mat$value.1 < 0.95] <- "Drift"
  
  null_mat %>% group_by(assembly) %>% summarise(n = n())
  #null_mat$assembly
  
  null_mat <-
    merge(null_mat, group[, c(1, 2)], by.x = "Var1", by.y = "#sample")
  null_mat <-
    merge(null_mat, group[, c(1, 2)], by.x = "Var2", by.y = "#sample")
  a <- null_mat[null_mat$type_merge.x == null_mat$type_merge.y,]
  
  ###
  null_mat[(str_detect(null_mat$Var2, "N") |
              str_detect(null_mat$Var1, "N")),] %>% group_by(assembly) %>% summarise(n =
                                                                                       n())
  
  null_mat$mean1 <- alpha(t_otu)[null_mat$Var2, 2]
  null_mat$mean2 <- alpha(t_otu)[null_mat$Var1, 2]
  null_mat$mean <- (null_mat$mean2 + null_mat$mean1) / 2
  null_mat$mean_minus <- abs(null_mat$mean1 - null_mat$mean2)
  
  library(robustbase)
  summary(lmrob(data = null_mat, value ~ mean_minus))
  
  
  

  
  
##correlations between prokaryotic and fungal biomass and diversity library(psych)

  tmp_targetsample <- colnames(otu_bac)[5:37]
  qpcr_all <- qpcr[tmp_targetsample, c(22, 21)]
  qpcr_all <-
    cbind(qpcr_all, alpha(t(otu_bac[, tmp_targetsample]))$Shannon, alpha(t(otu_fun[, tmp_targetsample]))$Shannon)
  colnames(qpcr_all)[1:4] <-
    c("Prokaryotic Biomass",
      "Fungal Biomass",
      "Prokaryotic Diversity",
      "Fungal Diveristy")
  corr.test(qpcr_all, adjust = "none")
  
  
  
  
###network robustness evaluation 
  #random attack could be repeated for some times and combined
  
  adj_matrix<-as_adjacency_matrix(igraph)
  
  nc <- function(adj_matrix) {
    #  0-1 matrix, 1 represent existence of edgeï¼0 represent nonexistence
    adj_matrix <- as.matrix(adj_matrix)
    adj_matrix[abs(adj_matrix) != 0] <- 1
    
    # Î»
    lambda <- eigen(adj_matrix, only.values = TRUE)$values
    lambda <- sort(lambda, decreasing = TRUE)
    
    # natural connectivity
    lambda_sum <- 0
    N = length(lambda)
    for (i in 1:N) lambda_sum = lambda_sum + exp(lambda[i])
    lambda_average <- log(lambda_sum/N, base = exp(1))
    lambda_average
  }
  
  ##natural_connectivity <- nc(adj_matrix)
  ##natural_connectivity
  
  net_attack<-function(adj_matrix){
    
    #
    g <- graph_from_adjacency_matrix(adj_matrix, mode = 'undirected', diag = FALSE)
    
    order_degree<-order(degree(g),decreasing = T)
    order_betwe<-order(betweenness(g),decreasing = T)
    order_random<-sample(1:dim(adj_matrix)[1],dim(adj_matrix)[1])
    
    
    degree_dist <- table(degree(g))
    degree_num <- as.numeric(names(degree_dist))
    degree_count <- as.numeric(degree_dist)
    names(degree_count) <- degree_num
    degs <- rep(degree_num, degree_count)
    
    adj_matrix_rand1 <- adj_matrix
    adj_matrix_rand2 <- adj_matrix
    adj_matrix_rand3 <- adj_matrix
    
    natural_connectivity_rand1 <- nc(adj_matrix_rand1)
    natural_connectivity_rand2 <- nc(adj_matrix_rand2)
    natural_connectivity_rand3 <- nc(adj_matrix_rand3)
    
    for (i in 1:(dim(adj_matrix)[1]-1)) {
      
      #specify the removing node 
      adj_matrix_rand1_remove <- adj_matrix_rand1[-order_degree[1:i],-order_degree[1:i]]
      adj_matrix_rand2_remove <- adj_matrix_rand2[-order_betwe[1:i],-order_betwe[1:i]]
      adj_matrix_rand3_remove <- adj_matrix_rand3[-order_random[1:i],-order_random[1:i]]
      #
      natural_connectivity_rand1 <- c(natural_connectivity_rand1, nc(adj_matrix_rand1_remove))
      natural_connectivity_rand2 <- c(natural_connectivity_rand2, nc(adj_matrix_rand2_remove))
      natural_connectivity_rand3 <- c(natural_connectivity_rand3, nc(adj_matrix_rand3_remove))
      if(i%%10==0) {print(i)}
    }
  
    dat <- data.frame(remove_node = rep(1:dim(adj_matrix)[1],3),
                      natural_connectivity = c(natural_connectivity_rand1, natural_connectivity_rand2, natural_connectivity_rand3),
                      network = c( rep('Degree Attack', dim(adj_matrix)[1]), rep('Betweenness Attack', dim(adj_matrix)[1]), rep('Random Attack', dim(adj_matrix)[1])))

    ###custom
    dat$nc_per<-dat$natural_connectivity/max(dat$natural_connectivity)
    dat$remove_node_per<-dat$remove_node/max(dat$remove_node)
    return(dat)
    
  }
  
