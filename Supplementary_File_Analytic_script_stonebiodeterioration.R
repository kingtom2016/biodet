
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



####FEAST

library(FEAST)
metadata_my <- metadata

FEAST_output <- data.frame(stringsAsFactors = F)
for (i in 1:dim(metadata_my)[1]) {
  metadata_new <- data.frame(stringsAsFactors = F)
  a <- (metadata_my[metadata_my$SourceSink == "Sink", ][i, ])
  a[, 3] <- i
  
  b <- metadata_my[metadata_my$SourceSink == "Source", ]
  b[, 3] <- i
  a <- rbind(a, b)
  a$Env <- paste0("sample", 1)###è¿ä¸ªæ²¡ç¨ï¼åªæ¯ä¸å®è¦æä¸å
  metadata_new <- rbind(metadata_new, a)
  
  FEAST_output1 <-
    FEAST(
      C = t(as.matrix(otu_animal)),
      metadata = metadata_new,
      different_sources_flag = 1,
      dir_path = getwd(),
      outfile = "demo"
    )
  FEAST_output <- rbind(FEAST_output, FEAST_output1)
}

FEAST_output$sample <- str_sub(rownames(FEAST_output), , -9)



### niche breadth and niche overlap
library(spaa)
t_otu <- t(otu_bac)
niche_breadth <- niche.width(t_otu, method = "levins")

t_otu_binary<-t_otu
t_otu_binary[t_otu_binary!=0]<-1
B_value_niche<-t_otu_binary
for (i in 1:dim(t_otu)[1]) { ###
  B_value_niche[i,]<- as.vector(as.matrix(B_value_niche[i,]))*as.vector(as.matrix(niche.width(t_otu,method="levins")))
}
rowMeans(B_value_niche,na.rm = T) ##Community level niche breadth (Bcom)


####Community level niche overlap (Ocom)
niche.overlap.comm <- function (otu_table, method = "morisita") {
  t_otu <- t(otu_table)
  a <- as.matrix(niche.overlap(t_otu, method = method))
  overlap_list <- list()
  
  for (i in 1:dim(t_otu)[1]) {
    p <- combn(colnames(t_otu)[t_otu[i, ] != 0], 2)
    overlap <- vector()
    for (j in 1:dim(p)[2]) {
      overlap[j] <- a[p[1, j], p[2, j]]
    }
    overlap_list[[i]] <- overlap
    print(i)
  }
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
  
