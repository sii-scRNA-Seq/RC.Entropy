
## does not work for more than sqrt(.Machine$integer.max) = 46340 points
## does not consider distances between 2 points, as long as their distance is less than dist

#' @importFrom magrittr %>%

count_coocurrences <- function(z, dist){

  data <- z

  ct <- as.vector(data$marks$marks.cat) # vector of the categories for all the points in partitions_sel
  I = length(table(ct)) # total number of categories observed in partitions_sel
  coocs_theo <- (I*I +I)/2 # maximum number of different possible co-occurences that could be observed for m=2

  # 1. calculate pseudo-adjacency matrix
  eucl <- spatstat.geom::pairdist(data)  # matrix of pair-wise euclidian distances between cells
  inradius <- (eucl <= dist) # matrix of pairs of cells at a distance less than threshold = pseudo-adjacency matrix (diagonal is 1)

  # 2. identify which pairs are in the same partition, and in which partition
  # creates n-by-n matrix where intersections are either 0 (not in the same partition) or the partition number e.g 1, 2...(ID)
  in_g <- data$marks$marks.partition
  sameg <- outer(in_g, in_g, "==") # matrix of pairs of cells that fall in the same partition

  # 3. Hadamard product of adjacency matrix and partition matrix
  # creates n-by-n matrix where intersections are either 0 (not in the same partition) or the partition number e.g 1, 2...(ID)
  ZLdg = inradius * sameg #ignore distances, gives diagonal of 1's
  #plot(as.im(ZLdg), clipwin=owin(xrange= c(0,100), yrange=c(0,100)))
  is_cooc <- ZLdg
  is_cooc[lower.tri(is_cooc)] <- 0 # matrix is symmetrical so co-occurences are double counted as ind1-ind2 and ind2-ind1, remove lower triangle
  is_cooc[diag(is_cooc)] <- 0 # remove diagonal to avoid counting self-links ind1-ind1
  #plot(as.im(is_cooc), clipwin=owin(xrange= c(0,100), yrange=c(0,100)))
  numing <- as.numeric(in_g)
  whichg <- matrix(numing, nrow=length(numing), ncol=length(numing), byrow=TRUE)
  cooc_byg <- is_cooc * whichg # puts the partition number (ID) where previously there was 1 indicating the 2 cells are in the same partition g
  #plot(as.im(cooc_byg), clipwin=owin(xrange= c(0,100), yrange=c(0,100)))

  # 4. find the points in the partitions of interest and their marks
  coocs <- which(cooc_byg != 0, arr.ind = T) #gives index position as row / col of all non-zero values in is_cooc
  ct1 <- ct[coocs[,1]] # position in row of all-vs-all matrix
  ct2 <- ct[coocs[,2]] # position in column of all-vs-all matrix
  ord_ct1 <- pmin(ct1, ct2)
  ord_ct2 <- pmax(ct1, ct2)
  ct_pairs <- paste(ord_ct1, ord_ct2, sep=" <-> ") # this is to order so that cat1-cat2 and cat2-cat1 both are counted as cat1-cat2
  #table(ct_pairs) #check

  # 5. make a vector of partition ids for each non zero values in coocs_byg
  cp <- as.vector(data$marks$marks.partition)
  g <- cp[coocs[,1]] # position in row of all-vs-all matrix # (sum(g == cp[coocs[,2]]) == length(ct_pairs)) # Returns TRUE, no need to order
  ct_pairs_byg <- data.frame(ct_pairs, ord_ct1, ord_ct2, g)
  #print(ct_pairs_byg)

  count_p_g <- ct_pairs_byg %>%
    dplyr::group_by(ct_pairs,ord_ct1, ord_ct2,g) %>% dplyr::summarize(n = n())

  return(count_p_g)
}

prep_b4_calc_H <- function(df, rep_sim, min_coocs){

  df["n"] <- df[rep_sim]
  df[rep_sim] <- NULL

  #remove pairs that have no occurences across all g's
  #start by identifying pairs with no occurences for safekeeping
  df_miss <- df %>% dplyr::filter(sum_n == 0)
  nb_pairs_miss <- length(unique(df_miss$ct_pairs))
  df <- df %>% dplyr::filter(sum_n >0)
  nb_pairs <- length(unique(df$ct_pairs))

  # replace NA's by 0 to prepare to add 1
  df$n[is.na(df$n)] <- 0

  #identify edge case of pairs not occuring in at least one partition
  undist <- df %>% dplyr::filter( n == 0)
  uneven_pairs_ls <- unique(undist$ct_pairs)

  # add 1 to occurence to each partition where there the pair does not occur in at least one partition
  df$n_mod <- df$n +  (df$ct_pairs %in% uneven_pairs_ls)
  rm(undist, uneven_pairs_ls)

  # recalculate total occurences per pair
  ncoocs_by_ctpairs <-  df[,c("ct_pairs", "n_mod")] %>% # Count the total number of cells from each co-occurence (across all partitions)
    dplyr::group_by(ct_pairs) %>% dplyr::summarise(sum_n_mod= sum(n_mod))
  df <- dplyr::inner_join(df, ncoocs_by_ctpairs, by="ct_pairs")

  # remove pairs where total occurences less than set minimum
  df <-  df %>% dplyr::filter( sum_n_mod > min_coocs)

  df["n"] <- df["n_mod"]
  df["sum_n"] <- df["sum_n_mod"]
  df["n_mod"] <- NULL
  df["sum_n_mod"] <- NULL

  return(df)
}

counts_to_H <- function(df, T_all){

  T_star <- T_all[[2]]
  T_tot <- T_all[[3]]

  df$p = df$n / df$sum_n
  df$plogtg_p = df$p * log(df$Tg / df$p)

  dec_H_r <- df %>%
    dplyr::group_by(ct_pairs) %>% dplyr::summarise(Hr = sum(plogtg_p))

  dec_H_r$Hr_rel <- (dec_H_r$Hr - log(T_star)) / (log(T_tot) - log(T_star))

  return(list(df, dec_H_r))
}


