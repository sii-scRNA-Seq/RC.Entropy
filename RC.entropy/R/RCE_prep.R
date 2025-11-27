# check correct formatting of input data and parameters
# data is a spatial point pattern (.ppp)
# partition_map is either a list (sites) or a pixel image (areas)

sanity_check_structs <- function(data, partition_map, partitions_sel, dist){

  spatstat.geom::is.ppp(data) # returns TRUE
  spatstat.geom::is.marked(data) # returns TRUE
  if(spatstat.geom::is.ppp(data) & !spatstat.geom::is.marked(data) & length(table(data$marks$cat))==1)
    stop("Data have only 1 category. Use Batty's entropy for this use case.")
  is.list(partition_map) # returns FALSE if "areas" and TRUE if "sites"
  spatstat.geom::is.im(partition_map) # returns TRUE if "areas" and FALSE if "sites"
  if(spatstat.geom::area.owin(data$window)==1)
    stop("Nemri's entropy cannot be computed on areas of total size = 1. You should change your measurement unit so that the window size is not 1.")
  #TO DO: check partitions_sel and dist in correct format
  # is.list(partition_sel)
  # for sites: dist == 0
  # for areas: dist >0
  # min_coocs >= 0 #add to function parameters
  }

# Collect information on partitions (T_g, T_star, T_tot)

partition_areas <- function(map, regions, map_type){
  # Get the area of each partition to table
  # if pixel image, we count the pixels in each partition
  if ( map_type == "areas") {
    T_g <- table(map$v) # Area of each partition
  }
  #if sites, we count the sites
  if ( map_type == "sites") {
    #print("map")
    T_g <- table(map) # "Area" of each partition
  }
  T_g_sub <- T_g[regions] # in partitions_sel
  T_star <- min(T_g_sub) # area of smallest partition
  T_tot <- sum(T_g_sub) # total area of partitions in partitions_sel

  return(list(T_g, T_star, T_tot))
}

# This function filters points that are in selected areas

select_points_in_partitions_areas <- function(z, morpho, partitions_sel){
  # Assign variables to function-to-be
  # z has just one column
  data <- z
  partition_map <- morpho
  # 1. Assign the cells to their partition based on the factor value for the pixel where they fall
  cutdata <- cut(spatstat.geom::unmark(data), z=partition_map )
  #table(cutdata$marks) #sanity check on the cell to partition attribution - before subsetting partition of interest
  partition <- cutdata$marks
  cutdata <- spatstat.geom::unmark(cutdata)
  cutdata$marks$partition <- partition
  cutdata$marks$cat <- data$marks$cat
  cutdata$marks$p_names <- data$marks$p_names
  cutdata <- subset(cutdata, partition %in% partitions_sel) # cells in other partitions and out of window are removed
  #table(cutdata$marks$marks.partition) #sanity check on the cell to partition attribution - after subsetting partition of interest

  return(cutdata)
}

# This function filters points that are in selected sites

select_points_in_partitions_sites <- function(z, partitions_sel){
  # z marks are 'p_names', 'cat' and 'partition
  data <- z
  cutdata <- subset(data, marks.partition %in% partitions_sel)
  return(cutdata)
}
