# QM helper functions

prep.packages = function(...)
{
  package.list = c(...)
  loaded = package.list %in% .packages()
  if ( all(loaded) ) return(invisible())
  
  package.list = package.list[!loaded]
  installed = package.list %in% .packages(TRUE)
  if ( !all(installed) ) install.packages(package.list[!installed])
  for ( p in package.list ) library(p,character.only=TRUE)
}

prep.packages('vegan','plotrix','MuMIn','compute.es','pwr','visreg','mvabund','ggplot2','DHARMa','r2glmm', 'dplyr')

# Define some defaults
d_grid_size = 120
d_quad_size = 10
d_quad_dist = 10
d_nb_dist = 5
d_should_plot = T
d_anim_interval = 0.02
d_reef_percent = 0.3
d_knoll_count = 50
d_knoll_size = 10

# Get a contingency table from a vector of samples
# This function also includes 0 counts in the table
# For instance, if we have quadrats with only 3 and 5 individuals
# The vector returned by this function will also contain 0 for 4 individuals
get_samp_table = function(samps,sum_below=0,min_samp=0,max_samp=max(samps))
{
  min_samp = min(min_samp,min(samps))
  max_samp = max(max_samp,max(samps))
  full_tbl = rep(0,max_samp-min_samp+1) # Create a 0-filled vector of the appropriate size
  names(full_tbl) = min_samp:max_samp # name the vector elements
  samp_tbl = table(samps) # calcuate the contingency table from the samples
  full_tbl[names(samp_tbl)] = samp_tbl # fill the 0 vector with higher than 0 counts
  low_vals = which(full_tbl < sum_below)
  # If the first bin (0) is below the threshold, we remove it
  # Otherwise all bins will be aggregated together
  low_vals = low_vals[low_vals>1]
  if ( length(low_vals) > 0 ) {
    first = min(low_vals)
    last = length(full_tbl)
    cur_sum = sum(full_tbl[first:last])
    if ( cur_sum < sum_below ) {
      first = max(first-1,1)
      cur_sum = sum(full_tbl[first:last])
    }
    full_tbl[first] = cur_sum
    full_tbl = full_tbl[1:first]
  }
  return(full_tbl) # return the full vector
}

# open a device window if needed
open_window = function()
{
  # setup the windows device, as the RStudio device doesn't work very well for animations
  os.dev = switch (.Platform$OS.type,
                   windows = list(fun=windows,devname="windows"),
                   unix = list(fun=X11,devname=c("X11","X11cairo")),
                   list(fun=quartz,devname=c("X11","quartz")))
  # first check if we already have an active window
  if ( names(dev.cur()) %in% os.dev$devname ) return(invisible())
  # The device is not active, look for other devices
  devs = dev.list() # get the device list, as maybe we have an inactive windows device open somewhere
  if ( is.null(dev.list) ) os.dev$fun()
  else {
    open.dev = intersect(os.dev$devname,names(devs))
    if ( length(open.dev) > 0 ) dev.set(devs[open.dev[1]]) # if we found a device, use it
    else os.dev$fun() # open a new device
  }
}

# Plot a square population matrix in black and white
plot_pop = function(pop,ind_col=1,add=FALSE,interval=d_anim_interval)
{
  # Define the plot axis using the grid size
  dim_range=1:nrow(pop)
  # open the device window if needed
  open_window()
  # Prevent flickers during the animation (this holds the drawing of graphics until everything is ready)
  dev.hold()

  # Draw the population matrix
  image(dim_range,dim_range,pop,col=c(0,ind_col),xlab="",ylab="",add=add)
  # Draw a nice box around the image
  if ( ! add ) box()
  # Draw the resulting image + box
  dev.flush()
  # wait a bit so the animation doesn't run too fast
  if ( interval > 0 ) Sys.sleep(interval)
}

# Generate a random population of a specified size in a square matrix
# This function creates a square matrix filled with zeros and ones
# The zeros are empty spaces and ones are occupied
# By default, the function will also plot the matrix as an image
get_random_pop = function (pop_size, grid_size=d_grid_size, should_plot=d_should_plot)
{
  # To avoid filling the same space twice, the function generates the population
  # in a single dimension vector (from 1 to grid_size^2).
  # Once the population data is generated, the final square matrix is filled.
    
  # Get the total size of the environment
  env_size = grid_size^2;
  
  # Create an empty population by using the function "rep" (short for "repeat")
  # This function will create a vector of the specified size
  # filled with a single value (0 in our case).
  pop_data = rep(0,env_size);
    
  # To create our population, we randomly choose the locations to occupy.
  # The function "sample" will select numbers from one to env_size without repetition.
  # The number of samples is specified by pop_size.
  occupied = sample(env_size, pop_size);
  
  # Now we have a list of occupied locations within the environment,
  # so we update our population data with it.
  pop_data[occupied] = 1;
  
  # Copy the single dimension population data into a 2-D matrix
  pop = matrix(pop_data,grid_size);
  
  # Plot the population if needed
  if ( should_plot ) plot_pop(pop)
  # Return the population matrix
  return(pop);
}

# Get the expected counts in a randomly distributed population
get_expected_counts = function(observed_table,observed_mean) {
  counts = as.numeric(names(observed_table))
  max_count = max(counts) # Find the highest count
  exp_probs = dpois(0:(max_count-1),observed_mean) # calculate the poisson probability of getting each count except the highest
  exp_probs[max_count+1] = 1-sum(exp_probs) # calculate the probability of the highest count and above
  names(exp_probs) = 0:max_count
  exp_probs * sum(observed_table) # multiply the probabilities by the sample size to get the expected counts
}

# A very simple function to get the settlement function
get_settle_func = function(max_dist,min_chance,max_chance)
{
  # Get the slope
  slope = (max_chance - min_chance) / max_dist;
  # If the distance is higher than max distance, the chance is max
  # Else, the chance is somewhere on the line defined by min_chance + slope * distance
  settle_func = function (dist) ifelse(dist >= max_dist, max_chance, min_chance + slope * dist);
  # Determine whether settlement should occur by comparing the settle function result with a random number
  # between 0 and 1
  should_settle = function (dist) ifelse(settle_func(dist) > runif(1),1,0);
  return (should_settle);
}

# Get the next generation of settlers, 
# where their settling chance depends on the settlement function
get_next_generation = function(current_pop, settler_count, max_dist, chance_at_zero_dist, chance_at_max_dist, death_chance, settle_func=get_settle_func(max_dist,chance_at_zero_dist,chance_at_max_dist), should_plot=d_should_plot)
{
  # Get the grid size:
  grid_size = nrow(current_pop)
  # Get the new settlers
  new_settlers = get_random_pop(settler_count,grid_size,F);
  
  # Find all settler locations
  settler_pos = which(new_settlers == 1, arr.ind=T);

  # Precalculate the neighborhood center
  orig_center = max_dist + 1;
  
  for ( i in 1:settler_count ) # Go over all settlers
  {
    # Get the current settler's position
    row = settler_pos[i,'row'];
    col = settler_pos[i,'col'];
    
    # Calculate the neighborhood area.
    # This will define a square area in which to look for neighbors.
    min_row = row - max_dist;
    max_row = row + max_dist;
    min_col = col - max_dist;
    max_col = col + max_dist;
    
    center_row = orig_center;
    center_col = orig_center;
    
    # If the neighborhood exceeds the map bounds, reset it and update the center if needed
    if ( min_row < 1 ) { min_row = 1 ; center_row = row }
    if ( min_col < 1 ) { min_col = 1 ; center_col = col }
    
    if ( max_row > grid_size ) max_row = grid_size;
    if ( max_col > grid_size ) max_col = grid_size;

    # Get the neighborhood matrix
    nbhood = current_pop[min_row:max_row,min_col:max_col];
    
    # Clear the target if it already exists, to avoid 0 distance
    nbhood[center_row,center_col] = 0;
    
    # Get the locations of all neighbors
    curr_nbs = which(nbhood == 1,arr.ind = T);
    
    # If no neighbors were found, use the max distance as nearest neighbor
    if ( nrow(curr_nbs) == 0 )
    {
      new_settlers[row,col] = settle_func(max_dist);
      next;
    }
    
    # Calculate the distance of each neighbor from the center
    row_dists = curr_nbs[,'row'] - center_row;
    col_dists = curr_nbs[,'col'] - center_col;
    
    # Find the nearest neighbor - sqrt is not very efficient, so we use it only once
    nearest_neighbor = sqrt(min(row_dists^2 + col_dists^2));
    
    # Get the settlement chance near that neighbor
    new_settlers[row,col] = settle_func(nearest_neighbor);
  }
  
  # Add the new settlers to the current population
  current_pop = current_pop + new_settlers;
  
  # Remove those who settled in an occupied spot
  current_pop[which(current_pop > 1)] = 1;
  
  # Before the next turn, we kill some of the settlers, it is sad.
  current_pop = kill_settlers(current_pop,death_chance);
  
  # Plot the population if needed
  if ( should_plot ) plot_pop(current_pop)
  
  return(current_pop);
}

# This function randomly samples 'sample_chance' percent of individuals from the population
sample_pop = function(pop,sample_chance)
{
  # Get all individual positions
  pop_positions = which(pop == 1);
  # Get a random number between 0 to 1 for each individual
  all_chances = runif(length(pop_positions));
  # Return all individuals for which the random number is smaller than the sample chance
  return(pop_positions[sample_chance >= all_chances]);
}

# Reproduce some individuals
# Each individual has a chance to create another individual
# If it does, the new individual is placed in a free spot, if available
reproduce = function (pop, birth_chance)
{
  # Get the number of reproducing settlers
  reproducing_count = length(sample_pop(pop,birth_chance));
  # Add settlers as the number of reproducing individuals
  return(settle(pop,reproducing_count))
}

# Settle individuals in empty positions, if available
settle = function(pop,settler_count)
{
  # Get all available empty positions
  empty_positions = which(pop == 0);
  empty_count = length(empty_positions);
  # If no position is available, no one settles
  if ( empty_count == 0 ) return(pop);
  # Do not allow more settlers than empty positions
  if ( settler_count > empty_count ) settler_count = empty_count;
  # Get the positions for the new settlers
  settle_pos = sample(empty_positions,settler_count);
  # Populate the new positions
  pop[settle_pos] = 1;
  return(pop);
}

# This function will kill a fraction of the population
kill_settlers = function(pop,death_chance)
{
  # Get all positions
  kill_positions = sample_pop(pop,death_chance);
  pop[kill_positions] = 0;
  return(pop);
}

# Plot a square environment matrix
plot_env = function(env,xlab="",ylab="")
{
  # Define the plot axis using the grid size
  dim_range=1:nrow(env)
  # open the device window if needed
  open_window()
  # Prevent flickers during the animation (this holds the drawing of graphics until everything is ready)
  dev.hold()
  nvals = length(unique(c(env)))
  col.levels = seq(1,0,length.out = nvals)
  cols = grey(col.levels)
  # Draw the population matrix
  image(dim_range,dim_range,env,col=cols,xlab=xlab,ylab=ylab)
  # Draw a nice box around the image
  box()
  # Draw the resulting image + box
  dev.flush()
  # wait a bit so the animation doesn't run too fast
  invisible()
}

# Get an environment with open water, knolls and one big reef
get_env = function(knoll_count=d_knoll_count,
                   knoll_size=d_knoll_size,
                   reef_percent=d_reef_percent,
                   grid_size=d_grid_size,
                   should_plot=d_should_plot) {
  
  env = matrix(0,ncol = grid_size, nrow=grid_size)
  
  cats = rbind(Knoll = c(knoll_count,knoll_size,knoll_size),
    Reef = c(1,round(grid_size*reef_percent),grid_size))
  
  for ( i in 1:nrow(cats) ) {
    quad_cnt = cats[i,1]
    quad_x = cats[i,2]
    quad_y = cats[i,3]
    quads = get_random_quad_positions(env,quad_cnt,quad_x = quad_x,quad_y = quad_y)
    for ( q in 1:nrow(quads) ) {
      x_start = quads[q,'min_x']
      x_end = quads[q,'max_x']
      x_overflow = quads[q,'x_overflow']
      
      y_start = quads[q,'min_y']
      y_end = quads[q,'max_y']
      y_overflow = quads[q,'y_overflow']
      
      corner = x_overflow > 0 & y_overflow > 0
      
      env[x_start:x_end,y_start:y_end] = i
      if ( x_overflow > 0 ) env[1:x_overflow,y_start:y_end] = i
      if ( y_overflow > 0 ) env[x_start:x_end,1:y_overflow] = i
      if ( corner ) env[1:x_overflow,1:y_overflow] = i
    }
  }
  if ( should_plot ) plot_env(env,ylab="Depth")
  attr(env,'habitats')=c('Open',rownames(cats))
  env
}

# plot the results of sample_env_comm
plot_sampled_env = function(env_list) {
  with(env_list,{
    plot_env(env,ylab="Depth")
    plot_community(comm,add=TRUE)
    plot_quads(quads,col='blue')
  })
}

# resample an existing env_comm using different quadrats
resample_env_comm = function(env_comm,
                             quad_size=d_quad_size,
                             dist_x=d_quad_dist,
                             dist_y=d_quad_dist,
                             should_plot=d_should_plot)
{
  env_comm$quads = get_sys_quad_positions(env_comm$env,quad_size,dist_x,dist_y)
  hab.names = attr(env_comm$env,'habitats')
  habitats = hab.names[apply(env_comm$quads,1,function(q){
    fx = q['min_x']
    tx = q['max_x']
    fy = q['min_y']
    ty = q['max_y']
    which.max(get_samp_table(env_comm$env[fx:tx,fy:ty],max_samp=max(env_comm$env)))
  })]
  env_comm$sites = data.frame(Habitat=habitats,
                     Depth=rowMeans(env_comm$quads[,c('min_y','max_y')]))
  
  env_comm$samps = sample_quads(env_comm$comm,env_comm$quads,FALSE)
  if ( should_plot ) plot_sampled_env(env_comm)
  env_comm
}
# Generate an environment, populate it with a community and sample that community
# Some species will prefer different habitats or depths
sample_env_comm = function(sps_cnt,coverage,
                            knoll_percent=0.5,
                            reef_percent=d_reef_percent,
                            quad_size=d_quad_size,
                            dist_x=d_quad_dist,
                            dist_y=d_quad_dist,
                            grid_size=d_grid_size,
                            should_plot=d_should_plot) {
  set.seed(1)
  # calculate the knoll count and size based on quad and grid size
  # knoll size should be about 25% larger than quadrate size
  knoll_size = quad_size * 1.25
  max_knolls = grid_size^2/knoll_size^2
  knoll_count = round(max_knolls*knoll_percent)
  env = get_env(knoll_count,knoll_size,reef_percent,grid_size,FALSE) # will be plotted later
  dim.vals = seq(-1,1,length.out = grid_size)
  pos = expand.grid(x=1:grid_size,y=1:grid_size)
  hab.names = attr(env,'habitats')
  get.cat = function(x,y) { hab.names[1+env[x,y]] }
  pos$cat = factor(mapply(get.cat,pos$x,pos$y))
  cat.vars = contr.sum(levels(pos$cat))
  colnames(cat.vars) = apply(cat.vars,2,function(x){
    n = names(x)
    paste(n[x < 0],n[x > 0],sep='-')
  })
  X = cbind(Depth=dim.vals[pos$y],cat.vars[pos$cat,])
  rownames(X) = NULL
  sp.names = paste0("sp",1:sps_cnt)
  coef.vals = seq(-1,1,by=0.01)
  coef.probs = exp(-0.05 * 1:100)
  coef.probs = c(coef.probs,0,rev(coef.probs))
  coefs = t(replicate(ncol(X),sample(coef.vals,sps_cnt,replace = T,prob = coef.probs)))
  dimnames(coefs) = list(colnames(X),sp.names)
  Y = X %*% coefs
  probs = 1/(1+exp(-Y))
  sps.range = 0:sps_cnt
  pos$sps = apply(probs,1,function(p){
    p = c(1-coverage,coverage*p/sum(p))
    sample(sps.range,1,replace=T,prob = p)
  })
  pos = pos[pos$sps > 0,]
  comm = array(0,c(dim(env),sps_cnt),list(NULL,NULL,paste0("sp",1:sps_cnt)))
  for ( i in 1:nrow(pos) ) {
    row = pos[i,]
    comm[row$x,row$y,row$sps] = 1
  }
  ret.list = list(comm=comm,env=env,
                  model=list(X=X,coefs=coefs,Y=Y,P=probs))
  ret.list = resample_env_comm(ret.list,quad_size,dist_x,dist_y,should_plot)
  set.seed(NULL)
  ret.list
}
# draw a scree plot and some stress plots
stress_test = function(samps,min_k,max_k,plot_stress=TRUE) {
  stress_values = sapply(min_k:max_k,function(k) {
    mds_res = metaMDS(samps,distance='bray',k=k)
    if ( plot_stress ) stressplot(mds_res,main=paste('Stress =',round(mds_res$stress,3),'k =',k))
    mds_res$stress
  })
  plot(min_k:max_k,stress_values,type='b',xlab='k',ylab='stress',main='Scree plot')
}

# settle a population over X generations, using interaction based settlement parameters
# This function can be used to create clumped, even, or random populations
settle_interaction_pop = function(generations, settler_count, max_dist, chance_at_zero_dist, chance_at_max_dist, death_chance, max_pop_size=Inf, grid_size=d_grid_size, should_plot=d_should_plot)
{
  # Get the settlement function which will determine the settlement chance
  settle_func = get_settle_func(max_dist,chance_at_zero_dist,chance_at_max_dist);
  #Initialize an empty field
  current_pop = matrix(0,grid_size,grid_size);
  for ( gen in 1:generations )
  {
    # Each turn we get a new generation of settlers
    current_pop = get_next_generation(current_pop,settler_count,max_dist,chance_at_zero_dist,chance_at_max_dist,death_chance,settle_func,should_plot);
    # Report our progress
    print(paste('Generation:',gen,'Size:',sum(current_pop)));
    if ( sum(current_pop) >= max_pop_size ) return(current_pop)
  }
  return(current_pop);
}

# Default values to generate a clumped population
settle_clumped_pop = function(generations,grid_size=d_grid_size,should_plot=d_should_plot)
{
  settle_interaction_pop(generations,
    settler_count=120,
    max_dist=4,
    chance_at_zero_dist=1.75,
    chance_at_max_dist=0.01,
    death_chance=0.1,
    grid_size = grid_size,
    should_plot = should_plot)
}

# Default values to generate an even population
settle_even_pop = function(generations,grid_size=d_grid_size,should_plot=d_should_plot)
{
  settle_interaction_pop(generations,
                         settler_count=100,
                         max_dist=6,
                         chance_at_zero_dist=-120,
                         chance_at_max_dist=0.999,
                         death_chance=0.01,
                         grid_size = grid_size,
                         should_plot = should_plot)
}

# Randomly settle either an even or clumped population using the above functions
settle_any_pop = function(generations,grid_size=d_grid_size,should_plot=d_should_plot)
{
  if(runif(1) < 0.5 ) settle_clumped_pop(generations,grid_size,should_plot)
  else settle_even_pop(generations,grid_size,should_plot)
}

# grow a population by adding settlers each turn
settle_random_pop = function(generations,min_settlers,max_settlers,death_chance,grid_size=d_grid_size,should_plot=d_should_plot)
{
  #Initialize an empty field
  current_pop = matrix(0,grid_size,grid_size);
  for ( gen in 1:generations )
  {
    settlers = sample(min_settlers:max_settlers,1)
    current_pop = settle(current_pop,settlers)
    current_pop = kill_settlers(current_pop,death_chance)
    if ( should_plot ) plot_pop(current_pop)
    print(paste('Generation:',gen,'Size:',sum(current_pop)));
  }
  return(current_pop)
}

get_random_quad_positions = function(grid,samp_count,quad_size=d_quad_size,quad_x=quad_size,quad_y=quad_size) {
  # Find the sampling limits (avoid too small quadrates)
  grid_size = nrow(grid)
  x_start = sample(grid_size,samp_count)
  x_end = x_start + quad_x - 1
  
  y_start = sample(grid_size,samp_count)
  y_end = y_start + quad_y - 1
  
  cbind(min_x=x_start,
        min_y=y_start,
        max_x=pmin(x_end,grid_size),
        max_y=pmin(y_end,grid_size),
        x_overflow=pmax(x_end-grid_size,0),
        y_overflow=pmax(y_end-grid_size,0))
}

# Calculate the quadrat positions for systematic sampling
get_sys_quad_positions = function(grid,quad_size=d_quad_size, dist_x=d_quad_dist, dist_y=d_quad_dist) {
  grid_size=nrow(grid)
  max_x = seq(quad_size,grid_size,dist_x + quad_size)
  max_y = seq(quad_size,grid_size,dist_y + quad_size)
  start_diff = quad_size - 1
  ret.mat = NULL
  for ( x in max_x ) {
    sx = x - start_diff
    for ( y in max_y ) {
      sy = y - start_diff
      ret.mat = rbind(ret.mat,c(min_x=sx,min_y=sy,max_x=x,max_y=y))
    }
  }
  cbind(ret.mat,x_overflow=0,y_overflow=0)
}

plot_quads = function(quads,col='red') {
  # we plot the rect around the actual quadrat - everything inside the border is counted
  sx = quads[,'min_x']-1
  sy = quads[,'min_y']-1
  ex = quads[,'max_x']+1
  ey = quads[,'max_y']+1
  ox = quads[,'x_overflow']+1
  oy = quads[,'y_overflow']+1
  
  # first draw the part which doesn't overflow
  rect(sx,sy,ex,ey,border=col,lwd=2)
  # we added 1 to the offsets, so we draw them only if they are 2 and up.
  # When the last line in the quadrat is the last line on the grid,
  # the grid border will be our limit
  x_of = ox > 1
  if ( any(x_of) ) rect(0,sy[x_of],ox[x_of],ey[x_of],border=col,lwd=2)
  
  y_of = oy > 1
  if ( any(y_of) ) rect(sx[y_of],0,ex[y_of],oy[y_of],border=col,lwd=2)
  
  xy_of = x_of & y_of
  if ( any(xy_of) ) rect(0,0,ox[xy_of],oy[xy_of],border=col,lwd=2)
}

# This function takes a population matrix and a set of quadrat coordinates and samples them
sample_quads = function(pop,quads,should_plot=d_should_plot) {
  # If we have more than 2 dimensions, we have a community
  pop_dim = dim(pop)
  if ( length(pop_dim) > 2 && pop_dim[3] == 1 ) pop = pop[,,1]
  is_single_pop = is.matrix(pop)

  # Draw the quadrats
  if ( should_plot )
  {
    if ( is_single_pop ) plot_pop(pop,interval=0)
    else plot_community(pop,interval=0)
    plot_quads(quads)
  }
  # This function samples a single quadrat
  single_sample = function(coords)
  {
    x_start = coords['min_x']
    x_end = coords['max_x']
    x_overflow = coords['x_overflow']
    
    y_start = coords['min_y']
    y_end = coords['max_y']
    y_overflow = coords['y_overflow']
    
    corner = x_overflow > 0 & y_overflow > 0
    
    sum_quad = function(cur_pop) {
      samp_sum = sum(cur_pop[x_start:x_end,y_start:y_end])
      if ( x_overflow > 0 ) samp_sum = samp_sum + sum(cur_pop[1:x_overflow,y_start:y_end])
      if ( y_overflow > 0 ) samp_sum = samp_sum + sum(cur_pop[x_start:x_end,1:y_overflow])
      if ( corner ) samp_sum = samp_sum + sum(cur_pop[1:x_overflow,1:y_overflow])
      samp_sum
    }
    if ( is_single_pop ) sum_quad(pop)
    else apply(pop,3,sum_quad)
  }
  result = apply(quads,1,single_sample)
  if ( is_single_pop ) result
  else t(result)
}

# Randomly sample a quadrate from the population or community
random_sampling = function(pop, samp_count, quad_size=d_quad_size,should_plot=d_should_plot)
{
  quads = get_random_quad_positions(pop,samp_count,quad_size)
  sample_quads(pop,quads,should_plot)
}

# Sample quadrats in equal distances from each other
systematic_sampling = function(pop, quad_size=d_quad_size, dist_x=d_quad_dist, dist_y=d_quad_dist, should_plot=d_should_plot)
{
  quads = get_sys_quad_positions(pop,quad_size,dist_x,dist_y)
  sample_quads(pop,quads,should_plot)
}

# Get the estimated population size from a sample
get_ests = function(samps,grid_size,quad_size)
{
  samps * grid_size^2/quad_size^2
}

# Simulate a random population using birth and death chances
# Sample each generation and return the vector of samples over all generations
get_pop_sample_sequence = function(generations,init_pop,birth_chance,death_chance,grid_size=d_grid_size,quad_size=d_quad_size,should_plot=d_should_plot)
{
  if ( is.matrix(init_pop) )
  {
    current_pop = init_pop
  }
  else
  {
    # Get the initial population
    current_pop = get_random_pop(init_pop,grid_size,should_plot)
  }
  # Make room for the samples
  samples = array(0,generations)
  for ( gen in 1:generations ) # Run for the requested number of generations
  {
    # Add more individuals with a replication chance for each existing individual
    # i.e. if we have 10 individuals and a birth chance of 0.5, then 5 individuals
    # will be added (if there is room for them)
    current_pop = reproduce(current_pop,birth_chance)
    # Kill some individuals, to keep things moving
    current_pop = kill_settlers(current_pop,death_chance)
        
    # Sample the population and store the result
    samples[gen] = random_sampling(current_pop,1,quad_size,should_plot)
    
    # Report our progress
    print(paste('Generation:',gen,'Size:',sum(current_pop), 'Sample:',samples[gen]))
    # If everyone is dead, we don't need to continue any longer
    # so we break from the loop
    if ( sum(current_pop) == 0 ) break
  }
  return(list(pop=current_pop,samples=samples))
}

# Functions for various patchiness indices
# lloyds index
lloyds = function(vals)
{
  m = mean(vals)
  v = var(vals)
  return(1 + ((v - m)/m^2))
}

# Coefficient of variation
cv = function(vals)
{
  sd(vals)/mean(vals)
}

# Index of dispersion (VMR - Variance to Mean Ratio)
vmr = function(vals)
{
  var(vals)/mean(vals)
}

get_state_data = function() {
  states = data.frame(state.x77)
  set.seed(1)
  states$Unemployment = states$Illiteracy * rnorm(nrow(states),mean=1,sd=0.5)
  set.seed(NULL)
  states
}

# Calculate the variance inflation factor
vif = function(lm.fit) {
  exp.vars = model.matrix(lm.fit)
  exp.vars = data.frame(exp.vars[,-1]) # the first column is the intercept
  ret.vals = NULL
  for ( v in colnames(exp.vars) ) {
    cur.form = formula(paste(v,'~ .'))
    cur.fit = summary(lm(cur.form,exp.vars))
    ret.vals = c(ret.vals,1/(1-cur.fit$r.squared))
  }
  names(ret.vals) = colnames(exp.vars)
  ret.vals
}

# Get samples from a discrete distribution with a given mean and variance
get_distributed_samps = function(samp_count,samp_mean,samp_var)
{
  if ( samp_var > samp_mean)
  {
    prob = samp_mean/samp_var
    size = samp_mean * prob / (1-prob)
    rnbinom(samp_count,
            size=size,
            prob=prob)
  }
  else if ( samp_var < samp_mean )
  {
    prob = 1-samp_var/samp_mean
    size = round(samp_mean/prob)
    rbinom(samp_count,size,prob)
  }
  else
  {
    rpois(samp_count,samp_mean)
  }
}

# Get a distribution index from a population with the variance as twice
get_distribution_index = function(last_mean,variance_ratio,index_function)
{
  samps = sapply(1:last_mean,function(x){get_distributed_samps(1000,x,x*variance_ratio)})
  apply(samps,2,index_function)
}

# A function to "mimic" the division of a random population into quadrats
get_random_quads = function(pop_size,quad_count)
{
  # Randomly assign a quadrat to each individual
  ind_quads = sample(quad_count,pop_size,replace=T)
  # Count the number of individuals in each quadrat
  return(tabulate(ind_quads,quad_count))
}

# Get the patchiness index of choice for several random populations
get_random_index = function(pop_size,quad_count,pop_count,index_function)
{
  replicate(pop_count,
  {
    quads = get_random_quads(pop_size,quad_count)
    index_function(quads)
  })
}

# Calculate an index of patchiness for a specific quadrate size
# This function will first divide the population into quadrates
# and count the number of individuals in each quadrate
# The mean and variance of those numbers will be used to
# find the patchiness level of the environment
get_pop_index = function(pop,quad_size,index_func)
{
  # First we need to divide the population matrix into quadrates
  # We find the number of quadrates per row and per column by dividing the
  # number of rows or columns by the quadrate size.
  # Please note, that this function supports non-square population matrices.
  quad_row_cnt = nrow(pop) %/% quad_size;
  quad_col_cnt = ncol(pop) %/% quad_size;
  
  # Now we use the "seq" function to generate the quadrate corner coordinates
  quad_rows = seq(from=1,by=quad_size,length.out=quad_row_cnt); # Get the quadrate row indices
  quad_cols = seq(from=1,by=quad_size,length.out=quad_col_cnt); # Get the quadrate column indices
  
  # We initialize an array Make room for the counts
  # The total number of quadrates is quad_row_cnt * quad_col_cnt
  individual_counts = array(NA,quad_row_cnt * quad_col_cnt);
  quad_count = 0; # Initialize a counter for the quadrates
  for ( first_row in quad_rows ) # Go over each quadrate row
  {
    # Calculate the last row of the current quadrate
    last_row = first_row + quad_size - 1;
    
    # Store the quadrate row range for use inside the column loop
    row_range = first_row:last_row;
    
    for ( first_col in quad_cols ) # For each row, go over each quadrate column
    {
      # Calculate the last column of the current quadrate
      last_col = first_col + quad_size - 1;
      
      # Get the quadrate column range
      col_range = first_col:last_col;
      
      # Get the quadrate matrix as a subset of the original population matrix
      curr_quad = pop[row_range,col_range];
      
      # Increment the quadrate count
      quad_count = quad_count + 1;
      
      # Store the number of individuals in the quadrate
      individual_counts[quad_count] = sum(curr_quad);
    }
  }
  
  # Return the relevant index:
  index_func(individual_counts);
}

get_pop_lloyds = function(pop,quad_size=d_quad_size) {
  get_pop_index(pop,quad_size,lloyds)  
}

get_pop_id = function(pop,quad_size=d_quad_size) {
  get_pop_index(pop,quad_size,vmr)  
}

get_pop_cv = function(pop,quad_size=d_quad_size) {
  get_pop_index(pop,quad_size,cv)  
}

# Get the average number of neighbors within a certain distance
get_neighbors = function(pop, dist=d_nb_dist)
{
  # Get the indices of all individuals in the population
  # By using arr.ind=T, we get both row and column indices.
  
  individuals = which(pop == 1, arr.ind=T);
  
  # Get the number of rows and columns in the matrix (needed to avoid the edge bias)
  rows = nrow(pop);
  cols = ncol(pop);
  
  # Prepare room for all neighbor counts
  individual_count = nrow(individuals);
  nbs = array(NA,individual_count);
  
  center = dist + 1; # Precalculate the neighborhood center coordinate for efficiency
  dist_sq = dist^2; # Precalculate the squared distance for efficiency
  
  for ( i in 1:individual_count ) # Go over all individuals in the population
  {
    # Get the current individual's position
    row = individuals[i,'row'];
    col = individuals[i,'col'];
    
    # Calculate the neighborhood area.
    # This will define a square area in which to look for neighbors.
    # We skip any indiviudals near the edges, to avoid an edge bias
    min_row = row - dist;
    if ( min_row < 1 ) next;
    max_row = row + dist;
    if ( max_row > rows ) next;
    min_col = col - dist;
    if ( min_col < 1 ) next;
    max_col = col + dist;
    if ( max_col > cols ) next;
    
    # Get the neighborhood matrix
    nbhood = pop[min_row:max_row,min_col:max_col];
    
    # Remove our target from the neighborhood, as we need to count only the neighbors
    nbhood[center,center] = 0;
    
    # Get the locations of all neighbors
    curr_nbs = which(nbhood == 1,arr.ind = T);
    
    # If no neighbors were found, set the count to 0 and move on
    if ( nrow(curr_nbs) == 0 )
    {
      nbs[i] = 0;
      next;
    }
    
    # Calculate the distance of each neighbor from the center
    row_dists = curr_nbs[,'row'] - center;
    col_dists = curr_nbs[,'col'] - center;
    nb_dists_sq = row_dists^2 + col_dists^2;
    
    # Include all neighbors within the requested distance
    nbs[i] = sum(dist_sq >= nb_dists_sq)
  }
  
  # As we skipped individuals near the edges, we need to remove NAs from the mean calculation
  return(mean(nbs,na.rm=T));
}

# Ex 4 - Community indices
# Plot a 3d array of community data
plot_community = function(comm,add=FALSE,interval=d_anim_interval)
{
  sps_cnt = dim(comm)[3]
  cols = rainbow(sps_cnt)
  plot_pop(comm[,,1],ind_col=cols[1],add=add,interval=interval)
  if ( sps_cnt > 1 )
  {
    sapply(2:sps_cnt,function(sps){
      plot_pop(comm[,,sps],ind_col=cols[sps],add=TRUE,interval=interval)
    })
  }
  invisible()
}

# Get a 3d array of community data
get_community = function(sps,meanlog,sdlog,grid_size=d_grid_size,should_plot=d_should_plot)
{
  # Generate initial species abundances (with a log-normal distribution)
  n = ceiling(rlnorm(sps,meanlog,sdlog))
  # Get the total number of occupied positions
  total_inds = sum(n)
  # Calculate the total number of positions in the environment
  env_size = grid_size^2
  # Stop if there are too many individuals
  if ( env_size < total_inds )
  {
    ratio = env_size / total_inds
    n = floor(n * ratio)
    total_inds = sum(n)
    warning(paste("Reduced abundances by",round(ratio,3),"to fit grid size"))
  }
  # Get all occupied indices
  occupied = sample(env_size,total_inds)
  # Divide the occupied indices to species
  last = 0
  comm = sapply(cumsum(n),function(cur_ind) {
    # The community structure is a 3d array of single-species population matrices
    # Generate the current species matrix
    pop = matrix(0,grid_size,grid_size)
    first = last + 1 # The first occupied index for the current species
    last <<- cur_ind # The last occupied index for the current species
    pop[occupied[first:last]] = 1 # Occupy the indices
    return(pop)
  },simplify='array')
  if ( should_plot ) plot_community(comm)
  dimnames(comm) = list(NULL,NULL,paste0("sp",1:sps))
  return(comm)
}

# Plot assemblage and sample species abundance histograms
# for log based distributions (e.g. log-normal)
plot_species_counts = function(comm,sample_data)
{
  # Get the relative abundances from the community
  rel_abundance = apply(comm,3,sum)
  rel_abundance = rel_abundance/sum(rel_abundance)
  
  # First remove species with 0 counts
  sample_data = sample_data[,colSums(sample_data)>0]
  # Get the relative abundances of the sample data
  sample_rel_abundance = colSums(sample_data)/sum(colSums(sample_data));
  
  # Structure the X axis of the histograms
  # First find the log of the smallest relative abndance
  min_log_ra = floor(log(min(c(rel_abundance,sample_rel_abundance))));
  # Use the smallest relative abundance as the leftmost index
  break_cnt = 1 + min_log_ra * -1;
  # The histogram breaks will be on log scale - from min_log_ra to 0
  breaks = seq(min_log_ra,0,length.out=break_cnt);

  # Get the assemblage histogram data (frequencies and counts)
  ra_freq = hist(log(rel_abundance),plot=F,breaks=breaks);
  # Get the sample histogram data (frequencies and counts)
  sample_ra_freq = hist(log(sample_rel_abundance),plot=F,breaks=breaks);
  
  # Store the histogram results in a barplot ready format
  bar_data = matrix(NA,2,length(ra_freq$counts));
  bar_data[1,] = ra_freq$counts;
  bar_data[2,] = sample_ra_freq$counts;
  # Create the X labels
  x_labels = paste("e^",breaks[-1],sep="")
  # Draw the histograms of both assemblage and sample data
  barplot(xlab='Relative abundance (% of total species count - Log scale)',
          ylab='Species Count',
          bar_data,
          beside=T,
          names.arg = parse(text=x_labels),
          legend.text=c('Assemblage','Sample'));
}

# Plot assemblage and sample species abundance histograms
# for log based distributions (e.g. log-normal)
hist_community = function(comm)
{
  # Get the relative abundances from the community
  abundances = apply(comm,3,sum)
  
  # Get the assemblage histogram data (frequencies and counts)
  ra_freq = hist(log(abundances),plot=F)
  
  breaks = round(exp(ra_freq$breaks))
  x_labels = paste0(breaks[-length(breaks)],"-",breaks[-1])
  # Create the X labels
  # Draw the histograms of both assemblage and sample data
  barplot(xlab='Number of individuals',las=3,cex.names=0.8,
          ylab='Species Count',
          ra_freq$counts,
          names.arg = x_labels)
  abline(h=0)
}

# Draw a box plot and return the outlier rows
box_plot = function(data,...)
{
  # First draw the boxplot as usual
  box_data = boxplot(data,...);
  # Get the outlier coordinates
  # Because each outlier may appear more than once
  # we work only with unique outliers
  outlier_info = unique(data.frame(Group=box_data$group,Value=box_data$out));
  # Count how many unique outliers we have
  outlier_cnt = nrow(outlier_info);
  # If there are no outliers plotted, return
  if ( outlier_cnt == 0 ) return(box_data);
  # Get the column names of the outliers (their X axis)
  outlier_info$GroupName = box_data$names[outlier_info$Group];
  # Create a function to find the label of each outlier
  get_label = function(outlier)
  {
    colname = outlier['GroupName'];
    value = as.numeric(outlier['Value']);
    # find all rows in the specified column with the outlier value
    # Note that because we may compare very small values, we need to use
    # Fuzzy comparison (all.equal) instead of exact '=='
    matching_rows = which(sapply(data[,colname],FUN=function(x) {isTRUE(all.equal(x,value,1e-06))}));
    return(matching_rows);
  }
  # Apply the get_label function to all outliers
  outlier_info$Rows = apply(outlier_info,1,get_label);
  # Add the outlier data to the return structure
  box_data$out_info = outlier_info;
  return(box_data);
}

pcacircle = function (pca.res)
{
  # Draws a circle of equilibrium contribution on a PCA plot 
  # generated from a vegan analysis.
  # vegan uses special constants for its outputs, hence 
  # the 'const' value below.
  p = length(pca.res$CA$eig)
  const = attr(summary(pca.res),"const")
  radius = const * (2/p)^0.5
  symbols(0, 0, circles=radius, inches=FALSE, add=TRUE, fg=2)
  # Get the species scores from the result (using the same scaling as the plot)
  species = scores(pca.res,scaling=1,display='sp')
  # Calculate the arrow lengths
  spec.dist = sqrt(species[,1]^2+species[,2]^2)
  # Get only the significant species:
  sig.spec = spec.dist[spec.dist > radius]
  # Sort from highest to lowest and return the sorted values:
  return(sort(sig.spec,decreasing=TRUE))
}

eigval_test = function(pca)
{
  #extract the eigenvalues and apply Kaiser-Guttman criterion for axes selection
  #***********************************
  ev=pca$CA$eig
  
  # Broken stick model
  n <- length(ev)
  bsm <- data.frame(j=seq(1:n), p=0)
  bsm$p[1] <- 1/n
  for (i in 2:n) {
    bsm$p[i] = bsm$p[i-1] + (1/(n + 1 - i))
  }
  bsm$p <- 100*bsm$p/n
  
  ### Plot eigenvalues and % of variance for each axis
  par(mfrow=c(2,1))
  barplot(ev, main="Eigenvalues", col="bisque", las=2)
  abline(h=mean(ev), col="red")  # average eigenvalue
  legend("topright", "Average eigenvalue", lwd=1, col=2, bty="n")
  barplot(t(cbind(100*ev/sum(ev),bsm$p[n:1])), beside=TRUE, 
          main="% variance", col=c("bisque",2), las=2)
  legend("topright", c("Variance Explained", "Broken stick model"), 
         pch=15, col=c("bisque",2), bty="n")
  par(mfrow=c(1,1))
}
