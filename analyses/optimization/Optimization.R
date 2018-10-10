############## Optimizing Window Settings #############

## Purpose: Given the number of frames and number of swath windows,
##          look at the coverage of library for each setting


############# Packages Required ######################

library(ggplot2)
library(stringr)
library(data.table)

############## Functions #############################

# Cleaning the data:
binColumn <- function(data, column_name, binwidth){
  # Purpose:
  #       Bin the data to shrink the data size
  # parameters:
  #       data          data.table  full maxquant output library
  #       column_name   str         columns to be binned
  cmi <- min(data[, get(column_name)])
  cma <- max(data[, get(column_name)])
  breaks <- seq(cmi, cma, by = binwidth)
  data[,bin:=findInterval(get(column_name), breaks, rightmost.closed = T)]
  data[, eval(column_name) := breaks[bin] +binwidth/2]
  return(data)
}

set_window_mz <- function(swath_windows){
  # Returns a list of swath windows, each element is a list of 
  # mz window start and end.
  
  # Parameters:
  # swath_windows   int     total number of swath windows
  mzmin <- 400
  mzmax <- 1200
  window_breaks <- seq(mzmin, mzmax, length.out = swath_windows + 1)
  window_mz <- lapply(1:(length(window_breaks)-1), function(i) window_breaks[i:(i+1)])
  return(window_mz)
}

getCoverage <- function(data_long, ow){
  # Parameters:
  # data_long   dataframe     filtered, binned data frame
  # ow          dataframe     window settings with mz and im windows
  # Returns the portions of precursors found in the windows
  total_precursors <- data_long[,.N]
  precursors <- numeric(0)
  for(i in 1:nrow(ow)){
    o <- ow[i,]
    precursors <- c(precursors, data_long[im >= o$imstart & im <= o$imend & mz >= o$mzstart & mz <= o$mzend, .N])
  }
  return(sum(precursors/total_precursors))
}

optimizeIMwindows <- function(data_long, window_mz, minim, maxim, swathpertims, im_peakwidth = 150,
                              variableWindows = F, considerIntensity =F, plot = F, fixed_ends = F){
  
  # Parameters:
  #     data_long      data table         binned data
  #     window_mz      list of lists      mz swath windows
  #     minim          int                min im scan number
  #     maxim          int                max im scan number
  #     swathpertims    int               number of mz windows per frame
  #     im_peakwith     int               im extraction peak width 
  #     fixed_end       boolean           fixed beginning and end IM extractions
  
  #assign frame groups
  groups <- rep(1:(length(window_mz)/swathpertims), swathpertims)

  plot_data <- copy(data_long)
  plot_data$window <- 0
  plot_data$group <- 1
  # get the window breaks
  allmz <- do.call("c", window_mz)
  plot_data$winmin <- 0
  plot_data$winmax <- 0
  plot_data <- plot_data[mz < max(allmz) & mz > min(allmz)] # filter values that do not fall into the windows
  # Assign the window and frame number to the data
  for(i in 1:length(window_mz)){
    plot_data[mz < max(window_mz[[i]]) & mz > min(window_mz[[i]])]$window <- i
    plot_data[mz < max(window_mz[[i]]) & mz > min(window_mz[[i]])]$winmin <- min(window_mz[[i]])
    plot_data[mz < max(window_mz[[i]]) & mz > min(window_mz[[i]])]$winmax <- max(window_mz[[i]])
    plot_data[mz < max(window_mz[[i]]) & mz > min(window_mz[[i]])]$group <- groups[i]
  }
  plot_data <- plot_data[window != 0]  
  if(considerIntensity){
    plot_data <- plot_data[, .(intensity = sum(intensity)), by = .(window, im, group, winmin, winmax)]
  }else{
    # sum up the number of data points in each window
    plot_data <- plot_data[, .(intensity = .N), by = .(window, im, group, winmin, winmax)]
  }
  if(variableWindows){
    print("Not yet implemented")
  }else{
    total_width <- im_peakwidth * swathpertims
    scan_range <- minim:(maxim - total_width)
    startIm <- numeric(0)
    
    for(g in unique(groups)){
      pd <- copy(plot_data[group == g])
      mobs <- data.frame(window = unique(sort(pd$window)))
      # Assign the window start and end im values
      mobs$imwindowmin <- minim + im_peakwidth * 0:(swathpertims-1)
      mobs$imwindowmax <- minim + im_peakwidth * 1:(swathpertims)
      pd <- merge(pd, mobs)
      ints <- numeric(0)
      # for each scan range, find the optimal number of data points
      for(s in scan_range){
        ints <-c(ints, sum(pd[im >= imwindowmin & im < imwindowmax]$intensity))
        pd$imwindowmin <- pd$imwindowmin + 1
        pd$imwindowmax <- pd$imwindowmax + 1
      }
      if(plot){
        p <-ggplot(pd, aes(x= im, y = intensity, color = window)) + geom_point() +
          geom_vline(xintercept = c(mobs$imwindowmin, mobs$imwindowmax) + which.max(ints)-1)
        print(p)
      }
      startIm <- c(startIm, which.max(ints)-1)
    }
    naiveims <- minim + ((maxim - minim) / swathpertims) * 0:swathpertims
    lower_bounds <- mobs$imwindowmin[rep(1:swathpertims, each = length(window_mz)/ swathpertims)] + rep(startIm, swathpertims)
    upper_bounds <- mobs$imwindowmax[rep(1:swathpertims, each = length(window_mz)/ swathpertims)] + rep(startIm, swathpertims) 
    cycle <- rep(1:swathpertims, each = length(window_mz)/swathpertims)
    res <- data.frame(window = 1:length(window_mz), group = groups, cycle = cycle, imstart = lower_bounds, imend = upper_bounds,
                      mzstart = sapply(window_mz, "[", 1), mzend = sapply(window_mz, "[", 2))
  }
  if(fixed_ends){
    res <- extendBoundaryWindows(res, minim, maxim)
  }
  return(res)
}


plot_window_settings <- function(ow, data_long, swathpertims){
  # Parameters:
  #   ow            list            list of window settings
  #   data_long     dataframe       binned ms data
  #   swathpertims  int             number of swath windows per frame
  #   plot.im_start colname         Specify colunm name for the im window start
  #   plot.im_end   colname         Specify colunm name for the im window end
  #   plot.mz_start colname         Specify column name for the mz window start
  #   plot.mz_end   colname         Specify column name for the mz window end
  mzmax <- 1200
  mzmin <- 400
  d <- data.frame(x1 = ow$mzstart, 
                  y1 = ow$imstart)
  d$x2 = ow$mzend
  d$y2 = ow$imend
  d$name <- rep(1:(swathes/swathpertims), swathpertims)
  cyclebounds <- mzmin + ((mzmax - mzmin) / swathpertims) * 0: swathpertims
  
  p <- ggplot(data_long) +
    geom_raster(aes(x = mz, y = im, fill = log(intensity)), interpolate=TRUE) +
    scale_fill_gradient2() +
    geom_rect(data = d, mapping= aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2), alpha = 0.6, color = "darkgrey") +
    geom_text(data=d, aes(x=x1+(x2-x1)/2, y=y1+(y2-y1)/2, label=name), size=4) +
    geom_vline(xintercept = cyclebounds) +
    theme_classic()
  p
}

optimize_window_overlaps <- function(ow, data_long){
  # Parameters:
  #   ow            list        list of optimized window frames from Optimize_IM_Window()
  #   data_long   dataframe     filtered and cleaned data of the im
  wins <- copy(ow)
  wins$imcutstart <- 0
  wins$imcutend <- 0
  wins[cycle==1,]$imcutstart <- 150
  for(i in unique(wins$group)){
    # isolate the windows by frame (group)
    frames <- wins[wins$group == i,] #assigns the list of windows for the first frame
    j <- 1
    while(j < length(frames$window)){
      mz1 <- df[mz > frames$mzstart[j] & mz < frames$mzend[j]]$im
      mz2 <- df[mz > frames$mzstart[j+1] & mz < frames$mzend[j+1]]$im
      d1 <- density(mz1, from=min(mz1), to=max(mz2), n=2048)
      d2 <- density(mz2, from=min(mz1), to=max(mz2), n=2048)
      inter <- d1$x[as.logical(abs(diff(d1$y < d2$y)))]
      # the maximum would be the last cutoff for intersection between 2 populations
      v <- max(inter)
      wins[wins$group==i,]$imcutend[j] <- v
      wins[wins$group==i,]$imcutstart[j+1] <- v
      j <- j+1
    }
    wins[wins$group == i,]$imcutend[j] <- 850
  }
  return(wins)
}


##############################
## Data import
##############################

# MaxQuant output
data <- fread("D:/TIMS-DIAPASEF/window_placement_setup/allPeptides.txt")
data <- data[Charge == 2 & !is.na(`MS/MS scan number`) & !is.na(Intensity)]
data <- data[,.(`m/z`, `Ion mobility index`, Intensity)]
names(data) <- c("mz", "im", "intensity")
data_long <- data[, isSC := FALSE]
data_long$im <-  max(data_long$im) -  data_long$im + min(data_long$im)
data_long <- binColumn(data_long, "im", 1)
data_long <- binColumn(data_long, "mz", 1)




swathpertims <- 4
ran <- seq((maxim-minim)/swathpertims - 60, (maxim-minim)/swathpertims, by =3)
win <- lapply(ran, function(x){
  optimizeIMwindows(data_long, window_mz, minim = minim, maxim = maxim, swathpertims = swathpertims, im_peakwidth = x)
})
coverage <- sapply(win, function(x) getCoverage(data_long, ow = x))
plot(ran, coverage)
im_pw <- ran[which.max(coverage)]
ow <- win[[which.max(coverage)]]

plot_window_settings(ow, data_long, swathpertims)
max(coverage) *100

im_span <- optimize_window_overlaps(ow, data_long)
im_new <- subset(im_span, select=-c(imstart,imend))
colnames(im_new)[6:7] <- c("imstart","imend")

plot_window_settings(im_new, data_long, swathpertims)

getCoverage(data_long, im_new)
# 0.869917
getCoverage(data_long,ow)
#0.865293

