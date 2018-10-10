library(data.table)
library(ggplot2)

############################
## Functions
############################

#' Divide the ion mobility range in windows of equal size and return a list with the windows
getWindowMobilities <- function(ion_mobilities, windows_per_cycle){
  m <- floor(length(ion_mobilities) / windows_per_cycle)
  res <- lapply(1:windows_per_cycle, function(w) ion_mobilities[((w-1)*m+1):(m*w)])
  return(res)
}

#' Plots a histogram of ions in different ion mobility ranges over the m/z dimension
plotMZhist <- function(data_long, window_mobilities, SConly = F){
  plot_data <- copy(data_long)
  plot_data$window <- as.character(1)
  allmob <- do.call("c", window_mobilities)
  plot_data <- plot_data[im < max(allmob) & im > min(allmob)]
  for(i in 1:length(window_mobilities)){
    plot_data[im < max(window_mobilities[[i]]) & im > min(window_mobilities[[i]])]$window <- paste0(min(window_mobilities[[i]]), " - ", max(window_mobilities[[i]]))
  }
  plot_data <- plot_data[window != "1"]
  plot_data <- plot_data[, .(intensity = sum(intensity)), by = .(window, binmz)]
  if(SConly) plot_data <- plot_data[isSC == T]
  ggplot(plot_data, aes(x=binmz, y = intensity, color = window)) +
    geom_point() +
    geom_smooth(span=0.2) #+
  # facet_wrap(~ window)
}


#' Plots a histogram of ions in different m/z ranges over the ion mobility dimension
plotIMhist <- function(data_long, window_mz, SConly = F, showPoints = T, groups = NULL, considerIntensity = T){
  plot_data <- copy(data_long)
  plot_data$window <- as.character(1)
  plot_data$group <- 1
  allmz <- do.call("c", window_mz)
  plot_data <- plot_data[mz < max(allmz) & mz > min(allmz)]
  for(i in 1:length(window_mz)){
    plot_data[mz < max(window_mz[[i]]) & mz > min(window_mz[[i]])]$window <- paste0(min(window_mz[[i]]), " - ", max(window_mz[[i]]))
    if(!is.null(groups)){
      plot_data[mz < max(window_mz[[i]]) & mz > min(window_mz[[i]])]$group <- groups[i]
    }
  }
  plot_data <- plot_data[window != "1"]
  if(considerIntensity){
    plot_data <- plot_data[, .(intensity = sum(intensity)), by = .(window, im, group)]
  }else{
    plot_data <- plot_data[, .(intensity = .N), by = .(window, im, group)]
  }
  if(SConly) plot_data <- plot_data[isSC == T]
  p <- ggplot(plot_data, aes(x=im, y = intensity, color = window)) +
    geom_smooth(span=0.2, se = F) #+
  # facet_wrap(~ window)
  if(showPoints) p <- p + geom_point()
  if(!is.null(groups)){
    p <- p + facet_wrap(~ group)
  }
  p
}

#' Determine window boundaries and plot m/z histograms for those windows
plotMobilitywindows <- function(data_long, windows_per_cycle, ion_mobilities, SConly = F, binwidth = 5){
  if(SConly) data_long <- data_long[isSC == T]
  breaks <- seq(min(data_long$mz), max(data_long$mz), by = binwidth)
  data_long[,bin:=findInterval(mz, breaks)]
  data_long[, binmz:= breaks[bin] +binwidth/2]
  window_mobilities <- getWindowMobilities(ion_mobilities, windows_per_cycle)
  plotMZhist(data_long = data_long, window_mobilities = window_mobilities)
}

#' Compute variable windows in m/z dimension given a distribution of ions
findWindows <- function(features, windows = 32, considerIntensity = F){
  
  if(!considerIntensity) features$intensity = 1
  totalFeatures <- sum(as.numeric(features$intensity), na.rm = T)
  binFeatures <- totalFeatures/windows
  
  features <- features[order(features$mz),]
  ints <- as.numeric(features$intensity)
  ind = 1
  boundaries <- data.frame(window = 1:windows, start = 0, end = 0)
  for (w in 1:windows) {
    bf = 0
    boundaries[w,"start"] <- features[ind, "mz"]
    while (bf < binFeatures & ind < nrow(features)) {
      bf = sum(bf, ints[ind], na.rm = T)
      ind = ind + 1
    }
    boundaries[w,"end"] <- features[ind, "mz"]
    # Dynamically adjust the remaining features so that errors don't propagate
    remainFeatures <- sum(ints[ind:length(ints)], na.rm = T)
    binFeatures <- remainFeatures / (windows - w)
  }
  return(boundaries)
}

#' 
getEfficiency <- function(windows_per_cycle, windows, boundaries, data_long, fun_2 = function(x) 0.8*x - 285, total_intensity){
  mzdist<- data_long[, .(intensity = sum(intensity)), by = .(mz, isSC)]
  boundaries <- lapply(windows, function(w) findWindows(mzdist[isSC == F], windows = w, considerIntensity = T))
  df <- copy(data_long[isSC == F])
  df$intensity <- as.numeric(df$intensity)
  efficiency <- rep(0, length(windows))
  for(w in windows){
    if(w == 1){
      efficiency[w] <- 1
    }else if(w %% windows_per_cycle !=0){
      efficiency[w] <- NA
    }else{
      mzmidpoints <- boundaries[[w]][(w/windows_per_cycle) * 1:(windows_per_cycle-1), "end"]
      mzmin <- boundaries[[w]][1,"start"]
      mzmax <- boundaries[[w]][nrow(boundaries[[w]]),"end"]
      imbreaks <- c(min(df$im), fun_2(mzmidpoints), max(df$im))
      
      mzbreaks <-sort(unique(c(boundaries[[w]][, "start"], boundaries[[w]][, "end"])))
      
      df <- df[, bin := findInterval(im, imbreaks, rightmost.closed = T)]
      df <- df[, mzbin := findInterval(mz, mzbreaks, rightmost.closed = T)]
      
      plot_data <- df[, .(intensity = sum(intensity)), by = .(bin, mz)]
      p <- ggplot(plot_data, aes(x=mz, y = intensity, color = factor(bin))) +
        geom_point() +
        geom_vline(xintercept = c(mzmin, mzmidpoints, mzmax))
      print(p)
      ints <- df[, .(intensity=sum(intensity)), by =.(bin, mzbin)]
      ints <- ints[(mzbin/(w/windows_per_cycle)) > bin & (mzbin/(w/windows_per_cycle)) <= bin+1]
      efficencies <- numeric(0)
      for(i in 1:(w/windows_per_cycle)) efficencies<- c(efficencies, sum(as.numeric(ints[seq(1, w, by = (w/windows_per_cycle)) + i -1,intensity])))
      efficiency[w] <- mean(efficencies)/total_intensity
    }
    
  }
  return(efficiency)
}

extendBoundaryWindows <- function(ow, minim, maxim){
  mincycle <- min(ow$cycle)
  maxcycle <- max(ow$cycle)
  owres <- ow
  owres[owres$cycle == as.numeric(mincycle),]$imstart <- minim
  owres[owres$cycle == as.numeric(maxcycle),]$imend <- maxim
  return(owres)
}


#' Optimize the placement of windows so that they cover as many precursors as possible
optimizeIMwindows <- function(data_long, window_mz, minim, maxim, swathpertims, im_peakwidth = 150,
                              variableWindows = F, considerIntensity =F, plot = F, fixed_ends = F){
  groups <- rep(1:(length(window_mz)/swathpertims), swathpertims)
  plot_data <- copy(data_long)
  plot_data$window <- 0
  plot_data$group <- 1
  allmz <- do.call("c", window_mz)
  plot_data$winmin <- 0
  plot_data$winmax <- 0
  plot_data <- plot_data[mz < max(allmz) & mz > min(allmz)]
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
      mobs$imwindowmin <- minim + im_peakwidth * 0:(swathpertims-1)
      mobs$imwindowmax <- minim + im_peakwidth * 1:(swathpertims)
      pd <- merge(pd, mobs)
      ints <- numeric(0)
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

#' Calculate the percentage of precursors covered bz a specific window setting
getCoverage <- function(data_long, ow){
  total_precursors <- data_long[,.N]
  precursors <- numeric(0)
  for(i in 1:nrow(ow)){
    o <- ow[i,]
    precursors <- c(precursors, data_long[im >= o$imstart & im <= o$imend & mz >= o$mzstart & mz <= o$mzend, .N])
  }
  return(sum(precursors/total_precursors))
}
getEff <- function(data_long, ow, considerIntensity = F){
  dl <- copy(data_long)
  if(!considerIntensity){
    dl[, intensity := 1]
  }

  gr <- unique(ow$group)
  ow$ncycle <- ceiling(ow$window/max(gr))
  efficiency <- numeric(0)
  for(c in unique(ow$ncycle)){
    cycleminmz <- min(ow[ow$ncycle == c,]$mzstart)
    cyclemaxmz <- max(ow[ow$ncycle == c,]$mzend)
    total_intensity <- sum(data_long[mz >= cycleminmz & mz <= cyclemaxmz, .N])
    # print(total_intensity)
    o <- ow[ow$ncycle == c,]
    winint <- apply(o, 1, function(x) sum(data_long[mz >= x["mzstart"] & mz <= x["mzend"] & im >= x["imstart"] & im <= x["imend"], .N]))
    # print(mean(winint))
    efficiency <- c(efficiency, mean(winint)/total_intensity)
  }
  # total_intensity <- sum(dl[,]$intensity)
  # winint <- apply(ow, 1, function(x) sum(dl[mz >= x["mzstart"] & mz <= x["mzend"] & im >= x["imstart"] & im <= x["imend"]]$intensity))
  # efficiency <- mean(winint)/total_intensity
  
  return(mean(efficiency))
}


#' Bin the values in a column of a data frame
binColumn <- function(data, column_name, binwidth){
  cmi <- min(data[, get(column_name)])
  cma <- max(data[, get(column_name)])
  breaks <- seq(cmi, cma, by = binwidth)
  data[,bin:=findInterval(get(column_name), breaks, rightmost.closed = T)]
  data[, eval(column_name) := breaks[bin] +binwidth/2]
  return(data)
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
data_long$im <-  max(data_long$im) -  data_long$im + min(data_long$im) @mfrank
data_long <- binColumn(data_long, "im", 1)
data_long <- binColumn(data_long, "mz", 1)



############################################
# Generate the window settings
############################################

## General Parameter settings
#setwd("~/roestlab/TIMS_DIA/results/2018-03-15/comparison_window_setup")
swathes <- 32
mzmin <- 400
mzmax <- 1200
window_breaks <- seq(mzmin, mzmax, length.out = swathes + 1)
window_mz <- lapply(1:(length(window_breaks)-1), function(i) window_breaks[i:(i+1)])
scanToIM <- function(scannr){
  scans <- 893
  im0 <- 0.575021
  imM <- 1.535878
  ims <- seq(from = im0, to = imM, length.out = scans) 
  return(ims[scannr])
}

ggplot(data_long[mz >= mzmin & mz <= mzmax], aes(x = im)) +
  geom_histogram(binwidth = 10) +
  geom_vline(xintercept = 850, color = "red") +
  geom_vline(xintercept = 150, color = "red")

ggsave("imcutoffs.png")

# Minimum ion mobility
scanToIM(150)
minim = 150
# Maximum ion mobility
scanToIM(850)
maxim = 850


## Schema 2
wins <- data.table(do.call(rbind, window_mz))
names(wins) <- c("start", "end")

win2 <- wins[seq(from =1, to = nrow(wins) -1, by = 2),]
win2$window_nr <- 1:nrow(win2)
win2[, width := end - start]
win2[, center := start + (end - start) / 2]
write.csv(win2, "windows_schema2.csv", quote = F, row.names = F)

win2_2 <- wins[seq(from =2, to = nrow(wins), by = 2),]
win2_2$window_nr <- 1:nrow(win2_2)
win2_2[, width := end - start]
win2_2[, center := start + (end - start) / 2]
write.csv(win2_2, "windows_schema2_2.csv", quote = F, row.names = F)

## Schema 3

  # Same as schema 2

## Schema 4

# Optimize the im peakwidth
swathpertims <- 4
ran <- seq((maxim-minim)/swathpertims - 60, (maxim-minim)/swathpertims, by =3)
win <- lapply(ran, function(x){
  optimizeIMwindows(data_long, window_mz, minim = minim, maxim = maxim, swathpertims = swathpertims, im_peakwidth = x)
})
coverage <- sapply(win, function(x) getCoverage(data_long, ow = x))
plot(ran, coverage)
im_pw <- ran[which.max(coverage)]


ow <- optimizeIMwindows(data_long, window_mz, minim = minim, maxim = maxim, swathpertims, im_peakwidth = im_pw, fixed_ends = T)

d <- data.frame(x1 = ow$mzstart, 
                y1 = ow$imstart)
d$x2 = ow$mzend
d$y2 = ow$imend
d$name <- rep(1:(swathes/swathpertims), swathpertims)
cyclebounds <- mzmin + ((mzmax - mzmin) / swathpertims) * 0: swathpertims

ggplot(data_long) +
  geom_raster(aes(x = mz, y = im, fill = log(intensity)), interpolate=TRUE) +
  scale_fill_gradient2() +
  geom_rect(data = d, mapping= aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2), alpha = 0.6, color = "darkgrey") +
  geom_text(data=d, aes(x=x1+(x2-x1)/2, y=y1+(y2-y1)/2, label=name), size=4) +
  geom_vline(xintercept = cyclebounds) +
  theme_classic()

ggsave("4windowsPerCycleVis.png")

win2 <- data.table(ow[seq(from =1, to = nrow(ow) -1, by = 2),])
win2[, mzwidth := mzend - mzstart]
win2[, mzcenter := mzstart + (mzend - mzstart) / 2]
setnames(win2, "imstart", "scanstart")
setnames(win2, "imend", "scanend")
win2[, imstart := scanToIM(scanstart)]
win2[, imend := scanToIM(scanend)]
win2[, imwidth := imend - imstart]
win2[, imcenter := imstart + (imend - imstart) / 2]
write.csv(win2, "windows_schema4.csv", quote = F, row.names = F)

win2_2 <- data.table(ow[seq(from =2, to = nrow(ow), by = 2),])
win2_2[, mzwidth := mzend - mzstart]
win2_2[, mzcenter := mzstart + (mzend - mzstart) / 2]
setnames(win2_2, "imstart", "scanstart")
setnames(win2_2, "imend", "scanend")
win2_2[, imstart := scanToIM(scanstart)]
win2_2[, imend := scanToIM(scanend)]
win2_2[, imwidth := imend - imstart]
win2_2[, imcenter := imstart + (imend - imstart) / 2]
write.csv(win2_2, "windows_schema4_2.csv", quote = F, row.names = F)

## Schema 5

# Optimize the im peakwidth

swathpertims <- 8

ran <- seq((maxim-minim)/swathpertims - 30, (maxim-minim)/swathpertims, by =3)
win <- lapply(ran, function(x){
  optimizeIMwindows(data_long, window_mz, minim = minim, maxim = maxim, swathpertims = swathpertims, im_peakwidth = x)
})
coverage <- sapply(win, function(x) getCoverage(data_long, ow = x))
plot(ran, coverage)
im_pw <- ran[which.max(coverage)]


ow <- optimizeIMwindows(data_long, window_mz, minim = minim, maxim = maxim, swathpertims = swathpertims, im_peakwidth = im_pw, fixed_ends = T)

d <- data.frame(x1 = ow$mzstart, 
                y1 = ow$imstart)
d$x2 = ow$mzend
d$y2 = ow$imend
d$name <- rep(1:(swathes/swathpertims), swathpertims)
cyclebounds <- mzmin + ((mzmax - mzmin) / swathpertims) * 0: swathpertims

ggplot(data_long) +
  geom_raster(aes(x = mz, y = im, fill = log(intensity)), interpolate=TRUE) +
  scale_fill_gradient2() +
  geom_rect(data = d, mapping= aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2), alpha = 0.6, color = "darkgrey") +
  geom_text(data=d, aes(x=x1+(x2-x1)/2, y=y1+(y2-y1)/2, label=name), size=4) +
  geom_vline(xintercept = cyclebounds) +
  theme_classic()

ggsave("8windowsPerCycleVis.png")


win2 <- data.table(ow[seq(from =1, to = nrow(ow) -1, by = 2),])
win2[, mzwidth := mzend - mzstart]
win2[, mzcenter := mzstart + (mzend - mzstart) / 2]
setnames(win2, "imstart", "scanstart")
setnames(win2, "imend", "scanend")
win2[, imstart := scanToIM(scanstart)]
win2[, imend := scanToIM(scanend)]
win2[, imwidth := imend - imstart]
win2[, imcenter := imstart + (imend - imstart) / 2]
write.csv(win2, "windows_schema5.csv", quote = F, row.names = F)

win2_2 <- data.table(ow[seq(from =2, to = nrow(ow), by = 2),])
win2_2[, mzwidth := mzend - mzstart]
win2_2[, mzcenter := mzstart + (mzend - mzstart) / 2]
setnames(win2_2, "imstart", "scanstart")
setnames(win2_2, "imend", "scanend")
win2_2[, imstart := scanToIM(scanstart)]
win2_2[, imend := scanToIM(scanend)]
win2_2[, imwidth := imend - imstart]
win2_2[, imcenter := imstart + (imend - imstart) / 2]
write.csv(win2_2, "windows_schema5_2.csv", quote = F, row.names = F)



