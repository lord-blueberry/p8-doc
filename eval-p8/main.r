library(rgdal)
library(rasterVis)
library(viridis)
library(RColorBrewer)
library(graphics)
library(ggplot2)
library(scales)

WriteMap2 <- function(x, at, scales, zlog=FALSE, xunits="arc seconds") {
  par(mar=c(0,0,0,0),plt=c(0,0,0,0),oma=c(0,0,0,0))
  levelplot(x,
            at = at,
            margin=FALSE,
            zscaleLog = zlog,
            col.regions= colorRampPalette(c("gray40","royalblue4", "dodgerblue4", "darkgreen", "darkolivegreen4" ,"darkgoldenrod", "orchid", "plum", "lightblue","lightcyan", "white"), bias = 1),
            #col.regions=colorRampPalette(brewer.pal(n =11, name="RdYlGn")),
            scales=list(x=scales, y=scales),
            xlab=xunits,
            ylab=""
  )
}

img.out <- function(name) {
  png(name,
      width = 5.0,
      height = 5.0,
      units = "in",
      res = 400)
}

read <- function(img, pixels, resolution) {
  axis <- round(0:(pixels-1) * resolution, digits=1)
  data <- read.table(img, header = FALSE, sep = ",")
  data = data.matrix(data)
  colnames(data) = axis
  rownames(data) = axis
  return(data)
}

calcline <- function(matrix, p0, p1, length.out=100) {
  a <- p1 - p0
  x0 <- p0 + a*((-p0[2]+1)/ a[2])
  if(x0[1] > nrow(matrix) | x0[1] < 1) {
    x0 <- p0 + a*((nrow(matrix)-p0[1])/ a[1])
  }
  y0 <- p0 + a*((-p0[1]+1)/ a[1])
  if(y0[2] > ncol(matrix) | y0[2] < 1 ) {
    y0 <- p0 + a*((ncol(matrix)-p0[2])/ a[2])
  }

  a <- y0 - x0
  line <-rep(c(0),length.out=length.out)
  i = 1
  for(t in seq(0, 1, length.out = length.out)) {
    pt <- x0+a*t
    pt.ceil <- ceiling(pt)
    pt.floor <- floor(pt)
    pt.factor <- pt.ceil - pt 
    
    #interpolate
    xt1 <- pt.factor[1] * matrix[pt.floor[1], pt.ceil[2]] + (1-pt.factor[1]) * matrix[pt.ceil[1], pt.ceil[2]]
    xt0 <- pt.factor[1] * matrix[pt.floor[1], pt.floor[2]] + (1-pt.factor[1]) * matrix[pt.ceil[1],pt.floor[2]]
    x <- pt.factor[2] * xt0 + (1-pt.factor[2])*xt1
    line[i] <- x
    
    i = i+ 1
  }
  return(line)
}

calcLineDF <- function(matrices, names, p0,p1, interpolation.length) {
  index <- 1:interpolation.length
  point <- c()
  value <- c()
  name <- c()
  for(i in 1:length(matrices)) {
    m <- matrices[[i]]
    n <- names[i]
    line <- calcline(m, p0, p1, interpolation.length)
    
    point <- c(point, index)
    value <- c(value, line)
    name <- c(name, rep(n, interpolation.length))
  }
  df <- data.frame(points=point, values=value, Legend=name)
  return(df)
}

asinh <- scales::trans_new(name = 'asinh', transform = function(x) asinh(x*1000), 
                           inverse = function(x) sinh(x)/1000)

folder <- "./sim01/"
tclean <- read(paste(folder, "tclean.csv", sep=""), 256, 0.5) - read(paste(folder,"tclean.residual.csv", sep=""), 256, 0.5)
cd <- read(paste(folder, "image0", sep=""), 256, 0.5)

skymodel <- t(read(paste(folder,"skymodel.csv", sep=""), 512, 0.5))
skymodel <- skymodel[129:384, 129:384]
model.axis <- round(0:(255) * 0.5, digits=1)
colnames(skymodel) = model.axis
rownames(skymodel) = model.axis
p0 <- c(128, 92)
p1 <- c(101, 159)
matrices <- list(skymodel, tclean, cd)
names <- c("Ground Truth", "CLEAN", "Coordinate Descent")

interpolation <- 10000
df <- calcLineDF(matrices, names, p0, p1, interpolation)
df$points <- df$points / interpolation* 0.5 * 256
png(paste("./points/contour_points", ".png",sep=""),
    width = 10.0,
    height = 4.0,
    units = "in",
    res = 400)
print(ggplot(data = df, aes(x=points, y=values, colour=Legend)) + 
        geom_line() +
        scale_y_continuous(trans=asinh, breaks=c(0, 0.001, 0.01, 0.1, 1, 1.4, 2.5)) +
        geom_polygon(aes(fill=Legend), alpha=0.1) +
        xlab("arc seconds") +
        ylab("Jansky/beam") +
        
        theme(legend.text=element_text(size=11), 
              legend.title=element_text(size=13)))
dev.off()

scales = list(at=c(1, 65, 129, 197, 255))
img.out("./points/tclean_points.png")
WriteMap2(tclean, at=seq(min(tclean), max(tclean), length.out=200), scales)
dev.off()
img.out("./points/cd_points.png")
cd_copy <- cd
cd_copy[cd_copy > 0.02] = 0.02 
colorbreaks <- seq(min(cd_copy), max(cd_copy), length.out=200)
WriteMap2(cd_copy, at=colorbreaks, scales)
dev.off()
png("./points/skymodel.png",
    width = 4.0,
    height = 4.0,
    units = "in",
    res = 200)
colorbreaks <- seq(min(skymodel), max(skymodel), length.out=200)
WriteMap2(skymodel, at=colorbreaks, scales)
dev.off()





folder <- "./sim00/"
resol <- 0.5/60
tclean <- read(paste(folder, "tclean.csv", sep=""), 1080, resol) - read(paste(folder,"tclean.residual.csv", sep=""), 1080, resol)
cd <- read(paste(folder, "image1", sep=""), 1080, resol)
cache <- read(paste(folder, "full_cache_debug", sep=""), 1080, resol)

skymodel <- t(read(paste(folder,"skymodel.csv", sep=""), 1080, resol))
model.axis <- round(0:(1079) * resol, digits=2)
colnames(skymodel) = model.axis
rownames(skymodel) = model.axis
pix = 1080
scales = list(at=c(1, pix/8+1, pix/4+1, pix/8*3+1 ,pix/2+1, pix/8*5+1, pix/4*3+1, pix/8*7+1, pix))
img.out("./mixed/mixed_clean.png")
colorbreaks <- seq(min(tclean), max(tclean), length.out=200)
print(WriteMap2(tclean, at=colorbreaks, scales, xunits="arc minutes"))
dev.off()
img.out("./mixed/mixed_cd.png")
cd_copy <- cd
cd_copy[cd_copy > 1] = 1
colorbreaks <- seq(min(cd_copy), max(cd_copy), length.out=200)
print(WriteMap2(cd_copy, at=colorbreaks, scales, xunits="arc minutes"))
dev.off()
img.out("./mixed/mixed_cache.png")
colorbreaks <- seq(min(cache), max(cache), length.out=200)
print(WriteMap2(cache, at=colorbreaks, scales, xunits="arc minutes"))
dev.off()

pix = 256
cut.row.names <- round(463:718*resol, digits=2)
cut.col.names <- round(round(752:1007*resol, digits=2))
cut.row <- 464:719
cut.col <- 464:1008
scales = list(at=c(1, pix/8+1, pix/4+1, pix/8*3+1 ,pix/2+1, pix/8*5+1, pix/4*3+1, pix/8*7+1, pix))
cut.model <- skymodel[cut.row, cut.col]
rownames(cut.model) =cut.row.names
colnames(cut.model) =cut.col.names
cut.tclean <- tclean[cut.row, cut.col]
rownames(cut.tclean) =cut.row.names
colnames(cut.tclean) =cut.col.names
cut.cd <- cd[cut.row, cut.col]
rownames(cut.cd) =cut.row.names
colnames(cut.cd) =cut.col.names
img.out("./mixed/mixed_cut_model.png")
cut.model_copy <- cut.model
cut.model_copy[cut.model_copy > 0.6] = 0.6
colorbreaks <- seq(min(cut.model_copy), max(cut.model_copy), length.out=200)
print(WriteMap2(cut.model_copy, at=colorbreaks, scales, xunits="arc minutes"))
dev.off()

matrices <- list(cut.model, cut.tclean, cut.cd)
names <- c("Ground Truth", "CLEAN", "Coordinate Descent")
interpolation <- 10000
p0 <- c(63+64, 60+64)
p1 <- c(67+64,65+64)
df <- calcLineDF(matrices, names, p0, p1, interpolation)
df$points <- resol*464 + df$points / interpolation* resol * 256
png("./mixed/mixed_cut0.png",
    width = 8.0,
    height = 4.0,
    units = "in",
    res = 400)
print(ggplot(data = df, aes(x=df$points, y=df$values, colour=Legend)) + 
        geom_line() +
        scale_y_continuous(trans=asinh, breaks=c(0, 0.001, 0.01, 0.1, 1, 10, 100,1000)) +
        geom_polygon(aes(fill=Legend), alpha=0.1) +
        xlab("arc seconds") +
        ylab("Jansky/beam") +
        theme(legend.text=element_text(size=11), 
              legend.title=element_text(size=13)))
dev.off()









pix = 512
cut.row.names <- round(271:782*resol, digits=2)
cut.col.names <- round(round(251:762*resol, digits=2))
cut.row <- 272:783
cut.col <- 252:763
scales = list(at=c(1, pix/8+1, pix/4+1, pix/8*3+1 ,pix/2+1, pix/8*5+1, pix/4*3+1, pix/8*7+1, pix))
cut.model <- skymodel[cut.row, cut.col]
rownames(cut.model) =cut.row.names
colnames(cut.model) =cut.col.names
cut.tclean <- tclean[cut.row, cut.col]
rownames(cut.tclean) =cut.row.names
colnames(cut.tclean) =cut.col.names
cut.cd <- cd[cut.row, cut.col]
rownames(cut.cd) =cut.row.names
colnames(cut.cd) =cut.col.names

matrices <- list(cut.model, cut.tclean, cut.cd)
names <- c("Ground Truth", "CLEAN", "Coordinate Descent")
#p1 <- c(531-271, 501-251)
#p0 <- c(546-271, 586-251)
p0 <- c(544-271, 586-251)
p1 <- c(531-271, 501-251)
#p1 <- c(513-271, 551-251)
df <- calcLineDF(matrices, names, p0, p1, interpolation)
png("./mixed/mixed_cut2.png",
    width = 8.0,
    height = 4.0,
    units = "in",
    res = 400)
print(ggplot(data = df, aes(x=df$points, y=df$values, colour=Legend)) + 
        geom_line() +
        scale_y_continuous(trans=asinh, breaks=c(0, 0.001, 0.01, 0.1, 1, 10, 100,1000)) +
        geom_polygon(aes(fill=Legend), alpha=0.1) +
        xlab("arc seconds") +
        ylab("Jansky/beam") +
        theme(legend.text=element_text(size=11), 
              legend.title=element_text(size=13)))
dev.off()
png("./mixed/mixed_cut_model2.png",
    width = 5.0,
    height = 5.0,
    units = "in",
    res = 400)
cut.model_copy <- cut.model
cut.model_copy[cut.model_copy > 0.7] = 0.7
colorbreaks <- seq(min(cut.model_copy), max(cut.model_copy), length.out=200)
print(WriteMap2(cut.model_copy, at=colorbreaks, scales, xunits="arc minutes"))
dev.off()





pix = 2048
scales = list(at=c(1, pix/8+1, pix/4+1, pix/8*3+1 ,pix/2+1, pix/8*5+1, pix/4*3+1, pix/8*7+1, pix))
arcmin <- 3.2/60
meerkat <- read("./sim_meerkat/meerkat.csv", pix, arcmin)

stopped <- meerkat
stopped[stopped < 0] = 0
stopped[stopped > 0.01] = 0.01
png("meerkat.png",
    width = 6.0,
    height = 6.0,
    units = "in",
    res = 600)
colorbreaks <- seq(min(stopped), max(stopped), length.out=800)
par(mar=c(0,0,0,0),plt=c(0,0,0,0),oma=c(0,0,0,0))
print(levelplot(stopped,
          at = colorbreaks,
          margin=FALSE,
          zscaleLog = FALSE,
          #col.regions=colorRampPalette(brewer.pal(n =11, name="RdYlGn")),
          col.regions= colorRampPalette(c("gray40","royalblue4", "dodgerblue4", "darkgreen", "darkolivegreen4" ,"darkgoldenrod","plum4", "orchid", "lightblue","lightcyan", "white"), bias = 1),  
          scales=list(x=scales, y=scales),
          xlab="arc minutes",
          ylab=""
))
dev.off()