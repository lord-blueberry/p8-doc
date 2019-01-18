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
            col.regions=colorRampPalette(brewer.pal(n =11, name="RdYlGn")),
            scales=list(x=scales, y=scales),
            xlab=xunits,
            ylab=""
  )
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
  df <- data.frame(points=point, values=value, names=name)
  return(df)
}

asinh <- scales::trans_new(name = 'asinh', 
                           transform = function(x) asinh(x*1000), 
                           inverse = function(x) sinh(x)/1000)

folder <- "./sim01/"
tclean <- read(paste(folder, "tclean.csv", sep=""), 256, 0.5) - read(paste(folder,"tclean.residual.csv", sep=""), 256, 0.5)
cd <- read(paste(folder, "image2", sep=""), 256, 0.5)

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

png(paste("./contour_points", ".png",sep=""),
    width = 8.0,
    height = 3.0,
    units = "in",
    res = 400)
print(ggplot(data = df, aes(x=df$points, y=df$values, colour=df$names)) + 
  geom_line() +
  scale_y_continuous(trans=asinh, breaks=c(0, 0.001, 0.01, 0.1, 1, 1.4, 2.5)) +
  xlab("arc seconds") +
  ylab("Jansky/beam") +
  labs(colour='Legend:') +
  scale_colour_brewer(palette = "Dark2") +
  theme(legend.text=element_text(size=11), 
        legend.title=element_text(size=13)))
dev.off()

scales = list(at=c(1, 65, 129, 197, 255))
png("tclean_points.png",
    width = 5.0,
    height = 5.0,
    units = "in",
    res = 400)
WriteMap2(tclean, at=seq(min(tclean), max(tclean), length.out=200), scales)
dev.off()
png("cd_points.png",
    width = 5.0,
    height = 5.0,
    units = "in",
    res = 400)
WriteMap2(cd, at=seq(min(cd), max(cd), length.out=200), scales)
dev.off()





folder <- "./sim00/"
resol <- 0.5/60
tclean <- read(paste(folder, "tclean.csv", sep=""), 1080, resol) - read(paste(folder,"tclean.residual.csv", sep=""), 1080, resol)
cd <- read(paste(folder, "image3", sep=""), 1080, resol)

skymodel <- t(read(paste(folder,"skymodel.csv", sep=""), 1080, resol))
model.axis <- round(0:(1079) * resol, digits=1)
colnames(skymodel) = model.axis
rownames(skymodel) = model.axis
pix = 1080
scales = list(at=c(1, pix/8+1, pix/4+1, pix/8*3+1 ,pix/2+1, pix/8*5+1, pix/4*3+1, pix/8*7+1, pix))
png("mixed_clean.png",
    width = 5.0,
    height = 5.0,
    units = "in",
    res = 400)
colorbreaks <- seq(min(tclean), max(tclean), length.out=200)
print(WriteMap2(tclean, at=colorbreaks, scales, xunits="arc minutes"))
dev.off()
png("mixed_cd.png",
    width = 5.0,
    height = 5.0,
    units = "in",
    res = 400)
colorbreaks <- seq(min(cd), 5, length.out=200)
colorbreaks <- c(colorbreaks, max(cd))
print(WriteMap2(cd, at=colorbreaks, scales, xunits="arc minutes"))
dev.off()


pix = 256
cut.row.names <- round(527:654*resol, digits=2)
cut.col.names <- round(round(816:943*resol, digits=2))
cut.row <- 528:655
cut.col <- 817:944
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
png("mixed_cut_model.png",
    width = 5.0,
    height = 5.0,
    units = "in",
    res = 400)
colorbreaks <- seq(min(cutout.model), 0.5, length.out=11)
colorbreaks <- c(colorbreaks, 0.58)
print(WriteMap2(cutout.model, at=colorbreaks, scales, xunits="arc minutes"))
dev.off()

matrices <- list(cut.model, cut.tclean, cut.cd)
names <- c("Ground Truth", "CLEAN", "Coordinate Descent")
interpolation <- 10000
p0 <- c(63+128, 60+128)
p1 <- c(67+128,65+128)
df <- calcLineDF(matrices, names, p0, p1, interpolation)
print(ggplot(data = df, aes(x=df$points, y=df$values, colour=df$names)) + 
        geom_line() +
        scale_y_continuous(trans=asinh, breaks=c(0, 0.001, 0.01, 0.1, 1, 1.4, 2.5)) +
        xlab("arc seconds") +
        ylab("Jansky/beam") +
        labs(colour='Legend:') +
        scale_colour_brewer(palette = "Dark2") +
        theme(legend.text=element_text(size=11), 
              legend.title=element_text(size=13)))