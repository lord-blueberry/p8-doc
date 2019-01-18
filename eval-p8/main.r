library(rgdal)
library(rasterVis)
library(viridis)
library(RColorBrewer)
library(graphics)
library(ggplot2)
library(scales)


WriteMap2 <- function(x, at, scales) {
  par(mar=c(0,0,0,0),plt=c(0,0,0,0),oma=c(0,0,0,0))
  levelplot(x,
            at = at,
            margin=FALSE,
            zscaleLog = NULL,
            col.regions=colorRampPalette(brewer.pal(n =11, name="RdYlGn")),
            scales=list(x=scales, y=scales),
            xlab="arc seconds",
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
  if(x0[1] > nrow(matrix)) {
    x0 <- p0 + a*((nrow(matrix)-p0[1])/ a[1])
  }
  y0 <- p0 + a*((-p0[1]+1)/ a[1])
  if(y0[2] > ncol(matrix)) {
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

png(paste("./contour_points", ".png",sep=""),
    width = 8.0,
    height = 3.0,
    units = "in",
    res = 400)
df <- calcLineDF(matrices, names, p0, p1, 10000)
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







lseq <- function(from=1, to=100000, length.out=6) {
  # logarithmic spaced sequence
  # blatantly stolen from library("emdbook"), because need only this
  exp(seq(log(from), log(to), length.out = length.out))
}

folder <- "./sim00/"
tclean <- read(paste(folder, "tclean.csv", sep=""), 1080, 0.5) - read(paste(folder,"tclean.residual.csv", sep=""), 1080, 0.5)
cd <- read(paste(folder, "image3", sep=""), 1080, 0.5)

skymodel <- t(read(paste(folder,"skymodel.csv", sep=""), 1080, 0.5))
model.axis <- round(0:(1079) * 0.5, digits=1)
colnames(skymodel) = model.axis
rownames(skymodel) = model.axis


scales = list(at=c(1, 32, 64, 96, 128, 256))
colorbreaks <- seq(min(tclean), max(tclean), length.out=200)

WriteMap2(tclean, at=seq(min(tclean), max(tclean), length.out=200), scales)

colorbreaks <- seq(min(cd), 10, length.out=200)
colorbreaks <- c(colorbreaks, max(cd))
WriteMap2(cd, at=colorbreaks, scales)