#!/usr/bin/env Rscript
### TODO 
### p-values
### labels, collapse

#############################################
# DEFINE configurable parameters
xmax <- 0.7
ymax <- 0.7
# do not label genes with low recurrence
min_recurrence_for_label <- 0.1
#min_pvalue_for_label <- 0.15

required.packages <- c("ggplot2", "binom", "plyr", "stats")

#############################################
# automatically install packages
new.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
invisible(lapply(required.packages, require, char=TRUE))

#############################################
# READ FIVE COLUMNS FROM STANDARD INPUT:
# 1) GENE NAME
# 2) NUMBER OF MUTATIONS IN CONDITION X (X-axis)
# 3) NUMBER OF MUTATIONS IN CONDITION Y (Y-axis)
# 4) NUMBER OF PATIENTS FOR CONDITION X (X-axis)
# 5) NUMBER OF PATIENTS FOR CONDITION Y (Y-axis)

d <- read.table(file("stdin"), header=TRUE)

# get the condition names from the header of the input file
conditionx <- names(d)[2]
conditiony <- names(d)[3]

# rename columns
names(d)[1:3] <- c("gene", "x", "y")

# get number of patients in each condition 
# from the first line only
nx <- d[1,4]
ny <- d[1,5]

# ALTERNATIVELY:
# d <- scan(file("stdin"), what=list(gene="",x=0,nx=0,y=0,ny=0), quiet=TRUE);
# d <- as.data.frame(d)

#############################################
# ARGUMENTS
args <- commandArgs(TRUE)
# conditionx <- args[[1]]
# nx <- as.integer(args[[2]])
# conditiony <- args[[3]]
# ny <- as.integer(args[[4]])

#############################################
# ROUND AXIS LIMITS TO WHOLE NUMBER OF PATIENTS
xmax <- ceiling(xmax*nx)/nx
ymax <- ceiling(ymax*ny)/ny

#############################################
# GROUP GENES WITH THE SAME RECURRENCE
d <- aggregate(data=d, gene ~ x + y, paste, collapse=", ")

#############################################
# CALCULATE CONFIDENCE INTERVAL USING BINOM PACKAGE
# 68% is one standard deviation
d <- adply(d, 1, transform, p.value = prop.test(c(x, y), c(nx, ny))[["p.value"]])

#############################################
# CALCULATE CONFIDENCE INTERVAL USING BINOM PACKAGE
confint <- function(x, nx, var) {
  # 0.68 is one standard deviation
  binom.confint(x, nx, conf.level=0.68, methods="exact")[,var] * nx
}
d <- adply(d, 1, transform, 
      x_lower = confint(x, nx, "lower"),
      x_upper = confint(x, nx, "upper"),
      y_lower = confint(y, ny, "lower"),
      y_upper = confint(y, ny, "upper")
)

#############################################
# PRINT TABLE WITH P-VALUES
options(width=200)
print(d[order(d$p.value),])

#############################################
# CALCULATE BACKGROUND P-VALUES
# nx <- d[1,]$nx
# ny <- d[1,]$ny
permutations <- expand.grid(x=0:(xmax*nx), nx=nx, y=0:(ymax*ny), ny=ny)
permutations <- adply(permutations, 1, transform, 
                      p.value = prop.test(c(x, y), c(nx, ny))[["p.value"]])

#############################################
# SET GENE LABELS
d$label <- d$gene
d[which(d$x/nx < min_recurrence_for_label & 
          d$y/ny < min_recurrence_for_label),]$label <- ""

#############################################
# SET DEFAULT PLOT SETTINGS
# white background (instead of grey default)
ggplot <- function(...) { ggplot2::ggplot(...) + theme_bw() }

#############################################
# BUILD PLOT

# gene coordinates
p <- ggplot(d, aes(x=x/nx, y=y/ny))

# axis scales, 'expand' removes margins
p <- p + scale_x_continuous(expand=c(0,0), limits=c(0,xmax)) 
p <- p + scale_y_continuous(expand=c(0,0), limits=c(0,ymax)) 

# p-value shading
# subtract 0.5 from x and y to center each box
p <- p + geom_tile(data=permutations, 
                   aes(x=(x-0.5)/nx, y=(y-0.5)/ny, 
                       fill=permutations$p.value)) 
p <- p + scale_fill_gradient("p-value", low="dark grey", high="white", trans="log10")

# center line
p <- p + geom_abline(a=0, b=1, col="dark grey")

# error bars
p <- p + geom_errorbarh(aes(xmin=d$x_lower/nx, xmax=d$x_upper/nx), col="dark grey", height=0)
p <- p + geom_errorbar(aes(ymin=d$y_lower/ny, ymax=d$y_upper/ny), col="dark grey", width=0) 

# labels for each point
p <- p + geom_text(aes(label=d$label), hjust=-0.2, vjust=0.5, angle = 45) 

# labels for the axes
p <- p + xlab(paste0("Recurrence of mutation in ", nx, " patients with ", conditionx))
p <- p + ylab(paste0("Recurrence of mutation in ", ny, " patients with ", conditiony))

# draw points
p <- p + geom_point()

# save to file
ggsave(p, file='RD.png')
