library("tess3r")
library("dplyr")

# get environment variable
input <- commandArgs(trailingOnly = TRUE)[1] # plink raw file
out <- commandArgs(trailingOnly = TRUE)[2] # output prefix
sp <- commandArgs(trailingOnly = TRUE)[3] # specoes

# prep output names
out_Rdat <- paste0(out, "pseudo.Rdata")
out_plot <- paste0(out, "pseudo-K2.png")
out_coords <- paste0(out, "pseudo_coords.txt")

cat("input file:", input, "\n")
cat("output prefix:", out, "\n")

###### convert to appropriate format
cat("\nreading genotype matrix\n")
dat2 <- read.table(gzfile(input), header = TRUE, na = "N")
# TRUE/FALSE whether same as major
df2 <- dat2$major == dat2[4:ncol(dat2)]
# FALSE = 0, TRUE = 2 bc no heterozygote
df2 <- df2*1 
# transpose it
df2 <- t(df2)
# write.table(df2, file = out_dat,
#             col.names = FALSE, row.names = FALSE)

##### prepare the coordinates
cat("reading sample info\n")
info <- read.table("data/sampleInfo.txt", sep = "\t", skip = 1, fill = TRUE)
# reorder the same as in input
bam <- read.csv("/dss/dsslegfs01/pr53da/pr53da-dss-0019/projects/2024__SexBias_GeneFlow/Ychromosome/lists/sexlist_CSLGSLSSL_nonrelatives.csv", header = TRUE, sep = "\t")
bam <- bam[bam$Species == sp, c("SampleID", "Population")]
info <- inner_join(bam, info, by = c("SampleID" = "V1"))

# convert to decimal degrees
coord <- info[, c("V8", "V7")]
# only S is negative
tmp <- grepl(" S", coord$V7)
coord$V7 <- as.numeric(gsub(" N| S", "", coord$V7))
coord$V7 <- ifelse(tmp == TRUE, -coord$V7, coord$V7)
# only W is negative
tmp <- grepl(" W", coord$V8)
coord$V8 <- as.numeric(gsub(" W| E", "", coord$V8))
coord$V8 <- ifelse(tmp == TRUE, -coord$V8, coord$V8)
coord <- as.matrix(coord)


###### run tess3
cat("running tess3\n")
cat("seed set to 31102025\n")
set.seed(31102025)
tess3 <- tess3(X = df2, coord = coord, K = 1:9, rep = 50, keep = "best",
                   method = "projected.ls", ploidy = 1) 
save(tess3, coord, info, file = out_Rdat)

# plot cross-validation scores
png(filename=out_plot)
plot(tess3, pch = 19, col = "green",
     xlab = "Number of ancestral populations",
     ylab = "Cross-validation score")
dev.off()
cat("TESS3 run complete\n")


## prepare map and ancestry coefficient plots for K=3
#load(out_Rdat)

# species raster
#sp <- gsub("^.*__|_.*$", "", input)
raster <- paste0(sp, "_raster.asc")
# species colours
if(sp == "GSL") {
        cols = c("#11819f", "#00c6c2") #"#d0d3fb", 
} else if(sp == "CSL") {
        cols = c("#ffc8b5", "#c55065", "#ff8178")
} else if(sp == "SSL") {
        cols = c("#f5f455", "#770076") #, "#e459b9"
}
my.palette <- CreatePalette(cols, 8)

# window size for map
winsize <- NULL # default window size
height <- 4680
# adjust for SSL 360 degree coords
if (unique(info$V3) == "SSL") {
        coord[,1] <- ifelse(coord[,1] < 0,
                          coord[,1] + 360,
                          coord[,1])
        winsize <- c(120, 240, 30, 70)
        height <- 2590
}


# ancestry coefficients 
cat("plotting ancestry coefficients for K=3\n")
q.matrix <- qmatrix(tess3, K = 2)
write.table(q.matrix, file = paste0(out, "-props_K2.txt"), 
            col.names = FALSE, row.names = FALSE, sep = "\t")
png(filename=paste0(out, "-props-pseudo-K2.png"))
barplot(q.matrix, border = NA, space = 0, 
        xlab = "Individuals", ylab = "Ancestry proportions for K = 2", 
        main = "Ancestry matrix", col.palette = my.palette, names.arg = rownames(q.matrix)) -> bp
#axis(1, at = 1:nrow(q.matrix), labels = bp$order, las = 3, cex.axis = .4) 
dev.off()

# map
cat("plotting map for K=3\n")
png(filename=paste0(out, "-map-pseudo_K2.png"), width = 5880, height = height, units = "px", res = 1200)
par(mar = c(0.2, 0.2, 0.2, 0.2))
plot(q.matrix, coord, method="map.max", cex = 1, 
     raster.filename = raster, interpol = FieldsKrigModel(10), 
     resolution = c(600,600), xaxt='n', ann=FALSE, yaxt='n',
     col.palette = my.palette, 
     window = winsize)
dev.off()

# get main cluster for clustering matrices
library("tidyr")
tmp <- cbind(bam, q.matrix)
tmp$main <- colnames(tmp)[apply(tmp, 1, which.max)]
tmp$mainL <- chartr("123", "ABC", tmp$main)
tmp$mainL <- gsub("^", paste0("X_", sp, "_"), tmp$mainL)
tmp <- unite(tmp, V1, mainL, col = "ID", sep = ";", )
tmp$ID
