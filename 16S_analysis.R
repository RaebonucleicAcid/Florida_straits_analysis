###ALL CODE FOR PROKARYOTES FIRST, THEN EUKARYOTES, THEN CHLOROPLASTS!

rm(list=ls())

##load packages
library("ggplot2");
library("vegan");
library("dplyr");
library("phyloseq")
library("data.table")

#set working directory
setwd('~/Desktop/R_stuff/R/FLORIDA_FILES/16S')

#confirm the necessary files are there
list.files()

##PROCESSING DATA:
#READ IN OTU + TAX + meta data dataframes & re-format 
##OTU table
otumat <- read.csv("OTU_table.tsv", header = TRUE, sep="\t");
row.names(otumat) <- otumat$OTU
otumat <- subset(otumat, select = 
                   -c(OTU, F028, F062, F103, F010, F012, F040, F051, F085));
otumat <- as.matrix(otumat)

#taxonomy table
taxmat <- read.csv("TAX_table.tsv", header = TRUE, sep = "\t");
row.names(taxmat) <- taxmat$OTU ;
taxmat <- subset(taxmat, select = -c(OTU));
taxmat <- as.matrix(taxmat)

#metadata, also adding gulf stream grouping 
metadata <- read.csv("Florida_subset_metadata.tsv", header = TRUE, sep="\t");
row.names(metadata) <- metadata$Sample_ID;
gulf_group <- c((rep(x="A", times=38)),(rep(x="B", times=11)))
gulf_group <- as.character(gulf_group)
metadata <- cbind(metadata, gulf_group)

##using PHYLOSEQ, it keeps everything organized 
##make phyloseq object
OTU <- otu_table(otumat, taxa_are_rows = TRUE);
TAX <- tax_table(taxmat);
SAMPLE_DATA <- sample_data(metadata)

physeq <- phyloseq(OTU, TAX, SAMPLE_DATA)
physeq
physeq_16S <- physeq
save(physeq_16S, file="physeq_object_FL_16S.rda")

head(OTU)
head(TAX)
head(SAMPLE_DATA)

#remove singletons + make separate phyloseq object with pruned OTU table
#so that both are available for subsequent analyses
pruned_otu <- prune_taxa(taxa_sums(OTU) > 1, OTU)
physeq_pruned <- phyloseq(pruned_otu, TAX, SAMPLE_DATA)
physeq_pruned

# Plot sequencing depth to make sure I removed all samples with low # of reads
seq_depth<- data.table(as(sample_data(physeq), "data.frame"),
                       TotalReads = sample_sums(physeq), keep.rownames = TRUE)
setnames(seq_depth, "rn", "SampleID")
plot_seq_depth <- ggplot(seq_depth, aes(TotalReads)) + geom_histogram() + 
  ggtitle("Sequencing Depth")
plot_seq_depth

#--------------------------------------------------------------------------------
##Alpha Diversity:
alpha <- estimate_richness(physeq, measures = 
                             c("Observed", "Shannon", "Simpson"))
plot_richness(physeq, measures = 
                c("Observed", "Shannon", "Simpson")) + geom_boxplot()
plot_richness(physeq, x="gulf_group", measures = 
                c("Observed", "Shannon", "Simpson")) + geom_boxplot()

?plot_richness
#------------------------------------------------------------------------------
##Exploring the data, abundance plots:
###FAMILY-LEVEL ABUNDANCE PLOTS
#family-level abundance
#transform to relative abundance
physeq_family <- transform_sample_counts(physeq, function(x) x / sum(x))
#transform to relative abundance
physeq_family <- filter_taxa(physeq_family, function(x) sum(x) > .12, TRUE) 
#filter out low abundance taxa for visualization purposes

physeq_family <- physeq_family %>%
  tax_glom(taxrank = "Family") %>%         # agglomerate at phylum level
  psmelt() %>%                             # Melt to long format
  arrange(Sample)                   # Sort data frame alphabetically by phylum   

##PLOTTING
# Set colors for plotting
family_colors <- c(
  "#CBD588", "#5F7FC7", "orange","#DA5724", "#508578", "#CD9BCD",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248", 
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861", "steelblue2","#618f5e",
  "#6e5e8f", "#eb102e", "#673770","#D14285", "#652926", "#C84248"
)

# Make abundance plot:
Prok_family_abundance_plot <- ggplot(data=physeq_family, 
                                     aes(x = Sample , y = Abundance, 
                                         fill = Family)) + 
  geom_bar(stat = "identity") +
  scale_fill_manual(values = family_colors) +
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +
  ylab("Relative Abundance ( > 5%) \n")

##Improve plot + throw it up there boom
Prok_family_abundance_plot <- Prok_family_abundance_plot +
  theme(legend.key.size = unit(0.5, "cm"),
        legend.title = element_text(size = 12, face = "bold"), 
        legend.text = element_text(size = 10)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

Prok_family_abundance_plot

head(otumat)

#------------------------------------------------------------------------------
###NMDS
otu.physeq <- otu_table(physeq)
otumat.transposed <- t(otumat)
bray <- vegdist(otumat.transposed, method = "bray")
distMat = as.dist(bray)

set.seed(23984)
NMDS1 = metaMDS(otumat.transposed, k = 2, trymax=999, autotransform = FALSE, 
                distance = "bray")
goodness(NMDS1)
stressplot(NMDS1)
summary(NMDS1)
str(NMDS1)

##PLOTTING NMDS
p1 <- plot_ordination(physeq, NMDS1, color = "gulf_group", label = "Sample_ID")
p1 + geom_point(size = 2) + theme_bw() + 
  scale_colour_brewer(type = "qual", palette = "Set1")

###conducting homogeneity of dispersal test + PERMANOVA
##see https://rpubs.com/DKCH2020/587758:
veganotu = function(physeq) {
  require("vegan")
  OTU = otu_table(physeq)
  if (taxa_are_rows(OTU)) {
    OTU = t(OTU)
  }
  return(as(OTU, "matrix"))
}

# export data from phyloseq to vegan-compatible object
physeq_vegan <- veganotu(physeq)

# make a data frame that can be used in vegan from the sample_data in the
# phyloseq object
sampledf <- data.frame(sample_data(physeq))
physeq_BRAY <- vegdist(wisconsin(sqrt(physeq_vegan)), method = "bray")

betadisp_physeq <- betadisper(physeq_BRAY, sampledf$gulf_group)
betadisp_physeq
par(mar = c(7, 4, 2, 2))
boxplot(betadisp_physeq, xlab = "", las = 2, cex.axis = 0.8)

#are the variances the same?
anova(betadisp_physeq)

###are the centroids of the clusters seen on the NMDS distinct?
sampledf <- data.frame(sample_data(physeq))
adonis(distMat ~ sampledf$gulf_group, data = sampledf)

#------------------------------------------------------------------------------
####DISTANCE DECAY FUNTION

#load a few more libraries 
library(FSA)
library(simba)

### function to plug in lat + lon to get geospatical distance back:

ReplaceLowerOrUpperTriangle <- function(m, triangle.to.replace){
  # If triangle.to.replace="lower", 
  #replaces the lower triangle of a square matrix with its upper triangle.
  # If triangle.to.replace="upper", 
  # replaces the upper triangle of a square matrix with its lower triangle.
  
  if (nrow(m) != ncol(m)) stop("Supplied matrix must be square.")
  if      (tolower(triangle.to.replace) == "lower") tri <- lower.tri(m)
  else if (tolower(triangle.to.replace) == "upper") tri <- upper.tri(m)
  else stop("triangle.to.replace must be set to 'lower' or 'upper'.")
  m[tri] <- t(m)[tri]
  return(m)
}

GeoDistanceInMetresMatrix <- function(df.geopoints){
  # Returns a matrix (M) of distances between geographic points.
  # M[i,j] = M[j,i] = Distance between (df.geopoints$lat[i], 
  # df.geopoints$lon[i]) and
  # (df.geopoints$lat[j], df.geopoints$lon[j]).
  # The row and column names are given by df.geopoints$name.
  
  GeoDistanceInMetres <- function(g1, g2){
    # Returns a vector of distances. (But if g1$index > g2$index, returns zero.)
    # The 1st value in the returned vector is the distance 
    #between g1[[1]] and g2[[1]].
    # The 2nd value in the returned vector is the distance between 
    #g1[[2]] and g2[[2]]. Etc.
    # Each g1[[x]] or g2[[x]] must be a list with 
    #named elements "index", "lat" and "lon".
    # E.g. g1 <- list(list("index"=1, "lat"=12.1, "lon"=10.1), 
    #list("index"=3, "lat"=12.1, "lon"=13.2))
    
    DistM <- function(g1, g2){
      require("Imap")
      return(ifelse(g1$index > g2$index, 0, gdist(lat.1=g1$lat, lon.1=g1$lon, 
                                                  lat.2=g2$lat, lon.2=g2$lon, 
                                                  units="m")))
    }
    return(mapply(DistM, g1, g2))
  }
  
  n.geopoints <- nrow(df.geopoints)
  
  # The index column is used to ensure we only do calculations for the 
  #upper triangle of points
  df.geopoints$index <- 1:n.geopoints
  
  # Create a list of lists
  list.geopoints <- by(df.geopoints[,c("index", "lat", "lon")], 
                       1:n.geopoints, function(x){return(list(x))})
  
  # Get a matrix of distances (in metres)
  mat.distances <- ReplaceLowerOrUpperTriangle(outer(list.geopoints, 
                                                     list.geopoints, 
                                                     GeoDistanceInMetres), 
                                               "lower")
  
  # Set the row and column names
  rownames(mat.distances) <- df.geopoints$name
  colnames(mat.distances) <- df.geopoints$name
  
  return(mat.distances)
}

#------------------------------------------------------------------------------
##making plot for distance decay
physeq_prop<- transform_sample_counts(physeq, function(x) x / sum(x))
metadata_pruned_physeq <-data.frame(sample_data(physeq_prop))
comm_pruned_physeq <- t((otu_table(physeq_prop)))
dist<-data.frame(row.names(comm_pruned_physeq))
geo_physeq <- data.frame("lat" = metadata_pruned_physeq$Latitude, 
                         "long" = metadata_pruned_physeq$Lontitude)

dist$lat<-geo_physeq$lat
dist$lon<-geo_physeq$long
colnames(dist)[colnames(dist)=="row.names.comm_pruned_physeq."] <- "name"

dist_m <-GeoDistanceInMetresMatrix(dist) 
#this does not round and keeps in meters
dist_m[upper.tri(dist_m)] <- NA

OTU1 = as(otu_table(physeq_prop), "matrix")
OTU1<-t(OTU1)
commdist <- vegdist(OTU1, method="bray",na.rm=TRUE,binary=FALSE)
#binary = TRUE makes it sorenson


###getting average and stdev of Bray Curtis values 
commdist_onerow <- as.matrix(commdist)
commdist_onerow
commdist_onerow <- commdist_onerow[1,]
commdist_onerow <- as.matrix(commdist_onerow)
commdist_onerow
commdist_onerow$dissimilarity <- 1-commdist_onerow
commdist_onerow
mean(commdist_onerow$dissimilarity)
sd(commdist_onerow$dissimilarity)


##this makes this a 3 column dataframe instead of matrix, 
#basically "melting" the data
comm.dist.ls <- liste(commdist, entry="comm")
comm.dist.ls
coord.dist.ls <- liste(dist_m, entry="dist")
coord.dist.ls <- na.omit(coord.dist.ls)
coord.dist.ls <- coord.dist.ls[!(coord.dist.ls$NBX==coord.dist.ls$NBY),]
coord.dist.ls <- as.data.frame(coord.dist.ls)

dim(as.data.frame(coord.dist.ls, rownames = FALSE))
colnames(coord.dist.ls)

###log model
model_log<-lm(log(1-comm.dist.ls$comm)~log(coord.dist.ls$dist))
summary(model_log)

### non-log model
model <- lm((1-comm.dist.ls$comm)~(coord.dist.ls$dist))
summary(model)

##plotting
df<-data.frame(sim=(1-comm.dist.ls$comm), dist= (coord.dist.ls$dist))

p<-ggplot(df, aes(x = dist, y =sim)) + geom_point()+
  geom_smooth(method = "lm", se=FALSE, color="red", aes(group=1))+
  labs(y = "Bray-Curtis Similarity")  +
  labs(x = "Spatial distance,m")

p<-p+ theme_bw() + theme(panel.border = element_blank(), 
                         panel.grid.major = element_blank(),
                         panel.grid.minor = element_blank(), 
                         axis.line = element_line(colour = "black"), 
                         axis.title.x = element_text(face="bold", size=16),
                         axis.title.y = element_text(face="bold", size=16))
p + ylim(0, 0.9) 


#------------------------------------------------------------------------------
###QUICK LITTLE SIDE EXPLORATION, DIDN'T MENTION IN PAPER:
##I TOOK SOME SAMPLES FROM A LAB MATE'S CRUISE WHICH CROSSED INTO THE GULF
# STREAM IN THE NORTHERN ATLANTIC TO SEE IF THOSE SAMPLES CLUSTERED WITH MY
# GULF STREAM SAMPLES:

####LOOKING AT GUULF GROUPS FROM JESSE'S CRUISE
##PROCESSING DATA:

#READ IN OTU + TAX + meta data dataframes & re-format 
##OTU table
otumat_AMT <- read.csv("OTU_table_FL_AMT.txt", header = TRUE, sep="\t");
otumat_AMT <-otumat_AMT[!duplicated(otumat_AMT$OTU), ]
row.names(otumat_AMT) <- otumat_AMT$OTU
otumat_AMT <- subset(otumat_AMT, select = 
                       -c(OTU, F028, F062, F103, F010, F012, F040, F051, F085, 
                          GA03_0134, GA03_0135, GA03_0052, GA03_0053));
otumat_AMT <- as.matrix(otumat_AMT)

#taxonomy table
taxmat_AMT <- read.csv("TAX_table_FL_AMT.txt", header = TRUE, sep = "\t");
dim(taxmat_AMT)
taxmat_AMT <- taxmat_AMT[!duplicated(taxmat_AMT$OTU), ]
row.names(taxmat_AMT) <- taxmat_AMT$OTU;
taxmat_AMT <- subset(taxmat_AMT, select = -c(OTU));
taxmat_AMT <- as.matrix(taxmat_AMT)

#metadata
metadata_AMT<- read.csv("Florida_AMT_subset_metadata.tsv", 
                        header = TRUE, sep="\t");
row.names(metadata_AMT) <- metadata_AMT$Sample_ID;

##using PHYLOSEQ to work through my exploratory stuff 
#because it keeps everything organized 

##make phyloseq object
OTU_AMT <- otu_table(otumat_AMT, taxa_are_rows = TRUE);
TAX_AMT <- tax_table(taxmat_AMT);
SAMPLE_DATA_AMT <- sample_data(metadata_AMT)

taxa_names(TAX_AMT)

physeq_AMT <- phyloseq(OTU_AMT, TAX_AMT, SAMPLE_DATA_AMT)
physeq_AMT

head(OTU)
head(TAX)
head(SAMPLE_DATA_AMT)

###NMDS FOR SIDE EXPLORATION
physeq_AMT_pruned = prune_taxa(taxa_sums(physeq_AMT) > 1, physeq_AMT)

otu.physeq.FL_AMT <- otu_table(physeq_AMT)

otumat_AMT.transposed <- t(otu.physeq.FL_AMT)
bray_AMT <- vegdist(otumat_AMT.transposed, method = "bray", na.rm = TRUE)
distMat_AMT = as.dist(bray_AMT)

set.seed(23984)
NMDS1_AMT = metaMDS(distMat_AMT, k = 2)
stressplot(NMDS1_AMT)
summary(NMDS1_AMT)
str(NMDS1_AMT)

##PLOTTING NMDS
p1 <- plot_ordination(physeq_AMT, NMDS1_AMT, 
                      color = "GROUP", label = "Sample_ID")
p1
p1 + geom_point(size = 2) + theme_bw() + 
  scale_colour_brewer(type = "qual", palette = "Set1")

####answer = nope!