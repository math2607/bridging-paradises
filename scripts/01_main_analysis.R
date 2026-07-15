###=###=###=###=###=###=###=###=###=###=###=###=###=###=###
# 1.1 Install and Load libraries ####
###=###=###=###=###=###=###=###=###=###=###=###=###=###=###
#=###=#=###=#=###=#=###=#=###=
# load packages
#=###=#=###=#=###=#=###=#=###=
# Data Manipulation and Strings
library("dplyr")
library("tidyr")
library("data.table")
library("reshape")
library("stringr")
library("stringi")
library("purrr")
library("rlang")
library("igraph")

# Phylogenetics and Evolution
library("ape")
library("phytools")
library("phangorn")
library("strap")
library("pegas")
library("splits")
library("delimtools")

# Bioinformatics and Sequencing
library("seqinr")
library("rentrez")
library("muscle")
library("msa")
library("Biostrings")
library("microseq")

# Biogeography and BioGeoBEARS
library("BioGeoBEARS")
library("cladoRcpp")
library("rexpokit")

# Spatial Data and Mapping
library("sf")
library("terra")
library("raster")
library("rnaturalearth")
library("rnaturalearthdata")
library("rgbif")
library("geodata")

# Plotting and Color Utilities
library("ggplot2")
library("colorRamps")
library("grDevices")
library("ggtree")
library("ggtreeExtra")
library("deeptime")
library("ggridges")
library("ggh4x")
library("patchwork")
library("ggspatial")

# Optimization, Statistics and Modeling
library("GenSA")
library("FD")
library("qpcR")
library("paran")
library("MASS")
library("MultinomialCI")

# Parallel Computing and Utilities
library("parallel")
library("snow")
library("devtools")
library("here")

# seed
set.seed(-236773-465642)  # replace seed!
#=###=#=###=#=###=#=###=#=###=
# load script
#=###=#=###=#=###=#=###=#=###=
source(here::here("00_functions.R"))

#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
# 2 GEOGRAPHIC INFORMATION ####
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=

# Searching squamata occurrences
# key <- name_backbone(name = "", class = "Squamata")$usageKey
# data <- occ_search(  
#   classKey = key,  
#   continent = c("north_america", "south_america"),
#   limit = 50000  
#   
# )  

# reading data downloaded
# ============================================================================
# NOTE:
# ============================================================================
# Run this line after running 00_download_gbif.R and retrieving the data from
# GBIF
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
data <- fread("../data/raw/gbif_occurrences.csv", sep = "\t")
names(data)

# Filter North and South America
sppdistrib <- data[,c("species", "decimalLongitude", 
                      "decimalLatitude","countryCode")]

# just a safe in case you lost the obj
safespdd <- sppdistrib

# check
nrow(sppdistrib)

# saving image

# Filter NA
sppdistrib <- sppdistrib[!is.na(sppdistrib$decimalLongitude) &
                           !is.na(sppdistrib$decimalLatitude), ]

# filtering species only in the land
vectA <- terra::vect("../data/geography/A.shp") # above
vectB <- terra::vect("../data/geography/B.shp") # on
vectC <- terra::vect("../data/geography/C.shp") # below
# rbinding vectors
vectALL <- rbind(vectA, vectB, vectC)

# convert to points
pts <- terra::vect(sppdistrib, geom=c("decimalLongitude", "decimalLatitude"), crs="EPSG:4326")

# divide into chuncks ro save memory
n_rows <- nrow(pts)
chunk_size <- 100000
n_chunks <- ceiling(n_rows / chunk_size)

# loop for testing interse (due memory)
results_list <- list()
for (i in 1:n_chunks) {
  # Define start and end indices
  start <- ((i - 1) * chunk_size) + 1
  end <- min(i * chunk_size, n_rows)
  
  # Slice the SpatVector
  pts_chunk <- pts[start:end, ]
  
  # Perform the intersection
  intersected_chunk <- terra::intersect(vectALL, pts_chunk)
  
  # Convert to dataframe immediately to save memory
  df_chunk <- as.data.frame(intersected_chunk, geom = "XY")
  
  # Keep only necessary columns to keep the list light
  cols_to_keep <- intersect(names(df_chunk), 
                            c("species", "x", "y", "countryCode"))
  results_list[[i]] <- df_chunk[, cols_to_keep]
  
  # print progress
  base::message(paste("Processed chunk", i, "of", n_chunks))
}
# combine
sppdistribCLEAN <- dplyr::bind_rows(results_list)
sppdistribCLEAN <- tibble::as_tibble(sppdistribCLEAN)

#check:
head(sppdistribCLEAN)
nrow(sppdistribCLEAN)

# check for all genera in the df
namesgenera <- str_extract(sppdistribCLEAN$species, "^[^ ]+")
namesgenera <- sort(unique(namesgenera)[!is.na(unique(namesgenera))])

world <- vectALL
worldPH <- terra::aggregate(world[3:11,])
world <- rbind(world[1:2,], worldPH)

# check for presence above, on and below isthmus
dfphTF <- data.frame("genera" = 1, "above" = 1, "on" = 1, "bellow" = 1)
dfphTF = dfphTF[!1,] # removing first row to make an empty df
for (j in 1:length(namesgenera)) {
  #grab all occurences for each genus
  dfphDIST <- sppdistribCLEAN[str_detect(sppdistribCLEAN$species, paste0("^", namesgenera[j], " ")),]
  dfphDIST <- terra::vect(dfphDIST, geom = c("x", "y")) # for intersecting using terra
  ints <- terra::relate(world, dfphDIST, "intersects") # check if occurs in each region
  dfdf <- c()
  for (i in 1:3){
    if (all(!ints[i,])){
      dfdf[i] <- 0 #  if not occuring, absence
    } else {
      dfdf[i] <- 1 # if occurs, presence
    }
  }
  
  # joing everything on the main df
  dfphTF[j,] <- c(namesgenera[j], dfdf)
  
  # print when done each gene
  base::cat(paste0("Remaining: ", nrow(dfphTF), "/", length(namesgenera), ". Last: ", namesgenera[j], "\n"))
  
}

# atrib motive and export df:
dfphTFFINAL <- data.frame("genera" = 1, "used" = 1, "reason" = 1)
dfphTFFINAL = dfphTFFINAL[!1,] # removing first row to make an empty df
for (i in 1:nrow(dfphTF)){
  # if only above:
  if (dfphTF[i,"above"] == 1 & dfphTF[i,"on"] == 0 & dfphTF[i,"bellow"] == 0){
    dfphTFFINAL[i,] <- c(dfphTF[i,1], "no", "only above")
  }
  # if only Isthmus:
  if (dfphTF[i,"above"] == 0 & dfphTF[i,"on"] == 1 & dfphTF[i,"bellow"] == 0){
    dfphTFFINAL[i,] <- c(dfphTF[i,1], "no", "only Isthmus")
  }
  # if only Below:
  if (dfphTF[i,"above"] == 0 & dfphTF[i,"on"] == 0 & dfphTF[i,"bellow"] == 1){
    dfphTFFINAL[i,] <- c(dfphTF[i,1], "no", "only below")
  }
  # if only above and Isthmus:
  if (dfphTF[i,"above"] == 1 & dfphTF[i,"on"] == 1 & dfphTF[i,"bellow"] == 0){
    dfphTFFINAL[i,] <- c(dfphTF[i,1], "no", "not below")
  }
  # if only above and below:
  if (dfphTF[i,"above"] == 1 & dfphTF[i,"on"] == 0 & dfphTF[i,"bellow"] == 1){
    dfphTFFINAL[i,] <- c(dfphTF[i,1], "no", "not on Isthmus")
  }
  # if only isthmus and below:
  if (dfphTF[i,"above"] == 0 & dfphTF[i,"on"] == 1 & dfphTF[i,"bellow"] == 1){
    dfphTFFINAL[i,] <- c(dfphTF[i,1], "no", "not above")
  }
  # if none:
  if (dfphTF[i,"above"] == 0 & dfphTF[i,"on"] == 0 & dfphTF[i,"bellow"] == 0){
    dfphTFFINAL[i,] <- c(dfphTF[i,1], "no", "not found/NA")
  }
  # if all (important):
  if (dfphTF[i,"above"] == 1 & dfphTF[i,"on"] == 1 & dfphTF[i,"bellow"] == 1){
    dfphTFFINAL[i,] <- c(dfphTF[i,1], "MAYBE", "CHECK")
  }
}

# write list of genera:
dir.create(here::here("..", "results", "supplementary"), recursive = TRUE, showWarnings = FALSE)
write.table(dfphTFFINAL, here::here("..","results", "supplementary", "S1_genera_screening.txt"),
            sep = ";", quote = F, row.names = F)


#plotting distributions (one-by-one, takes time)
worldsf <- st_as_sf(world)
worldsf$map_id <- c("Above Isthmus", "Isthmus", "Below Isthmus")

# genera which survived the filtering (above, on, below)
whowants <- c("Anomalepis", 
              "Amastridium", "Ameiva", "Amerotyphlops", "Anolis", "Atractus",
              "Basiliscus", "Boa", "Bothriechis", "Bothrops", "Brasiliscincus",
              "Chironius", "Clelia", "Cnemidophorus", "Coniophanes", "Conophis",
              "Corallus", "Corytophanes", "Crotalus", "Ctenosaura", "Dendrophidion",
              "Diploglossus", "Dipsas", "Drymarchon", "Drymobius", "Enuliophis",
              "Enulius", "Epicrates", "Epictia", "Erythrolamprus", "Geophis",
              "Gonatodes", "Gymnophthalmus", "Hemidactylus", "Holcosus", "Hydrophis",
              "Iguana", "Imantodes", "Indotyphlops", "Lachesis", "Lampropeltis",
              "Lepidoblepharis", "Lepidodactylus", "Leptodeira", "Leptophis", "Loxopholis",
              "Lygophis", "Mabuya", "Marisora", "Masticophis", "Mastigodryas",
              "Metlapilcoatlus", "Micrurus", "Ninia", "Nothopsis", "Oxybelis",
              "Oxyrhopus", "Phrynonax", "Phyllodactylus", "Pliocercus", "Polychrus",
              "Porthidium", "Rhadinaea", "Rhinobothryum", "Scaphiodontophis", "Sceloporus",
              "Sibon", "Sphaerodactylus", "Spilotes", "Stenorrhina", "Tantilla",
              "Tarentola", "Thamnophis", "Thecadactylus", "Tretanorhinus", "Tupinambis",
              "Urotheca", "Xenodon"
)

# check for n < 100
for (i in 1:length(whowants)){
  dfphDIST <- sppdistribCLEAN[str_detect(sppdistribCLEAN$species, paste0("^", whowants[i], " ")),]
  if (nrow(dfphDIST) < 100){
    base::message(paste0("Genus ", whowants[i], " has < 100 (n = ", nrow(dfphDIST), ")\n"))
  }else{
    base::cat(paste0("Genus ", whowants[i], " has > 100 (n = ", nrow(dfphDIST), ")\n"))
  }
}

# defining plot
whowants2 <- "Anomalepis" # <------ change this for each genus

# ploting and checking
ggplot() +
  geom_sf(data = worldsf, aes(fill = map_id)) +
  scale_fill_manual(values = c("grey85", "grey85", "red"))+
  coord_sf(
    xlim = c(-120, -40), # Longitude from -15 to 40
    ylim = c(-25, 55)   # Latitude from 35 to 65
  ) +
  geom_point(aes(x = sppdistribCLEAN[str_detect(sppdistribCLEAN$species, paste(whowants2, collapse = "|")),]$x,
                 y = sppdistribCLEAN[str_detect(sppdistribCLEAN$species, paste(whowants2, collapse = "|")),]$y,
                 colour = sppdistribCLEAN[str_detect(sppdistribCLEAN$species, whowants2),]$species)) +
  xlab("Longitude") + ylab("Lagitude") + theme_bw() +
  scale_color_discrete(name = "Species")


#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
# 3 GATHERING GENE INFORMATIONS ####
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=

#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
# 3.3 GENE RETRIEVAL ####
###=###=###=###=###=###=###=###=###=###=###=###=###=###=###
# title phylogeny #THIS IS JUST TO VIEW TIMES AND SISTER TAXA
# ============================================================================
# NOTE:
# ============================================================================
# Run this line after retrieving Title et al. (2024) phylogeny. See Main Document
# for more information. Save it in the data/raw folder with the name 
# main_titleetal_plain.tre
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
title <- read.tree("../data/raw/main_titleetal_plain.tre") # <----- don't forget to download first
# names of the genera you want to check
allsppused   <- "Anomalepis" # <----- change for the name of the genera you want
#names:
justgroup    <- allsppused
# grab node number for group of interest (see previous steps, this is MANUAL for better control)
node <- getMRCA(title, title$tip.label[str_detect(title$tip.label, paste0("^", justgroup, "_"))])
# check for null
node # if yes, maybe not in title, maybe only one species in title
# if, in Title, genus has only one species, use these lines instead:
# target_tip <- title$tip.label[str_detect(title$tip.label, paste0("^", justgroup, "_"))]
# tip_index <- which(title$tip.label == target_tip)
# node <- title$edge[which(title$edge[,2] == tip_index), 1]

# ploting only groups of interest
pnode <- Ancestors(title, node, type = "parent") # immediate ancestor
title2 <- extract.clade(title, node = pnode)
plotTree(title2) # check tree
# arrows to check
add.arrow(title2, tip = title2$tip.label[str_detect(title2$tip.label, paste0("^", justgroup, "_"))]) # group of interest
# genera in tree
unique(str_extract(title2$tip.label, ".*(?=_)"))[order(unique(str_extract(title2$tip.label, ".*(?=_)")))]

# NOTE: we decided to chose the main genera and one sister clade! folow next steps
allsppused <- c(allsppused, "Micruroiudes")           # choosing the sister-clade (CHANGE)
allsppused                                            # check
allsppusedph <- paste0(allsppused, collapse = "|")    # arranging for plotting
justsister   <- allsppused[2]                         # creating sister-clade object
add.arrow(title2, tip = title2$tip.label[str_detect(title2$tip.label, paste0("^", justsister, "_"))]) # group of interest
#nodeheight()

# FOR BEAST
nodelabels() # use the node label numbers here to check for dates for BEAST analysis if only one species in genus.
# all NOTE: these lines may not work if there is only one species for the genera used in Title
mrca_node <- ape::getMRCA(title2, title2$tip.label[str_detect(title2$tip.label, allsppusedph)])
round(nodeheight(title2, 1) - nodeheight(title2, mrca_node), 2)
# group of interest
mrca_node <- ape::getMRCA(title2, title2$tip.label[str_detect(title2$tip.label, justgroup)])
round(nodeheight(title2, 1) - nodeheight(title2, mrca_node), 2)
# sister-taxa
mrca_node <- ape::getMRCA(title2, title2$tip.label[str_detect(title2$tip.label, justsister)])
round(nodeheight(title2, 1) - nodeheight(title2, mrca_node), 2)

dev.off()

# unique(str_extract(extract.clade(title2, node = 1209)$tip.label, ".*(?=_)"))[order(unique(str_extract(extract.clade(title2, node = 1209)$tip.label, ".*(?=_)")))]
# plot(extract.clade(title2, node = 1209))
# nodeheight(extract.clade(title2, node = 1209), 1) - nodeheight(extract.clade(title2, node = 1209), 84)
# 
# tabletabletabletable <- do.call(rbind, ls2df[str_detect(ls2df, "^NA$") == FALSE])
# phph <- unique(tabletabletabletable$ORGANISM)[order(unique(tabletabletabletable$ORGANISM))]
# setdiff(gsub("_"," ",title2$tip.label), phph)

# setting genes for analysis:
# sppname = list()
# sppname[[1]]  =  c("Cnemidophorus", "Kentropyx")
# sppname[[2]]  =  c("Coniophanes", "Urotheca")
# sppname[[3]]  =  c("Bothrops asper", "Bothrops jararacussu")
# sppname[[4]]  =  c("Crotalus durissus", "Crotalus simus")
# sppname[[5]]  =  c("Dipsas", "Sibon", "Tropidodipsas")
# sppname[[6]]  =  c("Loxopholis", "Amapasaurus")
# sppname[[7]]  =  c("Marisora", "Varzea")
# sppname[[8]]  =  c("Micrurus", "Micruroides")
# sppname[[9]]  =  c("Oxybelis", "Leptophis")
# sppname[[10]] =  c("Sibon", "Tropidodipsas")
# sppname[[11]] = c("...", "...") ...
## NOTE: this is if you want to run as loop. We do not recommend it for better control
# You need to increase the list with all the genera if you want to retrieve molecular
# information in a loop
sppname = list()
sppname[[1]] = allsppused # define species names one by one
# without the loop, we will always rewrite the 1# element in list as our names vector
# deactivate if you want to run in loop. 

# check
sppname
# sppname2 = sppname

# input all names and variations into each object
genename = c()
gengen = list()
# cytb
gengen[[1]] = c("cytb", "cob", "cytochrome b", "CYB", "cyt-b", "cytochrome-b", "ubiquinol-cytochrome c reductase complex cytochrome b subunit", "MT-CYB")
genename = c(genename, gengen[[1]])
# # Cox1
gengen[[2]] = c("cox1", "COI", "CO1", "COX-1", "COXI", "cytochrome c oxidase subunit 1", "cytochrome c oxidase subunit I", "cytochrome oxidase subunit I")
genename = c(genename, gengen[[2]])
# nd2
gengen[[3]] = c("nd2", "NADH2", "ND-2", "NADH dehydrogenase subunit 2", "NADH dehydrogenase subunit II")
genename = c(genename, gengen[[3]])
# nd4
gengen[[4]] = c("nd4", "NADH4", "ND-4", "NADH dehydrogenase subunit 4", "NADH dehydrogenase subunit IV")
genename = c(genename, gengen[[4]])
#cmos
gengen[[5]] = c("c-mos", "mos", "oocyte maturation factor mos", "proto-oncogene serine/threonine-protein kinase mos", "v-mos") 
genename = c(genename, gengen[[5]])
#jun
gengen[[6]] = c("jun", "c-jun", "AP-1", "transcription factor AP-1", "jun proto-oncogene", "v-jun")
genename = c(genename, gengen[[6]])
# nt3
gengen[[7]] = c("nt3", "ntf3", "neurotrophin-3", "neurotrophic factor 3", "NT-3", "HDNF (nerve growth factor-2)")
genename = c(genename, gengen[[7]])
# rag1
gengen[[8]] = c("rag1", "recombination activating 1", "recombination activating gene 1", "RAG-1")
genename = c(genename, gengen[[8]])

# see all genes
genename 

# main loop for all species
#for (lol in 1:length(sppname)){ # THIS IS THE MAIN LOOP, reactivate this to run in loop.
# NOTE : if you desire to run in loop, you will need to reactivate line 893 as well

# deactivate this if you are running in loop
lol = 1 # we will always grab the first element in list without the loop

# Run searches for all genes
results <- search_by_gene(genes = genename, species = sppname[[lol]])

# Process all genes and combine results
all_results <- process_gene_results(gene_result =  results)
base::message("Retrieving genes completed")

# removing NAs
results_df <- all_results[!is.na(all_results$GENE),]
unique(results_df$GENE)

# renaming all genes to simplified version
for (la in 1:length(gengen)){
  
  #rename for first occurrence in vector (i.e. simplier name)
  phphph <- str_detect(results_df$GENE, regex(paste(gengen[[la]], collapse = "|"), ignore_case = T))
  if (identical(results_df[phphph,]$GENE, character(0))){
    next
  }else{
    results_df[phphph,]$GENE <- gengen[[la]][1]
  }
  
}
base::message("Rename completed")

# excluding genes that where not renamed (possibly errors)
uniGENE <- sapply(gengen, "[[", 1)
results_df <- results_df[str_detect(results_df$GENE, regex(paste(uniGENE, collapse = "|"), ignore_case = T)),]

# removing sequences without vouchers (NA)
results_df <- results_df[!is.na(results_df$VOUCHER),]

# check
head(results_df)
unique(results_df$GENE)

# removing West Indies occurrences from group of interest
west_indies_and_islands <- c(
  "Cuba", "Jamaica", "Haiti", "Dominican_Republic", "Puerto_Rico",
  "Antigua_and_Barbuda", "Saint_Kitts_and_Nevis", "Anguilla", 
  "Montserrat", "British_Virgin_Islands", "United_States_Virgin_Islands",
  "Virgin_Islands", "Dominica", "Saint_Lucia", "Saint_Vincent_and_the_Grenadines", 
  "Grenada", "Barbados", "Trinidad_and_Tobago", "Bahamas", "Turks_and_Caicos_Islands",
  "Guadeloupe", "Martinique", "Saint_Martin", "Sint_Maarten", "Curacao", "Aruba", 
  "Bonaire", "Saba", "Cayman", "Hawaii", "Palau", "Northern_Mariana", "Guam"
)

islands_space <- gsub("_", " ", west_indies_and_islands)
islands_all <- unique(c(west_indies_and_islands, islands_space))

islands_regex <- regex(paste0("\\b(", paste(islands_all, collapse = "|"), ")\\b"), ignore_case = TRUE)

results_df_ph <- results_df[!str_detect(results_df$GEO_LOC, islands_regex), ]

results_df_1 <- rbind(
  results_df_ph[str_detect(results_df_ph$ORGANISM, sppname[[lol]][1]), ],
  results_df[str_detect(results_df$ORGANISM, paste0(sppname[[lol]][-1], collapse = "|")), ]
)

# separate by voucher
results_split <- split(results_df_1, results_df_1$VOUCHER)

# should be equal
length(results_split)
length(unique(results_df$VOUCHER[results_df$VOUCHER %in% NA == FALSE]))

# grab genes by voucher
lsdf <- data.frame()
for (i in 1:length(results_split)){
  dfph <- results_split[[i]]
  for (j in 1:length(gengen)){
    if (nrow(dfph[str_detect(dfph$GENE, regex(gengen[[j]][1], ignore_case = T)),]) == 0){
      dfpheach <- data.frame(
        ORGANISM    = results_split[[i]]$ORGANISM[1], 
        VOUCHER     = results_split[[i]]$VOUCHER[1],
        GENE        = gengen[[j]][1],
        EXIST       = "NO",
        GEO_LOC     = NA,
        ACCESSION   = NA,
        TITLE       = NA)
    } else {
      dfpheach <- dfph[str_detect(dfph$GENE, regex(gengen[[j]][1], ignore_case = T)),]
      dfpheach <- cbind(dfpheach, EXIST = "YES")
      dfpheach <- dfpheach[,c(1,2,3,7,4,5,6)]
    }
    if (j == 1){
      dfphFINAL <- dfpheach
    } else {
      dfphFINAL <- rbind(dfphFINAL, dfpheach)
    }
  }
  if (i == 1){
    lsdf <- dfphFINAL
  } else {
    lsdf <- rbind(lsdf, dfphFINAL)
  }
}
base::message("lsdf completed")

head(lsdf)
nrow(lsdf)

# remove na from table (no organism)
lsdf <- lsdf[!is.na(lsdf$ORGANISM),]

# To attribute the same region for same vouchers without location in some genes.
# Just if authors submited different genes for same vouchers but forgot to 
# include the geo_loc in some of them
for (vch in unique(lsdf$VOUCHER)){
  
  #make placeholder table with matching  
  ph <- lsdf[lsdf$VOUCHER == vch, ]
  # ph <- ph[!is.na(ph$ORGANISM),]
  
  # check if ph has only NA in GEO_LOC and skip if true
  if (length(is.na(unique(ph$GEO_LOC))) == 1) next
  
  # check if more than one species have the same voucher (probably wrong, remove)
  if (any(duplicated(ph$GENE))){
    # grabing the line duplicated
    toremove <- ph[duplicated(ph$GENE),]
    ph <- ph[!duplicated(ph$GENE), ]
    
    # remove row from table
    lsdf <- anti_join(lsdf, toremove)
    
  }
  
  
  # matching the same location to vouchers with missing locations
  if (identical(unique(ph$GEO_LOC)[!is.na(unique(ph$GEO_LOC))], character(0))){
    ph$GEO_LOC <- ph$GEO_LOC #basically does nothing 
  }else{
    ph$GEO_LOC = unique(ph$GEO_LOC)[!is.na(unique(ph$GEO_LOC))]
  }
  lsdf[lsdf$VOUCHER == vch, ] <- ph
  
}
base::message("Removing duplicates and attrb same region completed")

#check
head(lsdf)
nrow(lsdf)

# cleaning voucher (for later)
lsdf$VOUCHER <- gsub("[^[:alnum:]]+", "_", lsdf$VOUCHER)
lsdf$VOUCHER <- gsub("^_|_$", "", lsdf$VOUCHER) # Removes leading or trailing underscores

# separating by gene
ls2df <- list()
for (i in 1:length(gengen)){
  genetxtph <- lsdf[str_detect(lsdf$GENE, regex(gengen[[i]][1], ignore_case = T)),]
  genetxtph <- genetxtph[genetxtph$EXIST %in% "NO" == FALSE,]
  genetxtph <- genetxtph[is.na(genetxtph$EXIST) == FALSE,]# take the ones that dont have place
  if (nrow(genetxtph) == 0){
    ls2df[[gengen[[i]][1]]] <- "NA"
  } else {
    genetxtph2 <- genetxtph[genetxtph$EXIST %in% "NO" == FALSE,]
    ls2df[[gengen[[i]][1]]] <- genetxtph2 
    # OLD 
    #rbind(genetxtph2, data.frame(
    #ORGANISM    = "chicken", 
    #VOUCHER     = "chicken",
    #GENE        = genename[i],
    #EXIST       = "YES",
    #GEO_LOC     = NA,
    #ACCESSION   = galus[i],
    #TITLE       = NA))
  }
}

length(ls2df)
base::message("ld2df completed")

# removing genes with no location (NA)
for (l in 1:length(ls2df)){
  ph <- ls2df[[l]]
  # skip NA 
  if (!is.data.frame(ph)){ next }
  ph <- ph[!is.na(ph$GEO_LOC),]
  if (nrow(ph) == 0){
    ls2df[[l]] = "NA"
  } else {
    ls2df[[l]] <- ph
  }
}

# IMPORTANT:
# check for only NAs here, if ONLY nas that means NO specimen have location, REMOVE
all(is.na(ls2df) | ls2df == "NA")

# making a df for all accession numbers, locations, etc
dfwithloc <- data.frame()
for (l in 1:length(ls2df)){
  group_name <- names(ls2df)[l]
  ph <- ls2df[[l]]
  if (!is.data.frame(ph)){if (ph == "NA") {next}}
  dfwithloc <- rbind(dfwithloc, ph)
}

# adding genus
phgenus   <- str_extract(dfwithloc$ORGANISM, "^[^ ]+")
dfwithloc <- cbind(dfwithloc, "GENUS" = phgenus)

# cleaning the table
dfwithloc <- dfwithloc[,c("GENE", "ACCESSION", "GENUS", "ORGANISM", "VOUCHER", "GEO_LOC", "TITLE")]
dfwithloc$GEO_LOC <- rename.locations(dfwithloc$GEO_LOC)
# removing ; from everything to save it as sep in table
for (l in 1:ncol(dfwithloc)){
  dfwithloc[,l] <- gsub(";", "_", dfwithloc[,l])
}

# just checking
head(dfwithloc)

# we will save the table later

# grabbing the genes from genbank by accession number
ls3df <- list()
for (i in seq_along(ls2df)) {
  x <- ls2df[[i]]
  group_name <- gengen[[i]][1]
  
  if (!is.data.frame(x) || nrow(x) == 0) {
    ls3df[[group_name]] <- NA 
  } else {
    # 1. Download all sequences at once (much faster than a loop)
    # read.GenBank can take a vector of accessions
    phgengen <- get_genbank_safe(x$ACCESSION, chunk_size = 100, genename = group_name)
    
    # 2. Assign the VOUCHER names directly
    names(phgengen) <- x$VOUCHER
    
    # 3. Store in your results list
    ls3df[[group_name]] <- phgengen
  }
}
base::message("ls3df completed")


# renaming directiories
sppnamePH <- sppname[[lol]][1]
sppnamePH <- word(sppnamePH, 1)
dir.create("../r_sequences", recursive = TRUE, showWarnings = FALSE)
dir.create(paste0("../r_sequences/", sppnamePH), recursive = TRUE, showWarnings = FALSE)

# write.fasta
for (i in 1:length(gengen)){
  write.dna(ls3df[[i]], paste0("../r_sequences/", sppnamePH, "/", sppnamePH, "_", gengen[[i]][1], ".fasta"), format = "fasta") 
}

fastafiles <- list.files(paste0("../r_sequences/", sppnamePH, "/"))
fastafiles
# cleaning fasta files to only fasta
fastafiles <- fastafiles[str_detect(fastafiles, ".fasta$")]

dir.create(paste0("../r_sequences/", sppnamePH, "/align"), recursive = TRUE, showWarnings = FALSE)
# fastafiles <- fastafiles[-1]

a = list()
for (i in 1:length(fastafiles)){
  a[[fastafiles[[i]]]] <- readBStringSet(paste0("../r_sequences/", sppnamePH, "/", fastafiles[i]),format="fasta") 
}

a

# write location table
write.table(dfwithloc, paste0("../r_sequences/", sppnamePH, "/1dfwithloc_", sppnamePH, ".txt"), sep = ";", quote = F, row.names = F)



#==#==#==#==#==##==#==#==#==#==##==#==#==#==#==##==#==#==#==#==#
## ## ## VIEW, ALIGN AND SAVE ####
#==#==#==#==#==##==#==#==#==#==##==#==#==#==#==##==#==#==#==#==#
# 2. Process all elements in the list
processed_list <- process_sequences(a) # this can take a while
base::message("Processing sequequences completed")

# Remove NULL elements (those with no valid sequences)
processed_list <- processed_list[!sapply(processed_list, is.null)]

dir.create(paste0("../r_sequences/", sppnamePH, "/align/genes"), recursive = TRUE, showWarnings = FALSE)
dir.create(paste0("../r_sequences/", sppnamePH, "/align/genes/correct"), recursive = TRUE, showWarnings = FALSE)

# write.fasta
for (gg in names(processed_list)){
  write.dna(processed_list[[gg]], paste0("../r_sequences/", sppnamePH, "/align/genes/", gg), format = "fasta") 
}

# changing names to chicken
#for (i in 1:length(processed_list)){
#  names(processed_list[[i]])[names(processed_list[[i]]) %in% galus] <- "chicken"
#}

#renaming for BEAST and testing for each gene (for possible astral)
melhorModelo = NULL
for (gg in names(processed_list)){
  
  # Optional: View the resulting alignment
  filefinal <- readDNAStringSet(paste0("../r_sequences/", sppnamePH, "/align/genes/", gg), format = "fasta")
  
  # OLD FILTERING, ignore
  # ph = data.frame()
  # # changing names
  # for (i in 1:length(ls2df)){
  #   if(str_detect(ls2df[[i]], "NA")[1] == TRUE){
  #     next
  #   }
  #   ph <- rbind(ph, ls2df[[i]])
  # }
  # 
  # ph <- ph[,c(1,2,5)]
  # names <- apply(ph, 1, paste0, collapse="|")
  # 
  # tiplabelfinal <- names(filefinal)
  # for (i in 1:length(tiplabelfinal)){
  #   escaped_pattern <- str_replace_all(
  #     tiplabelfinal[i], 
  #     "([.\\\\+*?\\[\\]\\{\\}\\(\\)\\^\\$\\|])", 
  #     "\\\\\\1"
  #   )
  #   newname <- unique(names[str_detect(names, paste0("\\|", fixed(escaped_pattern), "\\|"))])[1]
  #   cat(paste0("int n: ", i, " ||||"), tiplabelfinal[i], "||||", newname, "\n")
  #   tiplabelfinal[i] <- newname
  #   Sys.sleep(0.01)
  # }
  # 
  # # Clean names (if needed)
  # clean_names <- remove_colon_after_last_pipe(tiplabelfinal)
  # clean_names <- gsub("\\s+", "_", clean_names)  # Replace spaces with underscores
  # #clean_names <- gsub("[^a-zA-Z0-9|_]", "", clean_names)
  # 
  # clean_names[order(clean_names)]
  # all_locations <- all_locations_clean(clean_names)
  # 
  # #checkduplicates
  # dup1 <- names(filefinal)[duplicated(retrieve_names_wtlocation(names(filefinal)))]
  # dup2 <- tiplabelfinal[duplicated(retrieve_names_wtlocation(tiplabelfinal))]
  # dup3 <- clean_names[duplicated(retrieve_names_wtlocation(clean_names))]
  # dup1
  # dup2
  # dup3
  # 
  # match(dup1, names(filefinal))
  # match(dup2, tiplabelfinal)
  # match(dup3, clean_names)
  # 
  # #filefinal[c(868)]
  # 
  # all_locations
  # unique(all_locations)[order(unique(all_locations))]
  #
  #clean_names[!is.na(clean_names)]
  #names(filefinal)[c(128, 112, 113)]
  # clean_names <- rename.tips.align(clean_names)
  # names(filefinal) <- clean_names
  # 
  #names(filefinal) <- gsub(" ", "_", names(filefinal))
  # names(filefinal[grepl("^-+$", as.character(filefinal))])
  #filefinal <- filefinal[!grepl("^-+$", as.character(filefinal))]
  #names(filefinal)[c(128, 112, 113)]
  #filefinal <- filefinal[-c(128, 112, 113)]
  # filefinal <- filefinal[!is.na(names(filefinal))]
  # 
  # unique(str_extract(names(filefinal), ".*?(?=_)"))
  # 
  # names(filefinal)
  
  # remove duplicates
  names_dup <- names(which(table(base::labels(filefinal)) > 1))
  if (length(names_dup) > 0){
    pos <- which(names(filefinal) %in% names_dup)
    
    occurrence_index <- ave(
      seq_along(names(filefinal)[pos]), # Sequence of indices
      names(filefinal)[pos],            # Grouping variable (the string itself)
      FUN = seq_along         # Function applied to each group: return 1, 2, 3, ...
    )
    
    result_vector <- ifelse(
      occurrence_index == 2,
      sub("_(?!.*_)", "_2_", names(filefinal)[pos], perl = TRUE),
      names(filefinal)[pos]
    )
    
    names(filefinal)[pos] <- result_vector
    
  }
  
  # rename for convention SPECIES|VOUCHER
  jobname = c()
  for (nmb in 1:length(names(filefinal))){
    # grab line matching voucher ID from df
    phfromdfloc <- dfwithloc[str_detect(dfwithloc$VOUCHER, names(filefinal)[nmb]),]
    nameconv <- unique(phfromdfloc$ORGANISM)
    nameconv <- gsub(" ", "_", nameconv)
    # check for different species name, grab first
    if (length(nameconv) > 1){
      nameconv <- nameconv[1]
    }
    
    #changing name
    nameconv <- paste(nameconv, names(filefinal)[nmb], sep = "|")
    jobname[nmb] <- nameconv
    
  }
  # finally change
  names(filefinal) <- jobname
  
  #write
  dir.create(paste0("../r_sequences/", sppnamePH, "/align/genes/correct/final"), recursive = TRUE, showWarnings = FALSE)
  write.dna(filefinal, paste0("../r_sequences/", sppnamePH, "/align/genes/correct/", gg), format = "fasta") 
  write.dna(filefinal, paste0("../r_sequences/", sppnamePH, "/align/genes/correct/final/", gg), format = "fasta") 
  
  # write_phy
  # write.phyDat(msa_trimmed, paste0("../r_sequences/", sppnamePH, "/align/", sppnamePH, "_merged_trimmed_sequences.phy"), format = "phylip")
  
  # printing aligment
  final_dnabin <- readDNAStringSet(paste0("../r_sequences/", sppnamePH, "/align/genes/correct/final/", gg))
  final_dnabin <- as.DNAbin(final_dnabin)
  pdf(paste0("../r_sequences/", sppnamePH, "/", gg, ".pdf"), height = 20)
  par(mar = c(4,15,3,4))
  print(image.DNAbin(final_dnabin, legend = F , main = gg))
  dev.off() 
  base::message("Aligment completed")
  # length(unique(Biostrings::width(filefinal)))
  
  filefinal <- as.phyDat(final_dnabin) 
  
  # test
  if (length(filefinal) >= 3){
    tryCatch({
      mt <- modelTest(filefinal, model = c("JC", "F81", "HKY", "GTR"))
      
      # Get the best model based on AIC
      best_idx <- which.min(mt$AIC)
      melhorModeloPH <- as.character(mt$Model[best_idx])
      
      phphph <- data.frame(gene = gg, modelo = melhorModeloPH)
      melhorModelo <- rbind(melhorModelo, phphph)
    }, error = function(e) {
      # If the model fails, record it as NA instead of crashing the whole loop
      base::message(paste("Error in gene:", gg, "-", e$message))
      phphph <- data.frame(gene = gg, modelo = "FAILED")
      melhorModelo <<- rbind(melhorModelo, phphph)
    })
  }else{
    phphph <- data.frame(gene = gg, modelo = "NA")
    melhorModelo <- rbind(melhorModelo, phphph)
  }
  
}

# test result
print(melhorModelo)
base::message("Mt test completed")

# changing the list with the correct names
processed_list_final = list()
for (i in 1:length(processed_list)){
  processed_list_final[[i]] = readDNAStringSet(paste0("../r_sequences/", sppnamePH, "/align/genes/correct/final/", names(processed_list)[[i]]), format = "fasta")
  names(processed_list_final)[[i]] <- names(processed_list)[[i]]
}

# Extract all names from all sequences
all_names <- unique(unlist(lapply(processed_list_final, names)))

# Merge alignments (now with enforced equal lengths)
final_alignment <- merge_alignments(processed_list_final)

# Verify lengths are equal
stopifnot(length(unique(BiocGenerics::width(final_alignment))) == 1)

# write fasta
dnabin_sequences <- as.DNAbin(final_alignment)
write.FASTA(dnabin_sequences, paste0("../r_sequences/", sppnamePH, "/align/", sppnamePH, "_merged_sequences.fasta"))

# Read as text first
msa_file <- readFasta(paste0("../r_sequences/", sppnamePH, "/align/", sppnamePH, "_merged_sequences.fasta"))

msa_fileDNABIN <- readDNAStringSet(paste0("../r_sequences/", sppnamePH, "/align/", sppnamePH, "_merged_sequences.fasta"))
msa_fileDNABIN <- as.DNAbin(msa_fileDNABIN)

pdf(paste0("../r_sequences/", sppnamePH, "/", sppnamePH, "_alignBEFOREtrim.pdf"))
par(mar = c(4,15,3,4))
print(image.DNAbin(msa_fileDNABIN))
dev.off()

# n specimens
length(msa_fileDNABIN)

# removing rows with - in more than 50% of the rows
fastatrimmedPH <- readFasta(paste0("../r_sequences/", sppnamePH, "/align/", sppnamePH, "_merged_sequences.fasta"))

# # define gap end and gap mid parameters:
# gapendp <- 0.95
# gapmidp <- 0.95
# if (nrow(fastatrimmedPH) == 1){
#   fastatrimmedPH <- fastatrimmedPH
# } else {
#   fastatrimmedPH <- msaTrim(fastatrimmedPH, gap.end = gapendp, gap.mid = gapmidp)
# }
# base::message("Trim completed")

df <- as.data.frame(do.call(rbind, strsplit(fastatrimmedPH$Sequence, "")))
row.names(df) <- fastatrimmedPH$Header

# # test and remove rows with <50% but with columns also <50%
# # calculate - col
col_occupancy <- colMeans(df != "-") >= 0.5          # check col with >= 50%
dfPH <- df[,col_occupancy]                           # select cols with >= 50%
dfPH <- dfPH[rowSums(dfPH == "-") != ncol(dfPH), ]   # remove lines which failed

# transform it in vector with TRUE
nrow(dfPH)
dfPH <- row.names(dfPH)
phphphph <- rep(TRUE, length(dfPH))
names(phphphph) <- dfPH
length(phphphph)
col_tokeep <- phphphph

df_tokeep = c()
for (i in 1:nrow(df)){
  # calculate - row
  row_occupancy <- mean(df[i,] != "-")
  if (row_occupancy >= 0.30){ #check for row >30%
    df_tokeep[i] <- TRUE
    names(df_tokeep[i]) <- row.names(df[i,])
  } else { #check for row <50%
    if(!is.na(col_tokeep[row.names(df[i,])])){ #chech if, even tho that row is <30% you should maintain because column
      df_tokeep[i] <- TRUE
      names(df_tokeep[i]) <- row.names(df[i,])
    } else { #if not, then remove.
      df_tokeep[i] <- FALSE
      names(df_tokeep[i]) <- row.names(df[i,])
    }
  }
}

names(df_tokeep) <- row.names(df)
df <- cbind(keep = df_tokeep, df)
#selecting only rows to keep
df <- df[df$keep == TRUE, c(2:ncol(df))]
#remove black columns
df <- df[,colMeans(df != "-") >= 0.1]
# remove lines with only -
df <- df[!rowMeans(df == "-") == 1,]
# remove remaining lines with < 50%
#df <- df[!rowMeans(df == "-") > 0.5,]
# Extract the final result
sequences_collapsed <- apply(df, 1, paste, collapse = "")
names(sequences_collapsed) <- rownames(df)

# converting df into stringset
dnabintrimmed0.50 <- DNAStringSet(sequences_collapsed)

write.dna(dnabintrimmed0.50, paste0("../r_sequences/", sppnamePH, "/align/", sppnamePH, "_merged_trimmed_sequences_0.50.fasta"), format = "fasta")
# 
dnabintrimmed0.50 <- readFasta(paste0("../r_sequences/", sppnamePH, "/align/", sppnamePH, "_merged_trimmed_sequences_0.50.fasta"))

# trimming
gapendp <- 0.9 # <--------- change to your needs 
gapmidp <- 0.9 # <--------- change to your needs 
if (nrow(dnabintrimmed0.50) == 1){
  dnabintrimmed0.50 <- dnabintrimmed0.50
} else {
  dnabintrimmed0.50 <- msaTrim(dnabintrimmed0.50, gap.end = gapendp, gap.mid = gapmidp)
}
base::message("Trim completed")

writeFasta(dnabintrimmed0.50, paste0("../r_sequences/", sppnamePH, "/align/", sppnamePH, "_merged_trimmed_sequences_0.50.fasta"))

dnabintrimmed0.50 <- readDNAStringSet(paste0("../r_sequences/", sppnamePH, "/align/", sppnamePH, "_merged_trimmed_sequences_0.50.fasta"))

dnabintrimmed0.50 <- as.DNAbin(dnabintrimmed0.50)

pdf(paste0("../r_sequences/", sppnamePH, "/", sppnamePH, "_alignAFTERtrimAFTER0.50.pdf"))
par(mar = c(4,15,3,4))
print(image.DNAbin(dnabintrimmed0.50, legend = F , main = paste0("ALINED and TRIMMED for: ", sppnamePH)))
dev.off()


#n specimens
length(dnabintrimmed0.50)
# all species names:
jub <- sort(unique(str_extract(names(dnabintrimmed0.50), "^[^\\|]+")))
jub
# all genera
unique(str_extract(jub, "^[^_]+"))

filefinal <- readDNAStringSet(paste0("../r_sequences/", sppnamePH, "/align/", sppnamePH, "_merged_trimmed_sequences_0.50.fasta"))
filefinal <- as.phyDat(filefinal)

# test all
if (length(filefinal) >= 3){
  mt <- modelTest(filefinal, model = c("JC", "F81", "HKY", "GTR"))
  mt # veja o resultado
  
  subset(mt$Model,mt$AIC==min(mt$AIC)) # selecao do melhor modelo segundo o criterico de informacao de Akaike (AIC)
  melhorModeloPH = mt$Model[which.min(mt$AIC)] # outra forma de ver a mesma coisa
  phphph <- data.frame(gene = "ALL", modelo = melhorModeloPH)
  melhorModelo <- rbind(melhorModelo, phphph)
}else{
  phphph <- data.frame(gene = "ALL", modelo = "NA")
  melhorModelo <- rbind(melhorModelo, phphph)
}

print(melhorModelo)

write.table(melhorModelo, paste0("../r_sequences/", sppnamePH, "/2bestmodel_", sppnamePH, ".txt"), sep = ";")

# rewrite table with genbank code information
vchnamesff <- str_extract(names(filefinal), "(?<=\\|).*")
dfwithlocCLEAN<- dfwithloc[str_detect(dfwithloc$VOUCHER, paste0(vchnamesff, collapse = "|")), ]

# check for locations:
cat("\n\nYou have", 
    length(unique(dfwithlocCLEAN$VOUCHER[str_detect(dfwithlocCLEAN$ORGANISM, paste0("^", sppname[[lol]][1]))])), 
    "specimens from group of interest, (", sppnamePH, "), and",
    length(unique(dfwithlocCLEAN$VOUCHER[str_detect(dfwithlocCLEAN$ORGANISM, paste0("^", sppname[[lol]][2]))])),
    "from the sister-taxa, with", 
    length(unique(dfwithlocCLEAN$GEO_LOC)[unique(dfwithlocCLEAN$GEO_LOC) != "NA"]),
    "unique locations (including sister-taxa), them being:\n\n", paste(
      sort(unique(dfwithlocCLEAN$GEO_LOC[str_detect(dfwithlocCLEAN$ORGANISM, paste0("^", sppname[[lol]][1]))])), collapse = ", \n"),
    "\n\nfrom the main group, but also:\n\n",
    paste(
      sort(unique(dfwithlocCLEAN$GEO_LOC[str_detect(dfwithlocCLEAN$ORGANISM, paste0("^", sppname[[lol]][2]))])), collapse = ", \n"),
    "\n\nfrom the sister taxa.")
if (length(dfwithlocCLEAN$ORGANISM[str_detect(dfwithlocCLEAN$ORGANISM, paste0("^", sppname[[lol]][1]))]) < 30){
  base::message("Warning: You have less than 30 specimens from the group of interest. CAUTION!!!")
}

cat("\n\nYou have", 
    length(unique(dfwithlocCLEAN$VOUCHER[str_detect(dfwithlocCLEAN$ORGANISM, paste0("^", sppname[[lol]][1]))])), 
    "specimens from group of interest, (", sppnamePH, "), and",
    length(unique(dfwithlocCLEAN$VOUCHER[str_detect(dfwithlocCLEAN$ORGANISM, paste0("^", sppname[[lol]][2]))])),
    "from the sister-taxa, with", 
    length(unique(dfwithlocCLEAN$GEO_LOC)[unique(dfwithlocCLEAN$GEO_LOC) != "NA"]),
    "unique locations (including sister-taxa). Also, you have the number of genes:\n\n",
    "cytb:", nrow(dfwithlocCLEAN[str_detect(dfwithlocCLEAN$GENUS, sppname[[lol]][1]) & str_detect(dfwithlocCLEAN$GENE, "cytb"), ]), "\n",
    "cox1:", nrow(dfwithlocCLEAN[str_detect(dfwithlocCLEAN$GENUS, sppname[[lol]][1]) & str_detect(dfwithlocCLEAN$GENE, "cox1"), ]), "\n", 
    "nd2:",	 nrow(dfwithlocCLEAN[str_detect(dfwithlocCLEAN$GENUS, sppname[[lol]][1]) & str_detect(dfwithlocCLEAN$GENE, "nd2"), ]), "\n",
    "nd4:",	 nrow(dfwithlocCLEAN[str_detect(dfwithlocCLEAN$GENUS, sppname[[lol]][1]) & str_detect(dfwithlocCLEAN$GENE, "nd4"), ]), "\n",
    "c-mos:",nrow(dfwithlocCLEAN[str_detect(dfwithlocCLEAN$GENUS, sppname[[lol]][1]) & str_detect(dfwithlocCLEAN$GENE, "c-mos"), ]), "\n",	
    "jun:",	 nrow(dfwithlocCLEAN[str_detect(dfwithlocCLEAN$GENUS, sppname[[lol]][1]) & str_detect(dfwithlocCLEAN$GENE, "jun"), ]), "\n",
    "nt3:",	 nrow(dfwithlocCLEAN[str_detect(dfwithlocCLEAN$GENUS, sppname[[lol]][1]) & str_detect(dfwithlocCLEAN$GENE, "nt3"), ]), "\n",
    "rag1:", nrow(dfwithlocCLEAN[str_detect(dfwithlocCLEAN$GENUS, sppname[[lol]][1]) & str_detect(dfwithlocCLEAN$GENE, "rag1"), ]), "\n")
# rewrite location table
write.table(dfwithlocCLEAN, paste0("../r_sequences/", sppnamePH, "/1dfwithloc_", sppnamePH, ".txt"), sep = ";", quote = F, row.names = F)

#view n
filefinal

base::message("ALL completed")

# } # THE END OF THE MAIN LOOP, reactivate this if you want to run it!

#==#==#==#==#==##==#==#==#==#==##==#==#==#==#==##==#==#==#==#==#
## ## ## 4 MASTER SETUP FOR BIOGEOBEARS: IMPORTANT  ####
#==#==#==#==#==##==#==#==#==#==##==#==#==#==#==##==#==#==#==#==#

# setup master tables (only once)
ana_dis_ALL_Table_MASTER <- data.frame(matrix(ncol = 10, nrow = 0))
colnames(ana_dis_ALL_Table_MASTER) <- c(
  "node", "node.type", "current_rangetxt", "event_type", 
  "event_txt", "ana_dispersal_from", "dispersal_to", 
  "event_time", "abs_event_time", "nmap"
)

clado_dis_ALL_Table_MASTER <- data.frame(matrix(ncol = 11, nrow = 0))
colnames(clado_dis_ALL_Table_MASTER) <- c(
  "node", "node.type", "sampled_states_AT_nodes", "clado_event_type", 
  "clado_event_txt", "clado_dispersal_from", "clado_dispersal_to", 
  "time_bp", "left_desc_nodes", "right_desc_nodes", "nmap"
)

#DO NOT RUN THE CODE ABOVE MORE THAN ONCE, IT WILL OVERWRITE THE MASTER TABLE

# important for next steps
sppnameplotPH <- c("Anolis",
                   "Bothriechis",
                   "Bothrops",
                   "Dipsas",
                   "Marisora",
                   "Micrurus",
                   "Oxybelis",
                   "Sibon")

# if you want to run for a specific genera, uncomment this and dont run the loop,
## only each line:
# sppnameplot <- "Oxybelis"

## IMPORTANT: You need to setup the folders for each genus before running this
# loop.
## IMPORTANT2:Some genera requires special atention in the loop, so I suggest running
# each genera one at a time


#==#==#==#==#==##==#==#==#==#==##==#==#==#==#==##==#==#==#==#==#
## ## ## 5 GEODATA MANIPULATION  ####
#==#==#==#==#==##==#==#==#==#==##==#==#==#==#==##==#==#==#==#==#

# access biogeographical regions
regionA <- terra::vect("../data/geography/A.shp")
regionB <- terra::vect("../data/geography/B.shp")
regionC <- terra::vect("../data/geography/C.shp")
regionALL <- rbind(regionA, regionB, regionC)

# access data from Costa Rica and Panamá provinces
panama_provinces <- gadm(country = "PAN", level = 1, path = tempdir())
head(panama_provinces)
length(panama_provinces)

# access data from Costa Rica provinces
costarica_provinces <- gadm(country = "CRI", level = 1, path = tempdir())
head(costarica_provinces)
length(costarica_provinces)

# access data from Nicarágua provinces
nicaragua_provinces <- gadm(country = "NIC", level = 1, path = tempdir())
head(nicaragua_provinces)
length(nicaragua_provinces)

# ploting
plot(regionALL, col = c("white", "grey77", "white"), xlim = c(-90, -75), ylim = c(13, 5))
plot(panama_provinces, col = NULL, add = T, border = "red")
plot(costarica_provinces, col = NULL, add = T, border = "blue")
plot(nicaragua_provinces, col = NULL, add = T, border = "green")

# centroids
points(centroids(panama_provinces), col = "red", pch = 16)
points(centroids(costarica_provinces), col = "blue", pch = 16)
points(centroids(nicaragua_provinces), col = "green", pch = 16)

# filtering provinces in which the centroid intersect with region B
panama_provinces_int <- panama_provinces[terra::is.related(centroids(panama_provinces), regionB, relation = "intersects"),]
costarica_provinces_int <- costarica_provinces[terra::is.related(centroids(costarica_provinces), regionB, relation = "intersects"),]
nicaragua_provinces_int <- nicaragua_provinces[terra::is.related(centroids(nicaragua_provinces), regionB, relation = "intersects"),]

# ploting now only intersects
plot(regionALL, col = c("white", "grey77", "white"), xlim = c(-90, -75), ylim = c(13, 5))
plot(panama_provinces_int, col = NULL, add = T, border = "red")
plot(costarica_provinces_int, col = NULL, add = T, border = "blue")
plot(nicaragua_provinces_int, col = NULL, add = T, border = "green")

# centroids
points(centroids(panama_provinces_int), col = "red", pch = 16)
points(centroids(costarica_provinces_int), col = "blue", pch = 16)
points(centroids(nicaragua_provinces_int), col = "green", pch = 16)

# joining and checking for intersect (GEOG file)
provincesinregionB <- rbind(panama_provinces_int, costarica_provinces_int, nicaragua_provinces_int)
provincesinregionB <- provincesinregionB$NAME_1
provincesinregionB <- c(stri_trans_general(provincesinregionB, "latin-ascii"), provincesinregionB)
provincesinregionB <- str_replace_all(provincesinregionB, " ", "_")
provincesinregionB <- provincesinregionB[order(provincesinregionB)]
provincesinregionB <- unique(provincesinregionB)
provincesinregionB

# for filtering north and south america
north_america <- ne_countries(continent = "north america", returnclass = "sf")
north_america <- str_replace_all(north_america$name, " ", "_")
# remove countries to measure
north_america <- north_america[!str_detect(north_america, paste0(c("Panama", "Costa_Rica", "Nicaragua"), collapse = "|"))]
north_america <- c(north_america, "USA_", "United_States")
# south america
south_america <- ne_countries(continent = "south america", returnclass = "sf")
south_america <- str_replace_all(south_america$name, " ", "_")

#==#==#==#==#==##==#==#==#==#==##==#==#==#==#==##==#==#==#==#==#
## ## ## MAIN LOOP FOR GENERA  ####
#==#==#==#==#==##==#==#==#==#==##==#==#==#==#==##==#==#==#==#==#
for (sppnameplot in sppnameplotPH){
  
  # reading table
  dfforgeog <- fread(paste0("../r_sequences/", sppnameplot, "/1dfwithloc_", sppnameplot, ".txt"), sep = ";")
  head(dfforgeog)
  
  # removing duplicates (only one specimen)
  dfforgeog_nd <- dfforgeog %>% distinct(VOUCHER, .keep_all = TRUE)
  #dfforgeog_nd <- dfforgeog_nd[order(dfforgeog_nd$GEO_LOC),]
  
  dfforgeog_final         <- dfforgeog[0,]
  dfforgeog_final$GEOG    <- character(0)
  
  # checking islands
  introisland_regex <- regex(
    "(^|\\b|_)(Hawaii|Palau|Northern_Mariana|Guam|Puerto_Rico|Guadeloupe|Martinique|Saint_Martin|Virgin_Islands|Curacao|Aruba|Bonaire|Saba|Cayman)(\\b|_|$)", 
    ignore_case = TRUE
  )
  
  # checking
  for (i in 1:nrow(dfforgeog_nd)){
    ph <- dfforgeog_nd[i,]
    
    # if islands, CHECK
    if (str_detect(ph$GEO_LOC, introisland_regex)) {
      
      ph$GEOG <- base::as.character("CHECK")
      
      # if not:
    } else if (str_detect(ph$GEO_LOC, regex(paste(north_america, collapse = "|"), ignore_case = T))){
      ph$GEOG <- base::as.character("A")
      
      # if below panama
    } else if (str_detect(ph$GEO_LOC, regex(paste(south_america, collapse = "|"), ignore_case = T))){
      ph$GEOG <- base::as.character("C")
      
      # if in nic, costa rica or panama, check for intersect
    } else {
      # if intersect with panama unit
      if (str_detect(ph$GEO_LOC, regex(paste(provincesinregionB, collapse = "|"), ignore_case = T))){
        ph$GEOG <- base::as.character("B")
      } else {
        # if does not intersec (rest of Nicaragua)
        if (str_detect(ph$GEO_LOC, regex("Nicaragua", ignore_case = T))){
          ph$GEOG <- base::as.character("A")
          # if does not intersect (rest of Panama)
        } else if (str_detect(ph$GEO_LOC, regex("Panama", ignore_case = T))){
          ph$GEOG <- base::as.character("C")
          # If all failed, check for countries
        } else {
          if (str_detect(ph$GEO_LOC, regex("^Nicaragua$", ignore_case = T))) {
            ph$GEOG <- base::as.character("A")
          } else if (str_detect(ph$GEO_LOC, regex("^Panama$", ignore_case = T))){
            ph$GEOG <- base::as.character("B")
          } else if (str_detect(ph$GEO_LOC, regex("^Costa_Rica$", ignore_case = T))) {
            ph$GEOG <- base::as.character("B")
            # IF ALL FAILED, CHECK
          } else {
            ph$GEOG <- base::as.character("CHECK")
          }
        }
      }
    }
    
    dfforgeog_final <- rbind(dfforgeog_final, ph)
    
  }
  
  checks <- dfforgeog_final[dfforgeog_final$GEOG == "CHECK", ]
  if (nrow(checks) > 0) {
    base::message("\n[", sppnameplot, "] ", nrow(checks),
                  " especimens em CHECK. Serao REMOVIDOS antes do BioGeoBEARS:")
    print(unique(checks[, c("ORGANISM", "VOUCHER", "GEO_LOC")]))
    dfforgeog_final <- dfforgeog_final[dfforgeog_final$GEOG != "CHECK", ]
  } else {
    base::message("[", sppnameplot, "] Nenhum CHECK.")
  }
  
  # Todas as tres areas devem estar presentes
  stopifnot(all(c("A", "B", "C") %in% unique(dfforgeog_final$GEOG)))
  
  # CHECK ALL UNIQUE DISTRIBUTIONS. AT LEAST ALL REGIONS SHOULD BE PRESENT
  unique(dfforgeog_final$GEOG)
  
  # REMOVE? BE CAREFUL WITH THESE LINES!
  ### FOR SOME GROUPS YES!
  # dfforgeog_final <- dfforgeog_final[!dfforgeog_final$GEOG == "CHECK",]
  ### 
  
  # check
  dfforgeog_final$GEOG
  
  #==#==#==#==#==##==#==#==#==#==##==#==#==#==#==##==#==#==#==#==#
  ## ## ## 6 SPECIMENS COLLAPSING AND REORGANIZING  ####
  #==#==#==#==#==##==#==#==#==#==##==#==#==#==#==##==#==#==#==#==#
  # tree
  # ============================================================================
  # NOTE:
  # ============================================================================
  # Run this line after running the delimitation methods following the Main
  # Document. See Materials and Methods. Create a folder called "r_grouping" on 
  # root and create subfolders with the name of each genera. Inside, save the 
  # results for each delimm method and rename as "<genera>_result_<delim>.txt".
  # <delim> should be "mptp", "GYMC" or "asap". Pay attention because each
  # delimm method outputs different files, but you need to rename the main
  # one (only the txt file output with the OUT grouping) for the following
  # part. 
  #=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
  sppnamedir  <- paste0("../r_grouping/", sppnameplot, "/", sppnameplot)
  sppnamedir2 <- paste0("../r_grouping/", sppnameplot, "/")
  td_full <- treeio::read.beast(here::here("..", "data", "biogeobears", sppnameplot, paste0(sppnameplot, "_MCC.tre")))
  tree    <- treeio::as.phylo(td_full)
  plot(tree)
  # cleaning tips without name
  if(any(str_detect(tree$tip.label, "^\\|"))){
    toexcludeWNN <- tree$tip.label[str_detect(tree$tip.label, "^\\|")]
    tree <- drop.tip(tree, toexcludeWNN)
    plot(tree)
  }
  
  ###=###=###= ###=###=###= ###=###=###
  # OTU delimitation 
  ###=###=###= ###=###=###= ###=###=###
  # ACTIVATE THIS IF YOU WILL RUN THE ANALYSES, OR LEAVE FALSE IF YOU WANT TO 
  ## GRAB WHAT IS ALREADY ON THE FOLDER
  areyourunningOTU <- FALSE # REMEMBER: TRUE -> Runs (need the inputs), FALSE -> no run (need the output)
  
  # main loop for delimitation methods
  if (areyourunningOTU){
    
    ###=###=###
    ### mPTP
    ###=###=###
    
    #=#=#
    # run mPTP
    #=#=#
    
    mptp_path <- paste0(sppnamedir, "_result_mptp.txt")
    lines <- readLines(mptp_path)
    if (!any(grepl("^Species [0-9]+:", lines))) {
      stop(mptp_path, " does not look like raw mPTP output. ",
           "It may have been overwritten by an earlier run. Restore the original.")
    }
    
    # if you want the outgroup to be part of the delim
    lerelere <- unique(str_extract(tree$tip.label, ".+?(?=_)"))
    sppnameplotSIS <- lerelere[lerelere != sppnameplot]
    #Extract only the lines containing Species headers and sample names
    species_lines <- lines[grep(paste0("^Species [0-9]+:|", sppnameplot, "|", sppnameplotSIS, "_"), lines)]
    species_lines <- species_lines[-1] # removing header
    head(species_lines)
    
    # Loop through the lines to assign species IDs to each sample
    results <- list()
    current_species <- NA
    
    for (line in species_lines) {
      if (grepl("^Species", line)) {
        # Extract the number from "Species 1:"
        current_species <- as.numeric(str_extract(line, "\\d+"))
      } else {
        # It's a sample name, so pair it with the current species ID
        results[[length(results) + 1]] <- data.frame(
          mPTP_spec = current_species,
          sample_name = trimws(line),
          stringsAsFactors = FALSE
        )
      }
    }
    
    # Bind into a single table
    mptp_df <- bind_rows(results)
    head(mptp_df)
    write.csv(mptp_df, paste0(sppnamedir, "_result_mptp_parsed.csv"), row.names = FALSE)
    
    ###=###=###
    ### GMYC
    ###=###=###
    result <- gmyc(tree)
    summary(result)
    #plot(result)
    
    gmyc_res <- result
    1 - pchisq(2 * (max(gmyc_res$likelihood) - gmyc_res$likelihood[1]), 2)
    
    #See detailed results
    spec.list(result)
    detailed_result <- spec.list(result)
    write.csv(detailed_result, paste0(sppnamedir, "_result_gmyc_parsed.csv"), row.names = FALSE)
    
    gmyc_df <- read.csv(paste0(sppnamedir, "_result_gmyc_parsed.csv"))
    head(gmyc_df)
    
    ###=###=###
    ### ASAP 
    ###=###=###
    # writing to format accepted
    # ============================================================================
    # NOTE:
    # ============================================================================
    # Before running the ASAP, we recomend you running this line first. For this,
    # you will need the aligment also present on that r_grouping folder. We
    # aligned it in section 3, so grab the aligment "<genus>__merged_trimmed_
    # sequences_0.50.fasta" and put it in each genus' folder. For the ASAP,
    # use the output of the following lines!
    #=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
    mol_for_asap <- read.FASTA(paste0(sppnamedir, "_merged_trimmed_sequences_0.50.fasta"))
    
    # removing species not present in tree obj
    if(any(str_detect(names(mol_for_asap), "^\\|"))){
      toexcludeWNN <- names(mol_for_asap)[str_detect(names(mol_for_asap), "^\\|")]
      mol_for_asap <- mol_for_asap[names(mol_for_asap)%in% toexcludeWNN == FALSE]
    }
    names(mol_for_asap)
    
    # removing specimens with no common site
    dists <- dist.dna(mol_for_asap, model = "raw", pairwise.deletion = TRUE)
    dist_mat <- as.matrix(dists)
    na_counts <- apply(dist_mat, 1, function(x) sum(is.na(x))) # sum NAs
    to_keep <- names(na_counts[na_counts == 0]) # leave only non NAs
    cleaned_dna <- mol_for_asap[to_keep] # removing
    cat("Removed", 
        length(mol_for_asap) - length(cleaned_dna), 
        "sequences:\n", 
        paste0(names(na_counts)[na_counts > 0], "\n")) #printing for control
    
    # writing into new file for asap
    write.FASTA(cleaned_dna, 
                file = paste0(sppnamedir, 
                             "_merged_trimmed_sequences_0.50_NEW.fas")) # <- USE THIS FOR ASAP
    
    #=#=#
    # run ASAP
    #=#=#
    
    # reading
    asap_path <- paste0(sppnamedir, "_result_asap.txt")
    lines <- readLines(asap_path)
    if (!any(grepl("^Subset\\[ [0-9]", lines))) {
      stop(asap_path, " does not look like raw ASAP output. ",
           "It may have been overwritten by an earlier run. Restore the original.")
    }
    
    #Extract only the lines containing Species headers and sample names
    species_lines <- lines[grep("^Subset\\[ [0-9]", lines)]
    head(species_lines)
    
    # Loop through the lines to assign species IDs to each sample
    results <- list()
    current_species <- NA
    
    # loop in all lines
    for (line in species_lines) {
      # extract number related to grouping
      current_species <- as.numeric(str_extract(line, "\\d+"))
      # extract species
      sppinline <- str_extract(line, "(?<=id:\\s).*") %>% str_trim()
      sppinline <- str_split(sppinline, "\\s+")[[1]]
      # loop inside all species in grouping
      for (spp in sppinline){
        # create list for each spp
        results[[length(results) + 1]] <- data.frame(
          ASAP_spec = current_species,
          sample_name = spp,
          stringsAsFactors = FALSE)
      }
    }
    
    # Bind into a single table
    asap_df <- bind_rows(results)
    head(asap_df)
    write.csv(asap_df, paste0(sppnamedir, "_result_asap_parsed.csv"), row.names = FALSE)
    
  } else {
    
    mptp_parsed <- paste0(sppnamedir, "_result_mptp_parsed.csv")
    gmyc_parsed <- paste0(sppnamedir, "_result_gmyc_parsed.csv")
    asap_parsed <- paste0(sppnamedir, "_result_asap_parsed.csv")
    
    stopifnot(
      "mPTP parsed file not found -- run once with areyourunningOUT = TRUE" = file.exists(mptp_parsed),
      "GMYC parsed file not found -- run once with areyourunningOUT = TRUE" = file.exists(gmyc_parsed),
      "ASAP parsed file not found -- run once with areyourunningOUT = TRUE" = file.exists(asap_parsed)
    )
    
    mptp_df <- read.csv(mptp_parsed)
    head(mptp_df)
    
    gmyc_df <- read.csv(gmyc_parsed)
    head(gmyc_df)
    
    asap_df <- read.csv(asap_parsed)
    head(asap_df)
    
  }
  
  ###=###=###
  ### BIDING AND COLLAPSING 
  ###=###=###
  # joining everyone in the same table
  all_species_results <- full_join(gmyc_df, mptp_df, by = "sample_name") %>%
    full_join(asap_df, by = "sample_name")
  # reorganizing
  all_species_results <- all_species_results[,c("sample_name", "mPTP_spec", "GMYC_spec", "ASAP_spec")]
  head(all_species_results)
  
  # adjusting tables
  specimens <- all_species_results$sample_name
  n <- length(specimens)
  # initializing an empty adjacency matrix
  adj_matrix <- matrix(0, nrow = n, ncol = n)
  rownames(adj_matrix) <- specimens
  colnames(adj_matrix) <- specimens
  # Fill the matrix: Link specimens if at least 2 methods agree
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      
      # Custom function to check equality safely
      check_match <- function(a, b) {
        if (is.na(a) | is.na(b)) return(FALSE) 
        return(a == b)
      }
      
      matches <- sum(
        check_match(all_species_results$mPTP_spec[i], all_species_results$mPTP_spec[j]),
        check_match(all_species_results$GMYC_spec[i], all_species_results$GMYC_spec[j]),
        check_match(all_species_results$ASAP_spec[i], all_species_results$ASAP_spec[j])
      )
      
      # If 2 or 3 methods agree, they are the same OTU
      if (matches >= 2) {
        adj_matrix[i, j] <- 1
        adj_matrix[j, i] <- 1
      }
    }
  }
  
  # Find Connected Components
  g <- graph_from_adjacency_matrix(adj_matrix, mode = "undirected")
  clusters <- igraph::components(g)$membership
  
  # Final Join for final OTU
  all_species_results_final <- all_species_results %>%
    mutate(OTU_majority = clusters[sample_name])
  
  # checking
  print(all_species_results_final)
  
  # putting the geog information into this new file
  all_species_results_final_final <- all_species_results_final[0,]
  all_species_results_final_final$GEOG <- character(0)
  for (i in 1:nrow(all_species_results_final)){
    ph  <- all_species_results_final[i,]
    ph2 <- str_extract(ph$sample_name, "(?<=\\|).*")
    ph3 <- dfforgeog_final[str_detect(dfforgeog_final$VOUCHER, ph2),]
    # check if removed lines or specimens
    if (base::identical(ph3, character(0)) | nrow(ph3) == 0){
      ph4 <- "REMOVED"
    } else {
      ph4 <- ph3$GEOG 
    }
    # putting all together
    all_species_results_final_final <- rbind(all_species_results_final_final, cbind(ph, GEOG = ph4))
  }
  
  # Ordening by voucher
  dfforgeog_final <- dfforgeog_final[order(dfforgeog_final$VOUCHER),]
  all_species_results_final_final <- all_species_results_final_final[order(str_extract(all_species_results_final_final$sample_name, "(?<=\\|).*")),]
  
  # testing if identical
  geogin1 <- dfforgeog_final$GEOG
  geogin2 <- all_species_results_final_final$GEOG[!all_species_results_final_final$GEOG == "REMOVED"]
  setequal(geogin1, geogin2) # must be true!
  
  # creating a colum with the representative tip for collapsing:
  otu_labels <- all_species_results_final_final %>%
    # Extract name before the pipe
    mutate(species_name = str_extract(sample_name, "^[^|]+")) %>%  
    group_by(OTU_majority) %>%
    summarize(
      # Collapse unique species names
      new_label = paste(unique(species_name), collapse = "_"),
      # Collapse unique GEOG letters (e.g., if a group has B and C, result is "BC")
      GEOG_collapsed = paste(sort(unique(GEOG)), collapse = ""),
      # Keep one tip to "hold the place"
      representative_tip = dplyr::first(sample_name), 
      # List of all tips to be removed
      all_tips = list(sample_name)
    )
  # and changing for unique
  otu_labels$new_label <- make.unique(otu_labels$new_label, sep = "_")
  
  # collapsing tree
  collapsed_tree <- tree
  
  # loop for collapsing all tips
  for(i in 1:nrow(otu_labels)) {
    current_otu <- otu_labels[i, ]
    tips_to_remove <- setdiff(unlist(current_otu$all_tips), current_otu$representative_tip)
    
    # drop.tip if specimens were removed after cleaning at the beggining (missed for some reason, see steps above)
    if (current_otu$GEOG_collapsed == "REMOVED"){
      collapsed_tree <- drop.tip(collapsed_tree, unlist(current_otu$all_tips))
      next
    }
    
    # Remove the extra tips
    if(length(tips_to_remove) > 0) {
      collapsed_tree <- drop.tip(collapsed_tree, tips_to_remove)
    }
    
    # Rename the remaining representative tip
    collapsed_tree$tip.label[collapsed_tree$tip.label == current_otu$representative_tip] <- current_otu$new_label
  }
  
  # Clean up any singleton nodes created by dropping tips
  collapsed_tree <- multi2di(collapsed_tree)
  
  # This is the new tree with tips collapsed
  plotTree(ladderize(collapsed_tree))
  
  # removing outgroup
  otu_labels_NOOUT <- otu_labels[str_detect(otu_labels$new_label, sppnameplot),] # only in sppnameplot
  collapsed_tree <- drop.tip(collapsed_tree, collapsed_tree$tip.label[!str_detect(collapsed_tree$tip.label, sppnameplot)])
  
  # removing from tree the REMOVED labels
  otu_labelsPH <- otu_labels[str_detect(otu_labels$GEOG_collapsed, "^REMOVED$"),]
  otu_labelsPH <- base::unlist(otu_labelsPH$all_tips)
  tree <- drop.tip(tree, c(otu_labelsPH))
  
  # sometimes specimens mix REMOVED with some areas, lets leave the specimen with 1 representative without REMOVED
  otu_labelsPH <- otu_labels[!str_detect(otu_labels$GEOG_collapsed, "^REMOVED$"),] # removing all the "REMOVED" only
  otu_labelsPH$GEOG_collapsed   <- str_remove_all(otu_labelsPH$GEOG_collapsed, "REMOVED") # fixing remaining "REMOVED"
  otu_labels <- otu_labelsPH # now removing "REMOVED"
  
  # repeat for no out
  # sometimes specimens mix REMOVED with some areas, lets leave the specimen with 1 representative without REMOVED
  otu_labelsPH <- otu_labels_NOOUT[!str_detect(otu_labels_NOOUT$GEOG_collapsed, "^REMOVED$"),] # removing all the "REMOVED" only
  otu_labelsPH$GEOG_collapsed   <- str_remove_all(otu_labelsPH$GEOG_collapsed, "REMOVED") # fixing remaining "REMOVED"
  otu_labels_NOOUT <- otu_labelsPH # now removing "REMOVED"
  
  #==#==#==#==#==##==#==#==#==#==##==#==#==#==#==##==#==#==#==#==#
  ## ## ## 7 INPUT FILES FOR BIOGEOBEARS  ####
  #==#==#==#==#==##==#==#==#==#==##==#==#==#==#==##==#==#==#==#==#
  # now lets make the inputfiles for running the BioGeoBEARS
  dir.create(paste0("./BioGeoBEARS/", sppnameplot, "/"), recursive = TRUE, showWarnings = FALSE)
  
  ###=###=###
  ### TREE
  ###=###=###
  write.tree(ladderize(collapsed_tree), paste0("../data/biogeobears/", sppnameplot, "/", sppnameplot, "_tree_final.newick"))
  
  ###=###=###
  ### GEOG
  ###=###=###
  # define geog lookup
  geog_lookup <- c(
    "A"  = "100",
    "B"  = "010",
    "C"  = "001",
    "AB" = "110",
    "AC" = "101",
    "BC" = "011",
    "ABC" = "111"
  )
  
  # loop to create the geog file
  geogtext <- paste0(length(collapsed_tree$tip.label), "\t", "3\t(A\tB\tC)\n")
  for (i in 1:nrow(otu_labels_NOOUT)){
    
    #skip removed
    if (otu_labels_NOOUT[i,]$GEOG_collapsed == "REMOVED"){next}
    
    # starting
    otuforgeog <- otu_labels_NOOUT[i,]$new_label
    # testing for geog
    resultforgeog <- unname(geog_lookup[otu_labels_NOOUT[i,]$GEOG_collapsed])
    # check if NA
    if (is.na(resultforgeog)){
      stop(paste0("Loop terminated. Check NA in ", otu_labels_NOOUT[i,]$new_label))
      break
    }
    geogtextph <- paste0(otuforgeog, "\t", resultforgeog, "\n")
    
    # combining
    geogtext <- paste0(geogtext, geogtextph)
    
  }
  
  # checking:
  cat(geogtext)
  
  # saving:
  write(geogtext, paste0("../data/biogeobears/", sppnameplot, "/", sppnameplot, "_geog_final.txt"))
  
  #==#==#==#==#==##==#==#==#==#==##==#==#==#==#==##==#==#==#==#==#
  ## ## ## 8 PLOTING  ####
  #==#==#==#==#==##==#==#==#==#==##==#==#==#==#==##==#==#==#==#==#
  # organizing plot with outgroup
  otu_labels <- otu_labels %>%
    mutate(new_label = ifelse(str_detect(new_label, sppnameplot), 
                              new_label, 
                              "Outgroup"))
  otu_labels <- otu_labels %>%
    mutate(GEOG_collapsed = ifelse(new_label == "Outgroup", "", GEOG_collapsed))
  
  # organizing plot
  plot_data <- all_species_results_final %>%
    left_join(otu_labels %>% dplyr::select(OTU_majority, new_label), by = "OTU_majority")
  plot_data$new_label[is.na(plot_data$new_label)] <- "Outgroup"
  plot_data <- plot_data %>% dplyr::rename(label = sample_name)
  
  # Create the base tree
  p <- ggtree(tree) %<+% plot_data + 
    geom_tiplab(size = 5, align = TRUE, linetype = "dotted", linesize = 0.5) +
    hexpand(0.5, direction = 1) +# Add space for the bars
    vexpand(0.02, direction = 1)
  
  p <- revts(p)
  
  # grabing HPD 95%
  HPD_COL <- "height_0.95_HPD"
  
  if (HPD_COL %in% names(td_full@data)) {
    phy_full <- treeio::as.phylo(td_full)
    hpd_tbl  <- td_full@data
    hpd_tbl$node <- as.integer(hpd_tbl$node)
    
    # tips descendentes from each internal node from the FULL tree
    internal_full <- (ape::Ntip(phy_full) + 1):(ape::Ntip(phy_full) + phy_full$Nnode)
    desc_full <- lapply(internal_full, function(nd)
      sort(phy_full$tip.label[phangorn::Descendants(phy_full, nd, "tips")[[1]]]))
    names(desc_full) <- internal_full
    
    p_int <- p$data %>% dplyr::filter(!isTip)
    
    seg <- lapply(seq_len(nrow(p_int)), function(k) {
      nd_c <- p_int$node[k]
      tips_here <- sort(tree$tip.label[phangorn::Descendants(tree, nd_c, "tips")[[1]]])
      # no cheio cujo conjunto de tips contem exatamente estes tips focais
      ok <- which(vapply(desc_full, function(d) all(tips_here %in% d), logical(1)))
      if (length(ok) == 0) return(NULL)
      cand    <- as.integer(names(desc_full)[ok])
      nd_full <- cand[which.min(lengths(desc_full[as.character(cand)]))]  # o mais raso que contem
      hpd <- hpd_tbl[[HPD_COL]][hpd_tbl$node == nd_full]
      if (length(hpd) == 0 || is.null(hpd[[1]])) return(NULL)
      vals <- as.numeric(hpd[[1]])
      data.frame(y = p_int$y[k], x_min = -max(vals), x_max = -min(vals))
    })
    seg <- do.call(rbind, seg)
    
    if (!is.null(seg) && nrow(seg) > 0) {
      p <- p + geom_segment(
        data = seg, inherit.aes = FALSE,
        aes(x = x_min, xend = x_max, y = y, yend = y),
        color = "steelblue", alpha = 0.35, linewidth = 2.5, lineend = "round"
      )
    }
  } else {
    base::message("[", sppnameplot, "] HPD column '", HPD_COL,
                  "' not found. Run names(td_full@data) and set HPD_COL.")
  }
  
  # coordinates for ploting boxes
  tree_data <- p$data %>% 
    dplyr::filter(isTip) %>% 
    dplyr::select(label, y)
  
  # Merge coordinates into your plot_data
  plot_data_coord <- plot_data %>%
    left_join(tree_data, by = "label")
  
  # Define a function to get the boundaries of each group
  get_boundaries <- function(df, col_name) {
    df %>%
      dplyr::group_by(!!sym(col_name)) %>%
      dplyr::summarise(ymin = min(y) - 0.5, ymax = max(y) + 0.5, .groups = "drop")
  }
  
  # Create boundary data for each analysis
  mptp_bounds <- get_boundaries(plot_data_coord, "mPTP_spec")
  gmyc_bounds <- get_boundaries(plot_data_coord, "GMYC_spec")
  asap_bounds <- get_boundaries(plot_data_coord, "ASAP_spec")
  otu_bounds  <- get_boundaries(plot_data_coord, "OTU_majority")
  
  mptp_bounds_colored <- mptp_bounds %>%
    left_join(unique(plot_data_coord[, c("mPTP_spec", "OTU_majority")]), by = "mPTP_spec")
  gmyc_bounds_colored <- gmyc_bounds %>%
    left_join(unique(plot_data_coord[, c("GMYC_spec", "OTU_majority")]), by = "GMYC_spec")
  asap_bounds_colored <- asap_bounds %>%
    left_join(unique(plot_data_coord[, c("ASAP_spec", "OTU_majority")]), by = "ASAP_spec")
  
  #=#=#=# GEOMETRY scaled to tree depth =#=#=#
  xr <- diff(range(p$data$x, na.rm = TRUE))
  
  lab_w <- 0.42 * xr    # tip-label band (increase if tip names are clipped)
  bar_w <- 0.05  * xr   # width of each bar
  gap   <- 0.015 * xr   # gap between bars
  txt_w <- 0.60 * xr    # OTU-name band
  
  bar_start <- lab_w
  b1 <- bar_start
  b2 <- b1 + bar_w + gap
  b3 <- b2 + bar_w + gap
  b4 <- b3 + bar_w + gap
  x_txt <- b4 + bar_w * 0.4 + gap
  
  mptp_center  <- b1 + bar_w / 2
  gmyc_center  <- b2 + bar_w / 2
  asap_center  <- b3 + bar_w / 2
  final_center <- x_txt
  
  x_left   <- min(p$data$x, na.rm = TRUE)
  x_right  <- x_txt + txt_w
  header_y <- max(p$data$y, na.rm = TRUE) + 2
  
  #=#=#=#  OTU labels with line wrapping =#=#=#
  otu_text <- otu_bounds %>%
    dplyr::mutate(y_mid = (ymin + ymax) / 2) %>%
    dplyr::left_join(
      otu_labels %>% dplyr::select(OTU_majority, new_label, GEOG_collapsed) %>% dplyr::distinct(),
      by = "OTU_majority"
    ) %>%
    dplyr::filter(!is.na(new_label), new_label != "Outgroup") %>%
    dplyr::mutate(
      base = gsub("_", " ", new_label),
      base = ifelse(GEOG_collapsed == "" | is.na(GEOG_collapsed),
                    base, paste0(base, " (", GEOG_collapsed, ")")),
      label_wrapped = stringr::str_wrap(base, width = 30)
    )
  n_linhas_extra <- sum(stringr::str_count(otu_text$label_wrapped, "\n"))
  
  #=#=#=# p_clean =#=#=#
  rect_layer <- function(dat, x0) {
    geom_rect(data = dat,
              aes(xmin = x0, xmax = x0 + bar_w, ymin = ymin, ymax = ymax,
                  fill = as.factor(OTU_majority)),
              color = "black", linewidth = 0.2, inherit.aes = FALSE)
  }
  
  p_clean <- p +
    rect_layer(mptp_bounds_colored, b1) +
    rect_layer(gmyc_bounds_colored, b2) +
    rect_layer(asap_bounds_colored, b3) +
    geom_rect(data = otu_bounds,
              aes(xmin = b4, xmax = b4 + bar_w * 0.4,
                  ymin = ymin, ymax = ymax, fill = as.factor(OTU_majority)),
              color = "black", linewidth = 0.2, inherit.aes = FALSE) +
    
    geom_text(data = otu_text,
              aes(x = x_txt, y = y_mid, label = label_wrapped),
              hjust = 0, vjust = 0.5, size = 5, lineheight = 0.85,
              fontface = "italic", inherit.aes = FALSE) +
    
    annotate("text", x = mptp_center,  y = header_y, label = "mPTP",      size = 4, fontface = "bold") +
    annotate("text", x = gmyc_center,  y = header_y, label = "GMYC",      size = 4, fontface = "bold") +
    annotate("text", x = asap_center,  y = header_y, label = "ASAP",      size = 4, fontface = "bold") +
    annotate("text", x = final_center, y = header_y, label = "Final OTU", size = 4, fontface = "bold", hjust = 0) +
    
    theme_tree2() +
    scale_fill_viridis_d() +
    theme(legend.position = "none", plot.margin = margin(5, 5, 5, 5)) +
    
    coord_geo(dat = "periods", neg = TRUE, abbrv = TRUE,
              xlim = c(x_left, x_right),
              ylim = c(0, header_y + 1),
              expand = FALSE, clip = "off")
  
  #=#=#=# save, dimensions in INCHES =#=#=#
  n_tips <- ape::Ntip(tree)
  h_in <- max(6, n_tips * 0.16 + n_linhas_extra * 0.10)
  w_in <- 18
  
  # ============================================================================
  # NOTE:
  # ============================================================================
  # Remeber to create the r_grouping folder. 
  # Remeber to create the r_grouping folder. 
  # Remeber to create the r_grouping folder. 
  # Remeber to create the r_grouping folder. 
  # Remeber to create the r_grouping folder. 
  # Remeber to create the r_grouping folder. 
  # Remeber to create the r_grouping folder. 
  # (in root of course)
  #=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
  
  ggsave(paste0("../r_grouping/", sppnameplot, "/", sppnameplot, "_pdf.pdf"),
         plot = p_clean, width = w_in, height = h_in,
         device = "pdf", limitsize = FALSE)
  
  #==#==#==#==#==##==#==#==#==#==##==#==#==#==#==##==#==#==#==#==#
  ## ## ## 9 BIOGEOBEARS TIME  ####
  #==#==#==#==#==##==#==#==#==#==##==#==#==#==#==##==#==#==#==#==#
  # taking tree
  trfn = np(paste0("../data/biogeobears/", sppnameplot, "/", sppnameplot, "_tree_final.newick"))
  tr = read.tree(trfn)
  trfnPH = sppnameplot
  
  # taking geog
  geogfn = np(paste0("../data/biogeobears/", sppnameplot, "/", sppnameplot, "_geog_final.txt"))
  # Look at the raw geography text file:
  moref(geogfn)
  
  # set tipranges
  ranges = getranges_from_LagrangePHYLIP(lgdata_fn=geogfn)
  # write.csv(ranges@df, "ranges.txt", quote = F)
  any(is.na(ranges))
  
  # set tip ranges
  tipranges=order_tipranges_by_tr(ranges, tr)
  
  # set max range size
  max_range_size = 3
  
  # check number of sates
  numstates_from_numareas(numareas=length(tipranges@df), maxareas=max_range_size, include_null_range=TRUE)
  
  for (i in 1:6){
    # 12 is the number of models to run 
    # (base + with or without j + with or without TS + with or without trait + with or without m2+m3)
    
    # cleaning objects
    BioGeoBEARS_run_object = NULL
    resfna = NULL
    resfnb = NULL
    resfnc = NULL
    resfnd = NULL
    resfnPLOT = NULL
    resfn = NULL
    checkTS = FALSE
    checkJ = NULL
    checkTRAIT = FALSE
    checkM = FALSE
    
    # debug
    for (v in 1:10){
      tryCatch({
        dev.off()  # Turn off PDF
        dev.off()  # Turn off PDF
        dev.off()  # Turn off PDF
        dev.off()  # Turn off PDF
        dev.off()  # Turn off PDF
        dev.off()  # Turn off PDF
      }, error=function(e){})
    }
    
    ##=### DEFAULT CONFIG ##=##=##=##=##=##=##=##=##=##=#
    # Intitialize a default model (DEC model)
    BioGeoBEARS_run_object = NULL
    BioGeoBEARS_run_object = define_BioGeoBEARS_run()
    # Give BioGeoBEARS the location of the phylogeny Newick file
    BioGeoBEARS_run_object$trfn = trfn
    # Give BioGeoBEARS the location of the geography text file
    BioGeoBEARS_run_object$geogfn = geogfn
    # Input the maximum range size
    BioGeoBEARS_run_object$max_range_size = max_range_size
    # Min to treat tip as a direct ancestor (no speciation event)
    BioGeoBEARS_run_object$min_branchlength = 0.000001
    min_branchlength <- BioGeoBEARS_run_object$min_branchlength
    # set to FALSE for e.g. DEC* model, DEC*+J, etc.
    BioGeoBEARS_run_object$include_null_range = TRUE
    ##=##=##=##=##=##=##=##=##=##=##=##=##=##=##=##=##=##=##=##=#
    
    
    ##=### MORE DEFAULT CONFIG ##=##=##=##=##=##=##=##=##=##=#
    # Speed options and multicore processing if desired
    BioGeoBEARS_run_object$on_NaN_error = -1e50    
    BioGeoBEARS_run_object$speedup = TRUE
    
    BioGeoBEARS_run_object$use_optimx = "GenSA"    
    BioGeoBEARS_run_object$num_cores_to_use = parallel::detectCores()
    
    BioGeoBEARS_run_object$force_sparse = FALSE
    
    # This function loads the dispersal multiplier matrix etc. 
    # from the text files   into the model object. Required for these to work!
    # (It also runs some checks on these inputs for certain errors.)
    BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
    BioGeoBEARS_run_object$return_condlikes_table = TRUE
    BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
    BioGeoBEARS_run_object$calc_ancprobs = TRUE # get ancestral 
    ##=##=##=##=##=##=##=##=##=##=##=##=##=##=##=##=##=##=##=#
    
    ##=### DEFINING MODELS ##=##=##=##=##=##=##=##=##=##=
    # DEC models are defined as default, so just run the script! ;)
    ##=##=##=##=##=##=##=##=##=##=##=##=##=##=##=##=##=##
    
    if (i>=1 & i<=3){ ## IF IT DOES NOT HAVE +j 
      # (DEC, DIVALIKE AND BAYAREALIKE)
      
      checkJ = FALSE # check if model has j
      
      if (i%%3==1){  ## DEFINING DEC
        
        # use default options
        
      }else{
        
        if (i%%3==2){ ## DEFINING DIVALIKE
          
          # Set up DIVALIKE model
          # Remove subset-sympatry
          BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
          BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
          BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0
          
          BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "2-j"
          BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/2"
          BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "ysv*1/2"
          BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "ysv*1/2"
          
          # Allow classic, widespread vicariance; all events equiprobable
          BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","type"] = "fixed"
          BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","init"] = 0.5
          BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","est"] = 0.5
          
          
        }else{
          
          if (i%%3==0){ ## DEFINING BAYAREALIKE
            
            # Set up BAYAREALIKE model
            # No subset sympatry
            BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
            BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
            BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0
            
            # No vicariance
            BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "fixed"
            BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","init"] = 0.0
            BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","est"] = 0.0
            
            # Adjust linkage between parameters
            BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "1-j"
            BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/1"
            BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "1-j"
            
            # Only sympatric/range-copying (y) events allowed, and with 
            # exact copying (both descendants always the same size as the ancestor)
            BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","type"] = "fixed"
            BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","init"] = 0.9999
            BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","est"] = 0.9999
            
          }
          
        }
        
      }
      
    }else{ ## IF IT DOES HAVE +j (DEC+j, DIVALIKE+j AND BAYAREALIKE+j)
      
      checkJ = TRUE
      
      if (i%%3==1){  ## DEFINING DEC+j
        
        # Set up DEC+J model
        # Get the ML parameter values from the 2-parameter nested model
        # (this will ensure that the 3-parameter model always does at least as good)
        if (exists(ls(pattern = paste(trfnPH,  "resDEC$", sep = "_"))) == TRUE){
          resDECPH = lapply(ls(pattern = paste(trfnPH,  "resDEC$", sep = "_")), get)
        }else{
          errormessage <- "STOP ERROR: you are trying to run a +j model without running a base model. \nAs the +j models gets the start parameter values from the base models, \nthere is no way of running the model. RUN the base model first and THEN run the +j model"
          stop(errormessage)
        }
        
        dstart = resDECPH[[1]]$outputs@params_table["d","est"]
        estart = resDECPH[[1]]$outputs@params_table["e","est"]
        jstart = 0.009
        
        # Input starting values for d, e
        BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","init"] = dstart
        BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","est"] = dstart
        BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","init"] = estart
        BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","est"] = estart
        
        # Add j as a free parameter
        BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
        BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = jstart
        BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = jstart
        BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","min"] = 0.001
        BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","max"] = 2.999
        
        
      }else{
        
        if (i%%3==2){ ## DEFINING DIVALIKE+j
          
          # Set up DIVALIKE+J model
          # Get the ML parameter values from the 2-parameter nested model
          # (this will ensure that the 3-parameter model always does at least as good)
          if (exists(ls(pattern = paste(trfnPH, "resDIVALIKE$", sep = "_"))) == TRUE){
            resDIVALIKEPH = lapply(ls(pattern = paste(trfnPH, "resDIVALIKE$", sep = "_")), get)
          }else{
            errormessage <- "STOP ERROR: you are trying to run a +j model without running a base model. \nAs the +j models gets the start parameter values from the base models, \nthere is no way of running the model. RUN the base model first and THEN run the +j model"
            stop(errormessage)
          }
          dstart = resDIVALIKEPH[[1]]$outputs@params_table["d","est"]
          estart = resDIVALIKEPH[[1]]$outputs@params_table["e","est"]
          jstart = 0.009
          
          # Input starting values for d, e
          BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","init"] = dstart
          BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","est"] = dstart
          BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","init"] = 0.001
          BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","est"] = 0.001
          BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","min"] = 0.001
          BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","max"] = 4.999
          BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","min"] = 0.001
          BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","max"] = 4.999
          
          # Remove subset-sympatry
          BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
          BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
          BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0
          
          BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "2-j"
          BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/2"
          BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "ysv*1/2"
          BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "ysv*1/2"
          
          # Allow classic, widespread vicariance; all events equiprobable
          BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","type"] = "fixed"
          BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","init"] = 0.5
          BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","est"] = 0.5
          
          # Add jump dispersal/founder-event speciation
          BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
          BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = jstart
          BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = jstart
          
          # Under DIVALIKE+J, the max of "j" should be 2, not 3 (as is default in DEC+J)
          BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","min"] = 0.001
          BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","max"] = 1.999
          
          
        }else{
          
          if (i%%3==0){ ## DEFINING BAYAREALIKE+j
            
            # Set up BAYAREALIKE+J model
            # Get the ML parameter values from the 2-parameter nested model
            # (this will ensure that the 3-parameter model always does at least as good)
            if (exists(ls(pattern = paste(trfnPH, "resBAYAREALIKE$", sep = "_"))) == TRUE){
              resBAYAREALIKEPH = lapply(ls(pattern = paste(trfnPH, "resBAYAREALIKE$", sep = "_")), get)
            }else{
              errormessage <- "STOP ERROR: you are trying to run a +j model without running a base model. \nAs the +j models gets the start parameter values from the base models, \nthere is no way of running the model. RUN the base model first and THEN run the +j model"
              stop(errormessage)
            }
            dstart = resBAYAREALIKEPH[[1]]$outputs@params_table["d","est"]
            estart = resBAYAREALIKEPH[[1]]$outputs@params_table["e","est"]
            jstart = 0.009
            
            # Input starting values for d, e
            BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","init"] = dstart
            BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","est"] = dstart
            BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","init"] = estart
            BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","est"] = estart
            
            # No subset sympatry
            BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
            BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
            BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0
            
            # No vicariance
            BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "fixed"
            BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","init"] = 0.0
            BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","est"] = 0.0
            
            # *DO* allow jump dispersal/founder-event speciation (set the starting value close to 0)
            BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
            BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = jstart
            BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = jstart
            
            # Under BAYAREALIKE+J, the max of "j" should be 1, not 3 (as is default in DEC+J) or 2 (as in DIVALIKE+J)
            BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","max"] = 0.99999
            
            # Adjust linkage between parameters
            BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "1-j"
            BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/1"
            BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "1-j"
            
            # Only sympatric/range-copying (y) events allowed, and with
            # exact copying (both descendants always the same size as the ancestor)
            BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","type"] = "fixed"
            BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","init"] = 0.9999
            BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","est"] = 0.9999
            
            # NOTE (NJM, 2014-04): BAYAREALIKE+J seems to crash on some computers, usually Windows
            # machines. I can't replicate this on my Mac machines, but it is almost certainly
            # just some precision under-run issue, when optim/optimx tries some parameter value
            # just below zero.  The "min" and "max" options on each parameter are supposed to
            # prevent this, but apparently optim/optimx sometimes go slightly beyond
            # these limits.  Anyway, if you get a crash, try raising "min" and lowering "max"
            # slightly for each parameter:
            BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","min"] = 0.00001
            BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","max"] = 4.99999
            
            BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","min"] = 0.00001
            BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","max"] = 4.99999
            
            BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","min"] = 0.001
            BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","max"] = 0.999
            
            
            if (dstart < BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","min"]){
              BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","init"] = BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","min"]
              BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","est"] = BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","min"]
            }
            
            if (estart < BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","min"]){
              BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","init"] = BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","min"]
              BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","est"] = BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","min"]
            }
            
          }
          
        }
        
      }
      
    }
    
    ##=### RUNNING AFTER ALL ##=##=##=##=##=##=##=##=##=##=##=
    # LAST configs to run the models accordingly
    ##=##=##=##=##=##=##=##=##=##=##=##=##=##=##=##=##=##=##=#
    
    # check if its everything allright
    check_BioGeoBEARS_run(BioGeoBEARS_run_object)
    
    
    # check model base
    if (i == 0){resfna = "JUST"; a = 0}else{
      if (i%%3==1){resfna = "DEC"; a = 1}
      if (i%%3==2){resfna = "DIVALIKE"; a = 2}
      if (i%%3==0){resfna = "BAYAREALIKE"; a = 3}
    }
    
    # check TS
    if (checkTS == TRUE){resfnb = "TS"; b = 2}else{resfnb = ""; b = 1}
    
    # check j
    if (checkJ == TRUE){resfnc = "j"; c = 2}else{resfnc = ""; c = 1}
    
    # check trait
    if (checkTRAIT == TRUE){resfnd = "traitparam"; d = 2}else{resfnd = ""; d = 1}
    
    # check M
    resfne = ""; e = 1
    
    # configuring resfn 
    resfnPLOT = paste(resfna, resfnb, resfnc, resfnd, resfne, sep = "")
    resfn = paste("../results/",sppnameplot, "/", sppnameplot, "_", resfnPLOT,"_res_MPN.RData", sep = "")
    
    # debug
    for (v in 1:10){
      tryCatch({
        dev.off()  # Turn off PDF
        dev.off()  # Turn off PDF
        dev.off()  # Turn off PDF
        dev.off()  # Turn off PDF
        dev.off()  # Turn off PDF
        dev.off()  # Turn off PDF
      }, error=function(e){})
    }
    
    # running 
    runslow = TRUE
    
    if (runslow) {
      
      # trace(bears_optim_run, edit = T)
      res = bears_optim_run(BioGeoBEARS_run_object)
      save(res, file=resfn)
      assign(paste0(trfnPH, "_res", resfnPLOT), res)
      
    } else {
      # Loads to "res"
      load(resfn)
      assign(paste0(trfnPH, "_res", resfnPLOT), res)
      
    }
    
    ##=##=##=##=##=##=##=##=##=##
    # Plot ancestral states
    ##=##=##=##=##=##=##=##=##=##
    
    results_object =  get(paste0(trfnPH, "_res", resfnPLOT))
    
    # plot all of then
    print(paste("ploting", resfnPLOT, "using", trfnPH))
    pdffn = paste("../figures/",sppnameplot, "/", sppnameplot, "_", resfnPLOT, "_MPN.pdf", sep = "")
    pdf(pdffn, width=12, height=25)
    
    ##=##=##=##=##=##=##=##=##=##
    # Plot ancestral states
    ##=##=##=##=##=##=##=##=##=##
    analysis_titletxt =  paste("BioGeoBEARS", resfnPLOT, "on", trfnPH, sep = " ")
    
    # Setup
    scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))
    
    # plot geog and trait if model has trait
    print(plot_BioGeoBEARS_results(results_object, 
                                   analysis_titletxt, 
                                   addl_params=list("j"), 
                                   plotwhat="text", 
                                   label.offset=0.45, 
                                   tipcex=0.7, statecex=0.7, splitcex=0.6, 
                                   titlecex=0.8, plotsplits=TRUE, 
                                   cornercoords_loc=scriptdir, 
                                   include_null_range=TRUE, tr=tr, 
                                   tipranges=tipranges))
    
    # Pie chart
    tryCatch({
      print(plot_BioGeoBEARS_results(results_object, 
                                     analysis_titletxt, 
                                     addl_params=list("j"), 
                                     plotwhat="pie", 
                                     label.offset=0.45, 
                                     tipcex=0.7, statecex=0.7, splitcex=0.6, 
                                     titlecex=0.8, plotsplits=TRUE, 
                                     cornercoords_loc=scriptdir, 
                                     include_null_range=TRUE, tr=tr, 
                                     tipranges=tipranges))
    }, error=function(e){})
    
    # debug
    for (v in 1:10){
      tryCatch({
        dev.off()  # Turn off PDF
        dev.off()  # Turn off PDF
        dev.off()  # Turn off PDF
        dev.off()  # Turn off PDF
        dev.off()  # Turn off PDF
        dev.off()  # Turn off PDF
      }, error=function(e){})
    }
    
    ##=##=##=##=##=##=##=##=##=##=##=##=##=##=##=##=##=##=##=#
    # cleaning objects
    BioGeoBEARS_run_object = NULL
    resfna = NULL
    resfnb = NULL
    resfnc = NULL
    resfnd = NULL
    resfnPLOT = NULL
    resfn = NULL
    checkTS = FALSE
    checkJ = NULL
    checkTRAIT = FALSE
    checkM = FALSE
    
  }
  
  #==#==#==#==#==##==#==#==#==#==##==#==#==#==#==##==#==#==#==#==#
  ## ## ## 10 STATISTICS  ####
  #==#==#==#==#==##==#==#==#==#==##==#==#==#==#==##==#==#==#==#==#
  #=#=#=#=#
  # organizing models 
  #=#=#=#=#
  
  # Set up empty tables to hold the statistical results
  restable = NULL
  teststable = NULL
  
  # take all the res for tree
  trnames <- paste0(sppnameplot, "_res")
  
  # listing all results 
  resALLTREE <- mget(ls(pattern = trnames)) # getting res
  # resALL <- resALL[c("resDEC", "resDIVAL", "resBAYAREAL", "resDECTS", "resDIVALTS", "resBAYAREALTS", "resDECj", "resDIVALj", "resBAYAREALj", "resDECTSj", "resDIVALTSj", "resBAYAREALTSj")] 
  names(resALLTREE)
  resALLTREE <- resALLTREE[names(resALLTREE) %in% ls(pattern = "PH$") == FALSE]
  resALLTREE <- resALLTREE[names(resALLTREE) %in% ls(pattern = "table$") == FALSE]
  names(resALLTREE)
  # resALL2 <- mget(ls(pattern = "resDEC|resDIVA|resBAYAREA")) # getting res
  # View(resALL2)
  
  #=#=#=#=#=#=#=#=#=#=#
  # select models to compare
  #=#=#=#=#=#=#=#=#=#=#
  # select from bellow what you want to include in the statistics
  # this also will influence the plot bellow
  comparewhat = c(
    "ALL", ## <- name tag of the output files
    "DEC = TRUE", ## <- include DEC? TRUE or FALSE
    "DIVA = TRUE", ## <- include DIVALIKE? TRUE or FALSE
    "BAY = TRUE", ## <- include BAYAREALIKE? TRUE or FALSE
    "j = TRUE", ## <- include j? TRUE, FALSE or ONLY
    "TS = FALSE", ## <- include TS? TRUE, FALSE or ONLY
    "trait = FALSE", ## <- include TRAIT? TRUE, FALSE or ONLY
    "m = FALSE" ## <- include m2m3? TRUE, FALSE or ONLY
  )
  # choose from: DEC, DIVA and BAYAREA (input these strings), and if you want TS and j
  dirfortableresults <- paste0("../results/",sppnameplot)
  dir.create(dirfortableresults, recursive = TRUE, showWarnings = FALSE)
  
  #=#=#=#=#=#=#=#=#=#=#
  # separating models 
  #=#=#=#=#=#=#=#=#=#=#
  resALL <- resALLTREE
  
  # code
  for (j in 1:1){
    resCHOSEN = NULL
    NOPTEST = "FALSE"
    if (stringr::str_detect(comparewhat[2], 'TRUE') == TRUE){ # checking DEC
      resCHOSEN = c(resCHOSEN, resALL[names(resALL) %in% ls(pattern = "DEC") == TRUE])
    }
    if (stringr::str_detect(comparewhat[3], 'TRUE') == TRUE){ # checking DIVA
      resCHOSEN = c(resCHOSEN, resALL[names(resALL) %in% ls(pattern = "DIVA") == TRUE])
    }
    if (stringr::str_detect(comparewhat[4], 'TRUE') == TRUE){ # checking BAYAREA
      resCHOSEN = c(resCHOSEN, resALL[names(resALL) %in% ls(pattern = "BAY") == TRUE])
    }
    
    ##=##=##=##=##=#
    # checking j
    if (stringr::str_detect(comparewhat[5], 'TRUE') == TRUE){ 
      # does nothing
    }else{
      if (stringr::str_detect(comparewhat[5], 'FALSE') == TRUE){
        resCHOSEN = resCHOSEN[names(resCHOSEN) %in% ls(pattern = "_.*j") == FALSE]
      }else{
        if (stringr::str_detect(comparewhat[5], 'ONLY') == TRUE){
          resCHOSEN = resCHOSEN[names(resCHOSEN) %in% ls(pattern = "_.*j") == TRUE]
        }
      }
    }
    
    ##=##=##=##=##=#
    # checking TS
    if (stringr::str_detect(comparewhat[6], 'TRUE') == TRUE){
      # does nothing
    }else{
      if (stringr::str_detect(comparewhat[6], 'FALSE') == TRUE){
        resCHOSEN = resCHOSEN[names(resCHOSEN) %in% ls(pattern = "_.*TS") == FALSE]
      }else{
        if (stringr::str_detect(comparewhat[6], 'ONLY') == TRUE){
          resCHOSEN = resCHOSEN[names(resCHOSEN) %in% ls(pattern = "_.*TS") == TRUE]
        }
      }
    }
    
    ##=##=##=##=##=#
    # checking TRAIT
    if (stringr::str_detect(comparewhat[7], 'TRUE') == TRUE){
      # does nothing
    }else{
      if (stringr::str_detect(comparewhat[7], 'FALSE') == TRUE){
        resCHOSEN = resCHOSEN[names(resCHOSEN) %in% ls(pattern = "_.*trait") == FALSE]
      }else{
        if (stringr::str_detect(comparewhat[7], 'ONLY') == TRUE){
          resCHOSEN = resCHOSEN[names(resCHOSEN) %in% ls(pattern = "_.*trait") == TRUE]
        }
      }
    }
    
    ##=##=##=##=##=#
    # checking m
    if (stringr::str_detect(comparewhat[8], 'TRUE') == TRUE){
      # does nothing
    }else{
      if (stringr::str_detect(comparewhat[8], 'FALSE') == TRUE){
        resCHOSEN = resCHOSEN[names(resCHOSEN) %in% ls(pattern = "_.*m2") == FALSE] # check before if tree does not contain these characters
      }else{
        if (stringr::str_detect(comparewhat[8], 'ONLY') == TRUE){
          resCHOSEN = resCHOSEN[names(resCHOSEN) %in% ls(pattern = "_.*m2") == TRUE] # check before if tree does not contain these characters
        }
      }
    }
    
    ##=##=##=##=##=#
    # picking alt and null models
    NOPTEST = FALSE
    if (all(names(resCHOSEN) %in% ls(pattern = "_.*j")) == TRUE | any(names(resCHOSEN) %in% ls(pattern = "_.*j")) == FALSE){ # if all or none j
      if (all(names(resCHOSEN) %in% ls(pattern = "_.*trait")) == TRUE | any(names(resCHOSEN) %in% ls(pattern = "_.*trait")) == FALSE){ # if all or none trait
        if (all(names(resCHOSEN) %in% ls(pattern = "_.*m2")) == TRUE | any(names(resCHOSEN) %in% ls(pattern = "_.*m2")) == FALSE){ # if all or none m
          if (all(names(resCHOSEN) %in% ls(pattern = "_.*TS")) == TRUE | any(names(resCHOSEN) %in% ls(pattern = "_.*TS")) == FALSE){ # if all or none TS
            resALT = NULL
            resNULL = NULL
            base::message("WARNING: you DO NOT have an alternative model, as all your models\ninclude the j parameter and the TS matrix. For this, the\np value and the testable will not be generated.")
            NOPTEST = TRUE
          }else{ # if has some TS some not
            resALT = resCHOSEN[names(resCHOSEN) %in% ls(pattern = "_.*TS") == TRUE]
            resNULL = resCHOSEN[names(resCHOSEN) %in% ls(pattern = "_.*TS") == FALSE]
            base::message("See the models choosen for statistics")
            base::message(length(resCHOSEN), " models in total")
            print(names(resCHOSEN))
            base::message("See the models alternatives")
            base::message(length(resALT), " alternative models")
            print(names(resALT))
            base::message("See the models NULL")
            base::message(length(resNULL), " NULL models")
            print(names(resNULL))
            base::message("NOTE: you are comparing models with only or without +j+TRAIT+m.\nSo, for this, the alternative model considered here is the TS models")
          }
        }else{ # if has some m some not
          resALT = resCHOSEN[names(resCHOSEN) %in% ls(pattern = "_.*m2") == TRUE]
          resNULL = resCHOSEN[names(resCHOSEN) %in% ls(pattern = "_.*m2") == FALSE]
          base::message("See the models choosen for statistics")
          base::message(length(resCHOSEN), " models in total")
          print(names(resCHOSEN))
          base::message("See the models alternatives")
          base::message(length(resALT), " alternative models")
          print(names(resALT))
          base::message("See the models NULL")
          base::message(length(resNULL), " NULL models")
          print(names(resNULL))
          base::message("NOTE: you are comparing models with only or without +j+TRAIT.\nSo, for this, the alternative model considered here is the m models")
        }
      }else{ # if has some trait some not
        resALT = resCHOSEN[names(resCHOSEN) %in% ls(pattern = "_.*trait") == TRUE]
        resNULL = resCHOSEN[names(resCHOSEN) %in% ls(pattern = "_.*trait") == FALSE]
        base::message("See the models choosen for statistics")
        base::message(length(resCHOSEN), " models in total")
        print(names(resCHOSEN))
        base::message("See the models alternatives")
        base::message(length(resALT), " alternative models")
        print(names(resALT))
        base::message("See the models NULL")
        base::message(length(resNULL), " NULL models")
        print(names(resNULL))
        base::message("NOTE: you are comparing models with only or without +j.\nSo, for this, the alternative model considered here is the trait models")
      } 
    }else{ # if has some j some not
      resALT = resCHOSEN[names(resCHOSEN) %in% ls(pattern = "_.*j") == TRUE]
      resNULL = resCHOSEN[names(resCHOSEN) %in% ls(pattern = "_.*j") == FALSE]
      base::message("See the models choosen for statistics")
      base::message(length(resCHOSEN), " models in total")
      print(names(resCHOSEN))
      base::message("See the models alternatives")
      base::message(length(resALT), " alternative models")
      print(names(resALT))
      base::message("See the models NULL")
      base::message(length(resNULL), " NULL models")
      print(names(resNULL))
    }
    
    # testing p 
    restable = NULL
    teststable = NULL
    if (NOPTEST == FALSE){
      for (i in 1:(length(resCHOSEN)/2)){
        # i = 2
        tmp_tests = NULL
        
        LnL_2 = get_LnL_from_BioGeoBEARS_results_object(resALT[[i]])
        LnL_1 = get_LnL_from_BioGeoBEARS_results_object(resNULL[[i]])
        if (LnL_1 > 0){
          LnL_1 = -LnL_1
        }
        if (LnL_2 > 0){
          LnL_2 = -LnL_2
        }
        
        numparams1 = 3
        numparams2 = 2
        stats = AICstats_2models(LnL_2, LnL_1, numparams1, numparams2)
        stats
        
        
        res2 = NULL
        res1 = NULL
        # model for Likelihood Ratio Test (LRT)
        res2 = extract_params_from_BioGeoBEARS_results_object(results_object=resALT[[i]], returnwhat="table", addl_params=c("j"), paramsstr_digits=4)
        # alternative model for Likelihood Ratio Test (LRT)
        res1 = extract_params_from_BioGeoBEARS_results_object(results_object=resNULL[[i]], returnwhat="table", addl_params=c("j"), paramsstr_digits=4)
        
        # The null hypothesis for a Likelihood Ratio Test (LRT) is that two models
        # confer the same likelihood on the data. See: Brian O'Meara's webpage:
        # http://www.brianomeara.info/tutorials/aic
        # ...for an intro to LRT, AIC, and AICc
        
        rbind(res2, res1)
        tmp_tests = conditional_format_table(stats)
        
        restable = rbind(restable, res2, res1)
        teststable = rbind(teststable, tmp_tests)
      }
    }else{
      base::message("WARNING: AS STATED BEFORE, the P test and the testable will not be generated")
    }
  }
  
  #=#=#=#=#=#=#=#=#=#=#
  # namming models in tables
  #=#=#=#=#=#=#=#=#=#=#
  # naming the lines accordingly
  for (i in 1:(length(resCHOSEN)/2)){
    teststable$alt[[i]] = paste(names(resALT)[[i]], sep = "")
  }
  for (i in 1:(length(resCHOSEN)/2)){
    teststable$null[[i]] = paste(names(resNULL)[[i]], sep = "")
  }
  
  # naming the restble
  nameslnldf = NULL
  for (i in 1:nrow(restable)){
    #i = 1
    nameslnldfPH <- c(names(resCHOSEN)[i] , round(get_LnL_from_BioGeoBEARS_results_object(resCHOSEN[[i]]), 12))
    nameslnldf = rbind(nameslnldf, nameslnldfPH)
  }
  row.names(nameslnldf) <- NULL
  nameslnldf <- as.data.frame(nameslnldf)
  
  lnldf <- restable$LnL
  
  lnldf <- round(lnldf, 12)
  lnldf <- as.character(lnldf)
  
  restable$LnL <- lnldf
  
  class(restable$LnL)
  class(nameslnldf$V2)
  
  restable <- merge(restable, nameslnldf, by.x = "LnL", by.y = "V2")[,c(length(restable)+1, 1:length(restable))]
  colnames(restable)[1] <- "Model"
  restable$LnL <- -abs(as.numeric(restable$LnL))
  
  restable = put_jcol_after_ecol(restable)
  restable
  
  # Look at the results!!
  restable
  teststable
  
  #=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
  # Save the results tables for later -- check for e.g.
  # convergence issues
  #=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
  # ============================================================================
  # NOTE:
  # ============================================================================
  # These are the .RData for the tables. It is not necessary as we save the 
  # txt later. If you really want these, just uncomment it!
  #=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
  # # Loads to "restable"
  # save(restable, file=paste(dirfortableresults, "/", sppnameplot, "_restable_", comparewhat[1], "_MPN.Rdata", sep = ""))
  # load(file=paste(dirfortableresults, "/", sppnameplot, "_restable_", comparewhat[1], "_MPN.Rdata", sep = ""))
  # 
  # # Loads to "teststable"
  # save(teststable, file=paste(dirfortableresults, "/", sppnameplot, "_teststable_", comparewhat[1], "_MPN.Rdata", sep = ""))
  # load(file=paste(dirfortableresults, "/", sppnameplot, "_teststable_", comparewhat[1], "_MPN.Rdata", sep = ""))
  
  # Also save to text files
  write.table(restable, file=paste(dirfortableresults, "/", sppnameplot, "_restable_", comparewhat[1], "_MPN.txt", sep = ""), quote=FALSE, sep="\t", row.names = F)
  write.table(unlist_df(teststable), file=paste(dirfortableresults, "/", sppnameplot, "_teststable_", comparewhat[1], "_MPN.txt", sep = ""), quote=FALSE, sep="\t", row.names = F)
  
  
  #=#=#=#=#=#=#=#=#=#=#
  # AIC, AICc and weights
  #=#=#=#=#=#=#=#=#=#=#
  restable2 = restable
  
  # With AICs:
  AICtable = calc_AIC_column(LnL_vals=restable$LnL, nparam_vals=restable$numparams)
  restable = cbind(restable, AICtable)
  
  restable <- restable[order(restable$AIC, decreasing = FALSE),]
  
  restable_AIC_rellike = AkaikeWeights_on_summary_table(restable=restable, colname_to_use="AIC")
  restable_AIC_rellike = put_jcol_after_ecol(restable_AIC_rellike)
  restable_AIC_rellike
  
  # With AICcs -- factors in sample size
  samplesize = length(tr$tip.label)
  AICtable = calc_AICc_column(LnL_vals=restable$LnL, nparam_vals=restable$numparams, samplesize=samplesize)
  restable2 = cbind(restable2, AICtable)
  restable_AICc_rellike = AkaikeWeights_on_summary_table(restable=restable2, colname_to_use="AICc")
  restable_AICc_rellike = put_jcol_after_ecol(restable_AICc_rellike)
  restable_AICc_rellike
  
  # Also save to text files
  write.table(restable_AIC_rellike, file=paste(dirfortableresults, "/", sppnameplot, "_restable_AIC_rellike_", comparewhat[1], "_MPN.txt", sep = ""), quote=FALSE, sep="\t", row.names = F)
  write.table(restable_AICc_rellike, file=paste(dirfortableresults, "/", sppnameplot, "_restable_AICc_rellike_", comparewhat[1], "_MPN.txt", sep = ""), quote=FALSE, sep="\t", row.names = F)
  
  # Save with nice conditional formatting
  write.table(conditional_format_table(restable_AIC_rellike), file=paste(dirfortableresults, "/", sppnameplot, "_restable_AIC_rellike_formatted_", comparewhat[1], "_MPN.txt", sep = ""), quote=FALSE, sep="\t", row.names = F)
  write.table(conditional_format_table(restable_AICc_rellike), file=paste(dirfortableresults, "/", sppnameplot, "_restable_AICc_rellike_formatted_", comparewhat[1], "_MPN.txt", sep = ""), quote=FALSE, sep="\t", row.names = F)
  
  # code
  
  assign(paste0(sppnameplot,"_restable"), restable)
  assign(paste0(sppnameplot,"_teststable"), teststable)
  
  #==#==#==#==#==##==#==#==#==#==##==#==#==#==#==##==#==#==#==#==#
  ## ## ## 11 BSM  ####
  #==#==#==#==#==##==#==#==#==#==##==#==#==#==#==##==#==#==#==#==#
  dirforBSM <- paste0("../results/",sppnameplot, "/BSM")
  dir.create(dirforBSM, recursive = TRUE, showWarnings = FALSE)
  #=#=#=#=#=#=#=#=#=#=#
  # SETUP: Extension data directory
  #=#=#=#=#=#=#=#=#=#=#
  # When R packages contain extra files, they are stored in the "extdata" directory 
  # inside the installed package.
  #
  # BioGeoBEARS contains various example files and scripts in its extdata directory.
  # 
  # Each computer operating system might install BioGeoBEARS in a different place, 
  # depending on your OS and settings. 
  # 
  # However, you can find the extdata directory like this:
  extdata_dir = np(system.file("extdata", package="BioGeoBEARS"))
  extdata_dir
  list.files(extdata_dir)
  
  # Time-stratified files are here
  extdata_dir2 = np(slashslash(paste(extdata_dir, "examples/Psychotria_M3strat/BGB/", sep="/")))
  
  #=#=#=#=#=#=#=#=#=#=#
  # Phylogeny file
  #=#=#=#=#=#=#=#=#=#=#
  # Look at the raw Newick file:
  moref(trfn)
  
  # Look at your phylogeny:
  tr = read.tree(trfn)
  tr
  
  #=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#
  # Geography file
  #=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#
  # Look at the raw geography text file:
  moref(geogfn)
  
  # Look at your geographic range data:
  tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=geogfn)
  tipranges
  
  # Set the maximum number of areas any species may occupy; this cannot be larger 
  # than the number of areas you set up, but it can be smaller.
  max_range_size = 3
  
  #=#=#=#=#=#=#=#=#=#=#
  # Pick your model name:
  #=#=#=#=#=#=#=#=#=#=#
  # read least table with AIC results
  cond_restable_AICc_rellike<-read.csv(paste(dirfortableresults, "/", sppnameplot, "_restable_AIC_rellike_formatted_", comparewhat[1], "_MPN.txt", sep=""), sep="\t")
  
  # select best model
  rownames(cond_restable_AICc_rellike) <- cond_restable_AICc_rellike$Model
  best_model<-rownames(cond_restable_AICc_rellike)[cond_restable_AICc_rellike[,"AIC"]==min(cond_restable_AICc_rellike[,"AIC"])]
  
  res = get(ls(pattern = paste0(best_model, "$")))
  model_name = best_model
  model_name
  
  res$inputs$geogfn = geogfn
  res$inputs$trfn = trfn
  
  
  #=#=#=#=#=#=#=#=#=#=#
  # stochastic mapping
  #=#=#=#=#=#=#=#=#=#=#
  clado_events_tables = NULL
  ana_events_tables = NULL
  lnum = 0
  
  # initial config
  res$inputs$num_cores_to_use = detectCores()
  res$inputs$force_sparse = FALSE
  res$inputs$use_optimx = "GenSA"
  BSM_inputs_fn = paste0(dirforBSM, "/", sppnameplot, "_BSM_inputs_file.Rdata")
  BSM_inputs_fn = paste0(dirforBSM, "/", sppnameplot, "BSM_inputs_file_ERY.Rdata")
  closeAllConnections()
  runInputsSlow = TRUE
  if (runInputsSlow){
    stochastic_mapping_inputs_list = get_inputs_for_stochastic_mapping(res=res)
    save(stochastic_mapping_inputs_list, file=BSM_inputs_fn)
  } else {
    # Loads to "stochastic_mapping_inputs_list"
    load(BSM_inputs_fn)
  } # END if (runInputsSlow)
  
  # Check inputs (doesn't work the same on unconstr)
  # View(stochastic_mapping_inputs_list)
  summary(stochastic_mapping_inputs_list)
  names(stochastic_mapping_inputs_list[[1]])
  
  runBSMslow = TRUE
  
  if (runBSMslow == TRUE){
    BSM_output3 = runBSM(res, 
                         stochastic_mapping_inputs_list=stochastic_mapping_inputs_list, 
                         maxnum_maps_to_try=10000, nummaps_goal=1000, 
                         maxtries_per_branch=80000, save_after_every_try=TRUE, 
                         savedir=dirforBSM, seedval=12345, 
                         wait_before_save=0.01)
    
    RES_clado_events_tables = BSM_output3$RES_clado_events_tables
    RES_ana_events_tables = BSM_output3$RES_ana_events_tables
  } else {
    # Load previously saved...
    base::message("You need to setup for loading the data here")
  } 
  
  # Extract BSM output
  clado_events_tables = BSM_output3$RES_clado_events_tables
  ana_events_tables = BSM_output3$RES_ana_events_tables
  head(clado_events_tables[[1]])
  head(ana_events_tables[[1]])
  length(clado_events_tables)
  length(ana_events_tables)
  
  # if you want to sabe the .RDATA, uncomment this
  # save(BSM_output3, file = paste0(dirforBSM, "/", sppnameplot, "_BSM_output.RData"))
  
  #=#=#=#=#=#=#=#=#=#=#
  # Plot one stochastic map, manual metho
  #=#=#=#=#=#=#=#=#=#=#
  # (we have to convert the stochastic maps into event
  #  maps for plotting)
  
  #=#=#=#=#=#=#=#=#=#=#=#=
  # Get the color scheme
  #=#=#=#=#=#=#=#=#=#=#=#=
  
  include_null_range = TRUE
  areanames = names(tipranges@df)
  areas = areanames
  max_range_size = max_range_size
  
  # Note: If you did something to change the states_list from the default given the number of areas, you would
  # have to manually make that change here as well! (e.g., areas_allowed matrix, or manual reduction of the states_list)
  states_list_0based = rcpp_areas_list_to_states_list(areas=areas, maxareas=max_range_size, include_null_range=include_null_range)
  
  colors_list_for_states = get_colors_for_states_list_0based(areanames=areanames, states_list_0based=states_list_0based, max_range_size=max_range_size, plot_null_range=TRUE)
  
  #=#=#=#=#=#=#=#=#=#=#=#=
  # Setup for painting a single stochastic map
  #=#=#=#=#=#=#=#=#=#=#=#=
  scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))
  stratified=FALSE
  clado_events_table = clado_events_tables[[1]]
  ana_events_table = ana_events_tables[[1]]
  
  #=#=#=#=#=#=#=#=#=#=#=#=
  # Open a PDF
  #=#=#=#=#=#=#=#=#=#=#=#=
  pdffn = paste0(dirforBSM, "/", model_name, "_single_stochastic_map_n1.pdf")
  pdf(file=pdffn, width=6, height=6)
  
  # Convert the BSM into a modified res object
  master_table_cladogenetic_events = clado_events_tables[[1]]
  resmod = stochastic_map_states_into_res(res=res, master_table_cladogenetic_events=master_table_cladogenetic_events, stratified=stratified)
  
  plot_BioGeoBEARS_results(results_object=resmod, analysis_titletxt="Stochastic map", addl_params=list("j"), label.offset=0.5, plotwhat="text", cornercoords_loc=scriptdir, root.edge=TRUE, colors_list_for_states=colors_list_for_states, skiptree=FALSE, show.tip.label=TRUE)
  
  # Paint on the branch states
  paint_stochastic_map_branches(res=resmod, master_table_cladogenetic_events=master_table_cladogenetic_events, colors_list_for_states=colors_list_for_states, lwd=5, lty=par("lty"), root.edge=TRUE, stratified=stratified)
  
  plot_BioGeoBEARS_results(results_object=resmod, analysis_titletxt="Stochastic map", addl_params=list("j"), plotwhat="text", cornercoords_loc=scriptdir, root.edge=TRUE, colors_list_for_states=colors_list_for_states, skiptree=TRUE, show.tip.label=TRUE)
  
  #=#=#=#=#=#=#=#=#=#=#=#=
  # Close PDF
  #=#=#=#=#=#=#=#=#=#=#=#=
  dev.off()
  cmdstr = paste("open ", pdffn, sep="")
  system(cmdstr)
  
  #=#=#=#=#=#=#=#=#=#=#=#=
  # Plot all 1000 stochastic maps to PDF
  #=#=#=#=#=#=#=#=#=#=#=#=
  # Setup
  include_null_range = include_null_range
  areanames = areanames
  areas = areanames
  max_range_size = max_range_size
  states_list_0based = rcpp_areas_list_to_states_list(areas=areas, maxareas=max_range_size, include_null_range=include_null_range)
  colors_list_for_states = get_colors_for_states_list_0based(areanames=areanames, states_list_0based=states_list_0based, max_range_size=max_range_size, plot_null_range=TRUE)
  scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))
  stratified = stratified
  
  # Loop through the maps and plot to PDF
  pdffn = paste0(dirforBSM, "/", model_name, "_", length(clado_events_tables), "BSMs_v1.pdf")
  pdf(file=pdffn, width=6, height=6)
  
  nummaps_goal = length(clado_events_tables)
  for (i in 1:nummaps_goal){
    clado_events_table = clado_events_tables[[i]]
    analysis_titletxt = paste0(model_name, " - Stochastic Map #", i, "/", nummaps_goal)
    plot_BSM(results_object=res, clado_events_table=clado_events_table, stratified=stratified, analysis_titletxt=analysis_titletxt, addl_params=list("j"), label.offset=0.5, plotwhat="text", cornercoords_loc=scriptdir, root.edge=TRUE, colors_list_for_states=colors_list_for_states, show.tip.label=TRUE, include_null_range=include_null_range)
  } # END for (i in 1:nummaps_goal)
  
  dev.off()
  cmdstr = paste("open ", pdffn, sep="")
  system(cmdstr)
  
  #=#=#=#=#=#=#=#=#=#=#=#=
  # summarize stochastic map tables
  #=#=#=#=#=#=#=#=#=#=#=#=
  length(clado_events_tables)
  length(ana_events_tables)
  
  head(clado_events_tables[[1]][,-20])
  tail(clado_events_tables[[1]][,-20])
  
  head(ana_events_tables[[1]])
  tail(ana_events_tables[[1]])
  
  areanames = names(tipranges@df)
  actual_names = areanames
  actual_names
  
  # Get the dmat and times (if any)
  dmat_times = get_dmat_times_from_res(res=res, numstates=NULL)
  dmat_times
  
  # Extract BSM output
  BSM_output <- BSM_output3
  clado_events_tables = BSM_output$RES_clado_events_tables
  ana_events_tables = BSM_output$RES_ana_events_tables
  
  # Simulate the source areas
  BSMs_w_sourceAreas = simulate_source_areas_ana_clado(res, clado_events_tables, ana_events_tables, areanames)
  clado_events_tables = BSMs_w_sourceAreas$clado_events_tables
  ana_events_tables = BSMs_w_sourceAreas$ana_events_tables
  
  # printing to file all dispersal
  skipskip = FALSE
  ana_dis_ALL = list()
  for (i in 1:length(ana_events_tables)){ # i = 131
    if (!is.data.frame(ana_events_tables[[i]])){
      next
    }
    placeholder <- ana_events_tables[[i]][ana_events_tables[[i]]$event_type == "d",]
    if (nrow(placeholder) == 0){
      next
    }
    placeholder$nmap <- as.character(i)
    placeholder <- placeholder[,c("node", "node.type", "current_rangetxt", "event_type", "event_txt", "ana_dispersal_from", "dispersal_to", "event_time", "abs_event_time", "nmap")]
    ana_dis_ALL[[i]] <- placeholder
  }
  length(ana_dis_ALL)
  if(length(ana_dis_ALL) == 0){
    placeholder <- data.frame(node = NA, node.type = NA, current_rangetxt = NA, event_type = NA, event_txt = NA, ana_dispersal_from = NA, dispersal_to = NA, event_time = NA, abs_event_time = NA, nmap = NA)
    ana_dis_ALL_Table <- placeholder
    skipskip = TRUE
  }else{
    ana_dis_ALL_Table <- rbindlist(ana_dis_ALL)
    hist(ana_dis_ALL_Table$abs_event_time)
    skipskip = FALSE
  }
  
  # saving to file
  write.csv(ana_dis_ALL_Table, paste0(dirforBSM, "/", sppnameplot, "_ana_ALL_MASTER_table.txt"), quote = F)
  
  # adding to master table
  if (!skipskip){
    ana_dis_ALL_Table_MASTER <- rbind(ana_dis_ALL_Table_MASTER, ana_dis_ALL_Table) 
  }
  
  # printing to file all jump
  skipskip = FALSE
  clado_dis_ALL = list()
  for (i in 1:length(clado_events_tables)){ # i = 131
    if (!is.data.frame(clado_events_tables[[i]])){
      next
    }
    placeholder <- clado_events_tables[[i]][clado_events_tables[[i]]$clado_event_type == "founder (j)",]
    if (nrow(placeholder) == 0){
      next
    }
    placeholder$nmap <- as.character(i)
    placeholder <- placeholder[,c("node", "node.type", "sampled_states_AT_nodes", "clado_event_type", "clado_event_txt", "clado_dispersal_from", "clado_dispersal_to", "time_bp", "left_desc_nodes", "right_desc_nodes", "nmap")]
    clado_dis_ALL[[i]] <- placeholder
  }
  length(clado_dis_ALL)
  if(length(clado_dis_ALL) == 0){
    placeholder <- data.frame(node = NA, node.type = NA, sampled_states_AT_nodes = NA, clado_event_type = NA, clado_event_txt = NA, clado_dispersal_from = NA, clado_dispersal_to = NA, time_bp = NA, left_desc_nodes = NA, right_desc_nodes = NA, nmap = NA)
    clado_dis_ALL_Table <- placeholder
    skipskip = TRUE
  }else{
    clado_dis_ALL_Table <- rbindlist(clado_dis_ALL)
    hist(clado_dis_ALL_Table$time_bp)
    skipskip = FALSE
  }
  
  # saving to file
  write.csv(clado_dis_ALL_Table, paste0(dirforBSM, "/", sppnameplot, "_clado_ALL_MASTER_table.txt"), quote = F)
  
  # adding to master table
  if (!skipskip){
    clado_dis_ALL_Table_MASTER <- rbind(clado_dis_ALL_Table_MASTER, clado_dis_ALL_Table)
  }
  
  
  #==#==#==#==#==##==#==#==#==#==##==#==#==#==#==##==#==#==#==#==#
  ## ## ## 12 PLOTING BSM FOR EACH GENERA ####
  #==#==#==#==#==##==#==#==#==#==##==#==#==#==#==##==#==#==#==#==#
  # Prepare Anagenetic Data
  df_ana <- ana_dis_ALL_Table %>%
    dplyr::mutate(
      route = paste0(ana_dispersal_from, " -> ", dispersal_to),
      type = "Anagenetic",
      time = abs_event_time
    ) %>%
    dplyr::select(time, type, route)
  
  # Prepare Cladogenetic Data
  df_clado <- clado_dis_ALL_Table %>%
    dplyr::mutate(
      route = paste0(clado_dispersal_from, " -> ", clado_dispersal_to),
      type = "Cladogenetic",
      time = time_bp
    ) %>%
    dplyr::select(time, type, route)
  
  # Combine them and remove disp between A and B (not important)
  df_plot_routes <- rbind(df_ana, df_clado)
  df_plot_routes <- df_plot_routes[!df_plot_routes$route == "A -> B" & !df_plot_routes$route == "B -> A", ]
  
  # Calculate stats PER ROUTE
  route_stats <- df_plot_routes %>%
    dplyr::group_by(type, route) %>%
    dplyr::summarise(
      mu = mean(time, na.rm = TRUE),
      sd = sd(time, na.rm = TRUE),
      .groups = 'drop' 
    )
  
  # defining upper_limit
  upper_limit <- ceiling(max(df_plot_routes$time, na.rm = TRUE) / 5) * 5
  upper_limit <- max(upper_limit, 5)
  
  # change names of the faces in plot:
  facet_names <- c(
    "Anagenetic"   = "Anagenetic Dispersal (d)", 
    "Cladogenetic" = "Founder-Event Speciation (j)",
    "A -> C"       = "Above to Below (A > C)",
    "B -> C"       = "Panama to Bellow (B > C)",
    "C -> A"       = "Below to Above (C > A)",
    "C -> B"       = "Below to Panama (C > B)"
  )
  
  # Mid point for plotting
  route_mid <- df_plot_routes %>%
    dplyr::group_by(type, route) %>%
    do({
      h <- hist(.$time, breaks = seq(0, upper_limit + 1, by = 1), plot = FALSE)
      data.frame(y_mid = max(h$counts) / 2)
    }) %>%
    dplyr::left_join(route_stats, by = c("type", "route"))
  
  # Calculate total N for each route
  label_counts <- df_plot_routes %>%
    dplyr::group_by(type, route) %>%
    dplyr::summarise(total_n = dplyr::n(), .groups = 'drop')
  
  # filtering the epochs for plot
  clean_labels <- epochs
  if (upper_limit > 10){
    clean_labels$abbr[clean_labels$abbr %in% c("H", "Ple", "Pli")] <- "" 
  }
  
  # now the plot
  pgraph <- ggplot(df_plot_routes, aes(x = time, fill = type)) +
    # Geologic Time Scale
    coord_geo(dat = clean_labels, pos = "bottom", abbrv = TRUE, height = unit(1, "line"), xlim = c(upper_limit, 0)) +
    
    # Data Layers
    geom_histogram(aes(y = after_stat(count)), binwidth = 1, alpha = 0.6, color = "white") +
    geom_density(aes(y = after_stat(count)), alpha = 0.2) + 
    
    # 5 Ma Red Line
    geom_vline(xintercept = 5, color = "red", linetype = "dashed", linewidth = 0.8) +
    
    # Mean and SD indicators (using the dynamic y_mid)
    geom_errorbarh(data = route_mid, 
                   aes(xmin = mu - sd, xmax = mu + sd, y = y_mid, x = mu), 
                   height = 0.5, color = "black", inherit.aes = FALSE) +
    geom_point(data = route_mid, 
               aes(x = mu, y = y_mid), 
               color = "white", fill = "black", shape = 21, size = 2, inherit.aes = FALSE) +
    
    # Add the count box inside each panel
    geom_label(data = label_counts, 
               aes(label = paste0("n = ", total_n)),
               x = Inf, y = Inf,          # Position at the top-right corner
               hjust = 1.1, vjust = 1.1,  # Offset slightly so it doesn't touch the borders
               fill = "white", alpha = 0.8, 
               size = 3, fontface = "bold",
               inherit.aes = FALSE) +
    
    # Use facet_grid for better use of space with many routes
    facet_wrap(type ~ route, 
               scales = "free_y", 
               ncol = 4,             # Matches your 4 columns
               strip.position = "top", 
               labeller = as_labeller(facet_names)) +
    
    # Styling and Scales
    scale_x_continuous(breaks = seq(0, upper_limit, by = 5)) +
    scale_fill_viridis_d(option = "mako", begin = 0.3, end = 0.7) +
    theme_bw() +
    labs(
      title = paste0("Dispersal Frequency by Route: ", sppnameplot),
      subtitle = "Red line at 5 Ma. Points indicate Mean ± SD.",
      x = "Time (Ma)",
      y = "Frequency (Events per Route)"
    ) +
    theme(
      legend.position = "none",
      strip.background = element_rect(fill = "gray95"),
      strip.text = element_text(face = "bold", size = 9),
      panel.grid.minor = element_blank()
    )
  
  # saving to pdf:
  pdf(paste0("../figures/", sppnameplot, "_graph_results.pdf"), height = 7, width = 15)
  par(mar = c(0,0,0,0))
  print(pgraph)
  dev.off()
  
  #=#=#=#=#=#=#=#=#=#=#=#=
  # Other stuff
  #=#=#=#=#=#=#=#=#=#=#=#=
  
  # Count all anagenetic and cladogenetic events
  counts_list = count_ana_clado_events(clado_events_tables, ana_events_tables, areanames, actual_names)
  
  summary_counts_BSMs = counts_list$summary_counts_BSMs
  print(conditional_format_table(summary_counts_BSMs))
  
  
  dir.create(paste0(dirforBSM, "/out"), recursive = TRUE, showWarnings = FALSE)
  # Histogram of event counts
  hist_event_counts(counts_list, pdffn=paste0(dirforBSM, "/out/", model_name, "_histograms_of_event_counts.pdf"))
  
  #=#=#=#=#=#=#=#=#=#=#=#=
  # Print counts to files 
  #=#=#=#=#=#=#=#=#=#=#=#=
  tmpnames = names(counts_list)
  cat("\n\nWriting tables* of counts to tab-delimited text files:\n(* = Tables have dimension=2 (rows and columns). Cubes (dimension 3) and lists (dimension 1) will not be printed to text files.) \n\n")
  for (i in 1:length(tmpnames)){
    cmdtxt = paste0("item = counts_list$", tmpnames[i])
    eval(parse(text=cmdtxt))
    
    # Skip cubes
    if (length(dim(item)) != 2){
      next()
    }
    
    outfn = paste0(dirforBSM, "/out/", model_name, "_", tmpnames[i], ".txt")
    if (length(item) == 0){
      cat(outfn, " -- NOT written, *NO* events recorded of this type", sep="")
      cat("\n")
    } else {
      cat(outfn)
      cat("\n")
      write.table(conditional_format_table(item), file=outfn, quote=FALSE, sep="\t", col.names=TRUE, row.names=TRUE)
    } # END if (length(item) == 0)
  } # END for (i in 1:length(tmpnames))
  cat("...done.\n")
  
  #=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#
  # Check that ML ancestral state/range probabilities and
  # the mean of the BSMs approximately line up
  #=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#
  
  withr::with_dir(file.path(dirforBSM, "out"), {
    check_ML_vs_BSM(res, clado_events_tables, model_name, tr = NULL,
                    plot_each_node = FALSE, linreg_plot = TRUE, MultinomialCI = TRUE)
  })
  
  ##=# save ws ##=#
  # Uncomment this if you want to save the .RData
  # base::save.image(paste0("../results/", sppnameplot, "/", sppnameplot, ".RData"))
  #=#=#=#=#=#=#=#=#=# 
}
#==#==#==#==#==##==#==#==#==#==##==#==#==#==#==##==#==#==#==#==#
## ## ## 13 PLOTING BSM FOR ALL ####
#==#==#==#==#==##==#==#==#==#==##==#==#==#==#==##==#==#==#==#==#
# IMPORTANT: only run after the loop with all genera
# setup
sppnameplot <- "ALL"

# Prepare Anagenetic Data
df_ana <- ana_dis_ALL_Table_MASTER %>%
  dplyr::mutate(
    route = paste0(ana_dispersal_from, " -> ", dispersal_to),
    type = "Anagenetic",
    time = abs_event_time
  ) %>%
  dplyr::select(time, type, route)

# Prepare Cladogenetic Data
df_clado <- clado_dis_ALL_Table_MASTER %>%
  dplyr::mutate(
    route = paste0(clado_dispersal_from, " -> ", clado_dispersal_to),
    type = "Cladogenetic",
    time = time_bp
  ) %>%
  dplyr::select(time, type, route)

# Combine them and remove disp between A and B (not important)
df_plot_routes <- rbind(df_ana, df_clado)
df_plot_routes <- df_plot_routes[!df_plot_routes$route == "A -> B" & !df_plot_routes$route == "B -> A", ]

# Calculate stats PER ROUTE
route_stats <- df_plot_routes %>%
  dplyr::group_by(type, route) %>%
  dplyr::summarise(
    mu = mean(time, na.rm = TRUE),
    sd = sd(time, na.rm = TRUE),
    .groups = 'drop' 
  )

# change names of the faces in plot:
facet_names <- c(
  "Anagenetic"   = "Anagenetic Dispersal (d)", 
  "Cladogenetic" = "Founder-Event Speciation (j)",
  "A -> C"       = "Above to Below (A > C)",
  "B -> C"       = "Panama to Bellow (B > C)",
  "C -> A"       = "Below to Above (C > A)",
  "C -> B"       = "Below to Panama (C > B)"
)

upper_limit <- ceiling(max(df_plot_routes$time, na.rm = TRUE) / 5) * 5
upper_limit <- max(upper_limit, 5)

# Mid point for plotting
route_mid <- df_plot_routes %>%
  dplyr::group_by(type, route) %>%
  do({
    h <- hist(.$time, breaks = seq(0, upper_limit + 1, by = 1), plot = FALSE)
    data.frame(y_mid = max(h$counts) / 2)
  }) %>%
  dplyr::left_join(route_stats, by = c("type", "route"))

# Calculate total N for each route
label_counts <- df_plot_routes %>%
  dplyr::group_by(type, route) %>%
  dplyr::summarise(total_n = dplyr::n(), .groups = 'drop')

# filtering the period for plot
clean_labels <- epochs
if (upper_limit > 10){
  clean_labels$abbr[clean_labels$abbr %in% c("H", "Ple", "Pli")] <- "" 
}

# now the plot
pgraph <- ggplot(df_plot_routes, aes(x = time, fill = type)) +
  # Geologic Time Scale
  coord_geo(dat = clean_labels, pos = "bottom", abbrv = TRUE, height = unit(1, "line"), xlim = c(upper_limit, 0)) +
  
  # Data Layers
  geom_histogram(aes(y = after_stat(count)), binwidth = 1, alpha = 0.6, color = "white") +
  geom_density(aes(y = after_stat(count)), alpha = 0.2) + 
  
  # 5 Ma Red Line
  geom_vline(xintercept = 5, color = "red", linetype = "dashed", linewidth = 0.8) +
  
  # Mean and SD indicators (using the dynamic y_mid)
  geom_errorbarh(data = route_mid, 
                 aes(xmin = mu - sd, xmax = mu + sd, y = y_mid, x = mu), 
                 height = 0.5, color = "black", inherit.aes = FALSE) +
  geom_point(data = route_mid, 
             aes(x = mu, y = y_mid), 
             color = "white", fill = "black", shape = 21, size = 2, inherit.aes = FALSE) +
  
  # Add the count box inside each panel
  geom_label(data = label_counts, 
             aes(label = paste0("n = ", total_n)),
             x = Inf, y = Inf,          # Position at the top-right corner
             hjust = 1.1, vjust = 1.1,  # Offset slightly so it doesn't touch the borders
             fill = "white", alpha = 0.8, 
             size = 3, fontface = "bold",
             inherit.aes = FALSE) +
  
  # Use facet_grid for better use of space with many routes
  facet_wrap(type ~ route, 
             scales = "free_y", 
             ncol = 4,             # Matches your 4 columns
             strip.position = "top", 
             labeller = as_labeller(facet_names)) +
  
  # Styling and Scales
  scale_x_continuous(breaks = seq(0, upper_limit, by = 5)) +
  scale_fill_viridis_d(option = "mako", begin = 0.3, end = 0.7) +
  theme_bw() +
  labs(
    title = paste0("Dispersal Frequency by Route for All Genera"),
    subtitle = "Red line at 5 Ma. Points indicate Mean ± SD.",
    x = "Time (Ma)",
    y = "Frequency (Events per Route)"
  ) +
  theme(
    legend.position = "none",
    strip.background = element_rect(fill = "gray95"),
    strip.text = element_text(face = "bold", size = 9),
    panel.grid.minor = element_blank()
  )

# saving to pdf:
pdf(paste0("../figures/", sppnameplot, "_graph_results.pdf"), height = 7, width = 15)
par(mar = c(0,0,0,0))
print(pgraph)
dev.off()

#==#==#==#==#==##==#==#==#==#==##==#==#==#==#==##==#==#==#==#==#
## ## ## 14 Only gathering information for manuscript ####
#==#==#==#==#==##==#==#==#==#==##==#==#==#==#==##==#==#==#==#==#
treesALL <- list()
for (i in 1:length(sppnameplotPH)){
  sppspp <- sppnameplotPH[i]
  phtree <- read.tree(paste0("../data/biogeobears/", sppspp, "/", sppspp, "_tree_final.newick"))
  treesALL[[i]] <- phtree
}

# grabing n collapsed for each
tiplabelsALLPH <- sapply(treesALL, function (x) x$tip.label)
names(tiplabelsALLPH) <- sppnameplotPH

# unlisting for all
tiplabelsALL <- unlist(tiplabelsALLPH)

# checking numbers
sapply(tiplabelsALLPH, length) # for each
length(tiplabelsALL) # for all

#==#==#==#==#==##==#==#==#==#==##==#==#==#==#==##==#==#==#==#==#
## ## ## 15 Main Figure ####
#==#==#==#==#==##==#==#==#==#==##==#==#==#==#==##==#==#==#==#==#
# ------------------------------------------------------------------------------
# 1. USER INPUTS & SETTINGS
# ------------------------------------------------------------------------------
# This block only asserts it exists, so the section can still be run on its own.
stopifnot(exists("sppnameplotPH"), length(sppnameplotPH) == 8)

# Format geologic time scale labels (cleaning recent short epochs for aesthetics)
clean_labels <- deeptime::epochs
clean_labels$abbr[clean_labels$abbr %in% c("H", "Ple", "Pli")] <- ""
clean_labels2 <- deeptime::epochs # for timescales < 15 mya

# Define fixed colors to ensure consistency across all plots
event_colors <- c(
  "Anagenetic"   = "#D81B60", # Red-pink-ish
  "Cladogenetic" = "#004D40"  # Green
)

# Dictionary to map Route + Type to clean facet labels.
# This forces the facets to display only the route names.
combo_names <- c(
  "Anagenetic A -> C"   = "Above to Below (A > C)",
  "Anagenetic B -> C"   = "Panama to Bellow (B > C)",
  "Anagenetic C -> A"   = "Below to Above (C > A)",
  "Anagenetic C -> B"   = "Below to Panama (C > B)",
  "Cladogenetic A -> C" = "Above to Below (A > C)",
  "Cladogenetic B -> C" = "Panama to Bellow (B > C)",
  "Cladogenetic C -> A" = "Below to Above (C > A)",
  "Cladogenetic C -> B" = "Below to Panama (C > B)"
)

# ------------------------------------------------------------------------------
# 2. EXTRACT AND PROCESS DATA FROM BSM
# ------------------------------------------------------------------------------
general_list_df_clado <- list()
general_list_df_ana <- list()

for (spp in sppnameplotPH) {
  
  # Clear specific objects to prevent data bleeding between iterations
  rm(list = intersect(ls(), c("RES_clado_events_tables", "RES_ana_events_tables", "res")))
  
  path_bsm <- paste0("../results/", spp, "/BSM/")
  path_res <- paste0("../results/", spp, "/")
  
  # Load BSM results
  load(paste0(path_bsm, "RES_clado_events_tables.RData"))
  load(paste0(path_bsm, "RES_ana_events_tables.RData"))
  
  # Identify the best model based on AIC
  aic_path <- paste0(path_res, spp, "_restable_AIC_rellike_formatted_ALL_MPN.txt")
  cond_restable <- read.csv(aic_path, sep="\t")
  rownames(cond_restable) <- cond_restable$Model
  best_model_raw <- rownames(cond_restable)[cond_restable[,"AIC"] == min(cond_restable[,"AIC"])]
  best_model <- str_remove(best_model_raw[1], "res") 
  
  # Load ML results and simulate source areas
  load(paste0(path_res, best_model, "_res_MPN.RData"))
  areanames <- c("A", "B", "C")
  BSMs_w_sourceAreas <- simulate_source_areas_ana_clado(res, RES_clado_events_tables, RES_ana_events_tables, areanames)
  
  clado_list <- BSMs_w_sourceAreas$clado_events_tables
  ana_list   <- BSMs_w_sourceAreas$ana_events_tables
  
  # Process Cladogenetic (Founder) Events
  clado_processed <- imap(clado_list, function(x, idx) {
    if(!is.data.frame(x) || nrow(x) == 0) return(NULL)
    sub_df <- x[x$clado_event_type == "founder (j)", ]
    if(nrow(sub_df) > 0) {
      sub_df <- sub_df[, c("node", "clado_event_type", "clado_dispersal_from", "clado_dispersal_to", "time_bp")]
      sub_df$nmap   <- as.character(idx)
      sub_df$genera <- spp
      return(sub_df)
    }
    return(NULL)
  })
  general_list_df_clado[[spp]] <- bind_rows(clado_processed)
  
  # Process Anagenetic (Dispersal) Events
  ana_processed <- imap(ana_list, function(x, idx) {
    if(!is.data.frame(x) || nrow(x) == 0) return(NULL)
    sub_df <- x[x$event_type == "d", ]
    if(nrow(sub_df) > 0) {
      sub_df <- sub_df[, c("node", "event_type", "ana_dispersal_from", "dispersal_to", "abs_event_time")]
      sub_df$nmap   <- as.character(idx)
      sub_df$genera <- spp
      return(sub_df)
    }
    return(NULL)
  })
  general_list_df_ana[[spp]] <- bind_rows(ana_processed)
  
  base::message("Finished processing genus: ", spp)
}

# ------------------------------------------------------------------------------
# 3. PREPARE FINAL PLOTTING DATAFRAME
# ------------------------------------------------------------------------------
clado_dis_ALL_Table_MASTER <- bind_rows(general_list_df_clado)
ana_dis_ALL_Table_MASTER   <- bind_rows(general_list_df_ana)

# Format Anagenetic data
df_ana <- ana_dis_ALL_Table_MASTER %>%
  dplyr::mutate(
    route = paste0(ana_dispersal_from, " -> ", dispersal_to),
    type = "Anagenetic",
    time = abs_event_time
  ) %>%
  dplyr::select(time, type, route, genera)

# Format Cladogenetic data
df_clado <- clado_dis_ALL_Table_MASTER %>%
  dplyr::mutate(
    route = paste0(clado_dispersal_from, " -> ", clado_dispersal_to),
    type = "Cladogenetic",
    time = time_bp
  ) %>%
  dplyr::select(time, type, route, genera)

# Combine both event types
df_plot_routes <- rbind(df_ana, df_clado)

# Remove NA values and specific non-target routes (e.g., A -> B)
df_plot_routes <- df_plot_routes %>% 
  filter(!is.na(time) & !is.na(route)) %>%
  filter(!route %in% c("A -> B", "B -> A"))

# ------------------------------------------------------------------------------
# 4. DYNAMIC PLOTTING FUNCTION
# ------------------------------------------------------------------------------
create_dispersal_plot <- function(data_subset, plot_title, show_legend = FALSE) {
  
  # Calculate dynamic upper limit for the X-axis based on the maximum time in this subset
  # Rounds up to the nearest multiple of 5 (e.g., 12 Ma becomes 15 Ma)
  max_time <- max(data_subset$time, na.rm = TRUE)
  if(is.infinite(max_time) || is.na(max_time)) max_time <- 0
  dyn_limit <- max(ceiling(max_time / 5) * 5, 5) 
  
  clean_labels_what <- if (dyn_limit > 15) clean_labels else clean_labels2
  
  # Create a unique panel ID to fix the grid order and apply custom facet labels
  data_subset <- data_subset %>%
    dplyr::mutate(panel_id = factor(paste(type, route), levels = names(combo_names)))
  
  # Calculate route statistics (Mean and SD)
  route_stats <- data_subset %>%
    dplyr::group_by(type, route, panel_id) %>%
    dplyr::summarise(
      mu = mean(time, na.rm = TRUE),
      sd = sd(time, na.rm = TRUE),
      .groups = 'drop' 
    )
  
  # Calculate mid-point frequency for plotting the mean/SD indicators dynamically
  route_mid <- data_subset %>%
    dplyr::group_by(type, route, panel_id) %>%
    do({
      if(nrow(.) > 1) {
        h <- hist(.$time, breaks = seq(0, dyn_limit + 1, by = 1), plot = FALSE)
        data.frame(y_mid = max(h$counts) / 2)
      } else {
        data.frame(y_mid = 0.5) 
      }
    }) %>%
    dplyr::left_join(route_stats, by = c("type", "route", "panel_id"))
  
  # Calculate sample size (n) AND format the label text
  label_counts <- data_subset %>%
    dplyr::group_by(type, route, panel_id) %>%
    dplyr::summarise(
      total_n = dplyr::n(),
      mu_val  = mean(time, na.rm = TRUE),
      sd_val  = sd(time, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    dplyr::mutate(
      # Creates the formatted string, e.g.:
      # n = 100
      # 5.2 ± 1.1 Ma
      label_text = sprintf("n = %d\n%.1f \u00B1 %.1f Ma", total_n, mu_val, sd_val)
    )
  
  # Build the base plot
  p <- ggplot(data_subset, aes(x = time, fill = type)) +
    coord_geo(dat = clean_labels_what, pos = "bottom", abbrv = TRUE, height = unit(1, "line"), xlim = c(dyn_limit, 0)) +
    geom_histogram(aes(y = after_stat(count)), binwidth = 1, alpha = 0.6, color = "white") +
    geom_density(aes(y = after_stat(count)), alpha = 0.2) + 
    geom_vline(xintercept = 5, color = "red", linetype = "dashed", linewidth = 0.8) +
    
    # Statistical indicators (Points and bars)
    geom_errorbarh(data = route_mid, aes(xmin = mu - sd, xmax = mu + sd, y = y_mid, x = mu), height = 0.5, color = "black", inherit.aes = FALSE) +
    geom_point(data = route_mid, aes(x = mu, y = y_mid), color = "white", fill = "black", shape = 21, size = 2, inherit.aes = FALSE) +
    
    # Box with N, Mean, and SD
    geom_label(data = label_counts, aes(label = label_text), 
               x = Inf, y = Inf, hjust = 1.05, vjust = 1.05, 
               fill = "white", alpha = 0.85, size = 2.8, fontface = "bold", inherit.aes = FALSE) +
    
    # Faceting setup
    facet_wrap(~ panel_id, scales = "free_y", ncol = 4, strip.position = "top", labeller = as_labeller(combo_names)) +
    scale_x_continuous(breaks = seq(0, dyn_limit, by = 5)) +
    
    # Colors and theme
    scale_fill_manual(
      values = event_colors, 
      name = "Event Type:",
      labels = c("Anagenetic" = "Anagenetic Dispersal (d)", "Cladogenetic" = "Founder-Event Speciation (j)")
    ) +  
    theme_bw() +
    labs(
      title = plot_title,
      x = "Time (Ma)",
      y = "Frequency"
    ) +
    theme(
      legend.position = ifelse(show_legend, "top", "none"),
      legend.title = element_text(face = "bold"),
      strip.background = element_rect(fill = "gray95"),
      strip.text = element_text(face = "bold", size = 9),
      panel.grid.minor = element_blank()
    )
  
  return(p)
}

# ------------------------------------------------------------------------------
# 5. GENERATE INDIVIDUAL PANELS
# ------------------------------------------------------------------------------
# Panel A: All Genera Combined
panel_A <- create_dispersal_plot(df_plot_routes, "Panel A: All Genera Combined", show_legend = TRUE)

# Panel J: The Map
vectA <- sf::read_sf("../data/geography/A.shp") %>% 
  mutate(area_id = "Above Panama (A)") %>% 
  sf::st_transform(4326) %>%
  sf::st_make_valid() %>% 
  group_by(area_id) %>% 
  summarise()

sf::sf_use_s2(FALSE)

vectB <- sf::read_sf("../data/geography/B.shp") %>%
  sf::st_transform(4326) %>%
  sf::st_make_valid() %>%
  mutate(area_id = "On Panama (B)") %>%
  group_by(area_id) %>%
  summarise()

sf::sf_use_s2(TRUE) 

vectC <- sf::read_sf("../data/geography/C.shp") %>% 
  mutate(area_id = "Below Panama (C)") %>% 
  sf::st_transform(4326) %>%
  sf::st_make_valid() %>% 
  group_by(area_id) %>% 
  summarise()

# Combine
vectALL <- dplyr::bind_rows(vectA, vectB, vectC)

# for focusing
focusx <- c(-120, -50)
focusy <- c(-30, 50)

# main map
main_map <- ggplot(data = vectALL) +
  # Polygon
  geom_sf(aes(fill = area_id), color = "black", linewidth = 0.2) +
  
  # text
  annotate("text", x = -99, y = 37, label = "Above Panama (A)", size = 4, fontface = "bold") +
  annotate("text", x = -89, y = 7, label = "On Panama (B)", size = 4, fontface = "bold") +
  annotate("text", x = -63, y = -7, label = "Below Panama (C)", size = 4, fontface = "bold") +
  
  scale_fill_manual(values = c("Above Panama (A)" = "white", "On Panama (B)" = "grey80", "Below Panama (C)" = "white")) +
  coord_sf(xlim = focusx, ylim = focusy, expand = FALSE) + 
  
  # north star
  annotation_north_arrow(
    location = "tr", # "tr" = Top Right 
    which_north = "true",
    pad_x = unit(0.3, "in"), pad_y = unit(0.3, "in"),
    style = north_arrow_fancy_orienteering()
  ) +
  
  theme_minimal() +
  labs(title = "Panel J: Study Area") +
  theme(
    panel.background = element_rect(fill = "aliceblue"), 
    axis.title = element_blank(),
    legend.position = "none",
    axis.text = element_blank(),
    panel.grid = element_blank()
  )

# world map
world_map <- ne_countries(scale = "small", returnclass = "sf")

inset_map <- ggplot(data = world_map) +
  geom_sf(fill = "grey90", color = "white", linewidth = 0.1) +
  
  # rect
  annotate("rect", xmin = focusx[1], xmax = focusx[2], ymin = focusy[1], ymax = focusy[2],
           color = "red", fill = NA, linewidth = 1) +
  
  # focus on the MAP
  coord_sf(xlim = c(-140, -30), ylim = c(-60, 80), expand = FALSE) + 
  
  theme_void() +
  theme(
    panel.background = element_rect(fill = "white", color = "black", linewidth = 0.5) 
  )

#Join everything

panel_J <- main_map + 
  inset_element(inset_map, left = 0.02, bottom = 0.02, right = 0.35, top = 0.35)

# Generate Individual Genera Plots
# Extracting them manually to assign to specific grid locations
p_micrurus   <- create_dispersal_plot(df_plot_routes %>% filter(genera == "Micrurus"), "Genus: Micrurus")
p_marisora   <- create_dispersal_plot(df_plot_routes %>% filter(genera == "Marisora"), "Genus: Marisora")
p_oxybelis   <- create_dispersal_plot(df_plot_routes %>% filter(genera == "Oxybelis"), "Genus: Oxybelis")
p_sibon      <- create_dispersal_plot(df_plot_routes %>% filter(genera == "Sibon"), "Genus: Sibon")

p_anolis     <- create_dispersal_plot(df_plot_routes %>% filter(genera == "Anolis"), "Genus: Anolis")
p_bothrops   <- create_dispersal_plot(df_plot_routes %>% filter(genera == "Bothrops"), "Genus: Bothrops")

p_bothriechis <- create_dispersal_plot(df_plot_routes %>% filter(genera == "Bothriechis"), "Genus: Bothriechis")
p_dipsas     <- create_dispersal_plot(df_plot_routes %>% filter(genera == "Dipsas"), "Genus: Dipsas")

# ------------------------------------------------------------------------------
# 6. ASSEMBLE WITH DESIGN MATRIX & EXPORT
# ------------------------------------------------------------------------------
# The layout matrix controls exact placement. 
# 1 letter height = 1 row of facets. Genera with 2 facet rows get 2 letters.
# This prevents empty stretching for genera with only 1 dispersal type.

layout_matrix <- "
  AAAA
  AAAA
  AAAA
  AAAA
  BBMM
  CCMM
  DDMM
  EEMM
  FFHH
  FFHH
  GGII
  GGII
"

# Assemble using patchwork syntax
final_megaplot <- 
  panel_A +       # A
  p_micrurus +    # B
  p_marisora +    # C
  p_oxybelis +    # D
  p_sibon +       # E
  p_anolis +      # F
  p_bothrops +    # G
  p_bothriechis + # H
  p_dipsas +      # I
  panel_J +       # M
  plot_layout(design = layout_matrix)

# Export to a large PDF
pdf_filename <- "../figures/MEGA_graph_results_custom_layout.pdf"
pdf(pdf_filename, height = 25, width = 18)
print(final_megaplot)
dev.off()

base::message("Figure exported successfully to: ", pdf_filename)