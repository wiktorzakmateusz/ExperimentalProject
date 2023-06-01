

##### INSTALATION PART

install.packages("BiocManager")
install.packages("pacman")
install.packages('devtools')
install.packages("cli")
install.packages("rlang")

BiocManager::install(c("igraph","RColorBrewer","BH","qs","RSQLite","Cairo","RJSONIO","ggplot2",
                       "lpsymphony","data.table","metap","RandomWalkRestartMH","fitdistrplus",
                       "slam","dplyr","ggforce","graphlayouts","httr","plyr","pryr","stringr",
                       "rjson","tidyr"))

pacman::p_load(igraph, RColorBrewer, qs, rjson, RSQLite)

devtools::install_github("xia-lab/OmicsNetR", build = TRUE, build_opts = c("--no-resave-data", "--no-manual"), build_vignettes = TRUE, force = TRUE)
install.packages("https://www.xialab.ca/resources/OmicsNetR_1.0.0.tar.gz", repos = NULL, method = "libcurl");


library(devtools)
library(BiocManager)
library(rlang)
library(pacman)
library(cli)
library(OmicsNetR)

library(OmicsNetR)




##### list of genes

# Step 1. Initiate the dataSet object
dataSet<-Init.Data()

# Step 2. Map list of genes to the application
dataSet<-PrepareInputList(dataSet,"#Entrez  logFC
4495  61.12
4496  51.06
4499  23.79
6354  21.04
6369  19.76", "hsa", "gene", "entrez");

# Step 3. Identify interacting partners
dataSet<-QueryNet(dataSet, "gene", "innate")

# Step 4. Build interaction subnetwork
CreateGraph();

# Step 5. Prepare the network file to be used for visualization, the output will be in JSON format.
dataSet<-PrepareNetwork(dataSet, "subnetwork1", "omicsnet_list_of_genes.json")



##### Integration of genes and miRNA

library(OmicsNetR);
rm(list = ls()); # Clear .GlobalEnv of your R session

# Step 1. Initiate the dataSet object
dataSet<-Init.Data();

# Step 2. Map list of genes to the application
dataSet<-PrepareInputList(dataSet,
                          "#Entrez  logFC
4495  61.12
4496  51.06
4499  23.79
6354  21.04
6369  19.76", "hsa", "gene", "entrez");

# Step 2. Map list of miRNA to the application
dataSet<-PrepareInputList(dataSet,
                          "hsa-mir-101-3p
hsa-mir-133b
hsa-mir-147a
hsa-mir-3140-3p
hsa-mir-361-5p
hsa-mir-510-5p", "hsa", "mir", "mir_id");

# Step 3. Build PPI network from uploaded list of genes
dataSet<-QueryNet(dataSet, "gene", "innate", "gene");

# Step 4. Build miRNA-gene network from uploaded list of miRNA
dataSet<-QueryNet(dataSet, "mir", "mirtarbase", "mir");

# Step 5. Merge networks together through shared nodes
# and decompose into interconnected subnetworks
dataSet<-CreateGraph(dataSet);

# Step 6. Prepare the network file to be used for visualization, the output will be in JSON format.
dataSet<-PrepareNetwork(dataSet, "subnetwork1", "omicsnet_case2.json");



###### Integration of SNPs and Microbial Taxa

library(OmicsNetR);
rm(list = ls()); # Clear .GlobalEnv of your R session

# Step 1. Initiate the dataSet object
dataSet<-Init.Data();

# Step 2. Input the microbial taxa data
dataSet<-PrepareInputList(dataSet,
                          "Faecalibacterium_prausnitzii
Bacteroides_uniformis
Eubacterium_rectale
Alistipes_putredinis
Subdoligranulum_unclassified
Escherichia_coli
Bacteroides_vulgatus
Clostridium_clostridioforme
Klebsiella_pneumoniae
Clostridium_hathewayi
Alistipes_shahii
Ruminococcus_obeum", "microbiome", "mic", "species");

# Step 3. Input the SNP data
dataSet<-PrepareInputList(dataSet,
                          "rs1428554
rs41290504
rs3197999
rs516246
rs2228058", "hsa", "snp", "rsid");

# Step 4. SNP data annotation and network construction
SetSnp2GeneOpts("vep", "", "dis", 1, 5);
dataSet<-QueryNet(dataSet, "snp", "vep", "snp");

# Step 5. Protein-protein extension
SetPpiZero(FALSE);
dataSet<-QueryNet(dataSet, "gene", "huri", "snp");

# Step 6. Protein-metabolite extension based on recon3D database
dataSet<-QueryNet(dataSet, "met", "recon3D", "snp");

# Step 7. Microbial data annotation and network construction
SetOrganism("microbiome");
SetMetPotentialOpts("0.8", "TRUE", "TRUE", "TRUE");
dataSet<-QueryNet(dataSet, "mic", "agora", "mic");

# Step 8. Merge networks together through shared nodes
# and decompose into interconnected subnetworks
dataSet<-CreateGraph(dataSet);

# Step 9. Prepare the network file to be used for visualization, the output will be in JSON format.
dataSet<-PrepareNetwork(dataSet, "subnetwork1", "omicsnet_case4.json")

