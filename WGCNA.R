setwd("Brain_data/")

dataDirectory <- "."
list.files(dataDirectory, recursive=TRUE)

library(WGCNA);

options(stringsAsFactors = FALSE);
enableWGCNAThreads(16)

load(file = "5mc.and.5hmc.for.WGCNA.RData");
load(file = "Annotations.for.WGCNA.RData");

rownames(Meth) =  paste(Meth$chr, Meth$start, Meth$end, Meth$strand, sep = "_")
rownames(Hydroxy) =  paste(Hydroxy$chr, Hydroxy$start, Hydroxy$end, Hydroxy$strand, sep = "_")

#Remove the auxiliary data and transpose the expression data
datExpr0 = as.data.frame(t(Meth[, c(5:42)]));
names(datExpr0) = rownames(Meth);
rownames(datExpr0) = colnames(Meth)[c(5:42)];

datExpr0_h = as.data.frame(t(Hydroxy[, c(5:42)]));
names(datExpr0_h) = rownames(Hydroxy);
rownames(datExpr0_h) = colnames(Hydroxy)[c(5:42)];

#Check for genes and samples with too many missing values
gsg = goodSamplesGenes(datExpr0, verbose = 3);
gsg$allOK

gsg_h = goodSamplesGenes(datExpr0_h, verbose = 3);
gsg_h$allOK

## Use if the last command return FALSE
if (!gsg$allOK)
{
if (sum(!gsg$goodGenes)>0)
printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
if (sum(!gsg$goodSamples)>0)
printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}

if (!gsg_h$allOK)
{
if (sum(!gsg_h$goodGenes)>0)
printFlush(paste("Removing genes:", paste(names(datExpr0_h)[!gsg_h$goodGenes], collapse = ", ")));
if (sum(!gsg_h$goodSamples)>0)
printFlush(paste("Removing samples:", paste(rownames(datExpr0_h)[!gsg_h$goodSamples], collapse = ", ")));
datExpr0_h = datExpr0_h[gsg_h$goodSamples, gsg_h$goodGenes]
}

#Prepare the traits
traitData = meta

traitData$Sample.ID <- paste(traitData$Sample.ID, "-oxBS", sep="")
traitData$Sample.ID[2] <- paste(traitData$Sample.ID[2], "_R", sep="")
traitData$Sample.ID[20] <- paste(traitData$Sample.ID[20], "_R", sep="")
traitData$Sample.ID[26] <- paste(traitData$Sample.ID[26], "_R", sep="")

traitData$AgeDeath = as.numeric(traitData$AgeDeath)
traitData$PMI = as.numeric(traitData$PMI)
traitData$Race = as.numeric(traitData$Race) 
traitData$PTSD = as.numeric(traitData$PTSD) 
traitData$Opioid = as.numeric(traitData$Opioid) 
traitData$Smoking = as.numeric(traitData$Smoking) 

row.names(traitData)<-traitData$Sample.ID

traitData_h = meta

traitData_h$Sample.ID <- paste(traitData_h$Sample.ID, "-BS", sep="")
traitData_h$Sample.ID[2] <- paste(traitData_h$Sample.ID[2], "_R", sep="")
traitData_h$Sample.ID[20] <- paste(traitData_h$Sample.ID[20], "_R", sep="")
traitData_h$Sample.ID[26] <- paste(traitData_h$Sample.ID[26], "_R", sep="")

traitData_h$AgeDeath = as.numeric(traitData_h$AgeDeath)
traitData_h$PMI = as.numeric(traitData_h$PMI)
traitData_h$Race = as.numeric(traitData_h$Race) 
traitData_h$PTSD = as.numeric(traitData_h$PTSD) 
traitData_h$Opioid = as.numeric(traitData_h$Opioid) 
traitData_h$Smoking = as.numeric(traitData_h$Smoking) 

row.names(traitData_h)<-traitData_h$Sample.ID

#Select covariates
removedCovariates = traitData[,c(2,3,4,5,7)]
retainedCovariates = traitData[,c(6)]

removedCovariates_h = traitData_h[,c(2,3,4,5,7)]
retainedCovariates_h = traitData_h[,c(6)]

#Adjust for covariates
datEmp = empiricalBayesLM(datExpr0, removedCovariates, retainedCovariates, verbose = 3)

datEmp_h = empiricalBayesLM(datExpr0_h, removedCovariates_h, retainedCovariates_h, verbose = 3)

datEmp_adj = datEmp$adjustedData

datEmp_h_adj = datEmp_h$adjustedData

#cluster the samples to detect outliers.
sampleTree = hclust(dist(datEmp_adj), method = "average");
pdf(file = "CoMeth_postmortem_sampleClustering.pdf", width = 12, height = 9);
sizeGrWindow(12,9)
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
cex.axis = 1.5, cex.main = 2)
dev.off()

sampleTree_h = hclust(dist(datEmp_h_adj), method = "average");
pdf(file = "Hydroxy_postmortem_sampleClustering.pdf", width = 12, height = 9);
sizeGrWindow(12,9)
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree_h, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
cex.axis = 1.5, cex.main = 2)
dev.off()

# Update the tables names
datExpr = datEmp_adj
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

datExpr_h = datEmp_h_adj
nGenes_h = ncol(datExpr_h)
nSamples_h = nrow(datExpr_h)

# Data frame analogous to expression data
Samples = rownames(datExpr);
traitRows = match(Samples, traitData$Sample.ID);
datTraits = traitData[traitRows, -1];
rownames(datTraits) = traitData[traitRows, 1];

collectGarbage();

Samples_h = rownames(datExpr_h);
traitRows = match(Samples_h, traitData_h$Sample.ID);
datTraits_h = traitData_h[traitRows, -1];
rownames(datTraits_h) = traitData_h[traitRows, 1];

collectGarbage();

# Convert traits to a color representation
pdf(file = "CoMeth_postmortem_sampleClustering_Traits.pdf", width = 12, height = 9);
traitColors = numbers2colors(datTraits, signed = FALSE);
plotDendroAndColors(sampleTree, traitColors,
groupLabels = names(datTraits),
main = "Sample dendrogram and trait heatmap")
dev.off()
save(sampleTree, datExpr, datTraits, datEmp, file = "Postmortem-01-dataInput.RData")

pdf(file = "Hydroxy_postmortem_sampleClustering_Traits.pdf", width = 12, height = 9);
traitColors = numbers2colors(datTraits_h, signed = FALSE);
plotDendroAndColors(sampleTree_h, traitColors,
groupLabels = names(datTraits_h),
main = "Sample dendrogram and trait heatmap")
dev.off()
save(sampleTree_h, datExpr_h, datTraits_h, datEmp_h, file = "Hydroxy_Postmortem-01-dataInput.RData")

#Power analysis
pdf(file = "CoMeth_postmortem_Network_module.pdf", width = 9, height = 5);
powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5, blockSize = 500) ### ou aumentar
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
    xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
    main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
    labels=powers,cex=cex1,col="red");
abline(h=0.90,col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],
    xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
    main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

pdf(file = "Hydroxy_postmortem_Network_module.pdf", width = 9, height = 5);
powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft_h = pickSoftThreshold(datExpr_h, powerVector = powers, verbose = 5, blockSize = 500) ### ou aumentar
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
plot(sft_h$fitIndices[,1], -sign(sft_h$fitIndices[,3])*sft_h$fitIndices[,2],
    xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
    main = paste("Scale independence"));
text(sft_h$fitIndices[,1], -sign(sft_h$fitIndices[,3])*sft_h$fitIndices[,2],
    labels=powers,cex=cex1,col="red");
abline(h=0.90,col="red")
plot(sft_h$fitIndices[,1], sft_h$fitIndices[,5],
    xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
    main = paste("Mean connectivity"))
text(sft_h$fitIndices[,1], sft_h$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

#Calculate the networking 
net = blockwiseModules(datExpr, power = 5,  maxBlockSize = 50000,
TOMType = "signed", minModuleSize = 200,
    reassignThreshold = 0, mergeCutHeight = 0.25, 
    numericLabels = TRUE, pamRespectsDendro = FALSE,
    saveTOMs = TRUE,
    saveTOMFileBase = "postmortem_TOM",
    verbose = 3)

save(net, file = "CoMeth-02-net_int-auto.RData")

net_h = blockwiseModules(datExpr_h, power = 5,  maxBlockSize = 10000,
    TOMType = "signed", minModuleSize = 100,
    reassignThreshold = 0, mergeCutHeight = 0.25,
    numericLabels = TRUE, pamRespectsDendro = FALSE,
    saveTOMs = TRUE,
    saveTOMFileBase = "Hydroxy_TOM",
    saveConsensusTOMs = TRUE,
    consensusTOMFilePattern = "Hydroxy_consensusTOM-Block%b.RData",
    verbose = 3)

save(net_h, file = "Hydroxy-02-net_int-auto.RData")

#Modules identified
table(net$colors)

table(net_h$colors)

# hierarchical clustering dendrogram (tree)
pdf(file = "CoMeth_postmortem_Network_module_dendogram.pdf", width = 12, height = 9);
sizeGrWindow(12, 9)
mergedColors = labels2colors(net$colors)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
"Module colors",
dendroLabels = FALSE, hang = 0.03,
addGuide = TRUE, guideHang = 0.05)
dev.off()

pdf(file = "Hydroxy_Network_module_dendogram.pdf", width = 12, height = 9);
sizeGrWindow(12, 9)
mergedColors_h = labels2colors(net_h$colors)
plotDendroAndColors(net_h$dendrograms[[1]], mergedColors_h[net_h$blockGenes[[1]]],
"Module colors",
dendroLabels = FALSE, hang = 0.03,
addGuide = TRUE, guideHang = 0.05)
dev.off()

# Module assignment and module eigengene
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
save(MEs, moduleLabels, moduleColors, geneTree,
file = "CoMeth_Postmortem-02-networkConstruction-auto.RData")

moduleLabels_h = net_h$colors
moduleColors_h = labels2colors(net_h$colors)
MEs_h = net_h$MEs;
geneTree_h = net_h$dendrograms[[1]];
save(MEs_h, moduleLabels_h, moduleColors_h, geneTree_h,
file = "Hydroxy_Postmortem-02-networkConstruction-auto.RData")

datTraits = datTraits[,c(1,2,3,4,6,5)]
datTraits_h = datTraits_h[,c(1,2,3,4,6,5)]

# Quantifying module–trait associations
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

nGenes_h = ncol(datExpr_h);
nSamples_h = nrow(datExpr_h);
MEs0_h = moduleEigengenes(datExpr_h, moduleColors_h)$eigengenes
MEs_h = orderMEs(MEs0_h)
moduleTraitCor_h = cor(MEs_h, datTraits_h, use = "p");
moduleTraitPvalue_h = corPvalueStudent(moduleTraitCor_h, nSamples_h);

pdf(file = "CoMeth_postmortem_correlation_traits.pdf", width = 10, height = 140);
sizeGrWindow(10,140)
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
labeledHeatmap(Matrix = moduleTraitCor,
xLabels = names(datTraits),
yLabels = names(MEs),
ySymbols = names(MEs),
colorLabels = FALSE,
xLabelsAngle = 0,
xLabelsAdj = 0.5,  
colors = blueWhiteRed(50),
textMatrix = textMatrix,
setStdMargins = FALSE,
cex.text = 0.5,
zlim = c(-1,1),
main = paste("Module-trait relationships"))
dev.off()

pdf(file = "Hydroxy_postmortem_correlation_traits.pdf", width = 10, height = 140);
sizeGrWindow(10,140)
textMatrix = paste(signif(moduleTraitCor_h, 2), "\n(",
signif(moduleTraitPvalue_h, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor_h)
par(mar = c(6, 8.5, 3, 3));
labeledHeatmap(Matrix = moduleTraitCor_h,
xLabels = names(datTraits_h),
yLabels = names(MEs_h),
ySymbols = names(MEs_h),
colorLabels = FALSE,
xLabelsAngle = 0,
xLabelsAdj = 0.5,            
colors = blueWhiteRed(50),
textMatrix = textMatrix,
setStdMargins = FALSE,
cex.text = 0.5,
zlim = c(-1,1),
main = paste("Module-trait relationships"))
dev.off()

# Select best modules
moduleTraitCor2 = as.data.frame(moduleTraitCor)
moduleTraitCor2 <- moduleTraitCor2[((moduleTraitCor2$Opioid >= 0.40 |
                                     moduleTraitCor2$Opioid <= -0.40)),]
dim(moduleTraitCor2)

moduleTraitPvalue2 = as.data.frame(moduleTraitPvalue)
                                          
selected = rownames(moduleTraitCor2)
moduleTraitPvalue2 <- moduleTraitPvalue2[match(selected,rownames(moduleTraitPvalue2)),]

moduleTraitPvalue2 <- moduleTraitPvalue2[((moduleTraitPvalue2$Opioid <= 0.050)),]
dim(moduleTraitPvalue2) 

selected = rownames(moduleTraitPvalue2)
moduleTraitCor2 <- moduleTraitCor2[match(selected,rownames(moduleTraitCor2)),]

moduleTraitCor2 = data.matrix(moduleTraitCor2)
moduleTraitPvalue2 = data.matrix(moduleTraitPvalue2)

MEsColumn = match(selected, names(MEs));
MEs2 = MEs[, MEsColumn];

pdf(file = "CoMeth_postmortem_correlation_traits_cor0.40.pdf", width = 10, height = 5);
sizeGrWindow(10,5)
textMatrix = paste(signif(moduleTraitCor2, 2), "\n(",
signif(moduleTraitPvalue2, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor2)
par(mar = c(6, 8.5, 3, 3));
labeledHeatmap(Matrix = moduleTraitCor2,
xLabels = names(datTraits),
yLabels = names(MEs2),
ySymbols = names(MEs2),
colorLabels = FALSE,
xLabelsAngle = 0,
xLabelsAdj = 0.5,  
colors = blueWhiteRed(50),
textMatrix = textMatrix,
setStdMargins = FALSE,
cex.text = 0.5,
zlim = c(-1,1),
main = paste("Module-trait relationships"))
dev.off()

moduleTraitCor2_h = as.data.frame(moduleTraitCor_h)
moduleTraitCor2_h <- moduleTraitCor2_h[((moduleTraitCor2_h$Opioid >= 0.40 |
                                     moduleTraitCor2_h$Opioid <= -0.40)),]
dim(moduleTraitCor2_h) #4

moduleTraitPvalue2_h = as.data.frame(moduleTraitPvalue_h)
                                          
selected = rownames(moduleTraitCor2_h)
moduleTraitPvalue2_h <- moduleTraitPvalue2_h[match(selected,rownames(moduleTraitPvalue2_h)),]

moduleTraitPvalue2_h <- moduleTraitPvalue2_h[((moduleTraitPvalue2_h$Opioid <= 0.050)),]
dim(moduleTraitPvalue2_h) #35

selected = rownames(moduleTraitPvalue2_h)
moduleTraitCor2_h <- moduleTraitCor2_h[match(selected,rownames(moduleTraitCor2_h)),]

moduleTraitCor2_h = data.matrix(moduleTraitCor2_h)
moduleTraitPvalue2_h = data.matrix(moduleTraitPvalue2_h)

MEsColumn_h = match(selected, names(MEs_h));
MEs2_h = MEs_h[, MEsColumn_h];

pdf(file = "Hydroxy_postmortem_correlation_traits_cor0.45.pdf", width = 10, height = 5);
sizeGrWindow(10,5)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor2_h, 2), "\n(",
signif(moduleTraitPvalue2_h, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor2_h)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor2_h,
xLabels = names(datTraits_h),
yLabels = names(MEs2_h),
ySymbols = names(MEs2_h),
colorLabels = FALSE,
xLabelsAngle = 0,
xLabelsAdj = 0.5,            
colors = blueWhiteRed(50),
textMatrix = textMatrix,
setStdMargins = FALSE,
cex.text = 0.5,
zlim = c(-1,1),
main = paste("Module-trait relationships Hydroxymethylation"))
dev.off()

#Check moduleColors
table(moduleColors)
table(moduleColors_h)

# Define variable weight 
Opioid = as.data.frame(datTraits$Opioid);
names(Opioid) = "Opioid"
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(datExpr, Opioid, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(Opioid), sep="");
names(GSPvalue) = paste("p.GS.", names(Opioid), sep="");

Opioid_h = as.data.frame(datTraits_h$Opioid);
names(Opioid_h) = "Opioid"
modNames_h = substring(names(MEs_h), 3)
geneModuleMembership_h = as.data.frame(cor(datExpr_h, MEs_h, use = "p"));
MMPvalue_h = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership_h), nSamples_h));
names(geneModuleMembership_h) = paste("MM", modNames_h, sep="");
names(MMPvalue_h) = paste("p.MM", modNames_h, sep="");
geneTraitSignificance_h = as.data.frame(cor(datExpr_h, Opioid_h, use = "p"));
GSPvalue_h = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance_h), nSamples_h));
names(geneTraitSignificance_h) = paste("GS.", names(Opioid_h), sep="");
names(GSPvalue_h) = paste("p.GS.", names(Opioid_h), sep="");

#Module membership
modules_color = c("lightpink2","orchid4","steelblue1","cadetblue1","ghostwhite","wheat1","turquoise1",
                  "darkorchid2","turquoise3","aquamarine1")
#>= 0.4 of correlation
modules_color_c = c("lightpink2","orchid4","steelblue1","cadetblue1","ghostwhite","wheat1","turquoise1",
                  "darkorchid2","turquoise3","aquamarine1")

nmodules = length(modules_color)

pdf(file = "CoMeth_ModuleMembership_OPIOID.pdf", width = 7, height = 7);

for (mod2 in 1:nmodules)
{
module = modules_color[[mod2]]
col_mod = modules_color_c[[mod2]]
column = match(module, modNames);
moduleGenes = moduleColors==module;

sizeGrWindow(7, 7);    
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
    abs(geneTraitSignificance[moduleGenes, 1]),
    xlab = paste("Module Membership in", module, "module"),
    ylab = "Gene significance for Opioid",
    main = paste("Module membership vs. gene significance\n"),
    cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = col_mod)
}

dev.off()

modules_color_h = c("thistle1","rosybrown1","orange4","navajowhite3")
modules_color_c_h = c("thistle1","rosybrown1","orange4","navajowhite3")

nmodules_h = length(modules_color_h)

pdf(file = "Hydroxy_ModuleMembership_OPIOID.pdf", width = 7, height = 7);

for (mod2 in 1:nmodules_h)
{
module = modules_color_h[[mod2]]
col_mod = modules_color_c_h[[mod2]]
column = match(module, modNames_h);
moduleGenes = moduleColors_h==module;

sizeGrWindow(7, 7);    
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership_h[moduleGenes, column]),
    abs(geneTraitSignificance_h[moduleGenes, 1]),
    xlab = paste("Module Membership in", module, "module"),
    ylab = "Gene significance for Opioid",
    main = paste("Module membership vs. gene significance\n"),
    cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = col_mod)
}

dev.off()

