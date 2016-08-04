cancerTypes <- c("brca","coad","hnsc","kich","kirc","kirp","lihc","luad","lusc","prad","thca")
workingDir <- "/genomics/users/hector/smartas/results"
project <- "/projects_rg/TCGA/users/hector/SmartAS/"

colorPalette <- c("#CC79A7","#636363","#F0E442","#006D2C","#31A354","#74C476","#FC8D59","#08519C","#3182BD","#D55E00","#5E3C99","#000000","#969696","#EDF8FB","#B3CDE3","#8C96C6","#88419D","#D94701","#2171B5","#FD8D3C","#6BAED6","#FD8D3C","#33A02C","#E31A1C","#FF7F00","#6A3D9A","#B15928","#377EB8","#E41A1C")
names(colorPalette) <- c("brca","coad","hnsc","kich","kirc","kirp","lihc","luad","lusc","prad","thca","coad-hyper","coad-hypo","brca-basal","brca-her2","brca-luminala","brca-luminalb","expression-up","expression-down","psi-up","psi-down","odds","a3","a5","mx","ri","se","normal","tumor")

nPatients <- list()
nPatients[["brca"]] <- 1036
nPatients[["coad"]] <- 262
nPatients[["hnsc"]] <- 422
nPatients[["kich"]] <- 62
nPatients[["kirc"]] <- 505
nPatients[["kirp"]] <- 195
nPatients[["lihc"]] <- 197
nPatients[["luad"]] <- 488
nPatients[["lusc"]] <- 483
nPatients[["prad"]] <- 295
nPatients[["thca"]] <- 497
nPatients[["total"]] <- sum(1036,262,422,62,505,195,197,488,483,295,497)

nPatientsDf <- as.data.frame(do.call("rbind",nPatients))
nPatientsDf$Cancer <- rownames(nPatientsDf)
colnames(nPatientsDf) <- c("TotalPatients","Cancer")