source("~/smartas/notebook/data/env/variables.R")
source("~/wisdom/r/data_analysis_environment.R")
source("~/wisdom/r/clean_theme.R")
library(find.me)

switches.split <- "~/smartas/notebook/data/pancancer/candidateList_full.tumorSplit.tsv" %>%
  read_tsv
switches <- read_tsv("~/smartas/notebook/data/pancancer/candidateList_full.tsv") %>%
  filter(Reliable==1)

drivers.file <- "~/smartas/notebook/data/intogen_cancer_drivers-2014.12b/Mutational_drivers_per_tumor_type.tsv"
drivers <-  read_tsv(drivers.file, comment="#") %>%
  mutate(Tumor_type = ifelse(Tumor_type=="COREAD", "coad", Tumor_type),
         Tumor_type = ifelse(Tumor_type=="HC", "lihc", Tumor_type),
         Tumor_type = ifelse(Tumor_type=="RCCC", "kirc", Tumor_type),
         Tumor_type = tolower(Tumor_type) ) %>%
  set_colnames(c("Symbol","Tumor")) %>%
  mutate(Driver=TRUE)

proteome <- read_tsv("~/smartas/notebook/data/mutations/proteome_information.txt") %>%
  set_colnames(c("Tumor","GeneId","Symbol","Transcript","TPM","ProteinLength","asEvidence"))

wes <- read_tsv("~/smartas/notebook/data/mutations/wes_mutations.txt") %>%
  select(Tumor,GeneId,Symbol,Patient) %>%
  unique

# read interactions
ppi.file <- "~/smartas/notebook/data/eporta/raw_tables/Switched_interactions_consensus.txt"

## get max number of columns (necessary for reading)
no_col <- max(count.fields(ppi.file,sep = "\t"))
no_col.ppi <- (no_col-6)/2
ppi.cols <- paste(c("Origin","Interaction"), floor(seq(1,no_col.ppi,0.5)), sep="_")

## read table
ppi <- read.table(ppi.file,header=F,fill=T,col.names=1:no_col) %>%
  set_colnames(c("GeneId","Symbol","Normal_transcript","Tumor_transcript",
                 "partnerId","partnerSymbol",ppi.cols)) %>%
  # all Origin columns contail "DDI_match", so we can disregard them
  select(-starts_with("Origin_")) %>%
  # convert from wide to long table format
  melt(id.vars = c("GeneId","Symbol","Normal_transcript","Tumor_transcript",
                   "partnerId","partnerSymbol"),
       value.name = "Interaction") %>%
  select(-variable) %>%
  # remove cases with no interaction described
  filter(Interaction!="") %>%
  # split interaction information
  separate(Interaction, into=c("What","Transcript","Pfams"), sep="-") %>%
  # remove pfams columns (account for different domains for the same interaction)
  select(-Pfams) %>%
  # remove several instances of the same isoform
  unique %>%
  # annotate with switch info
  merge(switches.split) %>%
  # consider only the most abundant isoform as partner: 
  # one interaction per pair & only expressed genes
  merge(proteome,by=c("Tumor","Transcript"), suffixes = c(".switch",".partner"))

# filter out patients without RNAseq data available
rnaseqPatients <- strsplit(switches$Patients_affected,",") %>% unlist %>% unique
mutations <- wes %>%
  filter(Patient %in% rnaseqPatients) %>%
  mutate(Alteration2="MUT")

pannegative <- merge(wes,drivers,all.x=T) %>%
  mutate(Driver = ifelse(is.na(Driver), FALSE, TRUE)) %>%
  group_by(Patient) %>%
  summarise(Pannegative = ifelse(sum(Driver), "Mut+", "Mut-"))

for (x in unique(drivers$Symbol) ) {
  affectedSwitches <- ppi %>%
    filter(partnerSymbol==x & What!="Kept") %>%
    select(Tumor,Symbol.switch,Patients_affected,PatientNumber,What) %>%
    set_colnames(c("Tumor","Symbol","Patients_affected","PatientNumber","What"))
  
  # if present, add switches in the target gene
  affectedSwitches <- switches.split %>%
    filter(Symbol==x & Tumor %in% affectedSwitches$Tumor & IsFunctional==1) %>%
    select(Tumor,Symbol,Patients_affected,PatientNumber) %>%
    mutate(What="Switch") %>%
    rbind(affectedSwitches)
  
  suppressWarnings( nocols <- max(affectedSwitches$PatientNumber) )
  
  if (nocols <= 0)
    next
  
  suppressWarnings(
    affectedSwitches.long <- affectedSwitches %>%
      separate(Patients_affected, paste("Patient", 1:nocols, sep="_"), sep=",")
  )
  
  affectedSwitches.long <- affectedSwitches.long %>%
    select(Symbol,Tumor,What,starts_with("Patient_")) %>%
    melt(id.vars = c("Symbol","Tumor","What")) %>%
    select(-variable) %>%
    set_colnames(c("Symbol","Tumor","What","Patient")) %>%
    filter(!is.na(Patient)) %>%
    mutate(Alteration = ifelse(What=="Lost", "AMP", "GERMLINE"),
           Alteration = ifelse(What=="Gained", "HOMDEL", Alteration),
           Alteration = ifelse(What=="Switch", "SPLICING", Alteration))
  
  affectedMutations.long <- mutations %>%
    filter(Symbol==x)
  
  affected.long <- merge(affectedSwitches.long, affectedMutations.long,all=T) %>%
    mutate(Alteration=paste(Alteration,Alteration2,sep=";"))
  
  affected.wide <- affected.long %>%
    select(Symbol,Patient,Alteration) %>%
    spread(Patient, Alteration)
  
  affected.wide[is.na(affected.wide)] <- ""
  
  rownames(affected.wide) <- affected.wide$Symbol
  affected.wide <- as.matrix(affected.wide)
  affected.wide <- affected.wide[, !colnames(affected.wide) %in% c("Tumor","Symbol")]
  affected.wide <- affected.wide[,colSums(affected.wide=="") < nrow(affected.wide)]
  
  plot.colors <- c(colorPalette, "amp" = "firebrick", "del" = "blue", "up" = NA, 
                   "down" = NA, "splicing" = "forestgreen", "germline" = "purple",
                   "somatic" = "#36454F", "Mut+" = "firebrick", "Mut-" = "gray80", 
                   "MutUnknown"="white")
  
  patients <- affected.long %>% 
    select(Tumor,Patient) %>%
    unique %>%
    merge(pannegative,all.x=T) %>%
    mutate(Pannegative = ifelse(is.na(Pannegative), "MutUnknown", Pannegative))
  
  ngenes <- nrow(affected.wide)
  
  sorted.matrix <- getSortedMatrix(affected.wide)
  affected.wide <- affected.wide[c(x,setdiff(rownames(sorted.matrix$mutmat),x)),]
  
  p <- oncoprint(affected.wide, sortGenes=FALSE) + 
    geom_tile(data=patients, aes(x=Patient,y=ngenes+1,fill=Tumor), height=0.3) +
    geom_tile(data=patients, 
              aes(x=Patient, y=ngenes+.75, fill=as.character(Pannegative)), height=0.1) +
    labs(title=x, x="") + 
    clean_theme() +
    theme(axis.text.x=element_blank())
  
  suppressWarnings( p <- p + scale_fill_manual(values = plot.colors) )
  
  genesOncoprint <- unique(affected.long$Symbol)
  name <- paste(length(genesOncoprint),x,paste(setdiff(genesOncoprint,x),collapse="."),"png",sep=".")
  
  ggsave(paste0("~/smartas/notebook/results/oncoprints/",name),p, width = 10, height = 10)
}