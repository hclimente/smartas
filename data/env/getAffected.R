getAffected <- function(genes,swt,muts){
  
  # if present, add switches in the target gene
  affectedSwitches <- swt %>%
    filter(GeneId %in% genes & IsFunctional==1) %>%
    select(Tumor,Symbol,Patients_affected,PatientNumber) %>%
    mutate(What="Switch")
  
  suppressWarnings( nocols <- max(affectedSwitches$PatientNumber) )
  
  if (nocols <= 0)
    return(NA)
  
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
    mutate(Alteration = "SPLICING")
  
  affectedMutations.long <- muts %>%
    filter(GeneId %in% genes) %>%
    mutate(Alteration2 = "MUT")
  
  affected.long <- merge(affectedSwitches.long, affectedMutations.long,all=T) %>%
    mutate(Alteration=paste(Alteration,Alteration2,sep=";"))
  
  affected.wide <- affected.long %>%
    select(Symbol,Patient,Alteration) %>%
    unique %>%
    spread(Patient, Alteration)
  
  affected.wide[is.na(affected.wide)] <- ""
  
  rownames(affected.wide) <- affected.wide$Symbol
  affected.wide <- as.matrix(affected.wide)
  affected.wide <- affected.wide[, !colnames(affected.wide) %in% c("Tumor","Symbol")]
  affected.wide <- affected.wide[,colSums(affected.wide=="") < nrow(affected.wide)]
  
  list("wide" = affected.wide, "long" = affected.long)
}