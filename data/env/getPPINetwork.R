library(network)
library(ggnetwork)
library(igraph)
library(intergraph)

getEffect <- function(x,num){
  counts <- data.frame(x=x, num=num) %>%
    group_by(x) %>%
    summarise(num = sum(num)) %>%
    arrange(desc(num))
  
  if (nrow(counts) > 1 & length(unique(counts$num)) == 1)
    "Draw"
  else
    as.character(counts$x[1])
}

getPPINetwork <- function(geneset, ppi.network, ppi.changes){
  set.ppi.effect <- ppi.changes %>%
    filter(GeneId %in% geneset & partnerGeneId %in% geneset) %>%
    select(Symbol,partnerSymbol,Tag,PatientNumber) %>%
    unique
  
  set.ppi.effect$EdgeTag <- set.ppi.effect %>% 
    select(Symbol,partnerSymbol) %>%
    apply(1, function(x) sort(x) %>% paste(collapse="-") )
  
  set.ppi.effect <- set.ppi.effect %>%
    select(EdgeTag,Tag,PatientNumber) %>%
    set_colnames(c("EdgeTag","Effect","PatientNumber"))
  
  set.ppi <- ppi.network %>%
    filter(Gene1 %in% geneset & Gene2 %in% geneset)
  
  set.ppi$EdgeTag <- set.ppi %>% 
    select(Symbol1,Symbol2) %>%
    apply(1, function(x) sort(x) %>% paste(collapse="-") )
  
  set.ppi <- set.ppi %>%
    select(EdgeTag) %>%
    unique
  
  set.ppi.annotated <- merge(set.ppi, set.ppi.effect, all=T) %>%
    mutate(Effect = ifelse(is.na(Effect), "Unaffected", Effect)) %>%
    filter(Effect != "Unaffected") %>%
    group_by(EdgeTag) %>%
    summarise(Consensus = ifelse(length(unique(Effect))==1, "Yes", "No"),
              Effect = getEffect(Effect,PatientNumber)) %>%
    separate(EdgeTag, c("Gene1","Gene2"), "-")
  
  # Create network
  nw <- set.ppi.annotated %>%
    graph.data.frame(directed = FALSE)  %>%
    asNetwork
  
  return(nw)
}