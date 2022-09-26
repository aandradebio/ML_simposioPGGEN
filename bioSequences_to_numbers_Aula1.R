library(data.table)
library(dplyr)
library(kmer)
library(seqinr)
#library(taxonomizr) 

load("./data/bioSequences_to_numbers_Aula1.RData")

### Importar arquivos fasta
plant_vir <- read.fasta(file="./data/plant_vir.fasta", set.attributes = T)
mosquito_vir <- read.fasta(file="./data/mosquito_vir.fasta", set.attributes = T)

typeof(plant_vir)

### Remove os accession numbers do multifasta

AccessionNum_plant <- getAnnot(plant_vir)
AccessionNum_plant <- unname(unlist(AccessionNum_plant))
AccessionNum_plant <- gsub('>','',AccessionNum_plant)

AccessionNum_mosquito <- getAnnot(mosquito_vir)
AccessionNum_mosquito <- unname(unlist(AccessionNum_mosquito))
AccessionNum_mosquito <- gsub('>','',AccessionNum_mosquito)

### Adiciona informações de taxonomia a partir do header das sequências
#prepareDatabase('accessionTaxa.sql') Arquivo pesado. 19GB (não disponibilizado)
#taxaID_mosquito <-accessionToTaxa(AccessionNum_mosquito, "accessionTaxa.sql")
#taxaID_plant <-accessionToTaxa(AccessionNum_plant, "accessionTaxa.sql")
#tax_mosquito <- getTaxonomy(taxaID_mosquito,'accessionTaxa.sql')
#tax_plant <- getTaxonomy(taxaID_plant,'accessionTaxa.sql')

tax <- read.csv(file="./data/tax_example.csv")
tax_dup <- tax[duplicated(tax$taxaID), ] ###remover taxa duplicada
tax <- tax[!grepl("Togaviridae", tax$family),] ### remove herpesviridae
tax_subset <- subset(tax, select = c(superkingdom, phylum, class, order, family, genus, species))
tax_subset <- tax_subset[order(tax_subset$family),] ### ordenar
table(tax_subset$family) # N. de vírus por família
barplot(table(tax_subset$family)) # Porcentagem
prop.table(table(tax_subset$family)) # Porcentagem

### Escolher os genomas utilizados com base nos seus hospedeiros
virusHost <- read.csv(file="./data/virushostdb.tsv", header=T, sep="\t")

### Sumariza as informações biológicas em contagens das bases
### Chama os k-mers
plant_k2 <- kcount(plant_vir, k = 2, residues = NULL, gap = "-", named = TRUE,compress = TRUE, encode = FALSE) 
plant_vir = plant_k2[,colSums(plant_k2) > 2] # Filtrar separadamente 
plant_vir <- cbind(rownames(plant_vir), data.frame(plant_vir, row.names=NULL))
write.csv(plant_vir, file="./data/plant_vir_k2.csv", row.names=F)

mosquito_k2 <- kcount(mosquito_vir, k = 2, residues = NULL, gap = "-", named = TRUE,compress = TRUE, encode = FALSE) 
mosquito_vir = mosquito_k2[,colSums(mosquito_k2) > 2] # Filtrar separamente
mosquito_vir <- cbind(rownames(mosquito_vir), data.frame(mosquito_vir, row.names=NULL))
write.csv(mosquito_vir, file="./data/mosquito_vir_k2.csv", row.names=F)

save.image("./data/bioSequences_to_numbers_Aula1.RData")
