#Packages to load
library(biomformat) #1.7.0
library(dplyr) #1.1.1
library(tidyr) #1.1.3
library(ggplot2) #3.4.3
library(RColorBrewer) #1.1-3
library(xlsx) #0.6.5
library(forcats) #0.5.1
library(cowplot) #1.1.1
library(stringr) #1.5.0

#functions
se <- function(x) {sqrt(var(x)/length(x))}

#DADA2 stats
DADA2_stats = list()

for (i in list.files(pattern=".*DADA2.*.tsv")) {
  primer=gsub("CCA_", "", i)
  primer=gsub("_.*", "", primer)
  
  DADA2 <- read.delim(file = i, header = T)
  DADA2 <- DADA2[2:nrow(DADA2),]
  
  DADA2$primer <- rep(primer)
  
  DADA2_stats[[i]] <- DADA2
}

DADA2_table <- do.call(rbind, DADA2_stats)

DADA2_table_sum <- DADA2_table %>% select(-1) %>% mutate_at(c(1:8), as.numeric) 
DADA2_table_sum <- left_join(DADA2_table_sum[, c(1, 2, 4, 5, 7, 9)] %>% aggregate(.~primer, ., sum), DADA2_table_sum[, c(3, 6, 8, 9)] %>% aggregate(.~primer, ., mean)) %>% .[, c(1,2, 3, 7, 4, 5, 8, 6, 9)]

#Length of sequences
seq_len = list()

for (i in list.files(pattern="*.fasta$")) {
  primer=gsub("_.*", "", i)
  
  seq.ID <- read.delim(file = i, header = F) %>% filter(., grepl(">", V1)==T)
  seq.DNA <- read.delim(file = i, header = F) %>% filter(., grepl(">", V1)==F)
  
  seq.len <- cbind(seq.ID, seq.DNA)
  seq.len <- rename(seq.len, Feature.ID=1, Seq=2)
  seq.len$bp <- nchar(seq.len$Seq)
  seq.len$primer <- rep(primer)
  seq.len$Feature.ID <- gsub(">", "", seq.len$Feature.ID)
  
  seq_len[[i]] <- seq.len
}

seq.DADA2.output  <- do.call(rbind, seq_len)

seq.DADA2.output_sum <- seq.DADA2.output %>% select(primer, bp) %>% group_by(primer) %>% summarise(mean = mean(bp), max = max(bp), median = median(bp), min = min(bp), SE = se(bp), ASVs = length(bp))

DADA2_table_sum <- left_join(DADA2_table_sum %>% mutate(primer = if_else(primer == "23s", "23S", primer)), seq.DADA2.output_sum[, c(1, 2, 6, 7)])

#Load the sequence table - specially this is called an amplicon sequence variant table
#####BIOM TO CSV####
ASV_list = list()
for (i in list.files(pattern="*.biom")) {
  primer=gsub("_.*", "", i)
  ASV <- read_biom(biom_file = i)
  ASV <- as.data.frame(as.matrix(biom_data(ASV)))
  ASV$Feature.ID <- rownames(ASV)
  col_names <- colnames(ASV)
  col_excl <- col_names %>% grepl("^EB.*|Feature.ID", .)
  col_2_piv=col_names[c(col_excl)]
  ASV <- pivot_longer(ASV, cols = -all_of(x = col_2_piv), names_to = "CCA.ID", values_to = "Reads")
  ASV <- separate(ASV, col = "CCA.ID", into = c("CCA.ID", "Primer"), sep = "-")
  ASV$ext <- if_else(ASV$CCA.ID == "94EIA" | ASV$CCA.ID == "94EOA", "E1", "E2")
  colnames(ASV) <- sub("-.*", "", colnames(ASV))
  df_EB1 <- filter(ASV, ext == "E1")
  print(primer)
  try({df_EB1$Reads_nrm <- df_EB1$Reads
  df_EB1$Reads_nrm <- df_EB1$Reads - df_EB1$EB1})
  #df_EB1 can produce this error if cannot do subtraction as there are no negative reads. "Error: Assigned data `df_EB1$Reads - df_EB1[, grep(pattern = "^EB1.*", colnames(df_EB1))]` must be compatible with existing data." Using try to let the code try the line but does not break out of the loop if there is an error
  df_EB2 <- filter(ASV, ext == "E2")
  try({df_EB2$Reads_nrm <- df_EB2$Reads
  df_EB2$Reads_nrm <- df_EB2$Reads - df_EB2$EB2})
  
  ASV <- rbind(df_EB1, df_EB2)
  ASV <- filter(ASV, Reads_nrm > 10)
  ASV <- ASV %>% select(Feature.ID, CCA.ID, Primer, Reads_nrm)
  
  ASV_list[[i]] <- ASV
  
}

ASV_table <- do.call(rbind, ASV_list)

ASV_table$Primer <- gsub("23", "23S", ASV_table$Primer)
ASV_table$Primer <- gsub("16", "16S", ASV_table$Primer)

#Table S5
ASV_table %>% mutate(From = if_else(grepl("O", CCA.ID)==T, "Outer", "Inner"), CCA.no = if_else(grepl("94", CCA.ID)==T, "CCA694", "CCA643"), Preservation = if_else(grepl("[0-9][0-9]E", CCA.ID)==T, "Ethanol", "Silica")) %>% select(-CCA.ID) %>% aggregate(Reads_nrm~., ., sum) %>% group_by(Primer, Preservation, From, CCA.no) %>% summarise(Reads = sum(Reads_nrm), ASVs = length(Feature.ID)) %>% mutate(ASVs.Reads = paste(ASVs, " (", Reads, ")", sep = "")) %>% select(-Reads, -ASVs) %>% pivot_wider(., names_from = From, values_from = `ASVs.Reads`) %>% write.table(., "./TableS5.tsv", row.names = F, quote = F, sep = "\t")

DADA2_table_sum <- ASV_table %>% group_by(Primer) %>% select(Feature.ID, Reads_nrm) %>% summarise(ASV_n=length(unique(Feature.ID)), ReadsFiltered = sum(Reads_nrm)) %>% rename(., primer = Primer) %>% left_join(DADA2_table_sum, .)

#tablS4
#write.table(x = DADA2_table_sum, file = "./TableS4.tsv", quote = F, row.names = F)

#fasta seq for inner & outer
fasta <- inner_join(seq.DADA2.output, ASV_table)

for (i in c("16S", "18S", "23S", "COI", "ITS", "rbcL", "tufA")) {
  R.filter.fasta <- fasta %>% filter(., Primer==i) %>% filter(., grepl("^4", CCA.ID)==T) %>% select(Feature.ID, Seq) %>% unique()
  R.filter.fasta$Feature.ID <- gsub("^", ">", R.filter.fasta$Feature.ID) 
  R.filter.fasta <- data.frame(V1=c(do.call(rbind, lapply(R.filter.fasta, as.character))))
  
  write.table(x = R.filter.fasta, file = paste(i, "_43_fasta_R_filter.txt", sep=""), col.names = F, row.names = F, quote = F)
  
  R.filter.fasta <- fasta %>% filter(., Primer==i) %>% filter(., grepl("^9", CCA.ID)==T) %>% select(Feature.ID, Seq) %>% unique()
  R.filter.fasta$Feature.ID <- gsub("^", ">", R.filter.fasta$Feature.ID) 
  R.filter.fasta <- data.frame(V1=c(do.call(rbind, lapply(R.filter.fasta, as.character))))
  
  write.table(x = R.filter.fasta, file = paste(i, "_94_fasta_R_filter.txt", sep=""), col.names = F, row.names = F, quote = F)
}

#these fasta files were then used to create phylogenic trees in Qiime 2

#GenBank data
GB_list=list()
for (i in list.files(pattern="*blastn")) {
  primer=gsub("_.*", "", i)
  GB_tax <- read.delim(file = i, header = F)
  GB_tax <- rename(GB_tax, "Feature.ID" = "V1", "sskingdoms" = "V2", "qstart" = "V3", "qend" = "V4", "qlen" = "V5", "qseq" = "V6", "qcovs" = "V7", "pident" = "V8", "sseqid" = "V9", "sgi" = "V10", "sacc" = "V11", "sstart" = "V12", "send" = "V13", "staxids" = "V14", "sscinames" = "V15", "stitle" = "V16", "length" = "V17", "evalue" = "V18", "bitscore" = "V19", "qcovhsp" = "V20")
  print(primer)
  "The number of ASVs that have multiple matches to GenBank:"
  print(length(unique(GB_tax$Feature.ID)))
  
  dup <- GB_tax %>% group_by(Feature.ID) %>% filter(n() > 1) #no duplicates
  try({if (dup > 0) 
    {
    
    print(length(unique(dup$Feature.ID)))
    GB_tax <- GB_tax %>% group_by(Feature.ID) %>% arrange(desc(qcovhsp), desc(pident)) %>% slice_head(., n = 1)
    }
  })
  
  GB_tax <- as.data.frame(GB_tax)
  
  GB_tax <- filter(GB_tax, qcovhsp >= 90)
  GB_tax <- filter(GB_tax, pident >= 85)
  
  GB_tax$GB_Genus <- gsub(" .*", "", GB_tax$sscinames)
  GB_tax$Primer <- rep(primer)
  GB_list[[i]] <- GB_tax
}

GB_table <- do.call(rbind, GB_list)

#These are the number of ASVs that met filtering requirements to be kept as an assigned taxonomy. It is not the number of ASVs that used as it does not take into account e.g., >10 reads/ASV/sample
GB_table %>% filter(., grepl("uncultured", sscinames)) %>% select(Feature.ID, Primer) %>% group_by(Primer) %>% summarise(ASV_n=length(unique(Feature.ID)))

GB_table %>% filter(., grepl("uncultured", sscinames)==F) %>% select(Feature.ID, Primer) %>% group_by(Primer) %>% summarise(ASV_n=length(unique(Feature.ID)))

#WoRMs to get full classification. Select only those ASVs then based filtering. Remove 'uncultured' and then assign based on genus level. Assigning at genus level avoids issues with species
#GB_table %>% inner_join(ASV_table %>% select(Feature.ID), .) %>% select(sscinames) %>% mutate(ScientificName=if_else(grepl("uncultured", sscinames)==T, word(sscinames, 2), sscinames)) %>% mutate(ScientificName=gsub(" .*", "", ScientificName)) %>% select(ScientificName) %>% unique %>% write.csv(., "./CCA_WoRMs_genus_level.csv", row.names = F, quote = F)

worms <- read.xlsx2("./CCA_WoRMs_genus_level_matched.xlsx", sheetIndex = 1, header = T)

worms <- worms[,c("ScientificName", "Phylum", "Class", "Order", "Family", "Genus")] %>% filter(Phylum != "")

GB_worms <- left_join(GB_table %>% mutate(ScientificName = if_else(grepl("uncultured", sscinames)==T, word(sscinames, 2), sscinames)) %>% mutate(ScientificName=gsub(" .*", "", ScientificName)), worms)

ASV_tax <- left_join(ASV_table , GB_worms)

#not assigned by WoRMS
ASV_tax %>% filter(is.na(Phylum)) %>% select(Feature.ID, Primer) %>% group_by(Primer) %>% summarise(ASV_n=length(unique(Feature.ID)))

#assigned by WoRMS
ASV_tax %>% filter(Phylum != "") %>% select(Feature.ID, Primer) %>% group_by(Primer) %>% summarise(ASV_n=length(unique(Feature.ID)))

#reads assigned - dont think really need this
ASV_tax %>% mutate(Assigned = if_else(is.na(ASV_tax$qstart)==F, "T", "F")) %>% group_by(CCA.ID, Primer, Assigned) %>% summarise(total = sum(Reads_nrm)) %>% pivot_wider(., names_from = 2, values_from = 4, values_fill = 0) %>% View

#inner v outer
ASV_tax <- ASV_tax %>% mutate(ID.no=if_else(grepl("43", CCA.ID)==T, "43", "94"), From=if_else(grepl("O", CCA.ID)==T, "Outer", "Inner"), Pres=if_else(grepl("[1-9]E", CCA.ID)==T, "Ethanol", "Silica")) 

ASV_tax %>% select(Feature.ID, ID.no, From, Pres, CCA.ID, Primer) %>% unique() %>% group_by(Primer, CCA.ID, ID.no, From, Pres) %>% summarise(ASVs=length(Feature.ID)) %>% mutate(CCA.ID = gsub(".*[I,O]", "", CCA.ID)) %>%  pivot_wider(., names_from = From, values_from = ASVs) %>% mutate(IorO = Outer-Inner, Perc = Inner/Outer) %>% 
  #inner vs outer
  mutate(IvO = if_else(Perc > 1, "I", "O")) %>% ungroup() %>% select(IvO) %>% table()
  #silica v ethanol
  #filter(ID.no == "94") %>% mutate(total = Inner+Outer) %>% View

ASV_tax %>% select(Primer, Feature.ID, Phylum) %>% group_by(Primer, Phylum) %>% summarise(ASVs=length(unique(Feature.ID)))

#inner v outer by ASvs
#ASV_tax %>% 
  #filter(., Phylum %in% c("Rhodophyta", "Chlorophyta", "Ochrophyta")) %>% 
  #filter(Primer == "23S") %>%
algae %>%  
filter(grepl("^4", CCA.ID)) %>% group_by(Primer, From, Phylum) %>% summarise(tot_ASVs=length(unique(Feature.ID)), tot_Reads=sum(Reads_nrm)) %>% mutate(value=paste(tot_ASVs," (", tot_Reads, ")", sep = "")) %>% select(-tot_ASVs, -tot_Reads) %>% pivot_wider(., names_from = From, values_from = value, values_fill = "0 (0)") %>% mutate(sum_algae=paste(Inner, " | ", Outer, sep="")) %>% select(-Inner, -Outer)

#algae only
algae %>%  
  + group_by(Primer, From, Phylum) %>% summarise(tot_ASVs=length(unique(Feature.ID)), tot_Reads=sum(Reads_nrm)) %>% mutate(value=paste(tot_ASVs," (", tot_Reads, ")", sep = "")) %>% select(-tot_ASVs, -tot_Reads) %>% pivot_wider(., names_from = From, values_from = value, values_fill = "0 (0)") %>% mutate(sum_algae=paste(Inner, " | ", Outer, sep="")) %>% select(-Inner, -Outer) %>% pivot_wider(., names_from = Phylum, values_from = sum_algae, values_fill = "NA") %>% View


ASV_tax %>% mutate(GB_ID=if_else(grepl("[A-Z]", Phylum)==T, "Yes", "No")) %>% select(Primer, CCA.ID, Feature.ID, Reads_nrm, GB_ID) %>% aggregate(Reads_nrm~., ., sum) %>% group_by(Primer, CCA.ID, GB_ID) %>% summarise(ASV=length(Feature.ID), Reads=sum(Reads_nrm)) %>% ungroup() %>% group_by(Primer, CCA.ID) %>% mutate(total_ASVs=sum(ASV), ASVs_perc=ASV/total_ASVs*100) %>% View

ASV_tax %>% mutate(GB_ID=if_else(grepl("[A-Z]", Phylum)==T, "Yes", "No")) %>% select(Primer, CCA.ID, Feature.ID, Reads_nrm, GB_ID) %>% aggregate(Reads_nrm~., ., sum) %>% group_by(Primer, CCA.ID, GB_ID) %>% summarise(ASV=length(Feature.ID), Reads=sum(Reads_nrm)) %>% ungroup() %>% group_by(Primer, CCA.ID) %>% mutate(total_ASVs=sum(ASV), ASVs_perc=ASV/total_ASVs*100) %>% filter(., grepl("^94", CCA.ID)) %>% filter(., GB_ID == "No") %>% select(Primer, CCA.ID, ASVs_perc) %>% pivot_wider(., names_from = CCA.ID, values_from = ASVs_perc, values_fill = 0) %>% View

#ASV Venn diagram
library(ggVennDiagram)
test <- ASV_tax %>% filter(Primer == "ITS")
x <- list(
  EIA = test %>% filter(grepl("43EOA", CCA.ID)) %>% as.data.frame() %>% select(Class) %>% unique() %>% .[,1],
  EIB =   test %>% filter(grepl("43EOB", CCA.ID)) %>% as.data.frame() %>% select(Class) %>% unique() %>% .[,1],
  EIC =   test %>% filter(grepl("43EOC", CCA.ID)) %>% as.data.frame() %>% select(Class) %>% unique() %>% .[,1],
  EID =   test %>% filter(grepl("43EOD", CCA.ID)) %>% as.data.frame() %>% select(Class) %>% unique() %>% .[,1],
  EIE =   test %>% filter(grepl("43EOE", CCA.ID)) %>% as.data.frame() %>% select(Class) %>% unique() %>% .[,1]
  )

ggVennDiagram(x)

test <- ASV_tax %>% filter(Primer == "ITS")
x <- list(
  EIA = test %>% filter(grepl("94SIA", CCA.ID)) %>% as.data.frame() %>% select(Class) %>% unique() %>% .[,1],
  EIB =   test %>% filter(grepl("94SIB", CCA.ID)) %>% as.data.frame() %>% select(Class) %>% unique() %>% .[,1],
  EIC =   test %>% filter(grepl("94SIC", CCA.ID)) %>% as.data.frame() %>% select(Class) %>% unique() %>% .[,1]
)

ggVennDiagram(x)

#General assigned taxa
#How many Phyla per gene region? Minus 1 for 16S' "Bacteria incertae sedis"
ASV_tax %>% select(Primer, Phylum) %>% filter(is.na(Phylum)==F & Phylum != "" & Phylum != "Bacteria incertae sedis") %>% unique() %>% group_by(Primer) %>% summarise(PhylaNo = length(Phylum))
#micro and macro algae
algae <- rbind(ASV_tax %>% filter(grepl("phyta", Phylum)), ASV_tax %>% filter(grepl("Cyanobacteria", Phylum)), ASV_tax %>% filter(grepl("Dinophyceae", Class)))

algae %>% select(Primer, Phylum) %>% unique() %>% group_by(Primer) %>% summarise(algaeNo = length(Phylum))

algae %>% select(Primer, Class) %>% filter(., Class != "") %>% unique() %>% group_by(Primer) %>% summarise(algaeNo = length(Class))

algae %>% mutate(algae = rep("Yes")) %>% rbind(., ASV_tax %>% anti_join(., algae[,c("Feature.ID")]) %>% mutate(algae = rep("No"))) %>% group_by(Primer, algae) %>% select(Primer, Feature.ID, algae) %>% unique() %>% summarise(ASVs = length(Feature.ID)) %>% pivot_wider(., names_from = algae, values_from = ASVs) %>% mutate(perc = Yes/(Yes+No)*100)

algae %>% mutate(algae = rep("Yes")) %>% rbind(., ASV_tax %>% anti_join(., algae[,c("Feature.ID")]) %>% mutate(algae = rep("No"))) %>% group_by(Primer, algae) %>% select(Primer, Reads_nrm, algae) %>% aggregate(Reads_nrm~., ., sum) %>% pivot_wider(., names_from = algae, values_from = Reads_nrm) %>% mutate(perc = Yes/(Yes+No)*100)

algae %>% filter(., Primer == "23S" | Primer == "rbcL") %>% filter(Class != "") %>% aggregate(Reads_nrm~ID.no+From+Class+Primer, ., sum) %>% pivot_wider(., values_from = Reads_nrm, names_from = From, values_fill = 0) %>% mutate(moreX = Inner/Outer) %>% View

#TableS6
algae %>% mutate(algae = rep("Yes")) %>% rbind(., ASV_tax %>% anti_join(., algae[,c("Feature.ID")]) %>% mutate(algae = rep("No"))) %>% mutate(Identified = if_else(grepl("uncultured", sscinames)==T & grepl("[A-Z]", Phylum)==T, "Uncultured: assigned", if_else(grepl("uncultured", sscinames)==T & grepl("[A-Z]", Phylum)==F, "Uncultured: unassigned", if_else(is.na(ASV_tax$qstart)==F, "Yes", "No")))) %>% write.table(., "./TableS6.tsv", quote = F, row.names = F, sep = "\t")

ASV_tax %>% select(Primer, Class, Reads_nrm, Feature.ID, ID.no, From) %>% filter(is.na(Class)==F & Class != "") %>% aggregate(Reads_nrm~., ., sum) %>% group_by(Primer, ID.no, From) %>% summarise(Reads = sum(Reads_nrm), ASVs = length(Feature.ID)) %>% ungroup() %>% group_by(Primer, ID.no, From) %>% arrange(desc(Reads)) %>% slice_head(., n = 3) %>% View

#average results for GenBank categories for Figure 1
ASV_tax %>% mutate(Identified = if_else(grepl("uncultured", sscinames)==T & grepl("[A-Z]", Phylum)==T, "Uncultured: assigned", if_else(grepl("uncultured", sscinames)==T & grepl("[A-Z]", Phylum)==F, "Uncultured: unassigned", if_else(is.na(ASV_tax$qstart)==F, "Yes", "No")))) %>% group_by(CCA.ID, Primer, Identified) %>% summarise(total = sum(Reads_nrm), ASVs_n=length(unique(Feature.ID))) %>% ungroup() %>% group_by(CCA.ID, Primer) %>% mutate(Reads=total/sum(total)*100, ASVs=ASVs_n/sum(ASVs_n)*100) %>% pivot_longer(., cols = 6:7, names_to = "Perc") %>% select(2,3,6,7) %>% group_by(Primer, Identified, Perc) %>% summarise(mean = mean(value)) %>% filter(Perc == "ASVs")

FigureS1 <- 
ASV_tax %>% mutate(Identified = if_else(grepl("uncultured", sscinames)==T & grepl("[A-Z]", Phylum)==T, "Uncultured: assigned", if_else(grepl("uncultured", sscinames)==T & grepl("[A-Z]", Phylum)==F, "Uncultured: unassigned", if_else(is.na(ASV_tax$qstart)==F, "Yes", "No")))) %>% group_by(CCA.ID, Primer, Identified) %>% summarise(total = sum(Reads_nrm), ASVs_n=length(unique(Feature.ID))) %>% ungroup() %>% group_by(CCA.ID, Primer) %>% mutate(Reads=total/sum(total)*100, ASVs=ASVs_n/sum(ASVs_n)*100) %>% pivot_longer(., cols = 6:7, names_to = "Perc") %>%
  ggplot(.)+
  geom_bar(aes(x=Perc, y=value, fill=factor(Identified, levels=c("Yes", "Uncultured: assigned", "Uncultured: unassigned", "No"))), stat = "identity", color = "black", size = 0.1)+
  scale_fill_brewer(palette = "Paired")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 7, colour = "black"),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank(),
        strip.background = element_rect(fill = "white", color = "black"),
        panel.background = element_rect(fill = "white", color = "black"),
        axis.text.y = element_text(size = 7, colour = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 8, colour = "black"),
        strip.text.x = element_text(angle = 90, size = 7, colour = "black"),
        strip.text.y = element_text(angle = 0, size = 7, colour = "black"),
        legend.position = "bottom",
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 7, colour = "black"),
        legend.justification = "top",
        panel.spacing = unit(0.1, "lines"),
        legend.key.size = unit(0.3, "cm")
  )+
  facet_grid(Primer~CCA.ID, scales = "free", space = "free")+
  ylab("Relative abundance (%)")+
  labs(fill="GenBank Identification:")

#ggsave2(filename = "./FigureS1.jpg", plot = FigureS1, width = 17.5, height = 15, units = "cm",dpi = 300)

col_grad <- c((brewer.pal(8, "Dark2")[c(1:8)]), (brewer.pal(12, "Paired")[c(1:12)]), rep("black", times= 26))
col_grad <- c((brewer.pal(8, "Dark2")[c(1:8)]), (brewer.pal(12, "Paired")[c(1:12)]),  brewer.pal(8, "Pastel2")[c(1:8)],brewer.pal(8, "Dark2")[c(1:8)], (brewer.pal(12, "Paired")[c(1:12)]), "black")

all.bar <- ASV_tax %>% select(Primer, Feature.ID, Phylum, Reads_nrm, From) %>% group_by(Primer, From, Phylum) %>% summarise(tot_ASVs=length(unique(Feature.ID)), tot_Reads=sum(Reads_nrm)) %>% filter(., grepl("[A-Z]", Phylum)==T) %>% ungroup(Phylum) %>% mutate(Reads=tot_Reads/sum(tot_Reads)*100, ASVs=tot_ASVs/sum(tot_ASVs)*100) %>% select(-tot_ASVs,-tot_Reads) %>% pivot_longer(., cols = c("Reads", "ASVs"), names_to = "perc", values_to = "Relative abundance (%)") %>%
  ggplot(.)+
  geom_bar(stat = "identity", aes(x=From, y=`Relative abundance (%)`, fill=Phylum), color="black", size=0.1)+
  facet_grid(perc~Primer)+
  scale_fill_manual(values = col_grad)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 7, colour = "black"),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank(),
        strip.background = element_rect(fill = "white", color = "black"),
        panel.background = element_rect(fill = "white", color = "black"),
        axis.text.y = element_text(size = 7, colour = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 8, colour = "black"),
        strip.text.x = element_text(angle = 0, size = 7, colour = "black"),
        strip.text.y = element_text(angle = 0, size = 7, colour = "black"),
        legend.position = "right",
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 7, colour = "black"),
        legend.justification = "top",
        panel.spacing.x = unit(0.1, "lines"),
        panel.spacing.y = unit(0.4, "lines"),
        legend.key.size = unit(0.3, "cm")
  )+
  guides(fill=guide_legend(ncol=1,byrow=F))+
scale_y_continuous(expand = expansion(mult = c(0, 0)))

#ggsave2(filename = "./all.bar.pdf", plot = all.bar, width = 17.5, height = 15, units = "cm",dpi = 600)

col.pal <- read.xlsx("./col.pal.xlsx", sheetIndex = 2, header = T) %>% .[2:nrow(.),c("Colour")]

#algae only - all primers
FigureS2 <- algae %>% filter(., Class != "") %>% select(Feature.ID, Class, Reads_nrm, From, CCA.ID, Primer) %>% group_by(Primer, CCA.ID, From, Class) %>% summarise(tot_ASVs=length(unique(Feature.ID)), tot_Reads=sum(Reads_nrm)) %>% filter(., grepl("[A-Z]", Class)==T) %>% ungroup(Class) %>% mutate(Reads=tot_Reads/sum(tot_Reads)*100, ASVs=tot_ASVs/sum(tot_ASVs)*100) %>% select(-tot_ASVs,-tot_Reads) %>% pivot_longer(., cols = c("Reads", "ASVs"), names_to = "perc", values_to = "Relative abundance (%)") %>% mutate(Class = if_else(`Relative abundance (%)` > 5, Class, "Other")) %>% 
  #aggregate(`Relative abundance (%)`~., ., sum) %>%
  ggplot(.)+
  geom_bar(stat = "identity", aes(x=perc, y=`Relative abundance (%)`, fill=Class), color="black", size=0.1)+
  facet_grid(Primer~CCA.ID)+
  scale_fill_manual(values = col.pal)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 7, colour = "black"),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank(),
        strip.background = element_rect(fill = "white", color = "black"),
        panel.background = element_rect(fill = "white", color = "black"),
        axis.text.y = element_text(size = 7, colour = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 8, colour = "black"),
        strip.text.x = element_text(angle = 90, size = 7, colour = "black"),
        strip.text.y = element_text(angle = 0, size = 7, colour = "black"),
        legend.position = "bottom",
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 7, colour = "black"),
        legend.justification = "top",
        panel.spacing.x = unit(0.1, "lines"),
        panel.spacing.y = unit(0.4, "lines"),
        legend.key.size = unit(0.3, "cm")
  )+
  guides(fill=guide_legend(nrow=3,byrow=T))+
  scale_y_continuous(expand = expansion(mult = c(0, 0)))

ggsave2(filename = "./FigureS2.jpg", plot = FigureS2, width = 17, height = 17, units = "cm",dpi = 300)

#23S figures

GBid23S <- ASV_tax %>%  mutate(Identified = if_else(grepl("uncultured", sscinames)==T & grepl("[A-Z]", Phylum)==T, "Uncultured: assigned", if_else(grepl("uncultured", sscinames)==T & grepl("[A-Z]", Phylum)==F, "Uncultured: unassigned", if_else(is.na(ASV_tax$qstart)==F, "Yes", "No")))) %>% group_by(CCA.ID, Primer, Identified) %>% summarise(total = sum(Reads_nrm), ASVs_n=length(unique(Feature.ID))) %>% ungroup() %>% group_by(CCA.ID, Primer) %>% mutate(Reads=total/sum(total)*100, ASVs=ASVs_n/sum(ASVs_n)*100) %>% pivot_longer(., cols = 6:7, names_to = "Perc") %>% filter(., Primer == "23S") %>%
  ggplot(.)+
  geom_bar(aes(x=Perc, y=value, fill=factor(Identified, levels=c("Yes", "Uncultured: assigned", "Uncultured: unassigned", "No"))), stat = "identity", color = "black", size = 0.1)+
  scale_fill_brewer(palette = "Paired")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 7, colour = "black"),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank(),
        strip.background = element_rect(fill = "white", color = "black"),
        panel.background = element_rect(fill = "white", color = "black"),
        axis.text.y = element_text(size = 7, colour = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 8, colour = "black"),
        strip.text.x = element_text(angle = 90, size = 7, colour = "black"),
        strip.text.y = element_text(angle = 0, size = 7, colour = "black"),
        legend.position = "bottom",
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 7, colour = "black"),
        legend.justification = "top",
        panel.spacing = unit(0.1, "lines"),
        legend.key.size = unit(0.3, "cm")
  )+
  facet_grid(~CCA.ID, scales = "free", space = "free")+
  ylab("Relative abundance (%)")+
  labs(fill="GenBank Identification:")
  
ggsave2(filename = "./JoPhycology/GBid23S.pdf", plot = GBid23S, width = 17.5, height = 6, units = "cm",dpi = 600)

test_col <- c(brewer.pal(12, "Set3")[c(1:12)], brewer.pal(12, "Paired")[c(1:12)])

ASVGB23S <- algae %>% filter(., Primer == "23S") %>% filter(., Class != "") %>% select(Feature.ID, Class, Reads_nrm, From, CCA.ID) %>% group_by(CCA.ID, From, Class) %>% summarise(tot_ASVs=length(unique(Feature.ID)), tot_Reads=sum(Reads_nrm)) %>% filter(., grepl("[A-Z]", Class)==T) %>% ungroup(Class) %>% mutate(Reads=tot_Reads/sum(tot_Reads)*100, ASVs=tot_ASVs/sum(tot_ASVs)*100) %>% select(-tot_ASVs,-tot_Reads) %>% pivot_longer(., cols = c("Reads", "ASVs"), names_to = "perc", values_to = "Relative abundance (%)") %>% mutate(Class = if_else(`Relative abundance (%)` > 5, Class, "Other")) %>% 
  #aggregate(`Relative abundance (%)`~., ., sum) %>%
  ggplot(.)+
  geom_bar(stat = "identity", aes(x=perc, y=`Relative abundance (%)`, fill=Class), color="black", size=0.1)+
  facet_grid(~CCA.ID)+
  scale_fill_manual(values = test_col)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 7, colour = "black"),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank(),
        strip.background = element_rect(fill = "white", color = "black"),
        panel.background = element_rect(fill = "white", color = "black"),
        axis.text.y = element_text(size = 7, colour = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 8, colour = "black"),
        strip.text.x = element_text(angle = 90, size = 7, colour = "black"),
        strip.text.y = element_text(angle = 0, size = 7, colour = "black"),
        legend.position = "bottom",
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 7, colour = "black"),
        legend.justification = "top",
        panel.spacing = unit(0.1, "lines"),
        legend.key.size = unit(0.3, "cm")
  )+
  guides(fill=guide_legend(nrow=2,byrow=T))+
  scale_y_continuous(expand = expansion(mult = c(0, 0)))

ggsave2(filename = "./JoPhycology/ASVGB23S.pdf", plot = ASVGB23S, width = 17.5, height = 6, units = "cm",dpi = 600)


#tree 43 23S
tree <- read.tree("./43_unrooted_trees/23S_tree.nwk")
circ <- ggtree(tree, layout = "circular", branch.length = "branch.length", size = 0.01)

tip.label.order <- tree[["tip.label"]] %>% as.data.frame() %>% rename(., tip.label=1)

tree_meta <- ASV_tax[,c("Feature.ID", "Phylum", "Primer", "CCA.ID")] %>% filter(Primer == "23S") %>% mutate(Phylum=if_else(Phylum=="" | is.na(Phylum)==T | Phylum == "Bacteria incertae sedis", "Unassigned", Phylum)) %>% filter(., grepl("^4", CCA.ID)) %>% mutate(From=if_else(grepl("O", CCA.ID)==T, "Outer", "Inner")) %>% select(-CCA.ID) %>% unique() %>% aggregate(From~., ., toString) %>% unique() %>%  rename(., tip.label=Feature.ID) %>% unique() %>% left_join(tip.label.order, .) %>% mutate(From=if_else(grepl(",", From)==T, "Both", From))

rownames(tree_meta) <- tree_meta$tip.label

tree_meta$node <-  rep(1:nrow(tree_meta))

#options <- c(brewer.pal(12, "Paired")[c(1:12)],  brewer.pal(12, "Set3")[c(1:12)], brewer.pal(8, "Set2")[c(1:8)])
  #c("dark grey", "black", "light grey")

p1 <- gheatmap(circ, tree_meta %>% select(From), offset=0, width=.1, colnames = F) +
  scale_fill_manual(values = c("#BDBDBD", "#000000", "red"))+
  labs(fill="Layer:")+
  theme(legend.title = element_text(size = 7),
        legend.text = element_text(size = 7, colour = "black"),
        legend.key.size = unit(0.3, "cm"))

p2 <- p1+
  new_scale_fill()+
  geom_hilight(data = tree_meta, mapping=aes(node=node, fill=Phylum), type = "rect", align = "right", alpha = 0.5)+
  scale_fill_manual(values = c("#D9D9D9", #light grey #Acidobacteria
                    "#E6AB02", #light brown #Bacillariophyta
                    "#666666", #grey #Bacteroidetes
                    "#FFED6F", #light yellow #Cercozoa
                    "#66A61E",#dark green #Chlorophyta
                    "#FCCDE5", #blue #Cryptophyta
                    "#CCEBC5",#light green #Cyanobacteria
                    "#BEBADA", #light orange #Haptophyta
                    "#80B1D3", #light purple #Lentisphaerae
                    "#B3DE69", #purple #Myzozoa
                    "#A6761D", #brown #Ochrophyta
                    "#D95F02", #orange-red #Proteobacteria
                    "#FB8072",#red #Rhodophyta
                    "#1B9E77", #yellow-brown #Tracheophyta
                    "white",#Unassigned
                    "#8DD3C7" #grey #Verrucomicrobia
                    ))+
  geom_point(data = circ$data %>% mutate(bootstrap = if_else(grepl("0.9.*", label), "Bootstrap > 0.9", NA)) %>% filter(is.na(bootstrap)==F), aes(color=bootstrap), size=0.1)+
  scale_color_manual(values = brewer.pal(9, "YlOrRd")[c(9)])+
  labs(fill="Phylum:", color = "")+
  theme(legend.title = element_text(size = 7),
        legend.text = element_text(size = 7, colour = "black"),
        legend.key.size = unit(0.3, "cm")
  )

ggsave(filename = "./JoPhycology/CCA643_tree_23S.pdf", plot = p2, width = 24, height = 18, units = "cm",dpi = 600)

#tree 94 23S
tree <- read.tree("./94_unrooted_trees/23S_tree.nwk")
circ <- ggtree(tree, layout = "circular", branch.length = "branch.length", size = 0.01)

tip.label.order <- tree[["tip.label"]] %>% as.data.frame() %>% rename(., tip.label=1)

tree_meta <- ASV_tax[,c("Feature.ID", "Phylum", "Primer", "CCA.ID")] %>% filter(Primer == "23S") %>% mutate(Phylum=if_else(Phylum=="" | is.na(Phylum)==T | Phylum == "Bacteria incertae sedis", "Unassigned", Phylum)) %>% filter(., grepl("^9", CCA.ID)) %>% mutate(From=if_else(grepl("O", CCA.ID)==T, "Outer", "Inner")) %>% select(-CCA.ID) %>% unique() %>% aggregate(From~., ., toString) %>% unique() %>%  rename(., tip.label=Feature.ID) %>% unique() %>% left_join(tip.label.order, .) %>% mutate(From=if_else(grepl(",", From)==T, "Both", From))

rownames(tree_meta) <- tree_meta$tip.label

tree_meta$node <-  rep(1:nrow(tree_meta))

#options <- c(brewer.pal(12, "Paired")[c(1:12)],  brewer.pal(12, "Set3")[c(1:12)], brewer.pal(8, "Set2")[c(1:8)])
#c("dark grey", "black", "light grey")

p1 <- gheatmap(circ, tree_meta %>% select(From), offset=0, width=.1, colnames = F) +
  scale_fill_manual(values = c("#BDBDBD", "#000000", "red"))+
  labs(fill="Layer:")+
  theme(legend.title = element_text(size = 7),
        legend.text = element_text(size = 7, colour = "black"),
        legend.key.size = unit(0.3, "cm"))

p2 <- p1+
  new_scale_fill()+
  geom_hilight(data = tree_meta, mapping=aes(node=node, fill=Phylum), type = "rect", align = "right", alpha = 0.5)+
  scale_fill_manual(values = c(#"#D9D9D9", #light grey #Acidobacteria
                               "#E6AB02", #light brown #Bacillariophyta
                               "#666666", #grey #Bacteroidetes
                               "#FFED6F", #light yellow #Cercozoa
                               "#66A61E",#dark green #Chlorophyta
                               "#FCCDE5", #blue #Cryptophyta
                               "#CCEBC5",#light green #Cyanobacteria
                               "#BEBADA", #light orange #Haptophyta
                               #"#80B1D3", #light purple #Lentisphaerae
                               "#B3DE69", #purple #Myzozoa
                               "#A6761D", #brown #Ochrophyta
                               "#8DD3C7", #Plantomucetes
                               "#BEBADA", #Prasinodermatophyta
                               "#D95F02", #orange-red #Proteobacteria
                               "#FB8072",#red #Rhodophyta
                               #"#1B9E77", #yellow-brown #Tracheophyta
                               "white",#Unassigned
                               "#8DD3C7" #grey #Verrucomicrobia
  ))+
  geom_point(data = circ$data %>% mutate(bootstrap = if_else(grepl("0.9.*", label), "Bootstrap > 0.9", NA)) %>% filter(is.na(bootstrap)==F), aes(color=bootstrap), size=0.1)+
  scale_color_manual(values = brewer.pal(9, "YlOrRd")[c(9)])+
  labs(fill="Phylum:", color = "")+
  theme(legend.title = element_text(size = 7),
        legend.text = element_text(size = 7, colour = "black"),
        legend.key.size = unit(0.3, "cm")
  )

ggsave(filename = "./JoPhycology/CCA694_tree_23S.pdf", plot = p2, width = 24, height = 18, units = "cm",dpi = 600)

# p2 <- p1+
#   new_scale_fill()
# 
# p3 <- gheatmap(p2, tree_meta %>% select(Phylum) %>% mutate(Phylum = if_else(Phylum == "Unassigned" | Phylum == "Bacteria incertae sedis", "Unassigned", Phylum)), offset=2, width=.1, colnames = F) +
#   scale_fill_manual(values = inner_join(cols23S %>% arrange(Phylum), tree_meta %>% select(Phylum) %>% mutate(Phylum = if_else(Phylum == "Unassigned" | Phylum == "Bacteria incertae sedis", "Unassigned", Phylum)) %>% unique()) %>% .[,2])+
#   labs(fill="Phylum")
# 
# p3

#ggsave(filename = "./JoPhycology/CCA643_tree_23S.pdf", plot = p3, width = 24, height = 18, units = "cm",dpi = 600)


ASV_tax%>% select(Feature.ID, Reads_nrm, sscinames, Primer) %>% unique() %>% mutate(uncul = if_else(grepl("uncultured", sscinames)==T, "Uncultured", "species")) %>% group_by(Primer, uncul) %>% summarise(Reads=sum(Reads_nrm), ASVs=length(Feature.ID))
  
#order of interest
OOI <- data.frame(
  Order = c("Anaulales", "Bryopsidales", "Ceramiales", "Chlorachniida", "Corallinales", "Dictyochales", "Dictyotales", "Gelidiales", "Laminariales", "Naviculales", "Sarcinochrysidales", "Sphacelariales", "Thraustochytrida", "Ulvales", "Bryopsidales", "Corallinales", "Laminariales", "Naviculales", "Bacillariales", "Bryopsidales", "Chlorachniida", 'Corallinales', "Dictyotales", "Gymnodiniales", "Hapalidiales", "Laminariales", "Lyrellales", "Mamiellales", "Naviculales", "Peyssonneliales", "Sarcinochrysidales", "Sphaeropleales", "Surirellales", "Toxariales", "Ulvales", "Acrochaetiales", "Bacillariales", "Ceramiales", "Chlorachniida", "Dothideales", "Ectocarpales", "Gigartinales", "Nemaliales", "Saccharomycetales", "Ulvales", "Wallemiales", "Xylariales"),
  Group = c("Diatoms", "Green algae", "Red algae", "Cercozoa", "Coralline algae", "Silicoflagellate", "Brown algae", "Red algae", "Brown algae", "Diatoms", "Palmelloid or filamentous heterokonts", "Brown algae", "Stramenopiles", "Green algae", "Green algae", "Coralline algae", "Brown algae", "Diatoms", "Diatoms", "Green algae", "Cercozoa", "Coralline algae", "Silicoflagellates", 'Dinoflagellates', "Coralline algae", "Brown algae", "Diatoms", "Green algae", "Diatoms", "Red algae", "Heterokont algae", "Green algae", "Diatoms", "Diatoms", 'Green algae', "Red algae", "Diatoms", "Red algae", "Cercozoa", "Bitunicate fungi", "Brown algae", "Red algae", "Red algae", "Fungi", "Green algae", "Fungi", "Fungi")
)

OOI <- unique(OOI)

#rm
OOI <- OOI %>% filter(., Group != "Silicoflagellates") %>% filter(., Group != "Palmelloid or filamentous heterokonts")

inner_join(ASV_tax, OOI) %>% select(Phylum) %>% unique()

ASV_OOI %>% filter(., Primer == "rbcL") %>% aggregate(Reads_nrm~CCA.ID+Genus+Group, ., sum) %>% group_by(CCA.ID) %>% mutate(RRA=Reads_nrm/sum(Reads_nrm)*100) %>% filter(., Genus != "") %>%
  ggplot(.)+
  geom_point(stat= "identity", aes(x=CCA.ID, y= fct_rev(Genus), size = RRA), shape = 21)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 7, colour = "black"),
        axis.ticks.x = element_blank(),
        strip.background = element_rect(fill = "white", color = "black"),
        panel.background = element_rect(fill = "white", color = "black"),
        axis.text.y = element_text(size = 7, colour = "black", face = "italic"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 8, colour = "black"),
        strip.text.x = element_text(angle = 90, size = 7, colour = "black"),
        strip.text.y = element_text(angle = 0, size = 7, colour = "black"),
        legend.position = "bottom",
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 7, colour = "black"),
        legend.justification = "top",
        panel.spacing = unit(0.1, "lines"),
        legend.key.size = unit(0.3, "cm")
  )+
  facet_grid(Group~., space = "free", scales = "free")

Cyanobacteria.Genus <- 
ASV_tax %>% filter(., grepl("^[1-9].E", CCA.ID)) %>% select(Primer, CCA.ID, Reads_nrm, Phylum, Genus, ID.no, From) %>% mutate(ID_From=paste(ID.no, From, sep = "_")) %>% group_by(Primer, ID.no, From) %>% mutate(RRA=Reads_nrm/sum(Reads_nrm)*100) %>% filter(., grepl("[A-Z]", Genus)==T) %>% select(-Reads_nrm) %>% aggregate(RRA~., ., sum) %>% group_by(Primer, ID.no, From, Genus) %>% mutate(FOO=length(CCA.ID)) %>% select(-CCA.ID) %>% aggregate(RRA~., ., sum) %>% 
  #filter(., grepl("Ochrophyta", Phylum)) %>% left_join(., ASV_tax %>% filter(., grepl("Ochrophyta", Phylum)) %>% mutate(Cat=if_else(Class=="Bacillariophyceae", "Diatoms", "Other\nOchrophyta")) %>% select(Cat, Genus) %>% unique()) %>%
  filter(., grepl("Cyanobacteria", Phylum)) %>%
  ggplot(.)+
  geom_point(stat= "identity", aes(x=ID_From, y= fct_rev(Genus), size = FOO, fill= RRA), shape = 21)+
  scale_fill_gradientn(colours = c("dark blue", "light blue", "light yellow", "orange", "red"))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 7, colour = "black"),
        axis.ticks.x = element_blank(),
        strip.background = element_rect(fill = "white", color = "black"),
        panel.background = element_rect(fill = "white", color = "black"),
        axis.text.y = element_text(size = 7, colour = "black", face = "italic"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 8, colour = "black"),
        strip.text.x = element_text(angle = 0, size = 7, colour = "black"),
        strip.text.y = element_text(angle = 0, size = 7, colour = "black"),
        legend.position = "bottom",
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 7, colour = "black"),
        legend.justification = "top",
        panel.spacing = unit(0.1, "lines"),
        legend.key.size = unit(0.3, "cm")
  )+
  facet_grid(.~Primer, space = "free", scales = "free")+
  ylab("Genus")

Dinoflagellates.Genus <- 
  ASV_tax %>% filter(., grepl("^[1-9].E", CCA.ID)) %>% select(Primer, CCA.ID, Reads_nrm, Class, Genus, ID.no, From) %>% mutate(ID_From=paste(ID.no, From, sep = "_")) %>% group_by(Primer, ID.no, From) %>% mutate(RRA=Reads_nrm/sum(Reads_nrm)*100) %>% filter(., grepl("[A-Z]", Genus)==T) %>% select(-Reads_nrm) %>% aggregate(RRA~., ., sum) %>% group_by(Primer, ID.no, From, Genus) %>% mutate(FOO=length(CCA.ID)) %>% select(-CCA.ID) %>% aggregate(RRA~., ., sum) %>% 
  #filter(., grepl("Ochrophyta", Phylum)) %>% left_join(., ASV_tax %>% filter(., grepl("Ochrophyta", Phylum)) %>% mutate(Cat=if_else(Class=="Bacillariophyceae", "Diatoms", "Other\nOchrophyta")) %>% select(Cat, Genus) %>% unique()) %>%
  filter(., grepl("Dinophyceae", Class)) %>%
  ggplot(.)+
  geom_point(stat= "identity", aes(x=ID_From, y= fct_rev(Genus), size = FOO, fill= RRA), shape = 21)+
  scale_fill_gradientn(colours = c("dark blue", "light blue", "light yellow", "orange", "red"))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 7, colour = "black"),
        axis.ticks.x = element_blank(),
        strip.background = element_rect(fill = "white", color = "black"),
        panel.background = element_rect(fill = "white", color = "black"),
        axis.text.y = element_text(size = 7, colour = "black", face = "italic"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 8, colour = "black"),
        strip.text.x = element_text(angle = 0, size = 7, colour = "black"),
        strip.text.y = element_text(angle = 0, size = 7, colour = "black"),
        legend.position = "bottom",
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 7, colour = "black"),
        legend.justification = "top",
        panel.spacing = unit(0.1, "lines"),
        legend.key.size = unit(0.3, "cm")
  )+
  facet_grid(.~Primer, space = "free", scales = "free")+
  ylab("Genus")

ggsave2(filename = "./Dinoflag.Genus.pdf", plot = Dinoflagellates.Genus, width = 17.5, height = 20, units = "cm",dpi = 600)

#Figures
species_colours <- c((brewer.pal(8, "Set2")[c(1:8)]), (brewer.pal(12, "Paired")[c(1:12)]),  "black", "white", c("#238A8DFF", "#56C677FF", "#3F4758FF", "#FDE725FF"), (brewer.pal(12, "Set3")[c(1:12)]))


ASV_tax %>% select(CCA.ID, Primer, Order, Reads_nrm) %>% filter(., Order != "") %>% aggregate(Reads_nrm~., ., sum) %>% group_by(Primer, CCA.ID) %>% mutate(RRA = Reads_nrm/sum(Reads_nrm)*100) %>% 
  #filter(., Primer == "tufA"| Primer == "rbcL") %>%
  filter(., Primer == "COI"| Primer == "ITS") %>%
  #filter(., Primer == "23S" | Primer == "18S") %>%
  #filter(., Primer == "16S") %>%
  filter(., RRA > 5) %>%
  ggplot(.)+
  geom_bar(stat = "identity", aes(x = CCA.ID, y = RRA, fill = Order), color = "black")+
  theme_bw()+
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))+
  scale_fill_manual(values = species_colours)+
  facet_wrap(~Primer,nrow = 7)

#Supp data
colnames(ASV_tax)
write.xlsx2(ASV_tax[, c("Feature.ID", "Primer", "CCA.ID", "ID.no", "From", "Pres", "Reads_nrm", "qlen", "qcovhsp", "pident", "Phylum", "Class", "Order", "Family", "Genus", "sscinames")], "./Data_CCA_2022.xlsx")
