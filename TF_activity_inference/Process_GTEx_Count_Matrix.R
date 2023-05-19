### Generate count matrix for Muscle_Skeletal

counts <- CePa::read.gct("GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct") #Downloaded from the GTEx protal
counts <- as.data.frame(counts)
sample_attributes <- data.table::fread("GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt") #Downloaded from the GTEx protal

gene.info <- rtracklayer::import("gencode.v26.GRCh38.genes.gtf") #Downloaded from GENCODE
protein.coding.genes <- gene.info[gene.info$type == "gene" & gene.info$gene_type == "protein_coding",]

rownames(counts) <- counts$Name
counts <- counts[protein.coding.genes$gene_id,]

tissue <- "Muscle_Skeletal"
samples <- sample_attributes$SAMPID[sample_attributes$SMTSD == tissue & sample_attributes$SMAFRZE == "RNASEQ"]
counts_tissue <- counts[,colnames(counts) %in% samples]
colnames(counts_tissue) <- paste0("GTEX-", sapply(strsplit(colnames(counts_tissue), "-"), "[", 2))
counts_tissue$gene <- rownames(counts_tissue)


### Downsampling

sample_sums <- apply(counts_tissue[,samples], 2, sum)
size_min <- min(sample_sums)
  
downSampling <- function(i) { rmultinom(n = 1, size = size_min, prob = counts_tissue[,samples[i]]/sample_sums[i]) }
counts_down <- as.data.frame(sapply(1:length(samples), downSampling))
colnames(counts_down) <- samples
counts_down <- cbind(counts_tissue[,"gene"], counts_down)
  
write.csv(counts_down, paste0("example_data/", tissue, ".counts.downsampled.csv"), row.names = T, quote = F)
