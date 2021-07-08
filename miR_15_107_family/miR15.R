library(gplots)
library(Rtsne)
library(png)
library(sva)
library(plyr)
library(cowplot)
library(edgeR)
library(limma)
library(Glimma)
library(ggplot2)
library(RColorBrewer)
library(ggfortify)
library(RUVcorr)
library(reshape)
library(data.table)
library(EDASeq)
library(biomaRt)
library(ggvenn)
library(multiMiR)
library(org.Hs.eg.db)
library(clusterProfiler)
library(idmap3)
library(oligo)
library(survival)
library(survminer)
library(dplyr)
library(stringr)
library(htmlTable)
library(GGally)
library(network)
library(purrr)
library(tidyverse)
library(rsample)
library(pROC)

output <- 'output'
p.adj <- 'BH'
pvalue <- 0.05

mir.gene <- c(
  'hsa-mir-103a-1',
  'hsa-mir-103a-2',
  'hsa-mir-107',
  'hsa-mir-15a',
  'hsa-mir-15b', 
  'hsa-mir-16-1',
  'hsa-mir-195',
  'hsa-mir-497',
  'hsa-mir-503',
  'hsa-mir-424',
  'hsa-mir-646',
  'hsa-mir-6838'
  )
mature.mir <- c(
  'hsa-miR-103a-3p', 
  'hsa-miR-103a-3p', 
  'hsa-miR-107', 
  'hsa-miR-15a-5p', 
  'hsa-miR-15b-5p',
  'hsa-miR-16-5p',
  'hsa-miR-195-5p',
  'hsa-miR-497-5p',
  'hsa-miR-503-5p',
  'hsa-miR-424-5p',
  'hsa-miR-646',
  'hsa-miR-6838-5p')

mir.table <- data.frame(`miRNA`=mir.gene, `mature_miRNA`=mature.mir)

ss <- lapply(mir.table$miRNA, function(x){
	ss <- str_replace_all(x, '-', '_')
	return(ss)
})

tt <- lapply(mir.table$mature_miRNA, function(x){
	ss <- str_replace_all(x, '-', '_')
	return(ss)
})

mir.table$miRNA_gene <- unlist(ss)
mir.table$miRNA_ID <- unlist(tt)


data3 <- mirna.df

common.samples <- intersect(colnames(data3), survival.df$sample)

common.samples <- intersect(common.samples, pheno.df$`submitter_id.samples`)

data3 <- data3[,match(common.samples, colnames(data3))]

batch <- pheno.df[match(common.samples, pheno.df$`submitter_id.samples`),]$`sample_type.samples`


data3 <- normalizeQuantiles(data3, ties=TRUE)

data3 <- as.data.frame(data3)



mir.15 <- c('hsa-mir-103a-1','hsa-mir-103a-2','hsa-mir-15a','hsa-mir-15b','hsa-mir-107','hsa-mir-16-1','hsa-mir-195','hsa-mir-497','hsa-mir-503','hsa-mir-424','hsa-mir-646','hsa-mir-6838')

covariates <- c(mir.15)

mirna.df.sub <- data3[covariates,]

covariates <- lapply(covariates, function(x){
  return(str_replace_all(x, '-', '_'))
  })

covariates <- unlist(covariates)

rownames(mirna.df.sub) <- covariates

covariates <- mir.table$miRNA_gene

mirna.df.sub <- mirna.df.sub[covariates,]

mirna.df.sub <- 2^mirna.df.sub - 1

mirna.df.sub$Symbol <- mir.table$miRNA_ID
mirna.df.sub <- data.table(mirna.df.sub)

mirna.df.sub <- mirna.df.sub[, lapply(.SD,sum), by=.(Symbol)]
mirna.df.sub <- as.data.frame(mirna.df.sub)
rownames(mirna.df.sub) <- mirna.df.sub$Symbol
mirna.df.sub <- mirna.df.sub[,-1]

mirna.df.sub <- log1p(mirna.df.sub)

# Figure
for(sample in unique(pheno.df$TCGA_ID)){

  # sample <- 'TCGA_ACC'
  mirna.df.sub.t <- as.data.frame(t(mirna.df.sub))
  common.samples <- intersect(rownames(mirna.df.sub.t), pheno.df[which(pheno.df$TCGA_ID == sample),]$submitter_id.samples)

  common.samples <- intersect(common.samples, survival.df$sample)

  mirna.df.sub.t <- mirna.df.sub.t[match(common.samples, rownames(mirna.df.sub.t)),]
  mirna.df.sub.t <- cbind(mirna.df.sub.t, survival.df[match(common.samples, survival.df$sample),][,-1])

  if(dim(mirna.df.sub.t)[1]==0) next

  print('Finish combining sample')

  common.samples <- intersect(rownames(mirna.df.sub.t), pheno.df$`submitter_id.samples`)

  mirna.df.sub.t$Sample_type <- ''

  mirna.df.sub.t[match(common.samples, rownames(mirna.df.sub.t)),]$Sample_type <- pheno.df[match(common.samples, pheno.df$`submitter_id.samples`),]$`sample_type.samples`


  my_data <- mirna.df.sub.t
  names(my_data)[names(my_data) == "OS.time"] <- "OS_time"


  covariates <- unique(mir.table[which(mir.table$miRNA %in% mir.15),]$miRNA_ID)

  my_df <- my_data[,c(covariates, c('OS','_PATIENT','OS_time','TCGA_ID','Sample_type'))]
  if(dim(my_df)[1]==0)next

  # my_df <- my_df[which(my_df$Sample_type == 'Solid Tissue Normal'),]
   # my_df <- my_df[which(my_df$Sample_type == 'Primary Tumor'),]

  result <- NULL

  data <- my_df[,covariates]
  if(dim(data)[1]==0) next
  data <- melt(data)
  png(paste('output/miR15/TCGA.',sample, ".primary.tumor.miR15Family.boxplot.png", sep=""),height=5000, width=8000, res=600)
  p <- ggplot(data, aes(x=variable, y=value, color=variable, fill=variable)) + 
    # geom_boxplot() +
    geom_violin(trim = FALSE)
    # geom_dotplot(binaxis='y', stackdir='center', stackratio=1.5, dotsize=0.5, color='black', fill='black') +
  theme(legend.position="none", axis.text.x = element_text(angle = 0, vjust = 0.7, hjust=0.8)) +
  ylab("Gene Expression Level [log2(fpkm+1)]") +
  xlab("")
  print(p)
  dev.off()  
}

# Figure 
mirna.df.sub.t <- as.data.frame(t(mirna.df.sub))
# common.samples <- intersect(rownames(mirna.df.sub.t), pheno.df[which(pheno.df$TCGA_ID == sample),]$submitter_id.samples)

common.samples <- intersect(rownames(mirna.df.sub.t), survival.df$sample)

mirna.df.sub.t <- mirna.df.sub.t[match(common.samples, rownames(mirna.df.sub.t)),]
mirna.df.sub.t <- cbind(mirna.df.sub.t, survival.df[match(common.samples, survival.df$sample),][,-1])

print('Finish combining sample')

common.samples <- intersect(rownames(mirna.df.sub.t), pheno.df$`submitter_id.samples`)

mirna.df.sub.t$Sample_type <- ''

mirna.df.sub.t[match(common.samples, rownames(mirna.df.sub.t)),]$Sample_type <- pheno.df[match(common.samples, pheno.df$`submitter_id.samples`),]$`sample_type.samples`


my_data <- mirna.df.sub.t
names(my_data)[names(my_data) == "OS.time"] <- "OS_time"


covariates <- unique(mir.table[which(mir.table$miRNA %in% mir.15),]$miRNA_ID)

my_df <- my_data[,c(covariates, c('OS','_PATIENT','OS_time','TCGA_ID','Sample_type'))]


data <- NULL

for(sample in unique(pheno.df$TCGA_ID)){
  print(sample)
  tmp.data <- my_df[which(my_df$TCGA_ID == sample),]
  if(dim(tmp.data)[1]==0) next
  if(is.null(data)){
    data <- melt(tmp.data[,covariates])
    data$TCGA_ID <- sample
  }else{
    tmp <- melt(tmp.data[,covariates])
    tmp$TCGA_ID <- sample
    print(tmp)
    data <- rbind(data,tmp)
  }
}

png(paste('output/miR15/TCGA.all',".only.primary.tumor.miR15Family.boxplot.png", sep=""),height=5000, width=7000, res=400)
p <- ggplot(data, aes(x=variable, y=value, color=variable, fill=variable)) + 
  geom_boxplot() +
  # geom_violin(trim = FALSE) + 
   # facet_grid(.~TCGA_ID)+
   # geom_text(size=2) +
  # geom_dotplot(binaxis='y', stackdir='center', stackratio=1.5, dotsize=0.5, color='black', fill='black') +
theme(legend.position="none", axis.text.x = element_text(angle = 0, vjust = 0.7, hjust=0.8, size=12)) +
ylab("Gene Expression Level [log2(fpkm+1)]") +
xlab("")
print(p)
dev.off()  

# Figure

my_files <- list.files('output/Pan_Cancer/')

top.genes <- NULL
for(sample in unique(pheno.df$TCGA_ID)){
  name <- paste0('TCGA.',sample, '.contrast_Primary_Tumor_v_Solid_Tissue_Normal.miRNA.top.genes.csv')
  
  if(name %in% my_files){

    top.genes.tmp <- read.csv(paste0('output/Pan_Cancer/', name))
    top.genes.tmp$TCGA_ID <- sample

    if(is.null(top.genes)){
      top.genes <- top.genes.tmp
    }else{
      top.genes <- rbind(top.genes, top.genes.tmp)
    }

  }
}

covariates <- unique(mir.table[which(mir.table$miRNA %in% mir.15),]$miRNA_ID)

covariates <- lapply(covariates, function(x){
  return(str_replace_all(x, '_', '-'))
  })

covariates <- unlist(covariates)

top.genes.mir15 <- top.genes[which(top.genes$X %in% covariates),]

data <- NULL
comparisons <- list()
total_com <- list()
for(sample in unique(top.genes.mir15$TCGA_ID)){
  print(sample)

  # sample <- 'TCGA_BRCA'
  tmp.data <- my_df[which(my_df$TCGA_ID == sample),]

  if(dim(tmp.data)[1]==0) next

  covariates <- top.genes.mir15[which(top.genes.mir15$TCGA_ID == sample),]$X

  covariates <- lapply(covariates, function(x){
  return(str_replace_all(x, '-', '_'))
  })

  covariates <- unlist(covariates)
  print(unique(tmp.data$Sample_type))
  tmp.data.normal <- tmp.data[which(tmp.data$Sample_type == 'Solid Tissue Normal'),]
  tmp.data.tumor <- tmp.data[which(tmp.data$Sample_type == 'Primary Tumor'),]

  print(dim(tmp.data.normal))
  if(dim(tmp.data.normal)[1] == 0 || dim(tmp.data.tumor)[1] == 0) next

  tmp.normal <- melt(tmp.data.normal[,covariates,drop=FALSE])
  tmp.tumor <- melt(tmp.data.tumor[,covariates,drop=FALSE])

  tmp.normal$Sample_type <- 'N'
  tmp.tumor$Sample_type <- 'T'

  comparisons[[sample]] <-  list()
  # print(comparisons)

  print(tmp.normal[1:2,])
  print(tmp.tumor[1:2,])
  my_tmp.data <- rbind(tmp.normal,tmp.tumor)
  print(class(my_tmp.data))
  my_tmp.data$TCGA_ID <- sample
  my_tmp.data$p.value <- 0.05
  my_tmp.data$Group <- ''
  my_tmp.data.mir15 <- top.genes.mir15[which(top.genes.mir15$TCGA_ID == sample),]
  for(m in unique(my_tmp.data.mir15$X)){
    p.value <- my_tmp.data.mir15[which(my_tmp.data.mir15$X == m),]$FDR
    m <- str_replace_all(m, '-', '_')
    print(p.value)
    my_tmp.data[which(my_tmp.data$variable == m),]$Group <- paste0(m, '_',my_tmp.data[which(my_tmp.data$variable == m),]$Sample_type)
    my_tmp.data[which(my_tmp.data$variable == m),]$p.value <- p.value
    comparisons[[sample]][[m]] <- c(paste0(m, '_T'), paste0(m, '_N'))

    total_com[[paste0(sample,'_',m)]] <- c(paste0(m, '_T'), paste0(m, '_N'))
  }

  print(my_tmp.data[1:2,])
  if(is.null(data)){
    data <- my_tmp.data
  }else{
    print(data[1:2,])
    data <- rbind(data,my_tmp.data)
  }

}

my_plots_data <- c()
widths <- c()
for(sample in unique(data$TCGA_ID)){
  fig_data <- data[which(data$TCGA_ID == sample),]
  my_plots_data <- c(my_plots_data, sample)
  widths <- c(widths, length(unique(fig_data$variable)))
}

myplot <- function(sample){
    fig_data <- data[which(data$TCGA_ID == sample),]
    
    sum.data <- ddply(fig_data,.(variable),summarise,value=max(value)+0.3,p.value=max(p.value))
    annotation <- c()
    y_position <- c()
    xmin <- c()
    xmax <- c()
    print(dim(sum.data))
    for(idx in 1:dim(sum.data)[1]){
      print(sum.data[idx,'value'])
      y_position <- c(y_position, sum.data[idx,'value'])

      xmin = c(xmin, 1+(idx-1)*2)
      xmax = c(xmax, 2+(idx-1)*2)
      p <- sum.data[idx, 'p.value']
      if(p > 0.01){
        annotation <- c(annotation, "*")
      }else if(p > 0.001){
        annotation <- c(annotation, "**")
      }else{
        annotation <- c(annotation, "***")
      }
    }
    print(y_position)
    print(annotation)
    print(xmin)
    ggplot(fig_data, aes(x=Group, y=value, color=variable, fill=variable)) + 
    # geom_boxplot() +
    geom_violin(trim = FALSE) +
    # geom_signif(y_position=y_position, xmin=xmin, xmax=xmax,annotation = annotation) +
    facet_grid(~TCGA_ID)+
    geom_signif(test="wilcox.test",comparisons = comparisons[[sample]], 
                map_signif_level=TRUE) +
    # geom_signif(y_position=y_position, xmin=xmin, xmax=xmax,annotation = annotation) +
     # geom_text(size=2) +
    # geom_dotplot(binaxis='y', stackdir='center', stackratio=1.5, dotsize=0.5, color='black', fill='black') +
    theme(legend.position="none", axis.text.x = element_text(angle = 45, vjust = 0.7, hjust=0.8, size=10)) +
    ylab("Gene Expression Level [log2(fpkm+1)]") +
    xlab("")
}

my_plots <- map(my_plots_data, myplot)

png(paste('output/miR15/TCGA.all',".primary.tumor.vs.normal.miR15Family.v2.boxplot.png", sep=""),height=12000, width=14000, res=800)
ggarrange(plotlist=my_plots, nrow=4, ncol=4)
dev.off()



#Figure AUC/ROC


#Figure coxph multivariate model


#-------------------------miRNA----------

mir.15 <- c('hsa-mir-103a-1','hsa-mir-103a-2','hsa-mir-15a','hsa-mir-15b','hsa-mir-107','hsa-mir-16-1','hsa-mir-195','hsa-mir-497','hsa-mir-503','hsa-mir-424')

covariates <- c(mir.15)

mirna.df.sub <- mirna.df[covariates,]

covariates <- lapply(covariates, function(x){
  return(str_replace_all(x, '-', '_'))
  })
covariates <- unlist(covariates)

rownames(mirna.df.sub) <- covariates

covariates <- mir.table$miRNA_gene

mirna.df.sub <- mirna.df.sub[covariates,]

mirna.df.sub <- 2^mirna.df.sub - 1

mirna.df.sub$Symbol <- mir.table$miRNA_ID
mirna.df.sub <- data.table(mirna.df.sub)
mirna.df.sub <- mirna.df.sub[, lapply(.SD,sum), by=.(Symbol)]
mirna.df.sub <- as.data.frame(mirna.df.sub)
rownames(mirna.df.sub) <- mirna.df.sub$Symbol
mirna.df.sub <- mirna.df.sub[,-1]

mirna.df.sub <- log1p(mirna.df.sub)

mirna.df.sub.t <- as.data.frame(t(mirna.df.sub))

common.samples <- intersect(rownames(mirna.df.sub.t), survival.df$sample)

mirna.df.sub.t <- mirna.df.sub.t[match(common.samples, rownames(mirna.df.sub.t)),]
mirna.df.sub.t <- cbind(mirna.df.sub.t, survival.df[match(common.samples, survival.df$sample),][,-1])

print('Finish combining sample')

common.samples <- intersect(rownames(mirna.df.sub.t), pheno.df$`submitter_id.samples`)

mirna.df.sub.t$Sample_type <- ''

mirna.df.sub.t[match(common.samples, rownames(mirna.df.sub.t)),]$Sample_type <- pheno.df[match(common.samples, pheno.df$`submitter_id.samples`),]$`sample_type.samples`


my_data <- mirna.df.sub.t
names(my_data)[names(my_data) == "OS.time"] <- "OS_time"

covariates <- unique(mir.table$miRNA_ID)

my_df <- my_data[,c(covariates, c('OS','_PATIENT','OS_time','TCGA_ID','Sample_type'))]

my_df <- my_df[which(my_df$Sample_type == 'Primary Tumor'),]


covariates <- unique(mir.table[which(mir.table$miRNA %in% mir.15),]$miRNA_ID)

covariates <- lapply(covariates, function(x){
  return(str_replace_all(x, '-', '_'))
  })

covariates <- unlist(covariates)

result_mirna <- NULL
for(sample in unique(my_data$TCGA_ID)){
  # sample <- 'TCGA_BRCA'
  my_df.sub <- my_df[which(my_df$TCGA_ID == sample),]

  if(dim(my_df.sub)[1]==0) next

  my_means <- colMedians(as.matrix(my_df.sub[,covariates]))
  for(i in c(1:length(covariates))){
    v <- my_means[i]
    high <- my_df.sub[,covariates[i]] >= v
    my_df.sub[,covariates[i]][which(high==TRUE)] <- 'high'
    my_df.sub[,covariates[i]][which(high==FALSE)] <- 'low'
  }

  # print(my_df.sub)
  cv <- NULL
  for (covariate in covariates){
    if (is.null(cv)){
      cv <- covariate
    }else{
      cv <- paste(cv, covariate, sep='+')
    }
    
  }

  univ_formulas <- sapply(c(cv), function(x) as.formula(paste('Surv(OS_time, OS)~', x)))
  model_names <- names(univ_formulas)
    # print(model_names)
  i <- 1
  univ_models <- list()
  for(m in univ_formulas){
    skip_to_next <- FALSE
    p <- tryCatch({
        print(m)    
        fit <- coxph(m, data = my_df.sub)
        univ_models[[model_names[i]]] <- fit
        i <- i+1 
        
      },error = function(e){
        # i <- i+1
        print(e)
        print(sample)
        print('\n')
        skip_to_next <- TRUE
        }, 
      warning = function(cond) {
        message(cond)
        print(sample)
        print('\n')
        # i <- i+1
        skip_to_next <- TRUE
      })
    # print(p)
    if(skip_to_next == TRUE){next}
    
    } 

  if(length(univ_models) == 0) next
  x <- univ_models[[1]]
  fit.coxph <- x
  x <- summary(x)
  print(x)
  p.value<-signif(x$coefficients[,'Pr(>|z|)'], digits=2)

  wald.test<-signif(x$wald["test"], digits=2)
  beta<-signif(x$coef[,1], digits=2);#coeficient beta
  HR <-signif(x$coef[,2], digits=2);#exp(beta)
  total <- x$n
  nevent <- x$nevent
  HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
  HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
  HR.confint <- paste0(HR, " (",HR.confint.lower, "-", HR.confint.upper, ")")
  file_name <- rownames(x$conf.int)
  file_name <- unlist(lapply(file_name, function(x){str_split(x, 'l')[[1]][1]}))
  expr_gene <- NULL
  for(name in file_name){
    expr_gene.tmp <- table(my_df.sub[,name])
    expr_gene.tmp <- as.data.frame(expr_gene.tmp)
    expr_gene.df <- data.frame('high'=c(expr_gene.tmp[1,2]), 'low'=c(expr_gene.tmp[2,2]))
    if(is.null(expr_gene)){
      expr_gene <- expr_gene.df
    }else{
      expr_gene <- rbind(expr_gene, expr_gene.df)
    }
  }

  tmp <- data.frame("Symbol"=file_name,"high"=expr_gene[,1], "low"=expr_gene[,2],
    "Total"=total, "Event"=nevent, "beta"=beta, "HR"=HR,
     "HR_low"=HR.confint.lower, "HR_high"=HR.confint.upper,"HR (95% CI for HR)"=HR.confint, 
     "wald.test"=wald.test, "p.value"=p.value, "TCGA_ID"=sample)

  if(dim(tmp)[1] == 0) next

  tmp$Symbol <- unlist(lapply(tmp$Symbol, function(x){
      return(str_replace_all(x, '_', '-'))
      }))


  write.csv(tmp, paste0(output, '/miR15/TCGA.', sample,'.miR15Family.HR.csv'))

  tmp_sig <- tmp[which(tmp$p.value < 0.05),]

  if(dim(tmp_sig)[1] > 0){
    tmp_sig$logFC <- 0
    tmp_sig$logFC <- top.genes[match(tmp_sig$Symbol, rownames(top.genes)),]$logFC

    write.csv(tmp_sig, paste0(output, '/miR15/TCGA.', sample,'.miR15Family.sig.HR.csv'))    
  }


  if (is.null(result_mirna)){
    result_mirna <- tmp
  }else{
    result_mirna <- rbind(result_mirna, tmp)
  }

  rownames(result_mirna) <- c(1:dim(result_mirna)[1])
}

# rownames(result) <- c(1:dim(result)[1])
# result <- result[,c('miRNA',colnames(result)[-dim(result)[2]])]

write.csv(result_mirna, 'output/miR15/TCGA.miR15Family.multivariates.mature.HR.csv')


#Figure coxph univariate model


#-------------------------miRNA----------

mir.15 <- c('hsa-mir-103a-1','hsa-mir-103a-2','hsa-mir-15a','hsa-mir-15b','hsa-mir-107','hsa-mir-16-1','hsa-mir-195','hsa-mir-497','hsa-mir-503','hsa-mir-424','hsa-mir-646','hsa-mir-6838')


covariates <- c(mir.15)

mirna.df.sub <- mirna.df[covariates,]

covariates <- lapply(covariates, function(x){
  return(str_replace_all(x, '-', '_'))
  })
covariates <- unlist(covariates)

rownames(mirna.df.sub) <- covariates

covariates <- mir.table$miRNA_gene

mirna.df.sub <- mirna.df.sub[covariates,]

mirna.df.sub <- 2^mirna.df.sub - 1

mirna.df.sub$Symbol <- mir.table$miRNA_ID
mirna.df.sub <- data.table(mirna.df.sub)
mirna.df.sub <- mirna.df.sub[, lapply(.SD,sum), by=.(Symbol)]
mirna.df.sub <- as.data.frame(mirna.df.sub)
rownames(mirna.df.sub) <- mirna.df.sub$Symbol
mirna.df.sub <- mirna.df.sub[,-1]

mirna.df.sub <- log1p(mirna.df.sub)

mirna.df.sub.t <- as.data.frame(t(mirna.df.sub))

common.samples <- intersect(rownames(mirna.df.sub.t), survival.df$sample)

mirna.df.sub.t <- mirna.df.sub.t[match(common.samples, rownames(mirna.df.sub.t)),]
mirna.df.sub.t <- cbind(mirna.df.sub.t, survival.df[match(common.samples, survival.df$sample),][,-1])

print('Finish combining sample')

common.samples <- intersect(rownames(mirna.df.sub.t), pheno.df$`submitter_id.samples`)

mirna.df.sub.t$Sample_type <- ''

mirna.df.sub.t[match(common.samples, rownames(mirna.df.sub.t)),]$Sample_type <- pheno.df[match(common.samples, pheno.df$`submitter_id.samples`),]$`sample_type.samples`


my_data <- mirna.df.sub.t
names(my_data)[names(my_data) == "OS.time"] <- "OS_time"

covariates <- unique(mir.table$miRNA_ID)

my_df <- my_data[,c(covariates, c('OS','_PATIENT','OS_time','TCGA_ID','Sample_type'))]

my_df <- my_df[which(my_df$Sample_type == 'Primary Tumor'),]


covariates <- unique(mir.table[which(mir.table$miRNA %in% mir.15),]$miRNA_ID)

covariates <- lapply(covariates, function(x){
  return(str_replace_all(x, '-', '_'))
  })

covariates <- unlist(covariates)

result_mirna_all <- NULL
for(sample in unique(my_data$TCGA_ID)){
  # sample <- 'TCGA_ACC'
  my_df.sub <- my_df[which(my_df$TCGA_ID == sample),]

  if(dim(my_df.sub)[1]==0) next

  my_means <- colMedians(as.matrix(my_df.sub[,covariates]))
  for(i in c(1:length(covariates))){
    v <- my_means[i]
    high <- my_df.sub[,covariates[i]] >= v
    my_df.sub[,covariates[i]][which(high==TRUE)] <- 'high'
    my_df.sub[,covariates[i]][which(high==FALSE)] <- 'low'
  }

  univ_formulas <- sapply(covariates, function(x) as.formula(paste('Surv(OS_time, OS)~', x)))
  model_names <- names(univ_formulas)
    # print(model_names)
  i <- 1
  univ_models <- list()
  for(m in univ_formulas){
    skip_to_next <- FALSE
      p <- tryCatch({
          # print(m)    
          fit <- coxph(m, data = my_df.sub)
          univ_models[[model_names[i]]] <- fit
          i <- i+1 
          
        },error = function(e){
          # i <- i+1
          print(sample)
          print('\n')
          skip_to_next <- TRUE
          }, 
        warning = function(cond) {
          message(cond)
          print(sample)
          print('\n')
          # i <- i+1
          skip_to_next <- TRUE
        })
      # print(p)
      if(skip_to_next == TRUE){next}
    
    } 

  j <- 1
  km_models <- list()
  for(m in univ_formulas){
      skip_to_next <- FALSE
      p <- tryCatch({
          print(m)    
          print(model_names[j])
          fit <- surv_fit(m, data = my_df.sub)
          km_models[[model_names[j]]] <- fit
          j <- j+1 
          
        },error = function(e){
          # j <- j+1
          print(sample)
          print('\n')
          skip_to_next <- TRUE
          }, 
        warning = function(cond) {
          message(cond)
          print(sample)
          print('\n')
          # j <- j+1
          skip_to_next <- TRUE
        })
      # print(p)
      if(skip_to_next == TRUE){next}
    } 


    if(length(univ_models) == 0) next
    univ_results <- lapply(univ_models,
                           function(x){ 
                            fit.coxph <- x
                            x <- summary(x)
                            print(x)
                            # p.value<-signif(x$wald["pvalue"], digits=2)
                            p.value<-signif(x$coefficients[,'Pr(>|z|)'], digits=2)
                            wald.test<-signif(x$wald["test"], digits=2)
                            beta<-signif(x$coef[1], digits=2);#coeficient beta
                            HR <-signif(x$coef[2], digits=2);#exp(beta)
                            total <- x$n
                            nevent <- x$nevent
                            HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
                            HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
                            HR.confint <- paste0(HR, " (",HR.confint.lower, "-", HR.confint.upper, ")")
                            file_name <- rownames(x$conf.int)[1]
                            file_name <- str_split(file_name, 'l')[[1]][1]
                            expr_gene <- table(my_df.sub[,file_name])
                            expr_gene <- as.data.frame(expr_gene)
    tmp<-c(file_name, expr_gene[1,2],expr_gene[2,2],total, nevent, beta, HR, HR.confint.lower, HR.confint.upper, HR.confint, wald.test, p.value, sample)
                          # 

    names(tmp)<-c("Symbol",levels(expr_gene[1,1])[expr_gene[1,1]],levels(expr_gene[2,1])[expr_gene[2,1]],"Total", "Event", "beta", "HR", "HR_low", "HR_high","HR (95% CI for HR)", "wald.test", 
                  "p.value", "TCGA_ID")
    return(tmp)})

    univ_km <- lapply(km_models, function(x){
      print(x)
      my_fit <- x
      # print(fit)
      fit_pval <- surv_pvalue(my_fit)$pval
      x <- summary(x)
      print(x$table)
      file_name <- rownames(x$table)[1]
      print(file_name)
      if(!is.null(file_name) && fit_pval < 0.05){
        file_name <- str_split(file_name, '=')[[1]][1]
        sample_s <- str_split(sample, '_')[[1]][[2]]
        png(paste(output, '/miR15/TCGA.', sample,'.', file_name,'.KM.v2.png', sep=""), height=6000, width=6000, res=1000)
        print(my_fit)
        p <- ggsurvplot(
              my_fit,                     
              data = my_df.sub,            
              risk.table = TRUE,       
              risk.table.col = file_name, 
              pval = TRUE,         
              conf.int = TRUE,     
              palette = "Dark2",
              xlab = "Time in days", 
              ggtheme = theme_bw(), 

              conf.int.style = "step",
              legend.title = paste0(sample_s, ': \n', file_name),
              legend.labs = c("High", "Low"))
        print(p)
        dev.off()
      }
      

      return(file_name)
    })


    tmp <- t(as.data.frame(univ_results, check.names = FALSE))
    print(tmp)
    tmp <- as.data.frame(tmp)

    if(dim(tmp)[1] == 0) next

    tmp$Symbol <- unlist(lapply(tmp$Symbol, function(x){
        return(str_replace_all(x, '_', '-'))
        }))


    if (is.null(result_mirna_all)){
      result_mirna_all <- tmp
    }else{
      result_mirna_all <- rbind(result_mirna_all, tmp)
    }

    rownames(result_mirna_all) <- c(1:dim(result_mirna_all)[1])
}

write.csv(result_mirna_all, 'output/miR15/TCGA.miR15.Family.univariates.mature.HR.csv')



#Figure
result_mirna_all <- read.csv('output/miR15/TCGA.miR15.Family.univariates.mature.HR.csv')

my_result_mirna <- NULL
for(sample in unique(result_mirna$TCGA_ID)){
  # sample <- 'TCGA_ACC'
  tmp1 <- result_mirna[which(result_mirna$TCGA_ID == sample),]
  tmp2 <- result_mirna_all[which(result_mirna_all$TCGA_ID == sample),]
  tmp2 <- tmp2[match(tmp1$Symbol, tmp2$Symbol),]
  tmp2$Symbol <- tmp1$Symbol
  tmp2$TCGA_ID <- tmp1$TCGA_ID

  cols <- colnames(tmp2)

  cols <- lapply(cols, function(x){
    paste0(x,'_U')
    })
  cols <- unlist(cols)
  colnames(tmp2) <- cols
  tmp3 <- cbind(tmp1, tmp2[,-1])
  if(is.null(my_result_mirna)){
    my_result_mirna <- tmp3
  }else{
    my_result_mirna <- rbind(my_result_mirna, tmp3)
  }
}

for(sample in unique(my_result_mirna$TCGA_ID)){
  # sample <- 'TCGA_BRCA'
  my_plot_data <- my_result_mirna[which(my_result_mirna$TCGA_ID == sample),]

#  [1] "Symbol"               "high"                 "low"                  "Total"               
#  [5] "Event"                "beta"                 "HR"                   "HR_low"              
#  [9] "HR_high"              "HR..95..CI.for.HR."   "wald.test"            "p.value"             
# [13] "TCGA_ID"              "Symbol_U"             "high_U"               "low_U"               
# [17] "Total_U"              "Event_U"              "beta_U"               "HR_U"                
# [21] "HR_low_U"             "HR_high_U"            "HR..95..CI.for.HR._U" "wald.test_U"         
# [25] "p.value_U"            "TCGA_ID_U"
  my_plot_data$Index <- as.numeric(c(1:dim(my_plot_data)[1]))
  my_plot_data$HR <- as.numeric(my_plot_data$HR)
  my_plot_data$HR_low <- as.numeric(my_plot_data$HR_low)
  my_plot_data$HR_high <- as.numeric(my_plot_data$HR_high)


  ############################################
  ### CUSTOMIZE APPEARANCE WITH THESE     ####
  ############################################
  blankRows<-2    # blank rows under boxplot
  titleSize<-4
  dataSize<-4
  boxColor<-"pink"
  ############################################
  ############################################

  ## BASIC THEMES (SO TO PLOT BLANK GRID)
  theme_grid <- theme(
    axis.line = element_blank(), 
    axis.text.x = element_blank(), 
    axis.text.y = element_blank(),
    axis.ticks = element_blank(), 
    axis.title.x = element_blank(), 
    axis.title.y = element_blank(), 
    axis.ticks.length = unit(0.0001, "mm"),
    axis.ticks.margin = unit(c(0,0,0,0), "lines"), 
    legend.position = "none", 
    panel.background = element_rect(fill = "transparent"), 
    panel.border = element_blank(), 
    panel.grid.major = element_line(colour="grey"), 
    panel.grid.minor = element_line(colour="grey"), 
    panel.margin = unit(c(-0.1,-0.1,-0.1,-0.1), "mm"), 
    plot.margin = unit(c(5,0,5,0.01), "mm")
  )

  theme_bare <- theme_grid +
    theme(
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank()
    )
# 1,2,3,20,21,22,23,25,7,8,9,10,12,13
  group_low <- my_plot_data[,c(1,3,7,8,9,10,12)]
  group_high <- my_plot_data[,c(1,2,7,8,9,10,12)]
  colnames(group_low) <- c('miRNA', 'NoP', 'HR', 'HR_low', 'HR_high', 'HR (95% CI for HR)', 'p.value')
  colnames(group_high) <- c('miRNA', 'NoP', 'HR', 'HR_low', 'HR_high', 'HR (95% CI for HR)', 'p.value')
  group_low$Group <- 'Low'
  group_high$Group <- 'High'
  group_high[,c(3,4,5)] <- 1
  group_high[,c(6,7)] <- NA
  group_data <- rbind(group_low, group_high)
  group_data <- group_data[order(group_data$miRNA),]
  rownames(group_data) <- c(1:dim(group_data)[1])
  group_data$ID <- c(1:dim(group_data)[1])

  if(dim(group_data[which(group_data$p.value < 0.001),])[1] > 0){
    group_data[which(group_data$p.value < 0.001),]$`HR (95% CI for HR)` <- paste0(group_data[which(group_data$p.value < 0.001),]$`HR (95% CI for HR)`, ' ***')

  }
  if(dim(group_data[which(group_data$p.value > 0.001 & group_data$p.value < 0.01),])[1] > 0){
    group_data[which(group_data$p.value > 0.001 & group_data$p.value < 0.01),]$`HR (95% CI for HR)` <- paste0(group_data[which(group_data$p.value > 0.001 & group_data$p.value < 0.01),]$`HR (95% CI for HR)`, ' **')
  }
  if(dim(group_data[which(group_data$p.value > 0.01 & group_data$p.value < 0.05),])[1] > 0){
    group_data[which(group_data$p.value > 0.01 & group_data$p.value < 0.05),]$`HR (95% CI for HR)` <- paste0(group_data[which(group_data$p.value > 0.01 & group_data$p.value < 0.05),]$`HR (95% CI for HR)`, ' *')
  }

  group_u_low <- my_plot_data[,c(1,3,20,21,22,23,25)]
  group_u_high <- my_plot_data[,c(1,2,20,21,22,23,25)]
  colnames(group_u_low) <- c('miRNA', 'NoP', 'HR', 'HR_low', 'HR_high', 'HR (95% CI for HR)', 'p.value')
  colnames(group_u_high) <- c('miRNA', 'NoP', 'HR', 'HR_low', 'HR_high', 'HR (95% CI for HR)', 'p.value')
  group_u_low$Group <- 'Low'
  group_u_high$Group <- 'High'
  group_u_high[,c(3,4,5)] <- 1
  group_u_high[,c(6,7)] <- NA
  group_u_data <- rbind(group_u_low, group_u_high)
  group_u_data <- group_u_data[order(group_u_data$miRNA),]
  rownames(group_u_data) <- c(1:dim(group_u_data)[1])
  group_u_data$ID <- c(1:dim(group_u_data)[1])

  if(dim(group_u_data[which(group_u_data$p.value < 0.001),])[1] > 0){
    group_u_data[which(group_u_data$p.value < 0.001),]$`HR (95% CI for HR)` <- paste0(group_u_data[which(group_u_data$p.value < 0.001),]$`HR (95% CI for HR)`, ' ***')

  }
  if(dim(group_u_data[which(group_u_data$p.value > 0.001 & group_u_data$p.value < 0.01),])[1] > 0){
    group_u_data[which(group_u_data$p.value > 0.001 & group_u_data$p.value < 0.01),]$`HR (95% CI for HR)` <- paste0(group_u_data[which(group_u_data$p.value > 0.001 & group_u_data$p.value < 0.01),]$`HR (95% CI for HR)`, ' **')
  }
  if(dim(group_u_data[which(group_u_data$p.value > 0.01 & group_u_data$p.value < 0.05),])[1] > 0){
    group_u_data[which(group_u_data$p.value > 0.01 & group_u_data$p.value < 0.05),]$`HR (95% CI for HR)` <- paste0(group_u_data[which(group_u_data$p.value > 0.01 & group_u_data$p.value < 0.05),]$`HR (95% CI for HR)`, ' *')
  }
  

  hazard_data<-expand.grid(ID=1:nrow(group_data),HR=1)
  hazard_data$HR <- group_data$HR
  hazard_data<-rbind(hazard_data,ddply(group_data,.(miRNA),summarise,ID=max(ID)+0.1,HR=NA)[,2:3])

  hazard_data<-rbind(hazard_data,data.frame(ID=c(0,-1:(-2-blankRows),max(group_data$ID)+1,max(group_data$ID)+2),HR=NA))

  hazard_data$HR_low <- NA
  hazard_data$HR_high <- NA

  hazard_data[1:dim(group_data)[1],]$HR_low <- group_data$HR_low
  hazard_data[1:dim(group_data)[1],]$HR_high <- group_data$HR_high

  hazard_u_data<-expand.grid(ID=1:nrow(group_u_data),HR=1)
  hazard_u_data$HR <- group_u_data$HR
  hazard_u_data<-rbind(hazard_u_data,ddply(group_u_data,.(miRNA),summarise,ID=max(ID)+0.1,HR=NA)[,2:3])

  hazard_u_data<-rbind(hazard_u_data,data.frame(ID=c(0,-1:(-2-blankRows),max(group_u_data$ID)+1,max(group_u_data$ID)+2),HR=NA))

  hazard_u_data$HR_low <- NA
  hazard_u_data$HR_high <- NA

  hazard_u_data[1:dim(group_u_data)[1],]$HR_low <- group_u_data$HR_low
  hazard_u_data[1:dim(group_u_data)[1],]$HR_high <- group_u_data$HR_high


  hr_labels <- group_data[,c('ID', 'HR (95% CI for HR)')]
  colnames(hr_labels) <- c("ID", "lab")

  upper <- max(hazard_data[!is.na(hazard_data$HR_high),]$HR_high)
  lower <- min(hazard_data[!is.na(hazard_data$HR_low),]$HR_low)
  lower <- min(0.4, lower)
  upper <- max(2.8, upper)

  upper_int <- round(((upper - 1)/2),2)
  lower_int <- round(((1 - lower)/2),2)

  scale_data <- data.frame(ID=0,HR=c((1-2*lower_int-0.05),(1-lower_int), 1, (1 + upper_int), (1+2*upper_int+0.2)))


  hr_u_labels <- group_u_data[,c('ID', 'HR (95% CI for HR)')]
  colnames(hr_u_labels) <- c("ID", "lab")

  upper_u <- max(hazard_u_data[!is.na(hazard_u_data$HR_high),]$HR_high)
  lower_u <- min(hazard_u_data[!is.na(hazard_u_data$HR_low),]$HR_low)
  lower_u <- min(0.4, lower_u)
  upper_u <- max(2.8, upper_u)

  upper_u_int <- round(((upper_u - 1)/2),2)
  lower_u_int <- round(((1 - lower_u)/2),2)

  scale_u_data <- data.frame(ID=0,HR=c((1-2*lower_u_int-0.05),(1-lower_u_int), 1, (1 + upper_u_int), (1+2*upper_u_int+0.2)))

  group_mirna<-ddply(group_data,.(miRNA),summarise,y=max(ID)+0.1)


  hl_rows <- data.frame(ID=(1:floor(length(unique(hazard_data$ID[which(hazard_data$ID>0)]))/2))*2,col="lightgrey")
  hl_rows$ID <- hl_rows$ID+blankRows+1

  hl_rect <- function(col="white",alpha=0.5){
    rectGrob(   x = 0, y = 0, width = 1, height = 1, just = c("left","bottom"), gp=gpar(alpha=alpha, fill=col))
  }

  ## DATA FOR TEXT LABELS

  md_labels <- data.frame(x=c(rep(length(unique(hazard_data$ID))-0.2,times=1)),
                        y=c(1),
                        lab=c("Hazard Ratio"),drop=FALSE)

  rt_labels <- data.frame(x=c(rep(length(unique(hazard_data$ID))-0.2,times=2)),
                        y=c(1,4),
                        lab=c("Hazard Ratio\n(95% CI)","P Value"))

  lf_labels <- data.frame(x=c(rep(length(unique(hazard_data$ID))-0.2,times=2)),
                       y=c(0.5,4),
                       lab=c("MicroRNA","No. of\nPatients"))

  legend_labels <- data.frame(x=c(rep(1,times=2)),
                       y=c((1-lower_int),(1+upper_int)),
                       lab=c("Low","High"))
  variates_labels <- data.frame(x=c(length(unique(hazard_data$ID))-2),
                       y=c(4),
                       lab=c("Multivariates Model"))
  u_variates_labels <- data.frame(x=c(length(unique(hazard_data$ID))-2),
                       y=c(4),
                       lab=c("Univariate Model"))

  haz <- ggplot(hazard_data,aes(factor(ID),HR))+ labs(x=NULL, y=NULL) 

  ## MIDDLE PANEL WITH LOG SCALE
  middle_panel <- haz +
    apply(hl_rows,1,function(x)annotation_custom(hl_rect(x["col"],alpha=0.4),as.numeric(x["ID"])-0.5,as.numeric(x["ID"])+0.5,-20,20)) +
    # geom_segment(aes(x = 2, y = 1, xend = 1.5, yend = 1)) + 
    geom_hline(aes(yintercept=1),linetype=2, size=0.5)+
    geom_point() + 
    geom_errorbar(data=hazard_data, aes(ymin=HR_low, ymax=HR_high), width=.3)+
    # geom_boxplot(fill=boxColor,size=0.5, alpha=0.8)+ 
    scale_y_log10() + 
    # scale_y_continuous() +
    coord_flip() +
    geom_text(data=scale_data,aes(3,HR,label=HR),vjust=0.5, size=dataSize) +
    geom_text(data=md_labels,aes(x,y,label=lab, fontface="bold"), vjust=0.5, size=titleSize) +
    # geom_text(data=hr_labels,aes(factor(ID),0.4,label=lab),vjust=0.5, hjust=1, size=dataSize) +
    # geom_text(data=group_p,aes(factor(y),11,label=P, fontface="bold"),vjust=0.5, hjust=1, size=dataSize) +
    # geom_text(data=group_data,aes(factor(ID),5,label=p.value),vjust=0.5, hjust=1, size=dataSize) +
    # geom_text(data=legend_labels,aes(x,y,label=lab, fontface="bold"),hjust=0.5, vjust=1, size=titleSize) +
    geom_text(data=legend_labels,aes(x,y,label=lab, fontface="bold"),hjust=0.5, vjust=1, size=titleSize) +
    geom_point(data=scale_data,aes(2.5,HR),shape=3,size=3) + 
    geom_point(aes(2,6),shape=3,alpha=0,vjust=0) + 
    # geom_errorbar(width=.08)
    geom_segment(aes(x = 2.5, y = 0, xend = 2.5, yend = max(8,max(hazard_data[!is.na(hazard_data$HR_high),]$HR_high)))) + 
    geom_segment(aes(x = 2, y = 1, xend = 2, yend = max(hazard_data[!is.na(hazard_data$HR_high),]$HR_high)),arrow=arrow(),linetype=1,size=0.3) + 
    geom_segment(aes(x = 2, y = 1, xend = 2, yend = min(hazard_data[!is.na(hazard_data$HR_low),]$HR_low)),arrow=arrow(),linetype=1,size=0.3) + 
    theme_bare


  ## RIGHT PANEL WITH LOG SCALE
  right_panel <- haz + 
    apply(hl_rows,1,function(x)annotation_custom(hl_rect(x["col"],alpha=0.4),as.numeric(x["ID"])-0.5,as.numeric(x["ID"])+0.5,-20,20)) +
    # scale_y_log10() +
    coord_flip(ylim=c(0,5.5)) +
    geom_point(aes(x=factor(ID),y=1),shape=3,alpha=0,vjust=0) + 
    geom_text(data=variates_labels,aes(x,y,label=lab, fontface="bold"),vjust=0.5, hjust=1, size=titleSize) +
    geom_text(data=hr_labels,aes(factor(ID),3,label=lab),vjust=0.5, hjust=1, size=dataSize) +
    geom_text(data=rt_labels,aes(x,y,label=lab, fontface="bold"), vjust=0.5, size=titleSize) +
    # geom_text(data=group_p,aes(factor(y),11,label=P, fontface="bold"),vjust=0.5, hjust=1, size=dataSize) +
    geom_text(data=group_data,aes(factor(ID),5,label=p.value),vjust=0.5, hjust=1, size=dataSize) +
    geom_segment(aes(x = 2.5, y = 0, xend = 2.5, yend = 5.5)) + 
    # geom_segment(aes(x = 2, y = 1, xend = 2, yend = max(hazard_data[!is.na(hazard_data$HR_high),]$HR_high)),arrow=arrow(),linetype=1,size=0.3) + 
    # geom_segment(aes(x = 2, y = 1, xend = 2, yend = min(hazard_data[!is.na(hazard_data$HR_low),]$HR_low)),arrow=arrow(),linetype=1,size=0.3) + 
    theme_bare

  ## LEFT PANEL WITH NORMAL SCALE
  left_panel<-haz + 
    apply(hl_rows,1,function(x)annotation_custom(hl_rect(x["col"],alpha=0.4),as.numeric(x["ID"])-0.5,as.numeric(x["ID"])+0.5,-20,20)) +
    coord_flip(ylim=c(0,5.5)) +
    geom_point(aes(x=factor(ID),y=1),shape=3,alpha=0,vjust=0) + 
    geom_text(data=group_mirna,aes(factor(y),0.5,label=miRNA, fontface="bold"),vjust=0.5, hjust=0, size=dataSize) +
    geom_text(data=group_data,aes(factor(ID),1,label=Group),vjust=0.5, hjust=0, size=dataSize) +
    geom_text(data=group_data,aes(factor(ID),5,label=NoP),vjust=0.5, hjust=1, size=dataSize) +
    geom_text(data=lf_labels,aes(x,y,label=lab, fontface="bold"), vjust=0.5, hjust=0, size=4, size=titleSize) +
    geom_segment(aes(x = 2.5, y = 0, xend = 2.5, yend = 5.5)) + 
    theme_bare

  u_haz <- ggplot(hazard_u_data,aes(factor(ID),HR))+ labs(x=NULL, y=NULL) 

  ## MIDDLE PANEL WITH LOG SCALE
  middle_u_panel <- u_haz +
    apply(hl_rows,1,function(x)annotation_custom(hl_rect(x["col"],alpha=0.4),as.numeric(x["ID"])-0.5,as.numeric(x["ID"])+0.5,-20,20)) +
    # geom_segment(aes(x = 2, y = 1, xend = 1.5, yend = 1)) + 
    geom_hline(aes(yintercept=1),linetype=2, size=0.5)+
    geom_point() + 
    geom_errorbar(data=hazard_u_data, aes(ymin=HR_low, ymax=HR_high), width=.3)+
    # geom_boxplot(fill=boxColor,size=0.5, alpha=0.8)+ 
    scale_y_log10() + 
    # scale_y_continuous() +
    coord_flip() +
    geom_text(data=scale_u_data,aes(3,HR,label=HR),vjust=0.5, size=dataSize) +
    geom_text(data=md_labels,aes(x,y,label=lab, fontface="bold"), vjust=0.5, size=titleSize) +
    # geom_text(data=hr_labels,aes(factor(ID),0.4,label=lab),vjust=0.5, hjust=1, size=dataSize) +
    # geom_text(data=group_p,aes(factor(y),11,label=P, fontface="bold"),vjust=0.5, hjust=1, size=dataSize) +
    # geom_text(data=group_u_data,aes(factor(ID),5,label=p.value),vjust=0.5, hjust=1, size=dataSize) +
    geom_text(data=legend_labels,aes(x,y,label=lab, fontface="bold"),hjust=0.5, vjust=1, size=titleSize) +
    geom_point(data=scale_u_data,aes(2.5,HR),shape=3,size=3) + 
    geom_point(aes(2,6),shape=3,alpha=0,vjust=0) + 
    # geom_errorbar(width=.08)
    geom_segment(aes(x = 2.5, y = 0, xend = 2.5, yend = max(8,max(hazard_u_data[!is.na(hazard_u_data$HR_high),]$HR_high)))) + 
    geom_segment(aes(x = 2, y = 1, xend = 2, yend = max(hazard_u_data[!is.na(hazard_u_data$HR_high),]$HR_high)),arrow=arrow(),linetype=1,size=0.3) + 
    geom_segment(aes(x = 2, y = 1, xend = 2, yend = min(hazard_u_data[!is.na(hazard_u_data$HR_low),]$HR_low)),arrow=arrow(),linetype=1,size=0.3) + 
    theme_bare


  ## RIGHT PANEL WITH LOG SCALE
  right_u_panel <- u_haz + 
    apply(hl_rows,1,function(x)annotation_custom(hl_rect(x["col"],alpha=0.4),as.numeric(x["ID"])-0.5,as.numeric(x["ID"])+0.5,-20,20)) +
    # scale_y_log10() +
    coord_flip(ylim=c(0,5.5)) +
    geom_point(aes(x=factor(ID),y=1),shape=3,alpha=0,vjust=0) + 
    geom_text(data=hr_u_labels,aes(factor(ID),3,label=lab),vjust=0.5, hjust=1, size=dataSize) +
    geom_text(data=rt_labels,aes(x,y,label=lab, fontface="bold"), vjust=0.5, size=titleSize) +
    geom_text(data=u_variates_labels,aes(x,y,label=lab, fontface="bold"),vjust=0.5, hjust=1, size=titleSize) +
    # geom_text(data=group_p,aes(factor(y),11,label=P, fontface="bold"),vjust=0.5, hjust=1, size=dataSize) +
    geom_text(data=group_u_data,aes(factor(ID),5,label=p.value),vjust=0.5, hjust=1, size=dataSize) +
    geom_segment(aes(x = 2.5, y = 0, xend = 2.5, yend = 5.5)) + 
    # geom_segment(aes(x = 2, y = 1, xend = 2, yend = max(hazard_u_data[!is.na(hazard_u_data$HR_high),]$HR_high)),arrow=arrow(),linetype=1,size=0.3) + 
    # geom_segment(aes(x = 2, y = 1, xend = 2, yend = min(hazard_u_data[!is.na(hazard_u_data$HR_low),]$HR_low)),arrow=arrow(),linetype=1,size=0.3) + 
    theme_bare

  ## PLOT THEM BOTH IN A GRID SO THEY MATCH UP
  height <- 4000 * (floor((dim(my_plot_data)[1] - 1) / 5) + 1)
  resolution <- 400 + 100 * (floor((dim(my_plot_data)[1] - 1) / 5) + 1)
  print(hazard_data)
  png(paste(output, '/miR15/TCGA.', sample,'.miR15.Family.all.tmp.HR.png', sep=""), height=height, width=10000, res=resolution)
  p <- ggarrange(left_panel,middle_u_panel, right_u_panel,middle_panel,right_panel,  widths=c(3,4,3,4,3), ncol=5, nrow=1)
  print(p)
  dev.off() 
}


#Figure-----------------------miRNA cox-----------------
result_mirna_sig <- result_mirna_all[which(result_mirna_all$p.value < 0.05),]
result_mirna_sig <- result_mirna_sig[,-1]
for(sample in unique(result_mirna_sig$TCGA_ID)){
  my_plot_data <- result_mirna_sig[which(result_mirna_sig$TCGA_ID == sample),]

  my_plot_data$Index <- as.numeric(c(1:dim(my_plot_data)[1]))
  my_plot_data$HR <- as.numeric(my_plot_data$HR)
  my_plot_data$HR_low <- as.numeric(my_plot_data$HR_low)
  my_plot_data$HR_high <- as.numeric(my_plot_data$HR_high)


  ############################################
  ### CUSTOMIZE APPEARANCE WITH THESE     ####
  ############################################
  blankRows<-2    # blank rows under boxplot
  titleSize<-4
  dataSize<-4
  boxColor<-"pink"
  ############################################
  ############################################

  ## BASIC THEMES (SO TO PLOT BLANK GRID)
  theme_grid <- theme(
    axis.line = element_blank(), 
    axis.text.x = element_blank(), 
    axis.text.y = element_blank(),
    axis.ticks = element_blank(), 
    axis.title.x = element_blank(), 
    axis.title.y = element_blank(), 
    axis.ticks.length = unit(0.0001, "mm"),
    axis.ticks.margin = unit(c(0,0,0,0), "lines"), 
    legend.position = "none", 
    panel.background = element_rect(fill = "transparent"), 
    panel.border = element_blank(), 
    panel.grid.major = element_line(colour="grey"), 
    panel.grid.minor = element_line(colour="grey"), 
    panel.margin = unit(c(-0.1,-0.1,-0.1,-0.1), "mm"), 
    plot.margin = unit(c(5,0,5,0.01), "mm")
  )

  theme_bare <- theme_grid +
    theme(
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank()
    )
  group_low <- my_plot_data[,c(1,3,7,8,9,10,12)]
  group_high <- my_plot_data[,c(1,2,7,8,9,10,12)]
  colnames(group_low) <- c('miRNA', 'NoP', 'HR', 'HR_low', 'HR_high', 'HR (95% CI for HR)', 'p.value')
  colnames(group_high) <- c('miRNA', 'NoP', 'HR', 'HR_low', 'HR_high', 'HR (95% CI for HR)', 'p.value')
  group_low$Group <- 'Low'
  group_high$Group <- 'High'
  group_high[,c(3,4,5)] <- 1
  group_high[,c(6,7)] <- NA
  group_data <- rbind(group_low, group_high)
  group_data <- group_data[order(group_data$miRNA),]
  rownames(group_data) <- c(1:dim(group_data)[1])
  group_data$ID <- c(1:dim(group_data)[1])
  hazard_data<-expand.grid(ID=1:nrow(group_data),HR=1)
  hazard_data$HR <- group_data$HR
  hazard_data<-rbind(hazard_data,ddply(group_data,.(miRNA),summarise,ID=max(ID)+0.1,HR=NA)[,2:3])

  hazard_data<-rbind(hazard_data,data.frame(ID=c(0,-1:(-2-blankRows),max(group_data$ID)+1,max(group_data$ID)+2),HR=NA))

  hazard_data$HR_low <- NA
  hazard_data$HR_high <- NA

  hazard_data[1:dim(group_data)[1],]$HR_low <- group_data$HR_low
  hazard_data[1:dim(group_data)[1],]$HR_high <- group_data$HR_high

  hr_labels <- group_data[,c('ID', 'HR (95% CI for HR)')]
  colnames(hr_labels) <- c("ID", "lab")

  upper <- max(hazard_data[!is.na(hazard_data$HR_high),]$HR_high)
  lower <- min(hazard_data[!is.na(hazard_data$HR_low),]$HR_low)
  lower <- min(0.4, lower)
  upper <- max(2.8, upper)

  upper_int <- round(((upper - 1)/2),2)
  lower_int <- round(((1 - lower)/2),2)

  scale_data <- data.frame(ID=0,HR=c((1-2*lower_int-0.05),(1-lower_int), 1, (1 + upper_int), (1+2*upper_int+0.2)))

  group_mirna<-ddply(group_data,.(miRNA),summarise,y=max(ID)+0.1)


  hl_rows <- data.frame(ID=(1:floor(length(unique(hazard_data$ID[which(hazard_data$ID>0)]))/2))*2,col="lightgrey")
  hl_rows$ID <- hl_rows$ID+blankRows+1

  hl_rect <- function(col="white",alpha=0.5){
    rectGrob(   x = 0, y = 0, width = 1, height = 1, just = c("left","bottom"), gp=gpar(alpha=alpha, fill=col))
  }

  ## DATA FOR TEXT LABELS
  md_labels <- data.frame(x=c(rep(length(unique(hazard_data$ID))-0.2,times=1)),
                        y=c(1),
                        lab=c("Hazard Ratio"),drop=FALSE)

  rt_labels <- data.frame(x=c(rep(length(unique(hazard_data$ID))-0.2,times=2)),
                        y=c(1,4),
                        lab=c("Hazard Ratio\n(95% CI)","P Value"))

  lf_labels <- data.frame(x=c(rep(length(unique(hazard_data$ID))-0.2,times=2)),
                       y=c(0.5,4),
                       lab=c("MicroRNA","No. of\nPatients"))

  legend_labels <- data.frame(x=c(rep(1,times=2)),
                       y=c((1-lower_int),(1+upper_int)),
                       lab=c("Low Better","High Better"))

  haz <- ggplot(hazard_data,aes(factor(ID),HR))+ labs(x=NULL, y=NULL) 

  ## MIDDLE PANEL WITH LOG SCALE
  middle_panel <- haz +
    apply(hl_rows,1,function(x)annotation_custom(hl_rect(x["col"],alpha=0.4),as.numeric(x["ID"])-0.5,as.numeric(x["ID"])+0.5,-20,20)) +
    # geom_segment(aes(x = 2, y = 1, xend = 1.5, yend = 1)) + 
    geom_hline(aes(yintercept=1),linetype=2, size=0.5)+
    geom_point() + 
    geom_errorbar(data=hazard_data, aes(ymin=HR_low, ymax=HR_high), width=.3)+
    # geom_boxplot(fill=boxColor,size=0.5, alpha=0.8)+ 
    scale_y_log10() + 
    # scale_y_continuous() +
    coord_flip() +
    geom_text(data=scale_data,aes(3,HR,label=HR),vjust=0.5, size=dataSize) +
    geom_text(data=md_labels,aes(x,y,label=lab, fontface="bold"), vjust=0.5, size=titleSize) +
    # geom_text(data=hr_labels,aes(factor(ID),0.4,label=lab),vjust=0.5, hjust=1, size=dataSize) +
    # geom_text(data=group_p,aes(factor(y),11,label=P, fontface="bold"),vjust=0.5, hjust=1, size=dataSize) +
    # geom_text(data=group_data,aes(factor(ID),5,label=p.value),vjust=0.5, hjust=1, size=dataSize) +
    geom_text(data=legend_labels,aes(x,y,label=lab, fontface="bold"),hjust=0.5, vjust=1, size=titleSize) +
    geom_point(data=scale_data,aes(2.5,HR),shape=3,size=3) + 
    geom_point(aes(2,6),shape=3,alpha=0,vjust=0) + 
    # geom_errorbar(width=.08)
    geom_segment(aes(x = 2.5, y = 0, xend = 2.5, yend = max(8,max(hazard_data[!is.na(hazard_data$HR_high),]$HR_high)))) + 
    geom_segment(aes(x = 2, y = 1, xend = 2, yend = max(hazard_data[!is.na(hazard_data$HR_high),]$HR_high)),arrow=arrow(),linetype=1,size=0.3) + 
    geom_segment(aes(x = 2, y = 1, xend = 2, yend = min(hazard_data[!is.na(hazard_data$HR_low),]$HR_low)),arrow=arrow(),linetype=1,size=0.3) + 
    theme_bare


  ## RIGHT PANEL WITH LOG SCALE
  rightPanel <- haz + 
    apply(hl_rows,1,function(x)annotation_custom(hl_rect(x["col"],alpha=0.4),as.numeric(x["ID"])-0.5,as.numeric(x["ID"])+0.5,-20,20)) +
    # scale_y_log10() +
    coord_flip(ylim=c(0,5.5)) +
    geom_point(aes(x=factor(ID),y=1),shape=3,alpha=0,vjust=0) + 
    geom_text(data=hr_labels,aes(factor(ID),3,label=lab),vjust=0.5, hjust=1, size=dataSize) +
    geom_text(data=rt_labels,aes(x,y,label=lab, fontface="bold"), vjust=0.5, size=titleSize) +
    # geom_text(data=group_p,aes(factor(y),11,label=P, fontface="bold"),vjust=0.5, hjust=1, size=dataSize) +
    geom_text(data=group_data,aes(factor(ID),5,label=p.value),vjust=0.5, hjust=1, size=dataSize) +
    geom_segment(aes(x = 2.5, y = 0, xend = 2.5, yend = 5.5)) + 
    # geom_segment(aes(x = 2, y = 1, xend = 2, yend = max(hazard_data[!is.na(hazard_data$HR_high),]$HR_high)),arrow=arrow(),linetype=1,size=0.3) + 
    # geom_segment(aes(x = 2, y = 1, xend = 2, yend = min(hazard_data[!is.na(hazard_data$HR_low),]$HR_low)),arrow=arrow(),linetype=1,size=0.3) + 
    theme_bare

  ## LEFT PANEL WITH NORMAL SCALE
  leftPanel<-haz + 
    apply(hl_rows,1,function(x)annotation_custom(hl_rect(x["col"],alpha=0.4),as.numeric(x["ID"])-0.5,as.numeric(x["ID"])+0.5,-20,20)) +
    coord_flip(ylim=c(0,5.5)) +
    geom_point(aes(x=factor(ID),y=1),shape=3,alpha=0,vjust=0) + 
    geom_text(data=group_mirna,aes(factor(y),0.5,label=miRNA, fontface="bold"),vjust=0.5, hjust=0, size=dataSize) +
    geom_text(data=group_data,aes(factor(ID),1,label=Group),vjust=0.5, hjust=0, size=dataSize) +
    geom_text(data=group_data,aes(factor(ID),5,label=NoP),vjust=0.5, hjust=1, size=dataSize) +
    geom_text(data=lf_labels,aes(x,y,label=lab, fontface="bold"), vjust=0.5, hjust=0, size=4, size=titleSize) +
    geom_segment(aes(x = 2.5, y = 0, xend = 2.5, yend = 5.5)) + 
    theme_bare

  ## PLOT THEM BOTH IN A GRID SO THEY MATCH UP
  height <- 4000 * (floor((dim(my_plot_data)[1] - 1) / 5) + 1)
  resolution <- 400 + 100 * (floor((dim(my_plot_data)[1] - 1) / 5) + 1)
  print(hazard_data)
  png(paste(output, '/miR15/TCGA.', sample,'.miR15Family.all.HR.png', sep=""), height=height, width=6000, res=resolution)
  p <- ggarrange(leftPanel,middle_panel,rightPanel, widths=c(3,4,3), ncol=3, nrow=1)
  print(p)
  dev.off() 
}

my_result_mirna_sig <- my_result_mirna[which(my_result_mirna$p.value_U < 0.05 | my_result_mirna$p.value < 0.05),]

# out <- my_result_mirna_sig[,c(1,2,3,7,8,12,13)]
out <- my_result_mirna_sig[,c(1,2,3,10,12,23,25)]

colnames(out) <- c('miRNA', 'High', 'Low', 'HR (95% CI for HR)_u', 'p.value_u', 'HR (95% CI for HR)', 'p.value')
out$p.value_u <- format(out$p.value_u, scientific = TRUE)

out$p.value <- as.numeric(out$p.value)
out$p.value_u <- as.numeric(out$p.value_u)
if(dim(out[which(out$p.value < 0.001),])[1] > 0){
    out[which(out$p.value < 0.001),]$`HR (95% CI for HR)` <- paste0(out[which(out$p.value < 0.001),]$`HR (95% CI for HR)`, ' ***')

  }
if(dim(out[which(out$p.value > 0.001 & out$p.value < 0.01),])[1] > 0){
  out[which(out$p.value > 0.001 & out$p.value < 0.01),]$`HR (95% CI for HR)` <- paste0(out[which(out$p.value > 0.001 & out$p.value < 0.01),]$`HR (95% CI for HR)`, ' **')
}
if(dim(out[which(out$p.value > 0.01 & out$p.value < 0.05),])[1] > 0){
  out[which(out$p.value > 0.01 & out$p.value < 0.05),]$`HR (95% CI for HR)` <- paste0(out[which(out$p.value > 0.01 & out$p.value < 0.05),]$`HR (95% CI for HR)`, ' *')
}

if(dim(out[which(out$p.value_u < 0.001),])[1] > 0){
    out[which(out$p.value_u < 0.001),]$`HR (95% CI for HR)_u` <- paste0(out[which(out$p.value_u < 0.001),]$`HR (95% CI for HR)_u`, ' ***')

  }
if(dim(out[which(out$p.value_u > 0.001 & out$p.value_u < 0.01),])[1] > 0){
  out[which(out$p.value_u > 0.001 & out$p.value_u < 0.01),]$`HR (95% CI for HR)_u` <- paste0(out[which(out$p.value_u > 0.001 & out$p.value_u < 0.01),]$`HR (95% CI for HR)_u`, ' **')
}
if(dim(out[which(out$p.value_u > 0.01 & out$p.value_u < 0.05),])[1] > 0){
  out[which(out$p.value_u > 0.01 & out$p.value_u < 0.05),]$`HR (95% CI for HR)_u` <- paste0(out[which(out$p.value_u > 0.01 & out$p.value_u < 0.05),]$`HR (95% CI for HR)_u`, ' *')
}

out$p.value <- factor(out$p.value)
out$p.value_u <- factor(out$p.value_u)

out$High <- paste0(out$High, '(Ref)')
# out$miRNA <- paste0(my_result_mirna_sig$TCGA_ID, '_', out$miRNA)
out <- as.matrix(out)
rownames(out) <- out[,1]
cgroup <- c('NoP', 'Univariate Model', 'Multivariates Model<br />(miR15 Family)')
rgroup <- unique(my_result_mirna_sig$TCGA_ID)
n.rgroup <- as.numeric(as.matrix(table(my_result_mirna_sig$TCGA_ID))[,1])

out <- out[,-1]


colnames(out) <- c('High', 'Low', 'HR<br />(95% CI for HR)', '<i>p</i> value', 'HR<br />(95% CI for HR)', '<i>p</i> value')
colnames(out) <- sprintf('<b>%s</b>', colnames(out))

out %>% 
  addHtmlTableStyle(col.rgroup = c("none", "#F7F7F7"),
                    css.cell = "padding-left: .5em; padding-right: .5em; line-height: 1.8;color:black",
                    align="r") %>% 
  htmlTable(rowlabel = 'Tumor/miRNA', 
            rgroup=rgroup,
            n.rgroup=n.rgroup,
            ctable = TRUE, align = 'cccccc',
            ## number of columns that each cgroup label spans:
            n.cgroup = c(2, 2, 2), cgroup = cgroup,
            caption = "<font size=3 color=black>Table 1: Hazard ratios and <i>p</i> values of two cox models with miR15 family vs OS in different tumors</font>", 
            tfoot = '<font size=1><sup>&dagger;</sup>* <i>p</i><0.05,** <i>p</i><0.01, *** <i>p</i><0.001.</font>') %>% 
  save_kable(file = "output/TCGA.all.miR15Family.HR.table.png",zoom = 4)


#Figure
# miR15 Targets
mir15.target <- read.csv('output/TCGA.miR.cluster.pred.csv')

covariates <- unique(mir.table[which(mir.table$miRNA %in% mir.15),]$miRNA_ID)

covariates <- lapply(covariates, function(x){
  return(str_replace_all(x, '_', '-'))
  })

covariates <- unlist(covariates)

res_lt4 <- mir15.target[which(mir15.target$mirtarbase > 4),]

res_lt4.mir <- res_lt4[which(res_lt4$mature_mirna_id %in% covariates),]

gene <- unique(res_lt4.mir$target_symbol)

gene.df <- bitr(gene, fromType = "SYMBOL",
        toType = c("ENSEMBL","ENTREZID"),
        OrgDb = org.Hs.eg.db)

geneList <- gene.df$ENTREZID
geneList <- geneList[which(duplicated(geneList)==FALSE)]

ego <- enrichGO(gene          = geneList,
                OrgDb         = org.Hs.eg.db,
                ont           = "ALL",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                keyType = 'ENTREZID')

ekegg <- enrichKEGG(gene          = geneList,
                  organism = "hsa",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.01,
                  qvalueCutoff  = 0.05)

ego <- setReadable(ego, OrgDb = org.Hs.eg.db)
ekegg <- setReadable(ekegg, OrgDb = org.Hs.eg.db, keyType = 'ENTREZID')


tmp <- gene.df[match(geneList, gene.df$ENTREZID),]
tmp <- common.mrna[match(tmp$SYMBOL, common.mrna$Symbol),]
FClist <- tmp$logFC
names(FClist) <- gene.df[match(geneList, gene.df$ENTREZID),]$ENTREZID

# edo2 <- gseGO(FClist[order(-FClist)], OrgDb = org.Hs.eg.db, pvalueCutoff = 0.05)
png(paste(output, '/miR15/TCGA.all.miR15Family.predicted',".go.png", sep=""),height=2000, width=4000, res=400)
dotplot(ego, showCategory=50)
dev.off()


png(paste(output, '/miR15/TCGA.all.miR15Family.predicted',".cnet.go.png", sep=""), height=4000, width=4000, res=400)
cnetplot(ego,showCategory=50,circular=TRUE,colorEdge=FALSE)
dev.off()

png(paste(output, '/miR15/TCGA.all.miR15Family.predicted',".cnet.go.cluster.png", sep=""), height=4000, width=4000, res=400)
cnetplot(ego,colorEdge=FALSE)
dev.off()

write.csv(ego@result, 'output/miR15/TCGA.all.miR15Family.predicted.go.csv')

png(paste(output, '/miR15/TCGA.all.miR15Family.predicted',".kegg.png", sep=""),height=2000, width=4000, res=400)
dotplot(ekegg, showCategory=50)
dev.off()

png(paste(output, '/miR15/TCGA.all.miR15Family.predicted',".cnet.kegg.png", sep=""), height=4000, width=4000, res=400)
cnetplot(ekegg,showCategory=30,circular=TRUE,colorEdge=FALSE)
dev.off()

png(paste(output, '/miR15/TCGA.all.miR15Family.predicted',".cnet.kegg.cluster.png", sep=""), height=4000, width=4000, res=400)
cnetplot(ekegg,colorEdge=FALSE)
dev.off()

write.csv(ekegg@result, 'output/miR15/TCGA.all.miR15Family.predicted.kegg.csv')


#Figure mRNA expression in different tumors

mrna.top <- NULL
my_files <- list.files('output/Pan_Cancer')
my_files <- my_files[which(grepl('contrast_Primary_Tumor_v_Solid_Tissue_Normal.mRNA.top.genes', my_files))]
for(my_file in my_files){
  top.genes.tmp <- read.csv(paste0('output/Pan_Cancer/',my_file))
  top.genes.tmp$TCGA_ID <- str_split(my_file, '\\.')[[1]][2]
  if(is.null(mrna.top)){
    mrna.top <- top.genes.tmp
  }else{
    mrna.top <- rbind(mrna.top, top.genes.tmp)
  }
}

mrna.top.sig <- mrna.top[which(mrna.top$X %in% res_lt4.mir$target_symbol),]


# top.genes.mir15 <- mrna.top.sig


data <- NULL
gene.df <- NULL
comparisons <- list()
total_com <- list()
for(sample in unique(mrna.top.sig$TCGA_ID)){
  # sample <- 'TCGA_BRCA'
  print(sample)
  mrna_data <- fread(paste0('data/',sample, '.htseq_fpkm_uq.tsv.gz'))

  # mrna_data <- fread('TCGA_BRCA.htseq_fpkm_uq.tsv.gz')
  mrna_data <- as.data.frame(mrna_data)
  mrna_data$Ensembl_ID <- lapply(mrna_data$Ensembl_ID,  sub, pattern = "\\.\\d+$", replacement = "")

  rownames(mrna_data) <- mrna_data$Ensembl_ID
  mrna_data <- mrna_data[,-1]

  my_genes <- rownames(mrna_data)
  if(is.null(gene.df)){
      gene.df <- bitr(my_genes, fromType = "ENSEMBL",
      toType = c("ENSEMBL","ENTREZID", "SYMBOL"),
      OrgDb = org.Hs.eg.db)
  }

  mrna_data <- mrna_data[gene.df$ENSEMBL,]

  mrna_data <- 2^mrna_data - 1
  mrna_data$Symbol <- gene.df$SYMBOL

  mrna_data <- data.table(mrna_data)

  mrna_data <- mrna_data[, lapply(.SD,sum), by=.(Symbol)]
  mrna_data <- as.data.frame(mrna_data)
  rownames(mrna_data) <- mrna_data$Symbol
  mrna_data <- mrna_data[,-1]

  mrna_data <- log1p(mrna_data)

  tmp <- pheno.df[which(pheno.df$TCGA_ID == sample),]
  cc <- c('Primary Tumor', 'Solid Tissue Normal')
  # tmp <- tmp[match(colnames(tmp2), tmp[which(tmp$`sample_type.samples` %in% cc),]$submitter_id.samples),]
  tmp <- tmp[which(tmp$sample_type.samples %in% cc),]
  tmp1 <- intersect(tmp$submitter_id.samples, colnames(mrna_data))

  mirnas <- top.mir15[which(top.mir15$TCGA_ID == sample),]$X

  covariates <- mrna.top.sig[which(mrna.top.sig$TCGA_ID == sample),]$X

  covariates <- intersect(covariates, unique(res_lt4[which(res_lt4$mature_mirna_id %in% mirnas),]$target_symbol))
  print(covariates)

  if(length(covariates) == 0) next
  tmp.data <- mrna_data[covariates,][,tmp1]

  tmp.data <- as.data.frame(t(tmp.data))

  tmp.data$Sample_type <- ''
  tmp.data$Sample_type <- tmp[match(rownames(tmp.data),tmp$submitter_id.samples),]$sample_type.samples

  # sample <- 'TCGA_BRCA'
  # tmp.data <- my_df[which(my_df$TCGA_ID == sample),]

  if(dim(tmp.data)[1]==0) next
  # covariates <- unique(res_lt4.mir$target_symbol)

  tmp.data.normal <- tmp.data[which(tmp.data$Sample_type == 'Solid Tissue Normal'),]
  tmp.data.tumor <- tmp.data[which(tmp.data$Sample_type == 'Primary Tumor'),]

  if(dim(tmp.data.normal)[1] == 0 || dim(tmp.data.tumor)[1] == 0) next

  tmp.normal <- melt(tmp.data.normal[,covariates,drop=FALSE])
  tmp.tumor <- melt(tmp.data.tumor[,covariates,drop=FALSE])

  tmp.normal$Sample_type <- 'N'
  tmp.tumor$Sample_type <- 'T'

  # tmp.normal$Sample_type <- paste0(sample, '_N')
  # tmp.tumor$Sample_type <- paste0(sample, '_T')

  # comparisons[[sample]] <- c(paste0(sample, '_T'), paste0(sample, '_N'))
  comparisons[[sample]] <-  list()
  # print(comparisons)

  print(tmp.normal[1:2,])
  print(tmp.tumor[1:2,])
  my_tmp.data <- rbind(tmp.normal,tmp.tumor)
  print(class(my_tmp.data))
  my_tmp.data$TCGA_ID <- sample
  my_tmp.data$p.value <- 0.05
  my_tmp.data$Group <- ''
  my_tmp.data.mir15 <- mrna.top.sig[which(mrna.top.sig$TCGA_ID == sample),]
  for(m in unique(covariates)){
    p.value <- my_tmp.data.mir15[which(my_tmp.data.mir15$X == m),]$FDR
    # m <- str_replace_all(m, '-', '_')
    print(p.value)
    my_tmp.data[which(my_tmp.data$variable == m),]$Group <- paste0(m, '_',my_tmp.data[which(my_tmp.data$variable == m),]$Sample_type)
    my_tmp.data[which(my_tmp.data$variable == m),]$p.value <- p.value
    comparisons[[sample]][[m]] <- c(paste0(m, '_T'), paste0(m, '_N'))

    total_com[[paste0(sample,'_',m)]] <- c(paste0(m, '_T'), paste0(m, '_N'))
  }

  print(my_tmp.data[1:2,])
  if(is.null(data)){
    data <- my_tmp.data
  }else{
    print(data[1:2,])
    data <- rbind(data,my_tmp.data)
  }

}

# fig_data <- data[which(data$TCGA_ID == 'TCGA_BRCA'),]
# png(paste('output/miR15/TCGA.all',".primary.tumor.vs.normal.miR15Family.boxplot.png", sep=""),height=5000, width=20000, res=400)

my_plots_data <- c()
widths <- c()
for(sample in unique(data$TCGA_ID)){
  print(sample)
  fig_data <- data[which(data$TCGA_ID == sample),]
  my_plots_data <- c(my_plots_data, sample)
  widths <- c(widths, length(unique(fig_data$variable)))
}

myplot <- function(sample){
    # sample <- 'TCGA_BRCA'
    fig_data <- data[which(data$TCGA_ID == sample),]
    
    sum.data <- ddply(fig_data,.(variable),summarise,value=max(value)+0.3,p.value=max(p.value))
    annotation <- c()
    y_position <- c()
    xmin <- c()
    xmax <- c()
    print(dim(sum.data))
    for(idx in 1:dim(sum.data)[1]){
      print(idx)
      print(sum.data[idx,'value'])
      y_position <- c(y_position, sum.data[idx,'value'])

      xmin = c(xmin, 1+(idx-1)*2)
      xmax = c(xmax, 2+(idx-1)*2)
      p <- sum.data[idx, 'p.value']
      if(p > 0.01){
        annotation <- c(annotation, "*")
      }else if(p > 0.001){
        annotation <- c(annotation, "**")
      }else{
        annotation <- c(annotation, "***")
      }
    }
    print(y_position)
    print(annotation)
    print(xmin)
    ggplot(fig_data, aes(x=Group, y=value, color=variable)) + 
    # geom_boxplot() +
    geom_violin(trim = FALSE, aes(fill = factor(variable))) +
    # geom_signif(y_position=y_position, xmin=xmin, xmax=xmax,annotation = annotation) +
    facet_grid(~TCGA_ID)+
    geom_signif(test="wilcox.test",comparisons = comparisons[[sample]], 
                map_signif_level=TRUE) +
    # geom_signif(y_position=y_position, xmin=xmin, xmax=xmax,annotation = annotation) +
     # geom_text(size=2) +
    # geom_dotplot(binaxis='y', stackdir='center', stackratio=1.5, dotsize=0.5, color='black', fill='black') +
    theme(legend.position="none", axis.text.x = element_text(angle = 45, vjust = 0.7, hjust=0.8, size=12)) +
    ylab("Gene Expression Level [log2(fpkm+1)]") +
    xlab("")
}

my_plots <- map(my_plots_data, myplot)

png(paste('output/miR15/TCGA.all',".primary.tumor.vs.normal.miR15Family.targets.filtered.v2.boxplot.png", sep=""),height=12000, width=16000, res=800)
ggarrange(plotlist=my_plots, nrow=3, ncol=4)
dev.off()


# Figure mRNA cox univariate

mrna.top.mir.15 <- mrna.top[which(mrna.top$X %in% res_lt4.mir$target_symbol),]

result_mrna <- NULL
for(sample in unique(mrna.top.mir.15$TCGA_ID)){
  mrna.top.mir.15.sub <- mrna.top.mir.15[which(mrna.top.mir.15$TCGA_ID == sample),]
  my_genes <- mrna.top.mir.15[which(mrna.top.mir.15$TCGA_ID == sample),]$X
  mirnas <- top.mir15[which(top.mir15$TCGA_ID == sample),]$X
  my_genes <- intersect(my_genes, unique(res_lt4[which(res_lt4$mature_mirna_id %in% mirnas),]$target_symbol))

  # my_genes <- unique(res_lt5[which(res_lt5$mature_mirna_id %in% top.genes.mir15.sub$mature_miRNA),]$target_symbol)
  if(length(my_genes) < 1) next
  mrna_data <- fread(paste0('data/',sample,'.htseq_fpkm_uq.tsv.gz'))
  mrna_data <- as.data.frame(mrna_data)
  mrna_data$Ensembl_ID <- lapply(mrna_data$Ensembl_ID,  sub, pattern = "\\.\\d+$", replacement = "")
  rownames(mrna_data) <- mrna_data$Ensembl_ID
  mrna_data <- mrna_data[,-1]
  gene.df <- bitr(my_genes, fromType = "SYMBOL",
        toType = c("ENSEMBL","ENTREZID"),
        OrgDb = org.Hs.eg.db)
  gene.df <- as.data.frame(gene.df)
  common.genes <- intersect(gene.df$ENSEMBL, rownames(mrna_data))
  gene.df <- gene.df[match(common.genes, gene.df$ENSEMBL),]
  mrna_data.sub <- mrna_data[match(common.genes, rownames(mrna_data)),]

  mrna_data.sub <- 2^mrna_data.sub - 1

  mrna_data.sub$Symbol <- gene.df$SYMBOL
  mrna_data.sub <- data.table(mrna_data.sub)
  mrna_data.sub <- mrna_data.sub[, lapply(.SD,sum), by=.(Symbol)]
  mrna_data.sub <- as.data.frame(mrna_data.sub)
  rownames(mrna_data.sub) <- mrna_data.sub$Symbol
  mrna_data.sub <- mrna_data.sub[,-1]

  mrna_data.sub <- log1p(mrna_data.sub)

  mrna_data.sub.t <- as.data.frame(t(mrna_data.sub))

  common.samples <- intersect(rownames(mrna_data.sub.t), survival.df$sample)

  mrna_data.sub.t <- mrna_data.sub.t[common.samples,c(1:dim(mrna_data.sub.t)[2]),drop=FALSE]
  mrna_data.sub.t <- cbind(mrna_data.sub.t, survival.df[match(common.samples, survival.df$sample),][,-1])

  print('Finish combining sample')

  common.samples <- intersect(rownames(mrna_data.sub.t), pheno.df$`submitter_id.samples`)

  mrna_data.sub.t$Sample_type <- ''

  mrna_data.sub.t[match(common.samples, rownames(mrna_data.sub.t)),]$Sample_type <- pheno.df[match(common.samples, pheno.df$`submitter_id.samples`),]$`sample_type.samples`

  my_df.sub <- mrna_data.sub.t
  names(my_df.sub)[names(my_df.sub) == "OS.time"] <- "OS_time"

  my_df.sub <- my_df.sub[which(my_df.sub$Sample_type == 'Primary Tumor'),]
  covariates <- unique(rownames(mrna_data.sub))

  print(colnames(my_df.sub))
  print(covariates)
  my_means <- colMedians(as.matrix(my_df.sub[,covariates]))
  for(i in c(1:length(covariates))){
    v <- my_means[i]
    high <- my_df.sub[,covariates[i]] >= v
    my_df.sub[,covariates[i]][which(high==TRUE)] <- 'high'
    my_df.sub[,covariates[i]][which(high==FALSE)] <- 'low'
  }

  print('Start to analyze surv')
  # print(my_df.sub)
  univ_formulas <- sapply(covariates, function(x) as.formula(paste('Surv(OS_time, OS)~', x)))
  model_names <- names(univ_formulas)
  # print(model_names)
  i <- 1
  univ_models <- list()
  for(m in univ_formulas){
    skip_to_next <- FALSE
    p <- tryCatch({
        print(m)    
        print(model_names[i])
        fit <- coxph(m, data = my_df.sub)
        univ_models[[model_names[i]]] <- fit
        i <- i+1 
        
      },error = function(e){
        i <- i+1
        print(sample)
        print('\n')
        skip_to_next <- TRUE
        }, 
      warning = function(cond) {
        message(cond)
        print(sample)
        print('\n')
        i <- i+1
        skip_to_next <- TRUE
      })
    # print(p)
    if(skip_to_next == TRUE){next}
  
  } 
  j <- 1
  km_models <- list()
  for(m in univ_formulas){
    skip_to_next <- FALSE
    p <- tryCatch({
        print(m)    
        print(model_names[j])
        fit <- surv_fit(m, data = my_df.sub)
        km_models[[model_names[j]]] <- fit
        j <- j+1 
        
      },error = function(e){
        i <- i+1
        print(sample)
        print('\n')
        skip_to_next <- TRUE
        }, 
      warning = function(cond) {
        message(cond)
        print(sample)
        print('\n')
        j <- j+1
        skip_to_next <- TRUE
      })
    # print(p)
    if(skip_to_next == TRUE){next}
  } 


  print('Start to analyze coxph')
  univ_results <- lapply(univ_models,
                         function(x){ 
                          fit.coxph <- x
                          x <- summary(x)
                          print(x)
                          p.value<-signif(x$coefficients[,'Pr(>|z|)'], digits=2)
                          wald.test<-signif(x$wald["test"], digits=2)
                          beta<-signif(x$coef[1], digits=2);#coeficient beta
                          HR <-signif(x$coef[2], digits=2);#exp(beta)
                          total <- x$n
                          nevent <- x$nevent
                          HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
                          HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
                          HR.confint <- paste0(HR, " (", 
                                     HR.confint.lower, "-", HR.confint.upper, ")")
                          file_name <- rownames(x$conf.int)[1]
                          file_name <- str_split(file_name, 'l')[[1]][1]
                        
                          expr_gene <- table(my_df.sub[,file_name])
                          expr_gene <- as.data.frame(expr_gene)
                          tmp<-c(file_name, expr_gene[1,2],expr_gene[2,2],total, nevent, beta, HR, HR.confint.lower, HR.confint.upper, HR.confint, wald.test, p.value, sample)
                     
                          names(tmp)<-c("Symbol",levels(expr_gene[1,1])[expr_gene[1,1]],levels(expr_gene[2,1])[expr_gene[2,1]],"Total", "Event", "beta", "HR", "HR_low", "HR_high","HR (95% CI for HR)", "wald.test", 
                                                  "p.value", "TCGA_ID")
                          return(tmp)
                           })
  univ_km <- lapply(km_models, function(x){
    my_fit <- x
    # print(fit)
    x <- summary(x)
    print(x$table)
    file_name <- rownames(x$table)[1]
    file_name <- str_split(file_name, '=')[[1]][1]
    print(file_name)
    png(paste(output, '/miR15/TCGA.', sample,'.', file_name,'.KM.png', sep=""), height=6000, width=6000, res=1000)
    print(my_fit)
    sample_s <- str_split(sample,'_')[[1]][[2]]
    p <- ggsurvplot(
          my_fit,                     
          data = my_df.sub,            
          risk.table = TRUE,       
          risk.table.col = file_name, 
          pval = TRUE,         
          conf.int = TRUE,     
          palette = "Dark2",
          xlab = "Time in days", 
          ggtheme = theme_bw(), 

          conf.int.style = "step",
          legend.title = paste0(sample_s, ': \n',file_name),
          legend.labs = c("High", "Low"))
    print(p)
    dev.off()
    return(file_name)
    })
  tmp <- t(as.data.frame(univ_results, check.names = FALSE))
  print(tmp)
  tmp <- as.data.frame(tmp)

  if(dim(tmp)[1] == 0) next

  if (is.null(result_mrna)){
    result_mrna <- tmp
  }else{
    result_mrna <- rbind(result_mrna, tmp)
  }

  rownames(result_mrna) <- c(1:dim(result_mrna)[1])

}

write.csv(result_mrna, 'output/miR15/TCGA.mRNA.HR.miR15Family.csv')

result_mrna <- read.csv('output/miR15/TCGA.mRNA.HR.miR15Family.csv')[,-1]

result_mrna <- result_mrna[which(result_mrna$p.value < 0.05),]
for(sample in unique(result_mrna$TCGA_ID)){
  my_plot_data <- result_mrna[which(result_mrna$TCGA_ID == sample),]

  my_plot_data$Index <- as.numeric(c(1:dim(my_plot_data)[1]))
  my_plot_data$HR <- as.numeric(my_plot_data$HR)
  my_plot_data$HR_low <- as.numeric(my_plot_data$HR_low)
  my_plot_data$HR_high <- as.numeric(my_plot_data$HR_high)


  ############################################
  ### CUSTOMIZE APPEARANCE WITH THESE     ####
  ############################################
  blankRows<-2    # blank rows under boxplot
  titleSize<-4
  dataSize<-4
  boxColor<-"pink"
  ############################################
  ############################################

  ## BASIC THEMES (SO TO PLOT BLANK GRID)
  theme_grid <- theme(
    axis.line = element_blank(), 
    axis.text.x = element_blank(), 
    axis.text.y = element_blank(),
    axis.ticks = element_blank(), 
    axis.title.x = element_blank(), 
    axis.title.y = element_blank(), 
    axis.ticks.length = unit(0.0001, "mm"),
    axis.ticks.margin = unit(c(0,0,0,0), "lines"), 
    legend.position = "none", 
    panel.background = element_rect(fill = "transparent"), 
    panel.border = element_blank(), 
    panel.grid.major = element_line(colour="grey"), 
    panel.grid.minor = element_line(colour="grey"), 
    panel.margin = unit(c(-0.1,-0.1,-0.1,-0.1), "mm"), 
    plot.margin = unit(c(5,0,5,0.01), "mm")
  )

  theme_bare <- theme_grid +
    theme(
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank()
    )
  group_low <- my_plot_data[,c(1,3,7,8,9,10,12)]
  group_high <- my_plot_data[,c(1,2,7,8,9,10,12)]
  colnames(group_low) <- c('Gene', 'NoP', 'HR', 'HR_low', 'HR_high', 'HR (95% CI for HR)', 'p.value')
  colnames(group_high) <- c('Gene', 'NoP', 'HR', 'HR_low', 'HR_high', 'HR (95% CI for HR)', 'p.value')
  group_low$Group <- 'Low'
  group_high$Group <- 'High'
  group_high[,c(3,4,5)] <- 1
  group_high[,c(6,7)] <- NA
  group_data <- rbind(group_low, group_high)
  group_data <- group_data[order(group_data$Gene),]
  print(group_data)
  rownames(group_data) <- c(1:dim(group_data)[1])
  group_data$ID <- c(1:dim(group_data)[1])
  hazard_data<-expand.grid(ID=1:nrow(group_data),HR=1)
  hazard_data$HR <- group_data$HR
  hazard_data<-rbind(hazard_data,ddply(group_data,.(Gene),summarise,ID=max(ID)+0.1,HR=NA)[,2:3])

  hazard_data<-rbind(hazard_data,data.frame(ID=c(0,-1:(-2-blankRows),max(group_data$ID)+1,max(group_data$ID)+2),HR=NA))

  hazard_data$HR_low <- NA
  hazard_data$HR_high <- NA

  hazard_data[1:dim(group_data)[1],]$HR_low <- group_data$HR_low
  hazard_data[1:dim(group_data)[1],]$HR_high <- group_data$HR_high

  hr_labels <- group_data[,c('ID', 'HR (95% CI for HR)')]
  colnames(hr_labels) <- c("ID", "lab")

  upper <- max(hazard_data[!is.na(hazard_data$HR_high),]$HR_high)
  lower <- min(hazard_data[!is.na(hazard_data$HR_low),]$HR_low)
  lower <- min(0.4, lower)
  upper <- max(2.8, upper)

  upper_int <- round(((upper - 1)/2),2)
  lower_int <- round(((1 - lower)/2),2)

  scale_data <- data.frame(ID=0,HR=c((1-2*lower_int-0.05),(1-lower_int), 1, (1 + upper_int), (1+2*upper_int+0.2)))

  group_mirna<-ddply(group_data,.(Gene),summarise,y=max(ID)+0.1)


  hl_rows <- data.frame(ID=(1:floor(length(unique(hazard_data$ID[which(hazard_data$ID>0)]))/2))*2,col="lightgrey")
  hl_rows$ID <- hl_rows$ID+blankRows+1

  hl_rect <- function(col="white",alpha=0.5){
    rectGrob(   x = 0, y = 0, width = 1, height = 1, just = c("left","bottom"), gp=gpar(alpha=alpha, fill=col))
  }

  ## DATA FOR TEXT LABELS
  md_labels <- data.frame(x=c(rep(length(unique(hazard_data$ID))-0.2,times=1)),
                        y=c(1),
                        lab=c("Hazard Ratio"),drop=FALSE)

  rt_labels <- data.frame(x=c(rep(length(unique(hazard_data$ID))-0.2,times=2)),
                        y=c(1,4),
                        lab=c("Hazard Ratio\n(95% CI)","P Value"))

  lf_labels <- data.frame(x=c(rep(length(unique(hazard_data$ID))-0.2,times=2)),
                       y=c(0.5,4),
                       lab=c("Gene","No. of\nPatients"))

  legend_labels <- data.frame(x=c(rep(1,times=2)),
                       y=c((1-lower_int),(1+upper_int)),
                       lab=c("Low Better","High Better"))

  haz <- ggplot(hazard_data,aes(factor(ID),HR))+ labs(x=NULL, y=NULL) 

  ## MIDDLE PANEL WITH LOG SCALE
  middle_panel <- haz +
    apply(hl_rows,1,function(x)annotation_custom(hl_rect(x["col"],alpha=0.4),as.numeric(x["ID"])-0.5,as.numeric(x["ID"])+0.5,-20,20)) +
    # geom_segment(aes(x = 2, y = 1, xend = 1.5, yend = 1)) + 
    geom_hline(aes(yintercept=1),linetype=2, size=0.5)+
    geom_point() + 
    geom_errorbar(data=hazard_data, aes(ymin=HR_low, ymax=HR_high), width=.3)+
    # geom_boxplot(fill=boxColor,size=0.5, alpha=0.8)+ 
    scale_y_log10() + 
    # scale_y_continuous() +
    coord_flip() +
    geom_text(data=scale_data,aes(3,HR,label=HR),vjust=0.5, size=dataSize) +
    geom_text(data=md_labels,aes(x,y,label=lab, fontface="bold"), vjust=0.5, size=titleSize) +
    # geom_text(data=hr_labels,aes(factor(ID),0.4,label=lab),vjust=0.5, hjust=1, size=dataSize) +
    # geom_text(data=group_p,aes(factor(y),11,label=P, fontface="bold"),vjust=0.5, hjust=1, size=dataSize) +
    # geom_text(data=group_data,aes(factor(ID),5,label=p.value),vjust=0.5, hjust=1, size=dataSize) +
    geom_text(data=legend_labels,aes(x,y,label=lab, fontface="bold"),hjust=0.5, vjust=1, size=titleSize) +
    geom_point(data=scale_data,aes(2.5,HR),shape=3,size=3) + 
    geom_point(aes(2,6),shape=3,alpha=0,vjust=0) + 
    # geom_errorbar(width=.08)
    geom_segment(aes(x = 2.5, y = 0, xend = 2.5, yend = max(8,max(hazard_data[!is.na(hazard_data$HR_high),]$HR_high)))) + 
    geom_segment(aes(x = 2, y = 1, xend = 2, yend = max(hazard_data[!is.na(hazard_data$HR_high),]$HR_high)),arrow=arrow(),linetype=1,size=0.3) + 
    geom_segment(aes(x = 2, y = 1, xend = 2, yend = min(hazard_data[!is.na(hazard_data$HR_low),]$HR_low)),arrow=arrow(),linetype=1,size=0.3) + 
    theme_bare


  ## RIGHT PANEL WITH LOG SCALE
  rightPanel <- haz + 
    apply(hl_rows,1,function(x)annotation_custom(hl_rect(x["col"],alpha=0.4),as.numeric(x["ID"])-0.5,as.numeric(x["ID"])+0.5,-20,20)) +
    # scale_y_log10() +
    coord_flip(ylim=c(0,5.5)) +
    geom_point(aes(x=factor(ID),y=1),shape=3,alpha=0,vjust=0) + 
    geom_text(data=hr_labels,aes(factor(ID),3,label=lab),vjust=0.5, hjust=1, size=dataSize) +
    geom_text(data=rt_labels,aes(x,y,label=lab, fontface="bold"), vjust=0.5, size=titleSize) +
    # geom_text(data=group_p,aes(factor(y),11,label=P, fontface="bold"),vjust=0.5, hjust=1, size=dataSize) +
    geom_text(data=group_data,aes(factor(ID),5,label=p.value),vjust=0.5, hjust=1, size=dataSize) +
    geom_segment(aes(x = 2.5, y = 0, xend = 2.5, yend = 5.5)) + 
    # geom_segment(aes(x = 2, y = 1, xend = 2, yend = max(hazard_data[!is.na(hazard_data$HR_high),]$HR_high)),arrow=arrow(),linetype=1,size=0.3) + 
    # geom_segment(aes(x = 2, y = 1, xend = 2, yend = min(hazard_data[!is.na(hazard_data$HR_low),]$HR_low)),arrow=arrow(),linetype=1,size=0.3) + 
    theme_bare

  ## LEFT PANEL WITH NORMAL SCALE
  leftPanel<-haz + 
    apply(hl_rows,1,function(x)annotation_custom(hl_rect(x["col"],alpha=0.4),as.numeric(x["ID"])-0.5,as.numeric(x["ID"])+0.5,-20,20)) +
    coord_flip(ylim=c(0,5.5)) +
    geom_point(aes(x=factor(ID),y=1),shape=3,alpha=0,vjust=0) + 
    geom_text(data=group_mirna,aes(factor(y),0.5,label=Gene, fontface="bold"),vjust=0.5, hjust=0, size=dataSize) +
    geom_text(data=group_data,aes(factor(ID),1,label=Group),vjust=0.5, hjust=0, size=dataSize) +
    geom_text(data=group_data,aes(factor(ID),5,label=NoP),vjust=0.5, hjust=1, size=dataSize) +
    geom_text(data=lf_labels,aes(x,y,label=lab, fontface="bold"), vjust=0.5, hjust=0, size=4, size=titleSize) +
    geom_segment(aes(x = 2.5, y = 0, xend = 2.5, yend = 5.5)) + 
    theme_bare

  ## PLOT THEM BOTH IN A GRID SO THEY MATCH UP
  height <- 4000 * (floor((dim(my_plot_data)[1] - 1) / 5) + 1)
  resolution <- 400 + 100 * (floor((dim(my_plot_data)[1] - 1) / 5) + 1)
  print(hazard_data)
  png(paste(output, '/TCGA.', sample,'.miR15Family.targets.all.HR.png', sep=""), height=height, width=6000, res=resolution)
  p <- ggarrange(leftPanel,middle_panel,rightPanel, widths=c(3,4,3), ncol=3, nrow=1)
  print(p)
  dev.off() 
}


out <- result_mrna[,c(1,2,3,10,12)]

colnames(out) <- c('Gene', 'High', 'Low', 'HR (95% CI for HR)', 'p.value')
out$p.value <- as.numeric(out$p.value)
if(dim(out[which((out$p.value * 10) < 0.01),])[1] > 0){
    out[which(out$p.value < 0.001),]$`HR (95% CI for HR)` <- paste0(out[which(out$p.value < 0.001),]$`HR (95% CI for HR)`, ' ***')

  }
if(dim(out[which(out$p.value > 0.001 & out$p.value < 0.01),])[1] > 0){
  out[which(out$p.value > 0.001 & out$p.value < 0.01),]$`HR (95% CI for HR)` <- paste0(out[which(out$p.value > 0.001 & out$p.value < 0.01),]$`HR (95% CI for HR)`, ' **')
}
if(dim(out[which(out$p.value > 0.01 & out$p.value < 0.05),])[1] > 0){
  out[which(out$p.value > 0.01 & out$p.value < 0.05),]$`HR (95% CI for HR)` <- paste0(out[which(out$p.value > 0.01 & out$p.value < 0.05),]$`HR (95% CI for HR)`, ' *')
}

out$p.value <- factor(out$p.value)

out$High <- paste0(out$High, '(Ref)')
# out$Gene <- paste0(result_mrna$TCGA_ID, '_', out$Gene)
out <- as.matrix(out)
rownames(out) <- out[,1]
cgroup <- c('NoP', 'Univariate Model')
rgroup <- unique(result_mrna$TCGA_ID)
n.rgroup <- as.numeric(as.matrix(table(result_mrna$TCGA_ID))[,1])

out <- out[,-1]


colnames(out) <- c('High', 'Low', 'HR<br />(95% CI for HR)', '<i>p</i> value')
colnames(out) <- sprintf('<b>%s</b>', colnames(out))

out %>% 
  addHtmlTableStyle(col.rgroup = c("none", "#F7F7F7"),
                    css.cell = "padding-left: .5em; padding-right: .5em; line-height: 1.8;color:black",
                    align="r") %>% 
  htmlTable(rowlabel = 'Tumor/Gene', 
            rgroup=rgroup,
            n.rgroup=n.rgroup,
            ctable = TRUE, align = 'cccc',
            ## number of columns that each cgroup label spans:
            n.cgroup = c(2, 2), cgroup = cgroup,
            caption = "<font size=3 color=black>Table 2</font>", 
            tfoot = '<font size=1><sup>&dagger;</sup>* <i>p</i><0.05,** <i>p</i><0.01, *** <i>p</i><0.001.</font>') %>% 
  save_kable(file = "output/miR15/TCGA.all.miR15Family.HR.targets.table.png",zoom = 4)



#Figure
# nets <- list()
colors <- list()
nodes_list <- list()
edges_list <- list()

res_lt4.mir.net <- res_lt4.mir[which(res_lt4.mir$target_symbol %in% mrna.top.sig$X),]

for(mir in unique(res_lt4.mir.net$mature_mirna_id)){
  # mir <- 'hsa-miR-646'
  res_lt4.sub <- res_lt4.mir.net[which(res_lt4.mir.net$mature_mirna_id == mir),]
  tmp1.t <- res_lt4.sub %>% rename(source=mature_mirna_id, destination=target_symbol)

  source <- data.frame('label'=unique(res_lt4.sub$mature_mirna_id))
  miRNA_count <- dim(source)[1]

  destination <- data.frame('label'=unique(res_lt4.sub$target_symbol))
  mRNA_count <- dim(destination)[1]

  nodes <- full_join(source, destination, by='label')
  nodes$id <- rownames(nodes)

  per_route <- tmp1.t
  edges <- per_route %>% left_join(nodes, by = c("source" = "label")) %>% rename(from = id)

  edges <- edges %>% left_join(nodes, by = c("destination" = "label")) %>% rename(to = id)
  edges <- edges[,c("from", "to")]


  edges.t <- matrix(0,nrow=(miRNA_count+mRNA_count),ncol=(miRNA_count+mRNA_count))
  for(i in 1:dim(edges)[1]){
    edges.t[as.numeric(edges[i,1]),as.numeric(edges[i,2])] <- 1
  }

  # net <- network(edges.t)
  # network.vertex.names(net) <- nodes$label
  nodes_list[[mir]] <- nodes
  edges_list[[mir]] <- edges.t
  # print(nodes_list)
  # print(edges_list)
  color <- c(rep('steelblue',miRNA_count), rep('orange',mRNA_count))
  colors[[mir]] <- color
  # set.network.attribute(net, "color", color)
  # nets[[mir]] <- net
}

myplot <- function(mir){
  edges.t <- edges_list[[mir]]
  nodes <- nodes_list[[mir]]
  # print(nodes)
  net <- network(edges.t)
  network.vertex.names(net) <- nodes$label
  color <- colors[[mir]]
  set.network.attribute(net, "color", color)
  ggnet2(net, size=16,label = TRUE, label.size = 5,color=color)
}

my_plots <- map(unique(res_lt4.mir.net[-which(res_lt4.mir.net$mature_mirna_id %in% c('hsa-miR-6838-5p', 'hsa-miR-646')),]$mature_mirna_id), myplot)

png(paste('output/miR15/',"miR15.targets.network.png", sep=""),height=10000, width=15000, res=600)
ggarrange(plotlist=my_plots, nrow=3, ncol=3)
dev.off()

# Figure
pr.data.all <- NULL
roc.data.all <- NULL
for(sample in unique(pheno.df$TCGA_ID)){

  # sample <- 'TCGA_BRCA'
  print(sample)
  mirna.df.sub.t <- as.data.frame(t(mirna.df.sub))
  common.samples <- intersect(rownames(mirna.df.sub.t), pheno.df[which(pheno.df$TCGA_ID == sample),]$submitter_id.samples)

  common.samples <- intersect(common.samples, survival.df$sample)

  mirna.df.sub.t <- mirna.df.sub.t[match(common.samples, rownames(mirna.df.sub.t)),]
  mirna.df.sub.t <- cbind(mirna.df.sub.t, survival.df[match(common.samples, survival.df$sample),][,-1])

  if(dim(mirna.df.sub.t)[1]==0) next

  print('Finish combining sample')

  common.samples <- intersect(rownames(mirna.df.sub.t), pheno.df$`submitter_id.samples`)

  mirna.df.sub.t$Sample_type <- ''

  mirna.df.sub.t[match(common.samples, rownames(mirna.df.sub.t)),]$Sample_type <- pheno.df[match(common.samples, pheno.df$`submitter_id.samples`),]$`sample_type.samples`


  my_data <- mirna.df.sub.t
  names(my_data)[names(my_data) == "OS.time"] <- "OS_time"

  covariates <- unique(mir.table[which(mir.table$miRNA %in% mir.15),]$miRNA_ID)

  my_df <- my_data[,c(covariates, c('OS','_PATIENT','OS_time','TCGA_ID','Sample_type'))]

  tmp.data <- my_df[which(my_df$TCGA_ID == sample),]

  if(dim(tmp.data)[1]==0) next

  covariates <- top.mir15[which(top.mir15$TCGA_ID == sample),]$X
  covariates <- 

  covariates <- lapply(covariates, function(x){
  return(str_replace_all(x, '-', '_'))
  })

  covariates <- unlist(covariates)

  tmp.data.normal <- tmp.data[which(tmp.data$Sample_type == 'Solid Tissue Normal'),]
  tmp.data.tumor <- tmp.data[which(tmp.data$Sample_type == 'Primary Tumor'),]

  if(dim(tmp.data.normal)[1] == 0 || dim(tmp.data.tumor)[1] == 0) next

  pr.data <- NULL
  roc.data <- NULL
  for(covariate in covariates){
      print(covariate)
      # covariate <- 'hsa_miR_424_5p'
      tmp.normal <- tmp.data.normal[,covariate,drop=TRUE]
      tmp.tumor <- tmp.data.tumor[,covariate,drop=TRUE]
      if(length(tmp.normal) == 0 || length(tmp.tumor) == 0) next

      skip_to_next <- FALSE
      p <- tryCatch({
          # print(m)    
          roc <- roc.curve(scores.class0 = tmp.normal, scores.class1 = tmp.tumor , curve = TRUE)
          pr <- pr.curve(scores.class0 = tmp.normal, scores.class1 = tmp.tumor , curve = TRUE)

          print(paste0(sample, ' ', covariate, ' ', 'ROC: ', roc$auc, ', PR: ', pr$auc.integral))

          if(roc$auc > 0.75){
            tmp.roc <- as.data.frame(roc$curve)
            tmp.pr <- as.data.frame(pr$curve)
            tmp.roc$AUC <- paste0(covariate, ' AUC: ', format(roc$auc,digits=3))
            tmp.pr$AUC <- paste0(covariate, ' AUC: ', format(pr$auc.integral,digits=3))
            tmp.roc$auc.value <- roc$auc
            tmp.pr$auc.value <- pr$auc.integral
            tmp.roc$miRNA <- covariate
            tmp.pr$miRNA <- covariate
            tmp.roc$TCGA_ID <- sample
            tmp.pr$TCGA_ID <- sample

            if(is.null(roc.data)){
              roc.data <- tmp.roc
              pr.data <- tmp.pr
            }else{
              roc.data <- rbind(roc.data, tmp.roc)
              pr.data <- rbind(pr.data, tmp.pr)
            }
          }
          
        },error = function(e){
          print(covariate)
          print('\n')
          skip_to_next <- TRUE
          }, 
        warning = function(cond) {
          message(cond)
          print(covariate)
          print('\n')
          skip_to_next <- TRUE
        })
      # print(p)
      if(skip_to_next == TRUE){next}
  }

  if(!is.null(roc.data)){
    png(paste('output/miR15/TCGA.',sample,".primary.tumor.vs.normal.ROC.sig.png", sep=""),height=4000, width=4000, res=800)
    p <- ggplot(roc.data[which(roc.data$auc.value > 0.8),],aes(x=V1,y=V2,color=AUC)) + 
          geom_line() +
          labs(x="1-Specificity",y="Sensitivity",title=paste0(sample, ' ROC Curve'),colour="AUC") + 
          # labs(x="Recall",y="Precision",title='ROC Curve',colour="AUC") + 
          theme(legend.position=c(0.6,0.5),legend.direction = "vertical", plot.title = element_text(hjust = 0.5, size = 14))
    print(p)
    dev.off()

    png(paste('output/miR15/TCGA.',sample,".primary.tumor.vs.normal.PR.sig.png", sep=""),height=4000, width=4000, res=800)
    p <- ggplot(pr.data,aes(x=V1,y=V2,color=AUC)) + 
          geom_line() +
          # labs(x="1-Specificity",y="Sensitivity",title='PR Curve',colour="AUC") + 
          labs(x="Recall",y="Precision",title=paste0(sample, ' PR Curve'),colour="AUC") + 
          theme(legend.position=c(0.4,0.3),legend.direction = "vertical", plot.title = element_text(hjust = 0.5, size = 14))
    print(p)
    dev.off()

    if(is.null(roc.data.all)){
      roc.data.all <- roc.data
      pr.data.all <- pr.data
    }else{
      roc.data.all <- rbind(roc.data.all, roc.data)
      pr.data.all <- rbind(pr.data.all, pr.data)
    }
  }

}


my_roc_plot <- function(sample){
      roc.data <- roc.data.all[which(roc.data.all$TCGA_ID == sample),]
      ggplot(roc.data[which(roc.data$auc.value > 0.8),],aes(x=V1,y=V2,color=AUC)) + 
          geom_line() +
          labs(x="1-Specificity",y="Sensitivity",title=paste0(sample, ' ROC Curve'),colour="AUC") + 
          # labs(x="Recall",y="Precision",title='ROC Curve',colour="AUC") + 
          theme(legend.position=c(0.6,0.5),legend.direction = "vertical", plot.title = element_text(hjust = 0.5, size = 14))
}

my_plots <- map(unique(roc.data.all$TCGA_ID), my_roc_plot)

png(paste('output/miR15/',"TCGA.all.primary.tumor.vs.normal.ROC.png", sep=""),height=16000, width=12000, res=900)
ggarrange(plotlist=my_plots, nrow=4, ncol=3)
dev.off()


my_pr_plot <- function(sample){
      pr.data <- pr.data.all[which(pr.data.all$TCGA_ID == sample),]
      ggplot(pr.data,aes(x=V1,y=V2,color=AUC)) + 
          geom_line() +
          # labs(x="1-Specificity",y="Sensitivity",title='PR Curve',colour="AUC") + 
          labs(x="Recall",y="Precision",title=paste0(sample, ' PR Curve'),colour="AUC") + 
          theme(legend.position=c(0.4,0.3),legend.direction = "vertical", plot.title = element_text(hjust = 0.5, size = 14))
}

my_pr_plots <- map(unique(pr.data.all$TCGA_ID), my_pr_plot)

png(paste('output/miR15/',"TCGA.all.primary.tumor.vs.normal.PR.png", sep=""),height=16000, width=12000, res=900)
ggarrange(plotlist=my_pr_plots, nrow=4, ncol=3)
dev.off()


#Figure ROC V2
cvfun2 <- function(split,...){
  mod <- glm(Type ~ miRNA,  data=analysis(split), family=binomial)
  fit <- predict(mod, newdata=assessment(split), type="response")
  data.frame(fit = fit, y = model.response(model.frame(formula(mod), data=assessment(split))))
}

pr.data.all <- NULL
roc.data.all <- NULL
for(sample in unique(pheno.df$TCGA_ID)){

  # sample <- 'TCGA_BRCA'
  print(sample)
  mirna.df.sub.t <- as.data.frame(t(mirna.df.sub))
  common.samples <- intersect(rownames(mirna.df.sub.t), pheno.df[which(pheno.df$TCGA_ID == sample),]$submitter_id.samples)

  common.samples <- intersect(common.samples, survival.df$sample)

  mirna.df.sub.t <- mirna.df.sub.t[match(common.samples, rownames(mirna.df.sub.t)),]
  mirna.df.sub.t <- cbind(mirna.df.sub.t, survival.df[match(common.samples, survival.df$sample),][,-1])

  if(dim(mirna.df.sub.t)[1]==0) next

  print('Finish combining sample')

  common.samples <- intersect(rownames(mirna.df.sub.t), pheno.df$`submitter_id.samples`)

  mirna.df.sub.t$Sample_type <- ''

  mirna.df.sub.t[match(common.samples, rownames(mirna.df.sub.t)),]$Sample_type <- pheno.df[match(common.samples, pheno.df$`submitter_id.samples`),]$`sample_type.samples`


  my_data <- mirna.df.sub.t
  names(my_data)[names(my_data) == "OS.time"] <- "OS_time"

  covariates <- unique(mir.table[which(mir.table$miRNA %in% mir.15),]$miRNA_ID)

  my_df <- my_data[,c(covariates, c('OS','_PATIENT','OS_time','TCGA_ID','Sample_type'))]

  tmp.data <- my_df[which(my_df$TCGA_ID == sample),]

  if(dim(tmp.data)[1]==0) next

  covariates <- top.mir15[which(top.mir15$TCGA_ID == sample),]$X

  covariates <- lapply(covariates, function(x){
  return(str_replace_all(x, '-', '_'))
  })

  covariates <- unlist(covariates)

  tmp.data.normal <- tmp.data[which(tmp.data$Sample_type == 'Solid Tissue Normal'),]
  tmp.data.tumor <- tmp.data[which(tmp.data$Sample_type == 'Primary Tumor'),]

  if(dim(tmp.data.normal)[1] == 0 || dim(tmp.data.tumor)[1] == 0) next

  roc.data <- NULL
  for(covariate in unique(covariates)){
      print(covariate)
      # covariate <- 'hsa_miR_195_5p'
      tmp.normal <- tmp.data.normal[,covariate,drop=TRUE]
      tmp.tumor <- tmp.data.tumor[,covariate,drop=TRUE]
      if(length(tmp.normal) == 0 || length(tmp.tumor) == 0) next

      model.data <- data.frame('miRNA'=c(tmp.normal, tmp.tumor), 'Type'=c(rep(0,length(tmp.normal)),rep(1,length(tmp.tumor))))

      print(dim(model.data))
      cv_out_plot <- vfold_cv(model.data, v=10, repeats = 5, strata=Type) %>% 
        mutate(fit = map(splits, cvfun2)) %>% 
        unnest(fit) %>% 
        group_by(id) %>% 
        summarise(sens = roc(y, fit, plot=FALSE)$sensitivities, 
                 spec = roc(y, fit, plot=FALSE)$specificities, 
                 auc = rep(roc(y, fit, plot=FALSE)$auc, length(sens)),
                 reg = rep('', length(sens)),
                 obs = 1:length(sens))

      ave <- cv_out_plot %>% 
        ungroup %>% 
        group_by(obs) %>% 
        summarise(sens = mean(sens), 
                  spec = mean(spec), 
                  auc = mean(auc),
                  reg = '',
                  id = "Average")

      cv_out_plot <- bind_rows(cv_out_plot, ave) %>% 
        mutate(col = factor(ifelse(id == "Average", "Average", "Individual"), 
                            levels=c("Individual", "Average")))

      # auc_avg_mean <- mean(cv_out_plot[which(cv_out_plot$col == 'Average'),]$auc)
      # auc_individual_mean <- mean(cv_out_plot[which(cv_out_plot$col == 'Individual'),]$auc)

      # cv_out_plot[which(cv_out_plot$col == 'Average'),]$avg <- auc_avg_mean
      # cv_out_plot[which(cv_out_plot$col == 'Individual'),]$avg <- auc_individual_mean

      cv_out_plot <- cv_out_plot %>% 
        mutate(TCGA_ID = sample, Group = covariate, AUC = paste0(covariate, ' AUC: ', format(auc,digits=3)))

      if(is.null(roc.data)){
        roc.data <- cv_out_plot
      }else{
        roc.data <- bind_rows(roc.data, cv_out_plot)
      }
  }

  if(is.null(roc.data)) next


  if(is.null(covariates)) next
  roc.data.multi <- NULL

  panel <- unique(roc.data[which(roc.data$auc > 0.8 & roc.data$col == 'Average'),]$Group)
  # panel <- covariates
  cvfun <- function(split,...){ 
    mod <- glm(as.formula(paste("Type ~ ", paste(panel, collapse= "+"))),  data=analysis(split), family=binomial)
    fit <- predict(mod, newdata=assessment(split), type="response")
    df <- as.data.frame(mod$coefficients)
    df$`mod$coefficients` <- format(round(df$`mod$coefficients`, 2), nsmall = 2)
    rownames(df) <- c(1, panel)
    df$reg <- paste0('(',df$`mod$coefficients`, ')', '*', rownames(df))
    # print(df$reg)
    reg <- paste(df$reg, collapse="+")
    print(reg)
    data.frame(fit = fit, y = model.response(model.frame(formula(mod), data=assessment(split))), reg=rep(reg, length(fit)))
  }

  if(dim(tmp.data.normal[,panel,drop=FALSE])[1]==0 || dim(tmp.data.tumor[,panel,drop=FALSE])[1] == 0) next
  model.data <- rbind(tmp.data.normal[,panel,drop=FALSE], tmp.data.tumor[,panel,drop=FALSE])
  model.data$Type <- c(rep(0, dim(tmp.data.normal)[1]), rep(1,dim(tmp.data.tumor)[1]))
  # model.data <- data.frame('miRNA'=c(tmp.normal, tmp.tumor), 'Type'=c(rep(0,length(tmp.normal)),rep(1,length(tmp.tumor))))

  cv_out_plot <- vfold_cv(model.data, v=10, repeats = 5, strata=Type) %>% 
        mutate(fit = map(splits, cvfun)) %>% 
        unnest(fit) %>% 
        group_by(id) %>% 
        summarise(sens = roc(y, fit, plot=FALSE)$sensitivities, 
                 spec = roc(y, fit, plot=FALSE)$specificities, 
                 auc = rep(roc(y, fit, plot=FALSE)$auc, length(sens)),
                 reg = rep(reg[1],length(sens)),
                 obs = 1:length(sens))

  ave <- cv_out_plot %>% 
        ungroup %>% 
        group_by(obs) %>% 
        summarise(sens = mean(sens), 
                  spec = mean(spec), 
                  auc = mean(auc),
                  reg = '',
                  id = "Average")

  cv_out_plot <- bind_rows(cv_out_plot, ave) %>% 
        mutate(col = factor(ifelse(id == "Average", "Average", "Individual"), 
                            levels=c("Individual", "Average")))

  # auc_avg_mean <- mean(cv_out_plot[which(cv_out_plot$col == 'Average'),]$auc)
  # auc_individual_mean <- mean(cv_out_plot[which(cv_out_plot$col == 'Individual'),]$auc)

  # cv_out_plot[which(cv_out_plot$col == 'Average'),]$avg <- auc_avg_mean
  # cv_out_plot[which(cv_out_plot$col == 'Individual'),]$avg <- auc_individual_mean

  cv_out_plot <- cv_out_plot %>% 
        mutate(TCGA_ID = sample, Group = 'Panel', AUC = paste0('Panel', ' AUC: ', format(auc,digits=3)))

  roc.data.multi <- cv_out_plot



  if(!is.null(roc.data)){
    if(!is.null(roc.data.multi)){
      roc.data <- rbind(roc.data ,roc.data.multi)
    }

    if(is.null(roc.data.all)){
      roc.data.all <- roc.data
      # pr.data.all <- pr.data
    }else{
      roc.data.all <- rbind(roc.data.all, roc.data)
      # pr.data.all <- rbind(pr.data.all, pr.data)
    }
  }

}

tmp <- roc.data.all[which(grepl('\\+',roc.data.all$reg)),]
tmp <- tmp[,c(5,7,8,4)] %>% group_by(TCGA_ID) %>% filter(auc == max(auc))
tmp <- tmp %>% group_by(TCGA_ID) %>% arrange(desc(auc), .by_group = TRUE) %>%
    filter(row_number() == 1)
write.csv(tmp, 'output/miR15/miR15.glm.all.lt.0.8.csv')

my_roc_plot <- function(sample){
      roc.data <- roc.data.all[which(roc.data.all$TCGA_ID == sample),]
      # ggplot(roc.data[which(roc.data$col == 'Average'),],aes(x=1-spec,y=sens,color=Group)) + 
      ggplot(roc.data[which(roc.data$auc > 0.8 & roc.data$col == 'Average'),],aes(x=1-spec,y=sens,color=Group)) + 
          geom_line() +
          labs(x="1-Specificity",y="Sensitivity",title=paste0(sample, ' ROC Curve'), colour="Group") + 
          # labs(x="Recall",y="Precision",title='ROC Curve',colour="AUC") + 
          theme(legend.position=c(0.6,0.5),legend.direction = "vertical", plot.title = element_text(hjust = 0.5, size = 14))
}

my_plots <- map(unique(roc.data.all$TCGA_ID), my_roc_plot)

png(paste('output/miR15/',"TCGA.all.primary.tumor.vs.normal.ROC.v3.png", sep=""),height=16000, width=12000, res=900)
ggarrange(plotlist=my_plots, nrow=5, ncol=3)
dev.off()



#Figure 

corr.filtered <- NULL
for(i in 1:dim(res_lt4.mir)[1]){
  mir <- res_lt4.mir[i,3]
  mir <- str_replace_all(mir, '-', '_')
  target <- res_lt4.mir[i,4]

  tmp <- corr.data[which(corr.data$X == mir & corr.data$Y == target),]
  if(dim(tmp)[1] > 0){
    if(is.null(corr.filtered)){
      corr.filtered <- tmp
    }else{
      corr.filtered <- rbind(corr.filtered, tmp)
    }
  }
}


corr.final <- NULL
for(i in 1:dim(mrna.top.sig)[1]){
  gene <- mrna.top.sig[i,1]

  sample <- mrna.top.sig[i,7]

  tmp <- corr.filtered[which(corr.filtered$TCGA_ID == sample & corr.filtered$Y == gene),]
  if(dim(tmp)[1] > 0){
    if(is.null(corr.final)){
      corr.final <- tmp
    }else{
      corr.final <- rbind(corr.final, tmp)
    }
  }
}


corr.result <- NULL
for(i in 1:dim(top.genes.mir15)[1]){
  mir <- top.genes.mir15[i,1]
  mir <- str_replace_all(mir, '-', '_')

  sample <- top.genes.mir15[i,7]

  tmp <- corr.final[which(corr.final$TCGA_ID == sample & corr.final$X == mir),]
  if(dim(tmp)[1] > 0){
    if(is.null(corr.result)){
      corr.result <- tmp
    }else{
      corr.result <- rbind(corr.result, tmp)
    }
  }
}


corr.result$miRNA_logFC <- 0

for(i in 1:dim(top.genes.mir15)[1]){
  mir <- top.genes.mir15[i,1]
  mir <- str_replace_all(mir, '-', '_')

  sample <- top.genes.mir15[i,7]

  if(dim(corr.result[which(corr.result$TCGA_ID == sample & corr.result$X == mir),])[1]==0) next
  corr.result[which(corr.result$TCGA_ID == sample & corr.result$X == mir),]$miRNA_logFC <- top.genes.mir15[i,2]
  
}

# top.genes.mir15 <- mrna.top.sig

corr.result$Gene_logFC <- 0

for(i in 1:dim(mrna.top.sig)[1]){
  gene <- mrna.top.sig[i,1]

  sample <- mrna.top.sig[i,7]

  if(dim(corr.result[which(corr.result$TCGA_ID == sample & corr.result$Y == gene),])[1]==0) next

  corr.result[which(corr.result$TCGA_ID == sample & corr.result$Y == gene),]$Gene_logFC <- mrna.top.sig[i,2]
  
}

filter.set <- data.frame('TCGA_ID'=c('TCGA_KIRC', 'TCGA_KIRP', 'TCGA_KIRP','TCGA_LUAD','TCGA_LUAD'), 'Gene'=c('TFAP2A','TFAP2A','E2F7','CCNE1','E2F7'))

corr.result.tmp <- NULL
for(i in 1:dim(filter.set)[1]){
  gene <- filter.set[i,2]

  sample <- filter.set[i,1]

  tmp <- corr.result[which(corr.result$TCGA_ID == sample & corr.result$Y == gene),]

  if(dim(tmp)[1] > 0){
    if(is.null(corr.result.tmp)){
      corr.result.tmp <- tmp
    }else{
      corr.result.tmp <- rbind(corr.result.tmp, tmp)
    }
  }
  
}

corr.result <- corr.result.tmp

library(ggplot2)
library(dplyr)
library(tidyr)
# complete(All_Tissues_BP_Head, Tissue, GO.ID)

my_corr_plot <- function(sample){
      cur_corr.data <- corr.final[which(corr.final$TCGA_ID == sample),]
      ggplot(complete(cur_corr.data,X,Y), aes(X, Y, fill= r)) + 
        geom_tile() +
        geom_text(aes(label = round(r, 2))) +
        scale_fill_gradient2(low = "blue", mid="green", high = "yellow",na.value="white") +
        theme(legend.position='none',text = element_text(size=12), axis.text.x = element_text(angle = 45, vjust = 0.7, hjust=0.8)) +
        labs(title = sample) +
        xlab("") +
        ylab("")
}

my_plots <- map(unique(corr.final$TCGA_ID), my_corr_plot)

png(paste('output/miR15/',"TCGA.all.primary.tumor.vs.normal.corr.V2.png", sep=""),height=12000, width=16000, res=800)
ggarrange(plotlist=my_plots, nrow=3, ncol=4)
dev.off()

