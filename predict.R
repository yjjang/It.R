library(keras)
library(tidyverse)
library(GSVA)
library(GSVAdata)
library(GSEABase)
library(estimate)
library(MCPcounter)
library(caret)
library(pROC)
#devtools::install_github('oganm/MSigDB')
library(MSigDB)
library(randomForest)


#### SMC Lung Cancer data loading ##############################################

## TPM
lung.tpm <- read.csv("Data/SMC_lungCancer/New_RSEM_TPM_149_lung.txt",
                     stringsAsFactors = F, row.names = 1, sep = '\t')

## Clinical data
lung.cl <- read.csv("Data/SMC_lungCancer/RNA_list_to_EWHA_031419_process.txt",
                    stringsAsFactors = F, row.names = NULL, sep = '\t')
lung.cl <- lung.cl[!duplicated(lung.cl[, 1]), ]
lung.cl <- merge(data.frame(Sample = colnames(lung.tpm)), lung.cl, by ="Sample", all.x = T, all.y = F)
rownames(lung.cl) <- lung.cl[, 1]; lung.cl <- lung.cl[, -1]
table(lung.cl$Response)
lung.cl$binaryResponse <- NA
lung.cl$binaryResponse[grep('CR|PR', lung.cl$Response)] <- 'CR/PR'
lung.cl$binaryResponse[grep('SD|PD', lung.cl$Response)] <- 'SD/PD'
table(lung.cl$binaryResponse)

## Order of samples
lung.tpm <- lung.tpm[, rownames(lung.cl)]

## TruSeq V2
lung.cl.v2 <- lung.cl[lung.cl$Kit == "TruSeq RNA v2", ]
lung.tpm.v2 <- lung.tpm[, rownames(lung.cl.v2)]
lung.cl.v2.ad <- lung.cl.v2[lung.cl.v2$Pathology == "ADC", ]
lung.tpm.v2.ad <- lung.tpm[, rownames(lung.cl.v2.ad)]
lung.cl.v2.sc <- lung.cl.v2[lung.cl.v2$Pathology == "SQCC", ]
lung.tpm.v2.sc <- lung.tpm[, rownames(lung.cl.v2.sc)]

## TruSeq Access
lung.cl.ac <- lung.cl[lung.cl$Kit == "TruSeq RNA Access", ]
lung.tpm.ac <- lung.tpm[, rownames(lung.cl.ac)]
lung.cl.ac.ad <- lung.cl.ac[lung.cl.ac$Pathology == "ADC", ]
lung.tpm.ac.ad <- lung.tpm[, rownames(lung.cl.ac.ad)]
lung.cl.ac.sc <- lung.cl.ac[lung.cl.ac$Pathology == "SQCC", ]
lung.tpm.ac.sc <- lung.tpm[, rownames(lung.cl.ac.sc)]


#### ssGSEA scores for MSigDB, known genesets ##################################

## Calculate ssGSEA scores for MSigDB signatures using the given tpm values
ssgsea.msig <- function(tpm.dat) {
  ssgsea.dat <- rbind(gsva(as.matrix(log(tpm.dat+1)), MSigDB$HALLMARK,
                           min.sz = 2, max.sz = 500, method = "ssgsea"),
                      gsva(as.matrix(log(tpm.dat+1)), MSigDB$C1_POSITIONAL,
                           min.sz = 2, max.sz = 500, method = "ssgsea"),
                      gsva(as.matrix(log(tpm.dat+1)), MSigDB$C2_CURATED,
                           min.sz = 2, max.sz = 500, method = "ssgsea"),
                      gsva(as.matrix(log(tpm.dat+1)), MSigDB$C3_MOTIF,
                           min.sz = 2, max.sz = 500, method = "ssgsea"),
                      gsva(as.matrix(log(tpm.dat+1)), MSigDB$C4_COMPUTATIONAL,
                           min.sz = 2, max.sz = 500, method = "ssgsea"),
                      gsva(as.matrix(log(tpm.dat+1)), MSigDB$C5_GENE_ONTOLOGY,
                           min.sz = 2, max.sz = 500, method = "ssgsea"),
                      gsva(as.matrix(log(tpm.dat+1)), MSigDB$C6_ONCOGENIC_SIGNATURES,
                           min.sz = 2, max.sz = 500, method = "ssgsea"),
                      gsva(as.matrix(log(tpm.dat+1)), MSigDB$C7_IMMUNOLOGIC_SIGNATURES,
                           min.sz = 2, max.sz = 500, method = "ssgsea"))
  ssgsea.dat <- t(scale(t(ssgsea.dat)))
  return(ssgsea.dat)
}

## Calculate P-vaule & FDR
ssgsea.fdr <- function(ssgsea.dat, response.dat) {
  for(i in 1:nrow(ssgsea.dat)) {
    x <- t.test(as.numeric(ssgsea.dat[i, rownames(response.dat)[response.dat$binaryResponse == "CR/PR"]]),
                as.numeric(ssgsea.dat[i, rownames(response.dat)[response.dat$binaryResponse == "SD/PD"]]),
                paired = F, var.equal = T, conf.level = 0.9)$p.value
    if(i == 1) {ptest <- x} else {ptest <- c(ptest, x)}}
  names(ptest) <- rownames(ssgsea.dat)
  ftest <- p.adjust(ptest, method = "fdr")
  return(ftest)
}

## Calculate ssGSEA scores for MSigDB genesets (TruSeq V2)
#lung.msig.v2 <- ssgsea.msig(lung.tpm.v2)
#ftest <- ssgsea.fdr(lung.msig.v2, lung.cl.v2)
#head(sort(ftest))
## Filter by FDR=?
##nrow(lung.msig.v2[ftest<=0.3, ])
##lung.msig.v2 <- lung.msig.v2[ftest<=0.182, ]
## Filter top 50
#lung.msig.v2 <- lung.msig.v2[order(ftest)[1:50], ]
#write.csv(cbind(rownames(lung.msig.v2), lung.msig.v2), row.names=FALSE, "Data/MsigDB.lung.v2.Top50.txt")
lung.msig.v2 <- read.csv("Data/MsigDB.lung.v2.Top50.txt", stringsAsFactors = F, row.names = 1)

## Calculate ssGSEA scores for MSigDB genesets (TruSeq V2: Adenocarcinoma)
lung.msig.v2.ad <- ssgsea.msig(lung.tpm.v2.ad)
ftest <- ssgsea.fdr(lung.msig.v2.ad, lung.cl.v2.ad)
head(sort(ftest))
## Filter top 50
lung.msig.v2.ad <- lung.msig.v2.ad[order(ftest)[1:50], ]
write.csv(cbind(rownames(lung.msig.v2.ad), lung.msig.v2.ad), row.names=FALSE, "Data/MsigDB.lung.v2.ad.Top50.txt")
lung.msig.v2.ad <- read.csv("Data/MsigDB.lung.v2.ad.Top50.txt", stringsAsFactors = F, row.names = 1)

## Known Genesets (T cell inflamed GEP, IPRES, F-TBRS)
IPRESgs <- read.delim("Data/Known.Genesets.genes",
                      header=T, sep="\t", stringsAsFactors=FALSE, check.names=F, na.strings="")
# GSVA scores for known genesets
lung.gs.v2 <- gsva(as.matrix(log(lung.tpm.v2+1)), IPRESgs, min.sz=2, max.sz=500, method="ssgsea")
lung.gs.v2 <- t(scale(t(lung.gs.v2)))
lung.gs.v2.ad <- gsva(as.matrix(log(lung.tpm.v2.ad+1)), IPRESgs, min.sz=2, max.sz=500, method="ssgsea")
lung.gs.v2.ad <- t(scale(t(lung.gs.v2.ad)))


#### ESTIMATE scores ###########################################################

# V2
out.file <- tempfile(pattern="estimate", fileext=".gct")
outputGCT(as.data.frame(lung.tpm.v2), out.file)
in.file <- out.file
out.file <- tempfile(pattern="estimate", fileext=".txt")
estimateScore(in.file, out.file, platform=c("affymetrix", "agilent", "illumina")[3])
lung.est.v2 <- read.csv(out.file, stringsAsFactors=FALSE, row.names=1, sep='\t', skip=2)[, -1]
lung.est.v2 <- t(scale(t(lung.est.v2)))

# V2 & Adenocarcinoma
out.file <- tempfile(pattern="estimate", fileext=".gct")
outputGCT(as.data.frame(lung.tpm.v2.ad), out.file)
in.file <- out.file
out.file <- tempfile(pattern="estimate", fileext=".txt")
estimateScore(in.file, out.file, platform=c("affymetrix", "agilent", "illumina")[3])
lung.est.v2.ad <- read.csv(out.file, stringsAsFactors=FALSE, row.names=1, sep='\t', skip=2)[, -1]
lung.est.v2.ad <- t(scale(t(lung.est.v2.ad)))

#### MCPCounter scores #########################################################

# V2
lung.mcp.v2 <- MCPcounter.estimate(as.matrix(log2(lung.tpm.v2+1)), featuresType="HUGO_symbols")
lung.mcp.v2 <- t(scale(t(lung.mcp.v2)))

# V2 & Lungadenocarcinoma
lung.mcp.v2.ad <- MCPcounter.estimate(as.matrix(log2(lung.tpm.v2.ad+1)), featuresType="HUGO_symbols")
lung.mcp.v2.ad <- t(scale(t(lung.mcp.v2.ad)))


#### Tumor Mutation Burden #####################################################

## Tumor mutation burden for TruSeq V2
lung.tmb.v2 <- lung.cl.v2$TMB
lung.tmb.v2 <- as.vector(scale(as.numeric(lung.tmb.v2))); names(lung.tmb.v2) <- rownames(lung.cl.v2)
# Remove patients who do not have TMB score (NA);
lung.cl.v2 <- lung.cl.v2[!is.na(lung.tmb.v2), ]
lung.msig.v2 <- lung.msig.v2[, !is.na(lung.tmb.v2)]
lung.gs.v2 <- lung.gs.v2[, !is.na(lung.tmb.v2)]
lung.est.v2 <- lung.est.v2[, !is.na(lung.tmb.v2)]
lung.mcp.v2 <- lung.mcp.v2[, !is.na(lung.tmb.v2)]
lung.tmb.v2 <- lung.tmb.v2[!is.na(lung.tmb.v2)]

## Tumor mutation burden for TruSeq V2 & Lungadenocarcinoma
lung.tmb.v2.ad <- lung.cl.v2.ad$TMB
lung.tmb.v2.ad <- as.vector(scale(as.numeric(lung.tmb.v2.ad))); names(lung.tmb.v2.ad) <- rownames(lung.cl.v2.ad)
# Remove patients who do not have TMB score (NA);
lung.cl.v2.ad <- lung.cl.v2.ad[!is.na(lung.tmb.v2.ad), ]
lung.msig.v2.ad <- lung.msig.v2.ad[, !is.na(lung.tmb.v2.ad)]
lung.gs.v2.ad <- lung.gs.v2.ad[, !is.na(lung.tmb.v2.ad)]
lung.est.v2.ad <- lung.est.v2.ad[, !is.na(lung.tmb.v2.ad)]
lung.mcp.v2.ad <- lung.mcp.v2.ad[, !is.na(lung.tmb.v2.ad)]
lung.tmb.v2.ad <- lung.tmb.v2.ad[!is.na(lung.tmb.v2.ad)]


#### Save Signature Categories (All) ###########################################

sigcat.lung.v2.all <- list(TMB = c('TMB'),
                           MSigDB = rownames(lung.msig.v2),
                           KnownGenesets = rownames(lung.gs.v2),
                           ESTIMATE = rownames(lung.est.v2),
                           MCPCounter = rownames(lung.mcp.v2))
sigcat.lung.v2.all <- plyr::ldply(sigcat.lung.v2.all, cbind)[, 2:1]
colnames(sigcat.lung.v2.all) <- c('Signature', 'Category')
rownames(sigcat.lung.v2.all) <- sigcat.lung.v2.all$Signature
write.table(sigcat.lung.v2.all, 'Out/Lung.v2.allSig.txt', row.names = F)

sigcat.lung.v2.ad.all <- list(TMB = c('TMB'),
                              MSigDB = rownames(lung.msig.v2.ad),
                              KnownGenesets = rownames(lung.gs.v2.ad),
                              ESTIMATE = rownames(lung.est.v2.ad),
                              MCPCounter = rownames(lung.mcp.v2.ad))
sigcat.lung.v2.ad.all <- plyr::ldply(sigcat.lung.v2.ad.all, cbind)[, 2:1]
colnames(sigcat.lung.v2.ad.all) <- c('Signature', 'Category')
rownames(sigcat.lung.v2.ad.all) <- sigcat.lung.v2.ad.all$Signature
write.table(sigcat.lung.v2.ad.all, 'Out/Lung.v2.ad.allSig.txt', row.names = F)


#### Random sampling and feature selection #####################################

## Calculate AUCes of each features for discriminating responses
calc.auces <- function(train.dat, response.dat, N) {
  bootstrap_tags <- createResample(y=response.dat$Response, times=N)
  for (M in 1:N) {
    etrain2 <- train.dat[bootstrap_tags[[M]], ]
    ctrain2 <- response.dat[rownames(etrain2), ]
    if (sum(!duplicated(ctrain2$binaryResponse)) < 2) { next; }
    # Calculate AUC by each features
    fred <- c()
    for(i in 1:ncol(etrain2)){
      predicted <- as.numeric(etrain2[, i]); names(predicted) <- rownames(etrain2)
      label <- ctrain2$binaryResponse; label[label=="CR/PR"] <- 1; label[label=="SD/PD"] <- 0
      names(label) <- names(predicted)
      fred <- c(fred, pROC::roc(label, predicted)$auc)
    }; names(fred) <- colnames(etrain2)

    if(M==1){
      fred_re <- data.frame(fred, stringsAsFactors=FALSE); colnames(fred_re) <- M
    }else{
      fred_re <- cbind(fred_re, fred)
      #colnames(fred_re)[M] <- M
    }
  }
  return(fred_re)
}

## TruSeq V2
test.lung.cl.v2 <- lung.cl.v2
test.lung.v2 <- cbind(as.matrix(lung.tmb.v2),
                      t(as.matrix(lung.est.v2)),
                      t(as.matrix(lung.gs.v2)),
                      t(as.matrix(lung.msig.v2)),
                      t(as.matrix(lung.mcp.v2)));
colnames(test.lung.v2)[1] <- 'TMB'
ncol(test.lung.v2)
test.lung.v2 <- test.lung.v2[rownames(test.lung.cl.v2), ]
# Divide Whole data into traing and validation data set
#tag.lung.v2 <- createDataPartition(y=test.lung.cl.v2$Response, p=0.7, list=FALSE)
#write.csv(tag.lung.v2, 'Data/Random.Sampling.tag.lung.v2', row.names = F)
tag.lung.v2 <- as.matrix(read.csv("Data/Random.Sampling.tag.lung.v2", stringsAsFactors=FALSE))
# Training data
train.lung.v2 <- test.lung.v2[tag.lung.v2, ]
train.lung.cl.v2 <- test.lung.cl.v2[rownames(train.lung.v2), ]
nrow(train.lung.v2)
table(train.lung.cl.v2$binaryResponse)
# Validation data
valid.lung.v2 <- test.lung.v2[-tag.lung.v2, ]
valid.lung.cl.v2 <- test.lung.cl.v2[rownames(valid.lung.v2), ]
nrow(valid.lung.v2)
table(valid.lung.cl.v2$binaryResponse)
# Random sampling 10000 times and find constantly singnificant signatures
N <- 10000
#fred_re <- calc.auces(train.lung.v2, train.lung.cl.v2, N)
#fred_fi.lung.v2 <- fred_re > 0.6
#fred_fi.lung.v2 <- rowSums(fred_fi.lung.v2)/N
#write.csv(cbind("Name"=names(fred_fi.lung.v2), fred_fi.lung.v2), row.names=FALSE, "Data/Random.Sampling.auces.lung.v2")
# Find constantly singnificant signature that selected over 7500 count
fred_fi.lung.v2 <- read.csv("Data/Random.Sampling.auces.lung.v2", stringsAsFactors=FALSE)$fred_fi
summary(fred_fi.lung.v2)
names(fred_fi.lung.v2) <- colnames(train.lung.v2)
cutoff.lung.v2 <- 0.7
sum(fred_fi.lung.v2 > cutoff.lung.v2)

## V2 & Aenocarcinoma
test.lung.cl.v2.ad <- lung.cl.v2.ad
test.lung.v2.ad <- cbind(as.matrix(lung.tmb.v2.ad),
                         t(as.matrix(lung.est.v2.ad)),
                         t(as.matrix(lung.gs.v2.ad)),
                         t(as.matrix(lung.msig.v2.ad)),
                         t(as.matrix(lung.mcp.v2.ad)));
colnames(test.lung.v2.ad)[1] <- 'TMB'
ncol(test.lung.v2.ad)
test.lung.v2.ad <- test.lung.v2.ad[rownames(test.lung.cl.v2.ad), ]
# Divide Whole data into traing and validation data set
#tag.lung.v2.ad <- createDataPartition(y=test.lung.cl.v2.ad$Response, p=0.7, list=FALSE)
#write.csv(tag.lung.v2.ad, 'Data/Random.Sampling.tag.lung.v2.ad', row.names = F)
tag.lung.v2.ad <- as.matrix(read.csv("Data/Random.Sampling.tag.lung.v2.ad", stringsAsFactors=FALSE))
# Training data
train.lung.v2.ad <- test.lung.v2.ad[tag.lung.v2.ad, ]
train.lung.cl.v2.ad <- test.lung.cl.v2.ad[rownames(train.lung.v2.ad), ]
table(train.lung.cl.v2.ad$binaryResponse)
# Validation data
valid.lung.v2.ad <- test.lung.v2.ad[-tag.lung.v2.ad, ]
valid.lung.cl.v2.ad <- test.lung.cl.v2.ad[rownames(valid.lung.v2.ad), ]
table(valid.lung.cl.v2.ad$binaryResponse)
# Random sampling 10000 times and find constantly singnificant signatures
N <- 10000
fred_re <- calc.auces(train.lung.v2.ad, train.lung.cl.v2.ad, N)
fred_fi.lung.v2.ad <- fred_re > 0.6
fred_fi.lung.v2.ad <- rowSums(fred_fi.lung.v2.ad)/N
write.csv(cbind("Name"=names(fred_fi.lung.v2.ad), fred_fi.lung.v2.ad), row.names=FALSE, "Data/Random.Sampling.auces.lung.v2.ad")
## Find constantly singnificant signature that selected over 7500 count
fred_fi.lung.v2.ad <- read.csv("Data/Random.Sampling.auces.lung.v2.ad", stringsAsFactors=FALSE)$fred_fi
summary(fred_fi.lung.v2.ad)
names(fred_fi.lung.v2.ad) <- colnames(train.lung.v2.ad)
cutoff.lung.v2.ad <- 0.7
sum(fred_fi.lung.v2.ad > cutoff.lung.v2.ad)


### Plot: constantly singnificant signatures ###################################

signatures.plot <- function(fred_fi.it, cutoff.it,
                            nrow.est, nrow.gs, nrow.msig, nrow.mcp) {
  par(mar = c(6,25,4,1))
  col <- rep("gray", length(fred_fi.it))
  col[fred_fi.it >= cutoff.it] <- "limegreen"
  barx <- barplot(as.matrix(fred_fi.it), xlim=c(0, 1.2), xaxt="n", horiz=TRUE,
                  cex.main=1.8, border=FALSE, beside=TRUE, cex.names=1.5, cex.axis=1.5, cex.lab=1.5, col=col)
  mtext("Discriminating Features (Selected by AUC > 0.6)", side=1, line=3, cex=1.5, at=0.15)
  mtext("Selected Count", side=2, line=3, cex=1, las=2, at=length(fred_fi.it)+5) #1=bottom, 2=left, 3=top, 4=right
  axis(side=3, at=seq(0, 1, length.out=11), labels=seq(0, N, length.out=11), col="black", cex.axis=1, las=1)
  axis(side=2, at=seq(1, length(fred_fi.it), length.out=length(fred_fi.it))+0.5, labels=NA,
       col="black", cex.axis=1.3, las=2, srt=45)
  text(par("usr")[1]-0.02, barx, labels=names(fred_fi.it), adj=1, srt=0, cex=0.6, xpd=TRUE,
       col=rep(c("black", "forestgreen", "indianred1", "dodgerblue1", "goldenrod2"),
               c(1, nrow.est, nrow.gs, nrow.msig, nrow.mcp)))
  text(0.1, (1:length(fred_fi.it))+0.5, col="black", cex=0.4, labels=fred_fi.it*N, srt=0)
  abline(v=cutoff.it, lty=3, lwd=2, col="red")
  #legend("topleft", inset=c(0.7, 0.7), title="Data Source",
  #       fill=c("black", "forestgreen", "indianred1", "dodgerblue1", "goldenrod2"),
  #       text.col=c("black", "forestgreen", "indianred1", "dodgerblue1", "goldenrod2"),
  #       legend=c("TMB", "ESTIMATE", "Known Signatures", "MsigDB", "MCP-counter"),
  #       horiz=FALSE, border=NA, bty='n', cex=1, y.intersp = 0.6, title.adj = 0.23)
  #legend("topleft", inset=c(0.35, 0.98), sprintf("Cut-off = %s", cutoff.it*N), col="red", text.col="red",
  #       cex=1, bty="n", horiz=FALSE)
  par(mar = c(1,1,1,1))
}

signatures.plot(
  fred_fi.lung.v2, cutoff.lung.v2,
  nrow(lung.est.v2),
  nrow(lung.gs.v2),
  nrow(lung.msig.v2),
  nrow(lung.mcp.v2))

signatures.plot(
  fred_fi.lung.v2.ad, cutoff.lung.v2.ad,
  nrow(lung.est.v2.ad),
  nrow(lung.gs.v2.ad),
  nrow(lung.msig.v2.ad),
  nrow(lung.mcp.v2.ad))


#### Training and test data for model learning #################################

# Learning and evaluation data
x.train <- train.lung.v2[, fred_fi.lung.v2 > cutoff.lung.v2]
y.train <- as.matrix(as.numeric(train.lung.cl.v2[rownames(x.train), ]$binaryResponse=="CR/PR"))
#x.train <- train.lung.v2.ad[, fred_fi.lung.v2.ad > cutoff.lung.v2.ad]
#y.train <- as.matrix(as.numeric(train.lung.cl.v2.ad[rownames(x.train), ]$binaryResponse=="CR/PR"))
x.test <- valid.lung.v2[, fred_fi.lung.v2 > cutoff.lung.v2]
y.test <- as.matrix(as.numeric(valid.lung.cl.v2[rownames(x.test), ]$binaryResponse=="CR/PR"))
#x.test <- valid.lung.v2.ad[, fred_fi.lung.v2.ad > cutoff.lung.v2.ad]
#y.test <- as.matrix(as.numeric(valid.lung.cl.v2.ad[rownames(x.test), ]$binaryResponse=="CR/PR"))
nfeatures <- ncol(x.train)

# Prediction data
x.test <- test.lung.v2[, fred_fi.lung.v2 > cutoff.lung.v2]
y.test <- as.matrix(as.numeric(test.lung.cl.v2[rownames(x.test), ]$binaryResponse=="CR/PR"))
#x.test <- test.lung.v2.ad[, fred_fi.lung.v2.ad > cutoff.lung.v2.ad]
#y.test <- as.matrix(as.numeric(test.lung.cl.v2.ad[rownames(x.test), ]$binaryResponse=="CR/PR"))


#### Deep learning:: Sequential model ##########################################

model.seq <- keras_model_sequential() %>%
  layer_dense(units  = nfeatures, activation = 'relu',
              #kernel_initializer = 'he_normal',
              #kernel_regularizer = regularizer_l1_l2(l1 = 0.001, l2 = 0.001),
              input_shape = nfeatures) %>%
  layer_dropout(rate = 0.3) %>%
  layer_dense(units  = 20, activation  = 'relu',
              #kernel_initializer = 'he_normal',
              #kernel_regularizer = regularizer_l1_l2(l1 = 0.001, l2 = 0.001)
  ) %>%
  layer_dropout(rate = 0.2) %>%
  layer_dense(units  = 1, activation   = 'sigmoid')
#summary(model.seq)

## Compile
model.seq %>% compile(
  loss      = 'binary_crossentropy',
  optimizer = optimizer_adam(),
  metrics   = c('accuracy')
)

## Training and evaluation
history.seq = model.seq %>% fit(
  x.train, y.train,
  epochs = 150,
  #batch_size = 32,
  validation_split = 0.2
)
#plot(history.seq)

perf = model.seq %>% evaluate(x.test, y.test)

## ROC
y.pred = as.vector(model.seq %>% predict(x.test))
y.real = as.vector(y.test)
roc <- pROC::roc(y.real, y.pred)
roc$auc
plot.roc(roc, add=F, col="purple", lwd=2,
         print.auc=F, max.auc.polygon=F, print.thres=F)
legend(0.6, 0.35, col = 'purple',
       lty=1, lwd=3, bty='n', cex=1,
       paste0('FFN: AUC=', round(roc$auc, 2)))

## Save prediction results

save.descSigs <- function(sig.all, x.test, file.name) {
  write.table(sig.all[colnames(x.test), ], file.name, row.names = F)
}

save.prediction <- function(y.pred, y.real, x.test, file.name) {
  results <- cbind(Prediction = y.pred, Response = y.real, x.test)
  results <- cbind(Patient = rownames(results), results)
  write.table(results, file.name, row.names = F)
}

save.descSigs(sigcat.lung.v2.all, x.test, 'Out/Lung.v2.descSig.txt')
#save.descSigs(sigcat.lung.v2.ad.all, x.test, 'Out/Lung.v2.ad.descSig.txt')
save.prediction(y.pred, y.real, x.test, 'Out/Lung.v2.results.txt')
#save.prediction(y.pred, y.real, x.test, 'Out/Lung.v2.ad.results.txt')


#### Deep learning:: CNN model #################################################

x.train.3d <- array_reshape(x.train, c(nrow(x.train), ncol(x.train), 1))
x.test.3d <- array_reshape(x.test, c(nrow(x.test), ncol(x.test), 1))

model.cnn <- keras_model_sequential() %>%
  layer_conv_1d(filters = 64, kernel_size = 3, activation = 'relu',
                input_shape = c(nfeatures, 1)) %>%
  layer_conv_1d(64, 3, activation = 'relu') %>%
  layer_max_pooling_1d(3) %>%
  layer_conv_1d(128, 3, activation = 'relu') %>%
  layer_conv_1d(128, 3, activation = 'relu') %>%
  layer_global_average_pooling_1d() %>%
  layer_dropout(rate = 0.5) %>%
  layer_dense(units = 1, activation = 'sigmoid')
#summary(model.cnn)

## Compile
model.cnn %>% compile(
  loss      = 'binary_crossentropy',
  optimizer = optimizer_rmsprop(),
  metrics   = c('accuracy')
)

## Training and evaluation
history.cnn = model.cnn %>% fit(
  x.train.3d, y.train,
  epochs = 150,
  #batch_size = 16,
  validation_split = 0.2
)
#plot(history.seq)

perf = model.cnn %>% evaluate(x.test.3d, y.test)

## ROC
y.pred = as.vector(model.cnn %>% predict(x.test.3d))
y.real = as.vector(y.test)
roc <- pROC::roc(y.real, y.pred)
roc$auc
plot.roc(roc, add=T, col="forestgreen", lwd=2,
         print.auc=F, max.auc.polygon=F, print.thres=F)
legend(0.6, 0.3, col = 'forestgreen',
       lty=1, lwd=3, bty='n', cex=1,
       paste0('CNN: AUC=', round(roc$auc, 2)))


#### Random Forest ###################################################################

model.rf <- randomForest(x = x.train, y = y.train, ntree = 1000)

y.pred <- as.vector(predict(model.rf, x.test))
y.real = as.vector(y.test)
roc <- pROC::roc(y.real, y.pred)
roc$auc
plot.roc(roc, add=T, col="deepskyblue2", lwd=2,
         print.auc=F, max.auc.polygon=F, print.thres=F)
legend(0.6, 0.25, col = 'deepskyblue2',
       lty=1, lwd=3, bty='n', cex=1,
       paste0('Random Forest: AUC=', round(roc$auc, 2)))

