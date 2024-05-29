rm(list = ls())
gc()
setwd("~/if01/Analysis/Analysis_XL/Project02_LncRNA/Step03_DEGs")

library(ggplot2)
library(ggrepel)
library(ggthemes)
library(limma)

DEGs_affy <- function(data,group_list){
  
  group_list = factor(group_list)
  
  design <- model.matrix(~0+group_list)
  rownames(design) = colnames(data)
  colnames(design) <- levels(group_list)
  contrast.matrix<-makeContrasts(paste0(rev(unique(group_list)),collapse = "-"),levels = design)
  fit <- lmFit(data,design)
  fit2 <- contrasts.fit(fit, contrast.matrix) ##
  fit2 <- eBayes(fit2)  ## default no trend !!!
  tempOutput = topTable(fit2, coef=1, n=Inf)
  tempOutput = na.omit(tempOutput) 
  tempOutput$id<-rownames(tempOutput)
  return(tempOutput)
  
}

load("~/if01/Analysis/Analysis_XL/Project02_LncRNA/D06_KIRC_normal_count.RData")
load("~/if01/Analysis/Analysis_XL/Project02_LncRNA/D07_KIRC_normal_cli.RData")

quantile(KIRC_normal_cli$age_at_initial_pathologic_diagnosis)
KIRC_normal_cli$Group <- ifelse(KIRC_normal_cli$age_at_initial_pathologic_diagnosis >65,'Old','Young')
colnames(KIRC_normal_count) <- gsub("(.*?)\\.(.*?)\\.(.*?)\\.(.*?)\\..*","\\1\\-\\2\\-\\3",colnames(KIRC_normal_count))
KIRC_normal_cli <- KIRC_normal_cli[KIRC_normal_cli$bcr_patient_barcode %in% colnames(KIRC_normal_count),]
KIRC_normal_cli <- KIRC_normal_cli[!duplicated(KIRC_normal_cli), ]
KIRC_normal_count <- subset(KIRC_normal_count,select = KIRC_normal_cli$bcr_patient_barcode)
KIRC_normal_count[1:5,1:5]
group_list = KIRC_normal_cli$Group
combat_DEGs <- DEGs_affy(data = KIRC_normal_count,group_list = group_list)

load("~/if01/Analysis/Analysis_JLX/Project01_KIRC/1-RawData/TCGA_KIRC/AnnoTCGA.RData")

colnames(combat_DEGs)[7] <- 'ENSG'

res.fc <- merge(AnnoTCGA,combat_DEGs)

res.fc.plot <- res.fc[abs(res.fc$logFC) < 5, ]

DEGs_tt1 <- res.fc.plot
colnames(DEGs_tt1)[2] <- 'id'
DEGs_tt1$id <- as.character(DEGs_tt1$id)
rownames(DEGs_tt1) <- DEGs_tt1$id
DEGs_tt1<-na.omit(DEGs_tt1)
DEGs_tt1$log10pvalue<--log10(DEGs_tt1$P.Value )

library(ggplot2)
library(ggpubr)
library(ggthemes)
library(ggrepel)

Group <- rep('noSig',nrow(DEGs_tt1))
Group[which(DEGs_tt1$logFC >= 3)] <- 'Enrich in Old'
Group[which(DEGs_tt1$logFC <= -3)] <- 'Enrich in Young'
DEGs_tt1$Group <- Group
DEGs_tt1$Group <- factor(DEGs_tt1$Group,levels = c('Enrich in Old','noSig','Enrich in Young'),ordered = T)

mycol <- c('red','gray','red')

DEGs_tt1 <- DEGs_tt1[order(DEGs_tt1$logFC,decreasing = T),]
DEGs_tt1$Rank <- c(1:nrow(DEGs_tt1))

DEGs_tt1 <- DEGs_tt1[DEGs_tt1$logFC != 'Inf',]

library(ggrepel)

p <- ggplot(data = DEGs_tt1) + 
  
  geom_vline(xintercept = c(-2.5,2.5), linetype = c("dashed","dashed"), colour = c('gray','gray'),alpha = c(0.8,0.8)) +
  
  geom_point(aes(x = logFC,y = log10pvalue,color = Group)) +
  
  scale_color_manual(values = mycol) +
  
  xlab('log2 fold change') + ylab('-log10 adj P-value') +
  
  # xlim(c(min(degs_tissue$avg_log2FC),-min(degs_tissue$avg_log2FC))) +
  
  theme_classic() +
  
  # geom_segment(aes(x = -3.2, xend = -5.1,y = 0,yend = 0),arrow = arrow(length = unit(0.2,'cm')),size = 0.1) +
  # 
  # annotate('text', x=-4 , y=3 ,label='non-mPR' ) +
  # 
  # geom_segment(aes(x = 2.8, xend = 4.7,y = 0,yend = 0),arrow = arrow(length = unit(0.2,'cm')),size = 0.1) +
  # 
  # annotate('text', x=3.5 , y=3 ,label='mPR' ) +
  
  # scale_x_continuous(breaks = c(-3,-1,1,3)) +
  
geom_text_repel(aes(x = logFC, y = log10pvalue, 
                    label = ifelse(c(logFC > 2.5 & log10pvalue != 0) | c(logFC < -2.5 & log10pvalue != 0), id,"")),
                colour="darkred", size = 3, box.padding = unit(0.35, "lines"), 
                point.padding = unit(0.3, "lines")) +
  
  theme(
    panel.grid.major =element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    legend.position = 'none',
    text = element_text(family = 'serif')) +
  labs(title = 'KIRC Normal Tissue (Old vs. Young)')

p
ggsave(p,file="p.DEGs.KIRC.Normal.Tissue.pdf",width = 6,height = 6)

save(res.fc,file = 'D01.res.fc.RData')

######################

rm(list = ls())
gc()
setwd("~/if01/Analysis/Analysis_XL/Project02_LncRNA/Step03_DEGs")
load("~/if01/Analysis/Analysis_XL/Project02_LncRNA/D05_KIRC_normal_tpm.RData")
load("~/if01/Analysis/Analysis_XL/Project02_LncRNA/D07_KIRC_normal_cli.RData")

quantile(KIRC_normal_cli$age_at_initial_pathologic_diagnosis)
KIRC_normal_cli$Group <- ifelse(KIRC_normal_cli$age_at_initial_pathologic_diagnosis >65,'Old','Young')
colnames(KIRC_normal_tpm) <- gsub("(.*?)\\.(.*?)\\.(.*?)\\.(.*?)\\..*","\\1\\-\\2\\-\\3",colnames(KIRC_normal_tpm))
KIRC_normal_cli <- KIRC_normal_cli[KIRC_normal_cli$bcr_patient_barcode %in% colnames(KIRC_normal_tpm),]
KIRC_normal_cli <- KIRC_normal_cli[!duplicated(KIRC_normal_cli), ]
KIRC_normal_tpm <- subset(KIRC_normal_tpm,select = KIRC_normal_cli$bcr_patient_barcode)
KIRC_normal_tpm[1:5,1:5]
KIRC_normal_tpm <- data.frame(t(KIRC_normal_tpm))
KIRC_normal_tpm$bcr_patient_barcode <- rownames(KIRC_normal_tpm)
df <- merge(KIRC_normal_cli,KIRC_normal_tpm)
colnames(df)[1:5]
df <- df[,-c(1,2)]
colnames(df)[1] <- 'response'

save(df,file = 'D08_lasso_df.RData')

######
rm(list = ls())
gc()
setwd("~/if01/Analysis/Analysis_XL/Project02_LncRNA/Step03_DEGs")

load("~/if01/Analysis/Analysis_XL/Project02_LncRNA/Step03_DEGs/D08_lasso_df.RData")
load("~/if01/Analysis/Analysis_XL/Project02_LncRNA/Step03_DEGs/D01.res.fc.RData")

## Cor ##

intersect(colnames(df),c('SLC5A11','RP11.503C24.4'))

res.fc$GeneName <- as.character(res.fc$GeneName)
tmp <- res.fc$GeneName[which(abs(res.fc$logFC) > 3 & res.fc$P.Value < 0.05)]

exp <- subset(df,select = intersect(colnames(df),c(tmp,'SLC5A11')))

tmp <- colnames(exp)[!colnames(exp) %in% 'SLC5A11']

cor.res <- data.frame()

for (i in tmp) {
  
  # i <- tmp[1]
  
  res <- cor.test(exp[,'SLC5A11'],exp[,i],method = 'spearman')
  
  res_p <- res$p.value %>% round(.,2)
  res_r <- res$estimate %>% round(.,2)
  
  res_df <- data.frame(Gene = i,
                       pvalue = res_p,
                       rvalue = res_r)
  
  cor.res <- rbind(cor.res,res_df)
  
}

tmp <- cor.res$Gene[which(cor.res$pvalue < 0.05 & cor.res$rvalue > 0.4)]

##

df <- na.omit(df)
df$response <- ifelse(df$response=='Young',0,1)

colnames(df)[1:10]

table(df$response)

library(stringr)
library(gridExtra)
library(future)
library(sva)
library(e1071)
library(pROC)
library(ROCit)
library(caret)
library(doParallel)
library(cancerclass)
library(survival)
library(glmnet)
library(pbapply)
library(survcomp)
Sys.setenv(LANGUAGE = "en") #显示英文报错信息
options(stringsAsFactors = FALSE) #禁止chr转成factor

display.progress = function (index, totalN, breakN=20) {
  
  if ( index %% ceiling(totalN/breakN)  ==0  ) {
    cat(paste(round(index*100/totalN), "% ", sep=""))
  }
} 

df <- subset(df,select = c('response',tmp))

# 单变量cox模型筛选基因
Logoutput <- NULL
for(i in 2:ncol(df)){
  # i <- 2
  display.progress(index = i,totalN = ncol(df),breakN = 20)
  g <- colnames(df)[i]
  mod1 <- glm(response~get(colnames(df)[i]), family = binomial(link = 'logit'),data = df)
  fit <- summary(mod1)
  Logoutput=rbind(Logoutput,data.frame(gene=g,
                                       OR=as.numeric(fit$coefficients[,"Estimate"])[2],
                                       z=as.numeric(fit$coefficients[,"z value"])[2],
                                       pvalue=as.numeric(fit$coefficients[,"Pr(>|z|)"])[2],stringsAsFactors = F))
}

(log.res <- Logoutput[which(Logoutput$pvalue < 0.1),"gene"])

df <- subset(df,select = c('response',log.res))

df$response <- factor(df$response,levels = c(0,1),labels = c('Young','Old'))

#### Step 01 Lasso Logistic ####

# 运行1000次multivariate cox model with lasso penalty
iter.times <- 1000 # 设置迭代次数，速度非常慢请耐心，例文是1000次

lasso_fea_list <- list()
lambda_all <- c('lambda.1se','lambda.min')
lambda_choose <- lambda_all[2]
list.of.seed <- 1:iter.times
lasso_fea_list <- pblapply(list.of.seed, function(x){ # 大概运行2天
  set.seed(list.of.seed[x])
  
  outcome <- df$response
  
  xx <- as.matrix(df[,log.res])
  
  cvfit = cv.glmnet(xx,outcome,family="binomial")
  
  # 取出最优lambda
  fea <- rownames(coef(cvfit, s = lambda_choose))[coef(cvfit, s = lambda_choose)[,1]!= 0]
  if(is.element("(Intercept)", fea)) {
    lasso_fea <- sort(fea[-1]) # 去掉截距项并排序
  } else {
    lasso_fea <- sort(fea)
  }
  return(lasso_fea)
})

# 输出每次运行的基因集合
lasso_res <- NULL
for(i in 1:iter.times) {
  lasso_res <- rbind.data.frame(lasso_res,
                                data.frame(iteration = i,
                                           n.gene = length(lasso_fea_list[[i]]),
                                           genelist = paste0(lasso_fea_list[[i]], collapse = " | "),
                                           stringsAsFactors = F),
                                stringsAsFactors = F)
}

# 检查相同基因数目的集合成分是否一致
## 注意原文没有提到该做法，但的确存在相同基因数目的模型基因成分不完全一致的情况，已致邮给通讯作者未回复
## 这里可以看到32个基因的模型存在两种基因组分，虽不影响35基因被最终选中，但应考虑为不同模型
uniquelist <- unique(lasso_res$genelist)
uniquelab <- LETTERS[1:length(uniquelist)]
lasso_res$uniquelab <- NA
for (i in 1:length(uniquelist)) {
  lasso_res[which(lasso_res$genelist == uniquelist[i]),"uniquelab"] <- uniquelab[i]
}
lasso_res$label <- paste(lasso_res$n.gene,"genes",lasso_res$uniquelab,sep = "_") # 最终模型标签

table(lasso_res$label)

sel.iter <- lasso_res[which(lasso_res$label == "9_genes_A"),"iteration"][1] # 选中模型为35_genes_A的某一次迭代种子（随便算哪个，基因集合都是一样的）
set.seed(sel.iter) # 设置当前种子以复现该基因集

lasso_hubgenes <- strsplit(lasso_res$genelist[sel.iter],' | ')[[1]][which(strsplit(lasso_res$genelist[sel.iter],' | ')[[1]] != '|')]

pdf(file = 'Figure Lasso.pdf',height = 5,width = 5)

par(mfrow = c(1,1)) # 创建画图并分割成左右两块
## 图A
par(bty="o", mgp = c(1.9,.33,0), mar=c(5.1,3.1,2.1,2.1)+.1, las=1, tcl=-.25)
plotdata <- sort(table(lasso_res$label))
a <- barplot(plotdata,
             ylab = "Frequency",
             col = "#00CDCD",
             main = "Frequency of models",
             border = NA,
             xaxt = "n",
             ylim = c(0,max(plotdata) + 70))
axis(side = 1, at = a, names(plotdata), las = 2) # 添加x轴标签
text(a, as.numeric(plotdata) + 30, as.numeric(plotdata), xpd = T) # 添加基因数目
par(new = T, bty="o", mgp = c(1.9,.33,0), mar=c(5.1,3.1,2.1,2.1)+.1, las=1, tcl=-.25) # 补全黑色边框
plot(NULL, NULL,
     col = "white",
     xlim = range(a), ylim = range(plotdata),
     xlab = "", ylab = "",
     xaxt = "n", yaxt = "n")

dev.off()

save(lasso_hubgenes,file = 'D01_lasso_hubgenes.RData')

## LQV ##

# ensure results are repeatable
set.seed(1234)
# load the library
library(mlbench)
library(caret)
# # load the dataset
# data(PimaIndiansDiabetes)
# prepare training scheme
cl <- makePSOCKcluster(10)
clusterEvalQ(cl, .libPaths("/home/a01/R/x86_64-pc-linux-gnu-library/4.2"))
registerDoParallel(cl)
control <- trainControl(method="repeatedcv", number=10, repeats=3)
# train the model
model <- train(response~., data=df, preProcess="scale", trControl=control)
# estimate variable importance
importance <- varImp(model, scale=FALSE)
# summarize importance
# importance$importance <- importance$importance[1:10,]
print(importance)
# plot importance
stopCluster(cl)

cutoff_line <- 0.5
p <- ggplot(importance) + geom_hline(aes(yintercept=cutoff_line),col='red',lty=3) + theme_classic()

pdf(file = 'Figure LQV.pdf',height = 6,width = 4)
# plot(importance)
p
dev.off()

LQV_HubGenes <- rownames(importance$importance)[which(importance$importance$R > cutoff_line)]

rm(cutoff_line)

save(LQV_HubGenes,file = 'D02_LQV_HubGenes.RData')

### bagged trees ###

# ensure the results are repeatable
set.seed(1234)
# define the control using a random forest selection function
control <- rfeControl(functions=treebagFuncs, method="cv", number=30)
# run the RFE algorithm
cl <- makePSOCKcluster(10)
clusterEvalQ(cl, .libPaths("/home/a01/R/x86_64-pc-linux-gnu-library/4.2"))
registerDoParallel(cl)
treebag_results <- rfe(df[,2:ncol(df)], df[,1], sizes=c(2:ncol(df)), rfeControl=control)
# summarize the results
print(treebag_results)
pdf(file = 'Figure bagged trees.pdf',height = 6,width = 4)
plot(treebag_results, type=c("o")) 
dev.off()
stopCluster(cl)
# list the chosen features
(treebag_HubGenes <- predictors(treebag_results))
save(treebag_HubGenes,file = 'D03_treebag_HubGenes.RData')
write.csv(treebag_results$results,file = 'TableS_treebag_accuracy.csv',row.names = T)

### Boruta ####

library(Boruta)
set.seed(1234)
boruta.train <- Boruta(response~., data = df, doTrace = 3)
print(boruta.train)

pdf(file = 'Figure bpruta.pdf',height = 6,width = 8)
plot(boruta.train, xlab = "", xaxt = "n")
lz<-lapply(1:ncol(boruta.train$ImpHistory),function(i)
  boruta.train$ImpHistory[is.finite(boruta.train$ImpHistory[,i]),i])
names(lz) <- colnames(boruta.train$ImpHistory)
Labels <- sort(sapply(lz,median))
# text(x=seq(1,89,by=1),y=-5.5, srt = 45, adj = 1, labels = names(Labels),xpd = TRUE,cex = 0.5)
axis(side = 1,las=2,labels = names(Labels),
     at = 1:ncol(boruta.train$ImpHistory), cex.axis = 0.7)
dev.off()

boruta_HubGenes <- TentativeRoughFix(boruta.train)
print(boruta_HubGenes)

(boruta_HubGenes <- getSelectedAttributes(boruta_HubGenes, withTentative = F))

save(boruta_HubGenes,file = 'D04_boruta_HubGenes.RData')

## random forest ##

set.seed(1234)
cl <- makePSOCKcluster(10)
clusterEvalQ(cl, .libPaths("/home/a01/R/x86_64-pc-linux-gnu-library/4.2"))
registerDoParallel(cl)
# define the control using a random forest selection function
control <- rfeControl(functions=rfFuncs, method="cv", number=30)
# run the RFE algorithm
rfe_results <- rfe(df[,2:ncol(df)], df[,1], sizes=c(2:ncol(df)), rfeControl=control)
# summarize the results
print(rfe_results)
stopCluster(cl)

pdf(file = 'Figure random forest.pdf',height = 4,width = 6)
plot(rfe_results, type=c("o")) 
dev.off()

# list the chosen features
(REF_HubGenes <- predictors(rfe_results))
save(REF_HubGenes,file = 'D05_REF_HubGenes.RData')
write.csv(rfe_results$results,file = 'TableS_rfe_accuracy.csv',row.names = T)

### Bayesian ###

# install.packages("klaR")
library(klaR)

set.seed(1234)
# define the control using a random forest selection function
control <- rfeControl(functions=nbFuncs, method="cv", number=30)
cl <- makePSOCKcluster(10)
clusterEvalQ(cl, .libPaths("/home/a01/R/x86_64-pc-linux-gnu-library/4.2"))
registerDoParallel(cl)

# run the RFE algorithm
results <- rfe(df[,2:ncol(df)], df[,1], sizes=c(2:ncol(df)), rfeControl=control)

stopCluster(cl)

pdf(file = 'Figure Bayesian.pdf',height = 6,width = 4)
plot(results, type=c("o")) 
dev.off()

# list the chosen features
(Bayesian_HubGenes <- predictors(results))
save(Bayesian_HubGenes,file = 'D06_Bayesian_HubGenes.RData')
write.csv(results$results,file = 'TableS_nbe_accuracy.csv',row.names = T)

# ## SVM ##
# 
# library(tidyverse)
# library(glmnet)
# source('msvmRFE.R')   #文件夹内自带
# library(VennDiagram)
# library(sigFeature)
# library(e1071)
# library(caret)
# library(randomForest)
# 
# input <- df
# 
# #采用五折交叉验证 (k-fold crossValidation）
# svmRFE(input, k = 5, halve.above = 100) #分割数据，分配随机数
# 
# nfold = 5
# nrows = nrow(input)
# folds = rep(1:nfold, len=nrows)[sample(nrows)]
# folds = lapply(1:nfold, function(x) which(folds == x))
# results = lapply(folds, svmRFE.wrap, input, k=5, halve.above=100) #特征选择
# 
# top.features = WriteFeatures(results, input, save=F) #查看主要变量
# head(top.features)
# write.csv(top.features,"TableS_feature_svm.csv")
# 
# library(snowfall)
# 
# setwd("~/if01/Analysis/Analysis_JLX/Project07_scRNA&Bulk&ICI&PanCancer/Step6_Model_ICB")
# 
# # 运行时间主要取决于选择变量的个数，一般的电脑还是不要选择太多变量
# 
# featsweep = sfLapply(1:ncol(input), FeatSweep.wrap, results, input)
# featsweep
# 
# # 画图
# no.info = min(prop.table(table(input[,1])))
# errors = sapply(featsweep, function(x) ifelse(is.null(x), NA, x$error))
# 
# #dev.new(width=4, height=4, bg='white')
# pdf("B_svm-error.pdf",width = 5,height = 5)
# PlotErrors(errors, no.info=no.info) #查看错误率
# dev.off()
# 
# #dev.new(width=4, height=4, bg='white')
# pdf("B_svm-accuracy.pdf",width = 5,height = 5)
# Plotaccuracy(1-errors,no.info=no.info) #查看准确率
# dev.off()
# 
# ###

## Upset plot ##
rm(list = ls())
gc()
setwd("~/if01/Analysis/Analysis_XL/Project02_LncRNA/Step03_DEGs")

load("./D01_lasso_hubgenes.RData")
load("./D02_LQV_HubGenes.RData")
load("./D03_treebag_HubGenes.RData")
load("./D04_boruta_HubGenes.RData")
load("./D05_REF_HubGenes.RData")
load("./D06_Bayesian_HubGenes.RData")

# BiocManager::install('UpSetR')
library(UpSetR)

listInput <- list(Lasso = lasso_hubgenes, 
                  LQV = LQV_HubGenes, 
                  Boruta = boruta_HubGenes,
                  RandomForest = REF_HubGenes,
                  Bayesian = Bayesian_HubGenes,
                  Treebag = treebag_HubGenes)

pdf(file = 'Figure Upset.pdf',height = 5,width = 6,onefile = F)

upset(fromList(listInput), order.by = "freq",nsets = length(listInput))

dev.off()

Final_HubGenes <- c(lasso_hubgenes,LQV_HubGenes,
                    boruta_HubGenes,
                    REF_HubGenes,Bayesian_HubGenes,treebag_HubGenes)
table(Final_HubGenes)

Final_HubGenes_ICB <- which(table(Final_HubGenes) > 4) %>% names
save(Final_HubGenes_ICB,file = 'D07_Final_HubGenes_ICB.RData')

stopCluster(cl)






