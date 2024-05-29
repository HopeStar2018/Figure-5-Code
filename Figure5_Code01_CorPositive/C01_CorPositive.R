rm(list = ls())
gc()
setwd("~/if01/Analysis/Analysis_XL/Project02_LncRNA/Step02_CorPositive")

load("~/if01/Analysis/Analysis_XL/Project02_LncRNA/Step01_WGCNA_adjNormal/D01_Exp_Cli.RData")

exp <- data.frame(t(Exp_df))

intersect(colnames(exp),'SLC5A11')

exp <- data.frame(exp)

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

save(cor.res,file = 'D01.cor.res.RData')

#####
rm(list = ls())
gc()
setwd("~/if01/Analysis/Analysis_XL/Project02_LncRNA/Step02_CorPositive")

load("~/if01/Analysis/Analysis_XL/Project02_LncRNA/Step03_DEGs/D08_lasso_df.RData")

load("~/if01/Analysis/Analysis_XL/Project02_LncRNA/Step03_DEGs/D07_Final_HubGenes_ICB.RData")

exp <- subset(df,select = c('SLC5A11',Final_HubGenes_ICB))

dir.create('Figure_Cor_res')

setwd("./Figure_Cor_res")

for (j in Final_HubGenes_ICB) {
    
  # j <- HubGenes[1]
  
  res <- cor.test(exp[,'SLC5A11'],exp[,j],method = 'spearman')
  
  r_value <- res$estimate %>% unname() %>% round(.,2)
  p_value <- res$p.value %>% round(.,2)
  
  p <- ggplot(exp, aes(x=SLC5A11, y=get(j))) +
    geom_point() +
    geom_smooth(method=lm , color="red", fill="#69b3a2", se=TRUE) +
    ylab(j) + xlab('SLC5A11') +
    ggtitle(paste0("r = ", r_value,
                   "\nP = ", sprintf("%1.1e", p_value)))+
    theme_classic() +
    theme(text = element_text(family = 'serif'))
  
  pdf(file = paste0('./Figure_Corplot_','SLC5A11','_',j,'.pdf'),height = 3,width = 3)
  print(p)
  dev.off()
  
}

## Hub mRNA Boxplot ##
rm(list =ls())
gc()
setwd("~/if01/Analysis/Analysis_XL/Project02_LncRNA/Step02_CorPositive")

dir.create('Figure_HubBoxplot_res')

setwd("./Figure_HubBoxplot_res")

load("~/if01/Analysis/Analysis_XL/Project02_LncRNA/Step03_DEGs/D08_lasso_df.RData")

load("~/if01/Analysis/Analysis_XL/Project02_LncRNA/Step03_DEGs/D07_Final_HubGenes_ICB.RData")

df <- subset(df,select = c('response',Final_HubGenes_ICB))

library(ggpubr)

for (i in Final_HubGenes_ICB) {
  
  p <- ggboxplot(data = df,
                 x = 'response',
                 y = i,
                 fill = 'response',
                 xlab = '') + stat_compare_means() +
    # ylim(0,4) +
    theme(panel.grid.major =element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line.x = element_line(colour = "black"),
          axis.line.y = element_line(colour = "black"),
          text = element_text(family = 'serif'),
          legend.position = 'none')
  
  ggsave(filename = paste0('Figure Boxplot ',i,'.pdf'),height = 4,width = 3,plot = p)
  
}






















