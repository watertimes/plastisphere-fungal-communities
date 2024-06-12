###figure 7A the mantel test and correlation 

library(dplyr)
library(linkET)
library(ggplot2)
#Data Loading
speciese <- read.csv("asv.csv",header = T,row.names = 1)
env <- read.csv("env.csv",header = T,row.names = 1)

#mantel test
mantel01 <- mantel_test(speciese, env,
                        spec_select = list(MP = 1:19,
                                           soil = 20:38)) %>% 
  mutate(rd = cut(r, breaks = c(-Inf, 0.2, 0.4, Inf),
                  labels = c("< 0.2", "0.2 - 0.4", ">= 0.4")),
         pd = cut(p, breaks = c(-Inf, 0.01, 0.05, Inf),
                  labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))
mantel01

#plot
qcorrplot(correlate(env), type = "lower", diag = FALSE) +
  geom_square() +
  geom_couple(aes(colour = pd, size = rd), data = mantel01, curvature = 0.1) +
  scale_fill_gradientn(colours = RColorBrewer::brewer.pal(11, "RdBu")) +
  scale_size_manual(values = c(0.5, 1, 2)) +
  scale_colour_manual(values = color_pal(3)) +
  guides(size = guide_legend(title = "Mantel's r",
                             override.aes = list(colour = "grey35"), 
                             order = 2),
         colour = guide_legend(title = "Mantel's p", 
                               override.aes = list(size = 3), 
                               order = 1),
         fill = guide_colorbar(title = "Pearson's r", order = 3)) +
  labs(title = "Mantel test")



###figure 7b the random forest test
library(psych)
library(reshape2)
library(ggplot2) 
library(randomForest)
library(patchwork)

#Data Loading
myro <- as.data.frame((read.csv("myro.csv", header=TRUE)))
spearman <- corr.test(myro[,2:20], myro[,21:41], method = 'spearman', adjust = 'none')#取相应的列并计算相关系数
spearman #View the correlation coefficient matrix
#Organize the results for easier plotting
r <- data.frame(spearman$r)  #Obtain the correlation coefficient R square
r$myro <- rownames(r)
r <- melt(r, id = 'myro')
spearman <- cbind(r)
spearman#View the new data box reconstructed with R-values
#The above code processes the correlation coefficient matrix into three columns of data boxes
#plot
p1 <- ggplot() +
  geom_tile(data = spearman, aes(x = variable, y = myro, fill = value)) +
  scale_fill_gradientn(colors = c('#2D6DB1', 'white', '#DC1623'), limit = c(-1, 1)) +
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black'), legend.key = element_blank(), 
        axis.text.x = element_text(color = 'black', angle =45, hjust = 1, vjust = 1), axis.text.y = element_text(color = 'black'), axis.ticks = element_line(color = 'black')) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(y = '', x = '', fill = 'Correlation')
p1

#The important values of physical and chemical properties are calculated one by one
set.seed(123)
dis1_forest <- randomForest(dissimilarity~A253+E2+A254+A260+NO+NH+AP+TC+TN+ph+EC+long+lati+MAP+MAT+clay+silt+sand, data =myro, importance = TRUE, ntree = 500)

#...calculated one by one

#Use the function importance() to view the score indicating the importance of each predictor
dis1<- data.frame(importance(dis1_forest, scale = TRUE), check.names = FALSE)
IMportance_t <- data.frame(cbind(dis1$`%IncMSE`))
colnames(IMportance_t) <- c('dis_1')#name the rowname
rownames(IMportance_t) <- c('A253', 'E2', "A254") #Named column name,some names are omitted here 
IMportance_t[IMportance_t<0] <- 0#Filter all items of importance below 0
IMportance_t[is.na(IMportance_t)] <- 0#Assign a null value to 0
write.csv(IMportance_t,file = "IMportance_t1.CSV")

impdata <- read.csv("IMportance_t1.CSV",header = TRUE)#Write the reprocessing data to R

measure_name=setdiff(colnames(impdata),
                     c('Items'))#Assign one-to-many values by the first column
#prepare data to plot
data1=melt(impdata ,
           id.vars='Items', 
           measure.vars=measure_name,
           variable.name = "sample", 
           value.name = "expr")
#plot
p <- p1+geom_point(data = data1, aes(x = sample, y = Items, size = expr*10), shape = 1) +
  scale_size_continuous(range = c(0,6)) +#Control circle size
  labs(size = 'Importance (%)')#Add circles to the heatmap
p

exp <- read.csv("exp.CSV",header = TRUE)#Read the total interpretation value data box, here you need to manually count the relevant data
exp
#combine two plots
p2 <- ggplot(exp, aes(microbe, Values)) +
  geom_col(fill = "steelblue") +
  theme_bw() +
  theme_classic()+
  theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank())+coord_cartesian(ylim = c(-120, 50))+
  ylab("Explained variation(%)")
p2 +p + plot_layout(ncol= 1, widths = c(2, 1))

