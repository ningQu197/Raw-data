library(ggplot2)

df = subset(df, Score>0)
ggplot(df, aes(x=reorder(TF, -Score, median), y=Score))+
geom_boxplot(outlier.size=0.1, width=0.3, fill="steelblue1")+
ylab("Score of AUCell")+
theme_bw()+
theme(text=element_text(family="Times"), axis.text.x=element_text(angle=30, hjust=1, size=13, face="bold"), axis.text.y=element_text(size=13, face="bold"), axis.title.x=element_blank(), axis.title.y=element_text(size=13, face="bold"))




