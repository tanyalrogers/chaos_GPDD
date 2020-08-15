#Heatmaps

library("ggplot2")
library("dplyr")
library("tidyr")
library("rEDM")
library("purrr")
library(ComplexHeatmap)

load(file = "./data/gpdd_results.Rdata")
source("~/GRAD SCHOOL/R reference/ggplot themes rogers.R")

otherresults=read.csv("./data/ChaosMetaAnalysisEmpiricalResults.csv", stringsAsFactors = F)
gpdd_d=left_join(gpdd_d, otherresults, by="MainID")

colors1=c("tomato", "cornflowerblue", "darkolivegreen3")
gpdd_dp=arrange(gpdd_d, TaxonomicClass2, desc(gle_mo)) %>% filter(monotonicR2<0.5)

mat=as.matrix(select(gpdd_dp,Recurrence.Plot,Permutation.Entropy,Visibility.Algorithm))
rownames(mat)=gpdd_dp$CommonName
colnames(mat)=c("RP", "PE", "VA")
LEs=HeatmapAnnotation(GLE=anno_points(gpdd_dp$gle_mo, pch = 16, gp = gpar(col=ifelse(gpdd_dp$glebest>0,"tomato","black"))),
                      which="row", width = unit(3, "cm"), 
                      annotation_name_side = "top", annotation_name_rot = 0)
png("empirical_heatmap.png", width = 6, height = 12, res = 300, units = "in")
Heatmap(mat, col=colors1, name="Dynamics", width = unit(3, "cm"),
        row_split = gpdd_dp$TaxonomicClass2, right_annotation = LEs,
        column_names_side = "top", column_names_rot = 0, column_names_centered = T,
        row_names_gp = gpar(cex=0.5))
for(i in 1:n_distinct(gpdd_dp$TaxonomicClass2)) {
  decorate_annotation("GLE", slice=i,{grid.lines(unit(c(0, 0), "native"), unit(c(0, 1), "npc"), gp = gpar(lty=2))})
}
dev.off()

#CIs
gpdd_results$stabilitybest=cbind(select(gpdd_results, stability1:stability5),gpdd_d$bestmodel) %>% apply(1, function(x) {m=as.numeric(x["gpdd_d$bestmodel"]); x[m][[1]]})
gpdd_d$lle_ci90lower=map_dbl(gpdd_results$stabilitybest, function(x) {
  std=sd(x$lle, na.rm=T)
  n=length(which(!is.na(x$lle)))
  ci=x$lle_avg+std/sqrt(n)*qt(p=0.05, df=n-1)})
gpdd_d$lle_ci90upper=map_dbl(gpdd_results$stabilitybest, function(x) {
  std=sd(x$lle, na.rm=T)
  n=length(which(!is.na(x$lle)))
  ci=x$lle_avg+std/sqrt(n)*qt(p=0.95, df=n-1)})
gpdd_d$lle_sign=ifelse(gpdd_d$lle_ci90lower>0, "lle_positive", "lle_negative")
gpdd_d$gle_within=ifelse(gpdd_d$glebest>gpdd_d$lle_ci90lower & gpdd_d$glebest<gpdd_d$lle_ci90upper, T, F)
length(which(gpdd_d$gle_within))/length(which(!is.na(gpdd_d$gle_within)))

ggplot(gpdd_d, aes(y=factor(MainID))) +
  #facet_grid(TSlength~Classification) + 
  #geom_point(aes(y=lle_avgbest), color="black") +
  geom_errorbarh(aes(xmin=lle_ci90lower, xmax=lle_ci90upper, color=lle_sign)) +
  geom_point(aes(x=glebest, fill=glesign), pch=21) + 
  theme_bw() + geom_vline(xintercept = 0)

table(gpdd_d$lle_sign,gpdd_d$glesign)

#test datasets
# test=data.frame(Class=rep(c("Aves","Mammalia","Insecta"), each=4),
#                 Species=paste("species",1:(3*4)),
#                 gle=c(0.1,0.5,-0.02,-0.5,-0.25,-0.01, 0.2, 0.8, 0.1, 0.06, 0.7,-0.4),
#                 rp=c("chaotic","chaotic","periodic","stochastic","stochastic","periodic","stochastic","stochastic","stochastic","stochastic","stochastic","stochastic"),
#                 pe=c("periodic","chaotic","stochastic","stochastic","stochastic","stochastic","chaotic","chaotic","chaotic","stochastic","stochastic","periodic"),
#                 va=c("periodic","chaotic","stochastic","stochastic","periodic","stochastic","chaotic","stochastic","stochastic","stochastic","stochastic","stochastic"))


