library(ComplexHeatmap)
colors1=c("tomato", "cornflowerblue", "darkolivegreen3")
test=data.frame(Class=rep(c("Aves","Mammalia","Insecta"), each=4), 
                Species=paste("species",1:(3*4)),
                gle=c(0.1,0.5,-0.02,-0.5,-0.25,-0.01, 0.2, 0.8, 0.1, 0.06, 0.7,-0.4),
                rp=c("chaotic","chaotic","periodic","stochastic","stochastic","periodic","stochastic","stochastic","stochastic","stochastic","stochastic","stochastic"),
                pe=c("periodic","chaotic","stochastic","stochastic","stochastic","stochastic","chaotic","chaotic","chaotic","stochastic","stochastic","periodic"),
                va=c("periodic","chaotic","stochastic","stochastic","periodic","stochastic","chaotic","stochastic","stochastic","stochastic","stochastic","stochastic"))

mat=as.matrix(select(test,rp,pe,va))
rownames(mat)=test$Species
colnames(mat)=c("RP", "PE", "VA")
LEs=HeatmapAnnotation(GLE=anno_points(test$gle, pch = 16, gp = gpar(col=ifelse(test$gle>0,"tomato","black"))),
                      which="row", width = unit(3, "cm"), 
                      annotation_name_side = "top", annotation_name_rot = 0)
Heatmap(mat, col=colors1, name="Dynamics", 
        row_split = test$Class, right_annotation = LEs,
        column_names_side = "top", column_names_rot = 0, column_names_centered = T)
for(i in 1:n_distinct(test$Class)) {
  decorate_annotation("GLE", slice=i,{grid.lines(unit(c(0, 0), "native"), unit(c(0.5, 4.5), "native"), gp = gpar(lty=2))})
}

