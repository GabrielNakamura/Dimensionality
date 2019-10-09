library(vegan)
source(here::here("functions","simm_comm.R"))
source(here::here("functions", "simul_dimensionality.R"))
source(here::here("functions", "ImportanceVal_V2.R"))
source(here::here("functions", "dimensionality.R"))

######Scenario 1 Pd1, Fd1 and Riq1 - equal importance ####
limit=0
Ncomm= 50
Nsp= 20
filter= 0.9
gradual= TRUE
richness.equal= FALSE
runs<-999
method="max"
power<-3
matrix.M_S1_list<- vector(mode="list",length=runs)
ImportanceValues_S1<-vector(mode = "list",length = runs)
EE_S1<-matrix(NA,nrow = runs, ncol = 1,dimnames = list(paste("run",1:runs,sep = ""),"Evenness_eng"))
for(i in 1:runs){
 tree<- geiger::sim.bdtree(b = 0.001,n = 200,stop = "taxa")
 trait<- ape::rTraitCont(ape::compute.brlen(tree, power = power), model = "BM") #simulating the traits
 trait<- scales::rescale(trait, to = c(-1, 101)) #re-scaling the trait??
 comm.simul_S1<-sim.comm_2(tree = tree,limit = limit,Ncomm = Ncomm,Nsp = Nsp,richness.equal = richness.equal,filter = filter,gradual = gradual)
 rich<-rowSums(comm.simul_S1)
 PD<-picante::pd(samp = comm.simul_S1,tree = tree)$PD
 FD<-picante::pd(samp = comm.simul_S1,tree = ape::as.phylo(hclust(dist(trait))))$PD
 ses.PD<- picante::ses.pd(samp = comm.simul_S1,tree = tree,null.model = "taxa.labels",runs = 10,iterations = 10)$pd.obs.z
 ses.FD<- picante::ses.pd(samp = comm.simul_S1,tree = ape::as.phylo(hclust(dist(trait))),null.model = "taxa.labels",runs = 10,iterations = 10)$pd.obs.z
 matrix.Mses_S1<- cbind(rich,ses.PD,ses.FD) #matrix M with standardized effect size for functional and phylogenetic diversity
 matrix.M_S1<- cbind(rich,PD,FD) #matrix M with raw metrics for functional and phylogenetic diversity
 matrix.M_S1_list[[i]]<- matrix.M_S1
 dim_corr_M_S1<- dimensionality(matrix.M = matrix.M_S1,scale = TRUE,method = "standardize",evenness = "Camargo") #baixa dimensionalidade
 IV_M_S1<- ImportanceVal(matrix.M = matrix.M_S1,scale = TRUE,method = method,stopRule = FALSE)
 EE_S1[i,]<- dim_corr_M_S1
 ImportanceValues_S1[[i]]<- IV_M_S1
 print(paste("run",i,sep = ""))
}
scenario_S1<- matrix(NA,nrow= runs, ncol= ncol(ImportanceValues_S1[[1]]$IV.obs),dimnames= list(paste("run",1:runs,sep=""), colnames(ImportanceValues_S1[[1]]$IV.obs)))??

for(i in 1:length(ImportanceValues_S1)){
 scenario_S1[i,]<-colSums(ImportanceValues_S1[[i]]$IV.obs[c(1,2),])
}

scenario_S1<- cbind(scenario_S1,EE_S1) #join EE and IV values
nsuccess_S1<-length(which(scenario_S1[,2]>scenario_S1[,3])) #number of sucessess
mean_scenarioS1<- apply(scenario_S1, 2, mean) #mean values of IV and EE for scenario 1
sd_scenarioS1<-apply(scenario_S1, 2, sd)
mean(unlist(lapply(matrix.M_S1_list, function(i) cor(i[,2],i[,3])))) #media de correlacao entre PD e FD

#Gradient graphic for scenario 1
windows()
layout.show(layout(matrix(c(0,0,1,2,3,4,5),nrow=2,ncol=3,byrow=T)))
barplot(decostand(matrix.M_S1[,1],"max")[,1], main= "rich")
barplot(decostand(matrix.M_S1[,2],"max")[,1], main= "PD")
barplot(decostand(matrix.M_S1[,3],"max")[,1], main= "FD")
plot(tree, show.tip.label = FALSE)
barplot(trait[tree$tip.label], horiz=TRUE, space=0, names="")
plotrix::color2D.matplot(t(comm.simul_S1),ylab="",xlab="",main="", yrev=F,cs1=c(1,0),cs2=c(1,0),cs3=c(1,0),border="white",axes=FALSE)

#Figure 1 - setting the problem - Scenario With correlation#######
quartz()
windows()
plot(decostand(matrix.M_S1_list[[3]][,"PD"],"max"), decostand(matrix.M_S1_list[[3]][,"FD"],"max"), pch= 19, cex= 3,ylim=c(0,1))
points(decostand(matrix.M_S9_list[[3]][,"PD"],"max"), decostand(matrix.M_S9_list[[3]][,"FD"],"max"), pch= 2, cex= 3)

#######Scenario 2 - PD1, FD0, Ric0 - different importance Richness lower IV values than Phylo#####
limit=0
Ncomm= 50
Nsp= 20
filter= 0.9
gradual= TRUE
richness.equal= TRUE
runs<-999
method="max"
matrix.M_S2_list<- vector(mode="list",length=runs)
ImportanceValues_S2<-vector(mode = "list",length = runs)
EE_S2<-matrix(NA,nrow = runs, ncol = 1,dimnames = list(paste("run",1:runs,sep = ""),"Evenness_eng"))
for(i in 1:runs){
 tree<-geiger::sim.bdtree(b = 0.001,n = 200,stop = "taxa")
 comm.simul_S2<-sim.comm_2(tree = tree,limit = limit,Ncomm = Ncomm,Nsp = Nsp,richness.equal = richness.equal,filter = filter,gradual = gradual)
 rich<-c(rowSums(comm.simul_S2)+rnorm(n = nrow(comm.simul_S2),mean = 0.0001,sd = 0.001))
 newTrait<-runif(n = length(tree$tip.label),min = 100,max = 101)
 names(newTrait)<-tree$tip.label
 PD<-picante::pd(samp = comm.simul_S2,tree = tree)$PD
 FD<-picante::pd(samp = comm.simul_S2,tree = ape::as.phylo(hclust(dist(newTrait))))$PD #modificacao para FD0
 ses.PD<- picante::ses.pd(samp = comm.simul_S2,tree = tree,null.model = "taxa.labels",runs = 10,iterations = 10)$pd.obs.z
 ses.FD<- picante::ses.pd(samp = comm.simul_S2,tree = ape::as.phylo(hclust(dist(newTrait))),null.model = "taxa.labels",runs = 10,iterations = 10)$pd.obs.z
 matrix.Mses_S2<-cbind(rich,ses.PD,ses.FD) #matrix M with standardized effect size for functional and phylogenetic diversity
 matrix.M_S2<-cbind(rich,PD,FD) #matrix M with raw metrics for functional and phylogenetic diversity
 matrix.M_S2_list[[i]]<- matrix.M_S2
 dim_corr_M_S2<-dimensionality(matrix.M = matrix.M_S2,scale = TRUE,method = "standardize",evenness = "Camargo") #baixa dimensionalidade
 IV_M_S2<- ImportanceVal(matrix.M = matrix.M_S2,scale = TRUE,method = method,stopRule = FALSE)
 ImportanceValues_S2[[i]]<-IV_M_S2
 EE_S2[i,]<-dim_corr_M_S2
 print(paste("run",i,sep = ""))
}
scenario_S2<- matrix(NA,nrow= runs, ncol= ncol(ImportanceValues_S2[[1]]$IV.obs),dimnames= list(paste("run",1:runs,sep=""), colnames(ImportanceValues_S2[[1]]$IV.obs)))??

for(i in 1:length(ImportanceValues_S2)){
 scenario_S2[i,]<-colSums(ImportanceValues_S2[[i]]$IV.obs[c(1,2),])
}
scenario_S2<- cbind(scenario_S2,EE_S2) #join EE and IV values
nsuccess_S2<-length(which(scenario_S2[,2]>scenario_S2[,3]&scenario_S2[,2]>scenario_S2[,1]))
mean_scenarioS2<- apply(scenario_S2, 2, mean) #mean values of IV and EE for scenario S2
sd_scenarioS2<- apply(scenario_S2, 2, sd) #mean values of IV and EE for scenario S2

#Figure 1 - setting the problem
quartz()
windows()
plot(decostand(matrix.M_S2_list[[100]][,"PD"],"max"), decostand(matrix.M_S2_list[[100]][,"FD"],"max"), pch= 2, cex= 2)

#Gradient graphic for scenario 2
quartz()
windows()
layout.show(layout(matrix(c(0,0,1,2,3,4,5),nrow=2,ncol=3,byrow=T)))
barplot(decostand(matrix.M_S2[,1],"max")[,1], main= "rich")
barplot(decostand(matrix.M_S2[,2],"max")[,1], main= "PD")
barplot(decostand(matrix.M_S2[,3],"max")[,1], main= "FD")
plot(tree, show.tip.label = FALSE)
barplot(newTrait[tree$tip.label], horiz=TRUE, space=0, names="")
plotrix::color2D.matplot(t(comm.simul_S2),ylab="",xlab="",main="", yrev=F,cs1=c(1,0),cs2=c(1,0),cs3=c(1,0),border="white",axes=FALSE)


######Scenario 3 - PD1, FD1, Ric0 - different importance, richness low values of IV#####
limit=0
Ncomm= 50
Nsp= 20
filter= 0.9
gradual= TRUE
richness.equal= TRUE
runs<-999
power<-2
method="max"
matrix.M_S3_list<- vector(mode="list",length=runs)
ImportanceValues_S3<-vector(mode = "list",length = runs)
EE_S3<-matrix(NA,nrow = runs, ncol = 1,dimnames = list(paste("run",1:runs,sep = ""),"Evenness_eng"))
for(i in 1:runs){
 tree<-geiger::sim.bdtree(b = 0.001,n = 200,stop = "taxa")
 trait<-ape::rTraitCont(ape::compute.brlen(tree, power = power), model = "BM")
 trait<- scales::rescale(trait, to = c(-1, 101)) #re-scaling the trait??
 comm.simul_S3<-sim.comm_2(tree = tree,limit = limit,Ncomm = Ncomm,Nsp = Nsp,richness.equal = richness.equal,filter = filter,gradual = gradual)
 rich<-c(rowSums(comm.simul_S3)+rnorm(n = nrow(comm.simul_S3),mean = 0.0001,sd = 0.001))
 PD<-picante::pd(samp = comm.simul_S3,tree = tree)$PD
 FD<-picante::pd(samp = comm.simul_S3,tree = ape::as.phylo(hclust(dist(trait))))$PD
 matrix.M_S3<-cbind(rich,PD,FD) #matrix M with raw metrics for functional and phylogenetic diversity
 matrix.M_S3_list[[i]]<- matrix.M_S3
 dim_corr_M_S3<-dimensionality(matrix.M = matrix.M_S3,scale = TRUE,method = "standardize",evenness = "Camargo") #baixa dimensionalidade
 IV_M_S3<- ImportanceVal(matrix.M = matrix.M_S3,scale = TRUE,method = method, stopRule = FALSE)
 EE_S3[i,]<-dim_corr_M_S3
 ImportanceValues_S3[[i]]<- IV_M_S3
 print(paste("run",i,sep = ""))
}
scenario_S3<- matrix(NA,nrow= runs, ncol= ncol(ImportanceValues_S3[[1]]$IV.obs),dimnames= list(paste("run",1:runs,sep=""), colnames(ImportanceValues_S3[[1]]$IV.obs)))??

for(i in 1:length(ImportanceValues_S3)){
  scenario_S3[i,]<-colSums(ImportanceValues_S3[[i]]$IV.obs[c(1,2),])
}

scenario_S3<- cbind(scenario_S3,EE_S3) #join EE and IV values
nsuccess_S3<-length(which(scenario_S3[,2]>scenario_S3[,1]&scenario_S3[,3]>scenario_S3[,1]))
mean_scenarioS3<- apply(scenario_S3, 2, mean) #mean values of IV and EE for scenario 3
sd_scenarioS3<- apply(scenario_S3, 2, sd) #mean values of IV and EE for scenario 3

#Gradient graphic for scenario 3
windows()
layout.show(layout(matrix(c(0,0,1,2,3,4,5),nrow=2,ncol=3,byrow=T)))
barplot(decostand(matrix.M_S3[,1],"max")[,1], main= "rich")
barplot(decostand(matrix.M_S3[,2],"max")[,1], main= "PD")
barplot(decostand(matrix.M_S3[,3],"max")[,1], main= "FD")
plot(tree, show.tip.label = FALSE)
barplot(trait[tree$tip.label], horiz=TRUE, space=0, names="")
plotrix::color2D.matplot(t(comm.simul_S3),ylab="",xlab="",main="", yrev=F,cs1=c(1,0),cs2=c(1,0),cs3=c(1,0),border="white",axes=FALSE)

#####Scenario 6 - PD0, FD1, Ric0 - ok#####
limit=0
Ncomm= 50
Nsp= 20
filter= 0.9
gradual= TRUE
richness.equal= TRUE
runs<-999
power<-0.001
method="max"
matrix.M_S6_list<- vector(mode="list",length=runs)
ImportanceValues_S6<-vector(mode = "list",length = runs)
EE_S6<-matrix(NA,nrow = runs, ncol = 1,dimnames = list(paste("run",1:runs,sep = ""),"Evenness_eng"))
for(i in 1:runs){
 tree<-geiger::sim.bdtree(b = 0.001,n = 200,stop = "taxa")
 trait<-ape::rTraitCont(ape::compute.brlen(tree, power = power), model = "BM")
 trait<- scales::rescale(trait, to = c(-1, 101)) #re-scaling the trait??
 comm.simul_S6<-sim.comm_2(tree = tree,limit = limit,Ncomm = Ncomm,Nsp = Nsp,richness.equal = richness.equal,filter = filter,gradual = gradual)
 rich<-c(rowSums(comm.simul_S6)+rnorm(n = nrow(comm.simul_S6),mean = 0.0001,sd = 0.00001))
 PD<-picante::pd(samp = comm.simul_S6,tree = ape::compute.brlen(tree,power=0.0001))$PD #lengthen the terminal branches to reduce the variance
 FD<-picante::pd(samp = comm.simul_S6,tree = ape::as.phylo(hclust(dist(trait))))$PD
 matrix.M_S6<-cbind(rich,PD,FD) #matrix M with raw metrics for functional and phylogenetic diversity
 matrix.M_S6_list[[i]]<- matrix.M_S6
 dim_corr_M_S6<-dimensionality(matrix.M = matrix.M_S6,scale = TRUE,method = "standardize",evenness = "Camargo") #baixa dimensionalidade
 IV_M_S6<- ImportanceVal(matrix.M = matrix.M_S6,scale = TRUE,method = method,stopRule = FALSE)
 EE_S6[i,]<-dim_corr_M_S6
 ImportanceValues_S6[[i]]<- IV_M_S6
 print(paste("run",i,sep = ""))
}

scenario_S6<- matrix(NA,nrow= runs, ncol= ncol(ImportanceValues_S6[[1]]$IV.obs),dimnames= list(paste("run",1:runs,sep=""), colnames(ImportanceValues_S6[[1]]$IV.obs)))
for(i in 1:length(ImportanceValues_S6)){
 scenario_S6[i,]<-colSums(ImportanceValues_S6[[i]]$IV.obs[c(1,2),])
}

scenario_S6<- cbind(scenario_S6,EE_S6) #join EE and IV values
nsuccess_S6<-length(which(scenario_S6[,3]>scenario_S6[,1]&scenario_S6[,3]>scenario_S3[,2]))
mean_scenarioS6<- apply(scenario_S6, 2, mean) #mean values of IV and EE for scenario 6
sd_scenarioS6<- apply(scenario_S6, 2, sd) #sd values of IV and EE for scenario 6

#Figure 1 - setting the problem
quartz()
windows()
plot(decostand(matrix.M_S6_list[[10]][,"PD"],"max"), decostand(matrix.M_S6_list[[10]][,"FD"],"max"), pch= 2, cex= 2)


#Gradient graphic for scenario 6
windows()
layout.show(layout(matrix(c(0,0,1,2,3,4,5),nrow=2,ncol=3,byrow=T)))
barplot(decostand(matrix.M_S6[,1],"max")[,1], main= "rich")
barplot(decostand(matrix.M_S6[,2],"max")[,1], main= "PD")
barplot(decostand(matrix.M_S6[,3],"max")[,1], main= "FD")
plot(ape::compute.brlen(tree,power=0.000000001), show.tip.label = FALSE)
barplot(trait[tree$tip.label], horiz=TRUE, space=0, names="")
plotrix::color2D.matplot(t(comm.simul_S6),ylab="",xlab="",main="", yrev=F,cs1=c(1,0),cs2=c(1,0),cs3=c(1,0),border="white",axes=FALSE)


#########Scenario 8 - PD0, FD0, Ric0#####
limit=0
Ncomm= 50
Nsp= 20
filter= 0.9
gradual= FALSE
richness.equal= TRUE
power<-2
runs<-999
method="max"
matrix.M_S8_list<- vector(mode="list",length=runs)
ImportanceValues_S8<-vector(mode = "list",length = runs)
EE_S8<-matrix(NA,nrow = runs, ncol = 1,dimnames = list(paste("run",1:runs,sep = ""),"Evenness_eng"))
for(i in 1:runs){
 tree<-geiger::sim.bdtree(b = 0.001,n = 200,stop = "taxa")
 newTrait<-runif(n = length(tree$tip.label),min = 100,max = 101)
 names(newTrait)<-tree$tip.label
 comm.simul_S8<-sim.comm_2(tree = tree,limit = limit,Ncomm = Ncomm,Nsp = Nsp,richness.equal = richness.equal,filter = filter,gradual = gradual)
 rich<-c(rowSums(comm.simul_S8)+rnorm(n = nrow(comm.simul_S8),mean = 0,sd = 0.001))
 PD<-picante::pd(samp = comm.simul_S8,tree = ape::compute.brlen(tree,power=0.000000001))$PD+rnorm(n = nrow(comm.simul_S8),mean = 0,sd = 0.001) #null PD
 FD<-picante::pd(samp = comm.simul_S8,tree = ape::as.phylo(hclust(dist(newTrait))))$PD #null FD
 matrix.M_S8<-cbind(rich,PD,FD) #matrix M with raw metrics for functional and phylogenetic diversity
 matrix.M_S8_list[[i]]<- matrix.M_S8
 dim_corr_M_S8<-dimensionality(matrix.M = matrix.M_S8,scale = TRUE,method = "standardize",evenness = "Camargo") #baixa dimensionalidade
 IV_M_S8<- ImportanceVal(matrix.M = matrix.M_S8,scale = TRUE,method = method,stopRule = FALSE)
 EE_S8[i,]<-dim_corr_M_S8
 ImportanceValues_S8[[i]]<- IV_M_S8
 print(paste("run",i,sep = ""))
}
scenario_S8<- matrix(NA,nrow= runs, ncol= ncol(ImportanceValues_S8[[1]]$IV.obs),dimnames= list(paste("run",1:runs,sep=""), colnames(ImportanceValues_S8[[1]]$IV.obs)))
for(i in 1:length(ImportanceValues_S8)){
 scenario_S8[i,]<-colSums(ImportanceValues_S8[[i]]$IV.obs[c(1,2),])
}


scenario_S8<- cbind(scenario_S8,EE_S8) #join EE and IV values
apply(scenario_S8,2,mean) #mean IVs obtained from the two axis??
apply(rich_EqualScenarioOU,2,sd) #IVs standard deviation??
mean(unlist(lapply(matrix.M_S8_list, function(i) cor(i[,1],i[,2]))))
mean_scenarioS8<- apply(scenario_S8, 2, mean) #mean values of IV and EE for scenario 8

#####Scenario 9 - Rich 1, PD 1, FD 0 same correlation#####
limit=0
Ncomm= 50
Nsp= 20
filter= 0.9
gradual= TRUE
richness.equal= FALSE
runs<- 999
method="max"
alpha= 0.8
matrix.M_S9_list<- vector(mode="list",length=runs)
ImportanceValues_S9<-vector(mode = "list",length = runs)
EE_S9<-matrix(NA,nrow = runs, ncol = 1,dimnames = list(paste("run",1:runs,sep = ""),"Evenness_eng"))
for(i in 1:runs){
 tree<- geiger::sim.bdtree(b = 0.001,n = 200,stop = "taxa")
 trait<- ape::rTraitCont(tree, model = "OU", alpha = 0.8) #simulating the traits
 comm.simul_S9<-sim.comm_2(tree = tree,limit = limit,Ncomm = Ncomm,Nsp = Nsp,richness.equal = richness.equal,filter = filter,gradual = gradual)
 rich<-rowSums(comm.simul_S9)
 PD<-picante::pd(samp = comm.simul_S9,tree = tree)$PD
 FD<-picante::pd(samp = comm.simul_S9,tree = ape::as.phylo(hclust(dist(trait))))$PD
 ses.PD<- picante::ses.pd(samp = comm.simul_S9,tree = tree,null.model = "taxa.labels",runs = 10,iterations = 10)$pd.obs.z
 ses.FD<- picante::ses.pd(samp = comm.simul_S9,tree = ape::as.phylo(hclust(dist(trait))),null.model = "taxa.labels",runs = 10,iterations = 10)$pd.obs.z
 matrix.Mses_S9<- cbind(rich,ses.PD,ses.FD) #matrix M with standardized effect size for functional and phylogenetic diversity
 matrix.M_S9<- cbind(rich,PD,FD) #matrix M with raw metrics for functional and phylogenetic diversity
 matrix.M_S9_list[[i]]<- matrix.M_S9
 dim_corr_M_S9<- dimensionality(matrix.M = matrix.M_S9,scale = TRUE,method = "standardize",evenness = "Camargo") #baixa dimensionalidade
  IV_M_S9<- ImportanceVal(matrix.M = matrix.M_S9,scale = TRUE,method = method,stopRule = FALSE)
 EE_S9[i,]<- dim_corr_M_S9
 ImportanceValues_S9[[i]]<- IV_M_S9
 print(paste("run",i,sep = ""))
}

scenario_S9<- matrix(NA,nrow= runs, ncol= ncol(ImportanceValues_S9[[1]]$IV.obs),dimnames= list(paste("run",1:runs,sep=""), colnames(ImportanceValues_S9[[1]]$IV.obs)))
for(i in 1:length(ImportanceValues_S9)){
 scenario_S9[i,]<-colSums(ImportanceValues_S9[[i]]$IV.obs[c(1,2),])
}

scenario_S9<- cbind(scenario_S9,EE_S9) #join EE and IV values
mean(unlist(lapply(matrix.M_S9_list, function(i) cor(i[,1],i[,2]))))
mean_scenarioS9<- apply(scenario_S9, 2, mean) #mean values of IV and EE for scenario 8
sd_scenarioS9<- apply(scenario_S9, 2, sd)

#gradient graphic for scenario 9
windows()
layout.show(layout(matrix(c(0,0,1,2,3,4,5),nrow=2,ncol=3,byrow=T)))
barplot(decostand(matrix.M_S9[,1],"max")[,1], main= "rich")
barplot(decostand(matrix.M_S9[,2],"max")[,1], main= "PD")
barplot(decostand(matrix.M_S9[,3],"max")[,1], main= "FD")
plot(tree, show.tip.label = FALSE)
barplot(trait[tree$tip.label], horiz=TRUE, space=0, names="")
plotrix::color2D.matplot(t(comm.simul_S9),ylab="",xlab="",main="", yrev=F,cs1=c(1,0),cs2=c(1,0),cs3=c(1,0),border="white",axes=FALSE)

#Fig 1 - setting the problem
windows()
quartz()
plot(decostand(matrix.M_S9_list[[1]][,"PD"], method= "max"),decostand(matrix.M_S9_list[[1]][,"FD"], method= "max"), pch=2, cex=2)

######Scenario 10 - rich 1, PD 1 and FD1 no correlation########
limit=0
Ncomm= 50
Nsp= 20
filter= 0.9
gradual= TRUE
richness.equal= TRUE
runs<-999
method="max"
power<-0.000000000001
matrix.M_S10_list<- vector(mode="list",length=runs)
ImportanceValues_S10<-vector(mode = "list",length = runs)
EE_S10<-matrix(NA,nrow = runs, ncol = 1,dimnames = list(paste("run",1:runs,sep = ""),"Evenness_eng"))
for(i in 1:runs){
 tree<- geiger::sim.bdtree(b = 0.001,n = 200,stop = "taxa")
 trait<- ape::rTraitCont(ape::compute.brlen(tree, power = power), model = "BM") #simulating the traits
 comm.simul_S10<-sim.comm(tree,SimClus=0.9,SimOver=0,limit=0,Ncomm=Ncomm,Nind=30,Nsp= Nsp,filter=0.99,stop="richness",LNDist=FALSE,pattern="random",gradual=TRUE) #assembly of communities by phylogeny
 comm.simul_S10<-ifelse(comm.simul_S10$ResSim>=1,1,0)
 rich<-c(rowSums(comm.simul_S10)+rnorm(n = nrow(comm.simul_S10),mean = 0.0001,sd = 0.001))
 PD<-picante::pd(samp = comm.simul_S10,tree = tree)$PD
 FD<-picante::pd(samp = comm.simul_S10,tree = ape::as.phylo(hclust(dist(trait))))$PD
 matrix.M_S10<- cbind(rich,PD,FD) #matrix M with raw metrics for functional and phylogenetic diversity
 matrix.M_S10_list[[i]]<- matrix.M_S10
 dim_corr_M_S10<- dimensionality(matrix.M = matrix.M_S10,scale = TRUE,method = "standardize",evenness = "Camargo") #baixa dimensionalidade
 IV_M_S10<- ImportanceVal(matrix.M = matrix.M_S10,scale = TRUE,method = method,stopRule = FALSE)
 #IV_Mses_S1<- ImportanceVal(matrix.M = matrix.Mses_S1[c(-1,-2,-3),],IV.bootstrap = FALSE,scale = TRUE,method = "standardize",stopRule = TRUE)
 EE_S10[i,]<- dim_corr_M_S10
 ImportanceValues_S10[[i]]<- IV_M_S10
 print(paste("run",i,sep = ""))
}
scenario_S10<- matrix(NA,nrow= runs, ncol= ncol(ImportanceValues_S10[[1]]$IV.obs),dimnames= list(paste("run",1:runs,sep=""), colnames(ImportanceValues_S10[[1]]$IV.obs)))??

for(i in 1:length(ImportanceValues_S1)){
 scenario_S10[i,]<-colSums(ImportanceValues_S10[[i]]$IV.obs[c(1,2),])
}

scenario_S10<- cbind(scenario_S10,EE_S10) #join EE and IV values
nsuccess_S10<-length(which(scenario_S10[,2]>scenario_S10[,3])) #number of sucessess
mean_scenarioS10<- apply(scenario_S10, 2, mean) # mean values of IV and EE for scenario 10
sd_scenarioS10<- apply(scenario_S10, 2, sd) # sd values of IV and EE for scenario 10
mean(unlist(lapply(matrix.M_S10_list, function(i) cor(i[,2],i[,3])))) #media de correlacao entre PD e FD

#Fig 1 - setting the problem - Scenario without correlation (Scenario S2 diff variation)#########
windows()
quartz()
plot(decostand(matrix.M_S10_list[[10]][,"PD"], method= "max"),decostand(matrix.M_S10_list[[10]][,"FD"], method= "max"), pch=19, cex=3, ylim=c(0.6,1))
points(decostand(matrix.M_S6_list[[10]][,"PD"], method= "max"),decostand(matrix.M_S6_list[[10]][,"FD"], method= "max"), pch=2, cex=3)
which(scenario_S10[,3]==min(scenario_S10[,3]))



######Scenario 11 - PD1, FD1, Ric0 - different importance, low correlation - similar scenario 10#####
limit=0
Ncomm= 50
Nsp= 20
filter= 0.9
gradual= TRUE
richness.equal= TRUE
runs<-999
power<-0.001
method="max"
matrix.M_S11_list<- vector(mode="list",length=runs)
ImportanceValues_S11<-vector(mode = "list",length = runs)
EE_S11<-matrix(NA,nrow = runs, ncol = 1,dimnames = list(paste("run",1:runs,sep = ""),"Evenness_eng"))
for(i in 1:runs){
 tree<-geiger::sim.bdtree(b = 0.001,n = 200,stop = "taxa")
 trait<-ape::rTraitCont(ape::compute.brlen(tree, power = power), model = "BM")
 trait<- scales::rescale(trait, to = c(-1, 101)) #re-scaling the trait?
 comm.simul_S11<-sim.comm_2(tree = tree,limit = limit,Ncomm = Ncomm,Nsp = Nsp,richness.equal = richness.equal,filter = filter,gradual = gradual)
 rich<-c(rowSums(comm.simul_S11)+rnorm(n = nrow(comm.simul_S11),mean = 0.0001,sd = 0.001))
 PD<-picante::pd(samp = comm.simul_S11,tree = tree)$PD
 FD<-picante::pd(samp = comm.simul_S11,tree = ape::as.phylo(hclust(dist(trait))))$PD
 matrix.M_S11<-cbind(rich,PD,FD) #matrix M with raw metrics for functional and phylogenetic diversity
 matrix.M_S11_list[[i]]<- matrix.M_S11
 dim_corr_M_S11<-dimensionality(matrix.M = matrix.M_S11,scale = TRUE,method = "standardize",evenness = "Camargo") #baixa dimensionalidade
 IV_M_S11<- ImportanceVal(matrix.M = matrix.M_S11,scale = TRUE,method = method, stopRule = FALSE)
 EE_S11[i,]<-dim_corr_M_S11
 ImportanceValues_S11[[i]]<- IV_M_S11
 print(paste("run",i,sep = ""))
}
scenario_S11<- matrix(NA,nrow= runs, ncol= ncol(ImportanceValues_S11[[1]]$IV.obs),dimnames= list(paste("run",1:runs,sep=""), colnames(ImportanceValues_S11[[1]]$IV.obs)))??

for(i in 1:length(ImportanceValues_S11)){
 scenario_S11[i,]<-colSums(ImportanceValues_S11[[i]]$IV.obs[c(1,2),])
}


scenario_S11<- cbind(scenario_S11,EE_S11) #join EE and IV values
nsuccess_S11<-length(which(scenario_S11[,2]>scenario_S11[,1]&scenario_S11[,3]>scenario_S11[,1]))
mean_scenarioS11<- apply(scenario_S11, 2, mean) #mean values of IV and EE for scenario 3







