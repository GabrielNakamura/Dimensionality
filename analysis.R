################################
#IV calculation on small mammals 
################################


####libraries and functions####
library(ape)
library(vegan)
library(picante)
library(EcoSimR)
library(ade4)
library(FD)
source(here::here("functions", "ImportanceVal_V2.R"))
source(here::here("functions ", "dimensionality_modific.R"))
source(here::here("functions","test_TryCatch_metricsFunc_27-7-2019.R")) #TryCatch function with Functional metrics separated
source(here::here("functions", "test_TryCatch_metricsPhylo_27-7-2019.R")) #TryCatch function with Functional metrics separated
source(here::here("functions","test_TryCatch_metrics_27-7-2019.R")) #TryCatch function for all metrics
source(here::here("functions", "Camargo_function.R"))


###reading data####
comm<-read.table(here::here("data","comm.txt"),header=T)
phy<-read.tree(here::here("data","tree.txt"))
par(mar=c(1,1,1,1))
plot.phylo(phy, type = "fan", cex= 0.3)
axisPhylo(1)
traits_complete<- read.csv(here::here("data","traits_completo.csv"), header = T, sep = ";", row.names = 1)
traits<-read.table(here::here("data", "traits.txt"),header=T)
composition_coords<- read.csv(here::here("data", "points_composition.csv"), sep = ";",header=T)
coords<- composition_coords[,c(1,2)]
abs(coords[1,1]- coords[-1,1])

###preparing functional traits#####
quant_traits<- data.frame(t(traits_complete[1:3,]))
bin_traits<- t(traits_complete[4:nrow(traits_complete),])
diet<- prep.binary(data.frame(bin_traits[,5:ncol(bin_traits)]), col.blocks= 8, label="dieta")
locomotion<- prep.binary(data.frame(bin_traits[,1:4]), col.blocks= 4, label="locomotion")
ktab1 <- ktab.list.df(list(quant_traits, diet, locomotion))
distrait <- dist.ktab(ktab1, c("Q", "B", "B"), c("scaledBYrange"))
funct_dendro<- as.phylo(hclust(distrait, method = "average"))
distrait_orig<- as.matrix(distrait) 
rownames(distrait_orig)<- colnames(comm)
colnames(distrait_orig)<- colnames(comm)
distrait_orig<- as.dist(distrait_orig) #trait distance matrix


###calculating metrics####
#phylo metrics
PDfaith<-pd(comm,phy)$PD #phylo diversity
mntd<- mntd(samp = comm, dis = cophenetic(phy))
PSV<- psv(samp = comm,tree = phy,compute.var=TRUE,scale.vcv=TRUE)$PSVs
#PSE<- pse(samp = comm, tree = phy, scale.vcv = TRUE)$PSEs #helmus pse
match(rownames(cophenetic(phy)), colnames(comm))
DBPhylo<- dbFD(x = cophenetic(phy), a = comm[,match(rownames(cophenetic(phy)), colnames(comm))], calc.FRic = T, w.abun = FALSE, calc.FDiv = TRUE, calc.CWM = FALSE, calc.FGR = FALSE) #phylogenetic db measures (Vill?ger)
Peve<- DBPhylo$FEve
Peve[which(is.na(Peve))]<- 0 #Phylogenetic evenness
#functional metrics
FDfaith<-pd(comm,FD_tree)$PD #func diversity
DBFunc<- dbFD(x = distrait, a = comm, calc.FRic = T, w.abun = FALSE, calc.FDiv = TRUE, calc.CWM = FALSE, calc.FGR = FALSE)
FEve<- DBFunc$FEve #Functional evenness
FEve[which(is.na(FEve))]<- 0
FDiv<- DBFunc$FDiv #Functional divergence
FDiv[which(is.na(FDiv))]<- 0
#taxonomic metric
rich<- rowSums(comm)


##Matrix M with different configurations#####
matrix_Mobs<- as.matrix(data.frame(PDfaith, mntd, PSV, Peve, FDfaith, FEve, FDiv, richness= rich))
matrix_Mobs_Func<- matrix_Mobs[,c("FDfaith","FEve","FDiv","richness")]
matrix_Mobs_Phylo<- matrix_Mobs[,c("PDfaith","mntd","PSV","Peve","richness")]
matrix_Mobs_PhyloFunc<- matrix_Mobs[,c("PDfaith","mntd","PSV","Peve","FDfaith","FEve","FDiv")]

#calculating IV metric
IVobs_mammals<-ImportanceVal(matrix.M = matrix_Mobs, scale= TRUE, method= "max", stopRule= TRUE) #observed IV
IV_total<- sort(apply(IVobs_mammals$IV.obs_stopRule, 2, sum))
IV_percentage<- sort(IV_total/sum(IV_total))
IV_boot_list<- vector(mode= "list", length= 999) # list with bootstrap IVs for all metrics
IV_boot_listFunc<- vector(mode= "list", length= 999) # list with bootstrap IVs for functional metrics
IV_boot_listPhylo<- vector(mode= "list", length= 999) # list with bootstrap IVs for phylogenetic metrics
IV_boot_listPhyloFunc<- vector(mode= "list", length= 999) # list with bootstrap IVs for phylo and functional metrics
IVobs_mammals_Func<- ImportanceVal(matrix.M = matrix_Mobs_Func, scale= TRUE, method= "max", stopRule= TRUE) #observed IV for functional matrix M
IV_total_Func<- sort(apply(IVobs_mammals_Func$IV.obs_stopRule,2,sum))
IV_percentage_Func<- sort(IV_total_Func/sum(IV_total_Func))
IVobs_mammals_Phylo<- ImportanceVal(matrix.M = matrix_Mobs_Phylo, scale= TRUE, method= "max", stopRule= TRUE) #observed IV for functional matrix M
IV_total_Phylo<- sort(apply(IVobs_mammals_Phylo$IV.obs_stopRule,2,sum))
IV_percentage_Phylo<- sort(IV_total_Phylo/sum(IV_total_Phylo))
IVobs_mammals_PhyloFunc<- ImportanceVal(matrix.M = matrix_Mobs_PhyloFunc, scale= TRUE, method= "max", stopRule= TRUE) #observed IV for functional matrix M
IV_total_PhyloFunc<- sort(apply(IVobs_mammals_PhyloFunc$IV.obs_stopRule,2,sum))
IV_percentage_PhyloFunc<- sort(IV_total_PhyloFunc/sum(IV_total_PhyloFunc))

#bootstrap IVs for all metrics
for(i in 1:999){
  IV_boot_list[[i]]<- ImportanceVal(matrix_Mobs[sample(1:nrow(matrix_Mobs), size = nrow(matrix_Mobs), replace = TRUE), ], scale=T, method= "max", 
                stopRule = TRUE)
}
IV_boot_all<- matrix(unlist(lapply(IV_boot_list, function(i) {
  if(is.matrix(i$IV.obs_stopRule)){
    colSums(i$IV.obs_stopRule)
  } else{
    i$IV.obs_stopRule
  }
})), nrow= 999, ncol= ncol(matrix_Mobs), byrow = T, dimnames= list(paste("boot", 1:999, sep= ""), colnames(matrix_Mobs))) #matrix with bootstrap IVs

IV_confint_all<- apply(IV_boot_all,MARGIN = 2, function(i) quantile(sort(i), c(0.025,0.975))) #confidence interval for bootstrap IVs
IV_meanboot_all<- apply(IV_boot_all,MARGIN = 2, function(i) mean(i)) #mean IV for bootstrap IVs
(sort(IV_meanboot_all)/sum(IV_meanboot_all))
#bootstrap Functional IVs
for(i in 1:999){
  IV_boot_listFunc[[i]]<- ImportanceVal(matrix_Mobs_Func[sample(1:nrow(matrix_Mobs_Func), size = nrow(matrix_Mobs_Func), replace = TRUE), ], scale=T, method= "max", 
                                        stopRule = TRUE)
}
IV_boot_Func<- matrix(unlist(lapply(IV_boot_listFunc, function(i) {
  if(is.matrix(i$IV.obs_stopRule)){
    colSums(i$IV.obs_stopRule)
  } else{
    i$IV.obs_stopRule
  }
  })), nrow= 999, ncol= ncol(matrix_Mobs_Func), byrow = T, dimnames= list(paste("boot", 1:999, sep= ""), colnames(matrix_Mobs_Func))) #matrix with bootstrap IVs for functional metrics
IV_confint_Func<- apply(IV_boot_Func,MARGIN = 2, function(i) quantile(sort(i), c(0.025,0.975))) #confidence interval for bootstrap IVs
IV_meanboot_Func<- apply(IV_boot_Func,MARGIN = 2, function(i) mean(i)) #mean IV for bootstrap IVs
sort(IV_meanboot_Func)/sum(IV_meanboot_Func)

#bootstrap Phylo IVs
for(i in 1:999){
  IV_boot_listPhylo[[i]]<- ImportanceVal(matrix_Mobs_Phylo[sample(1:nrow(matrix_Mobs_Phylo), size = nrow(matrix_Mobs_Phylo), replace = TRUE), ], scale=T, method= "max", 
                                        stopRule = TRUE)
}
IV_boot_Phylo<- matrix(unlist(lapply(IV_boot_listPhylo, function(i) {
  if(is.matrix(i$IV.obs_stopRule)){
    colSums(i$IV.obs_stopRule)
  } else{
    i$IV.obs_stopRule
  }
})), nrow= 999, ncol= ncol(matrix_Mobs_Phylo), byrow = T, dimnames= list(paste("boot", 1:999, sep= ""), colnames(matrix_Mobs_Phylo))) #matrix with bootstrap IVs for functional metrics
IV_confint_Phylo<- apply(IV_boot_Phylo,MARGIN = 2, function(i) quantile(sort(i), c(0.025,0.975))) #confidence interval for bootstrap IVs
IV_meanboot_Phylo<- apply(IV_boot_Phylo,MARGIN = 2, function(i) mean(i)) #mean IV for bootstrap IVs
sort(IV_meanboot_Phylo)/sum(IV_meanboot_Phylo)
#bootstrap Phylo and func IVs
for(i in 1:999){
  IV_boot_listPhyloFunc[[i]]<- ImportanceVal(matrix_Mobs_PhyloFunc[sample(1:nrow(matrix_Mobs_PhyloFunc), size = nrow(matrix_Mobs_PhyloFunc), replace = TRUE), ], scale=T, method= "max", 
                                         stopRule = TRUE)
}
IV_boot_PhyloFunc<- matrix(unlist(lapply(IV_boot_listPhyloFunc, function(i) {
  if(is.matrix(i$IV.obs_stopRule)){
    colSums(i$IV.obs_stopRule)
  } else{
    i$IV.obs_stopRule
  }
})), nrow= 999, ncol= ncol(matrix_Mobs_PhyloFunc), byrow = T, dimnames= list(paste("boot", 1:999, sep= ""), colnames(matrix_Mobs_PhyloFunc))) #matrix with bootstrap IVs for functional metrics
IV_confint_PhyloFunc<- apply(IV_boot_PhyloFunc,MARGIN = 2, function(i) quantile(sort(i), c(0.025,0.975))) #confidence interval for bootstrap IVs
IV_meanboot_PhyloFunc<- apply(IV_boot_PhyloFunc,MARGIN = 2, function(i) mean(i)) #mean IV for bootstrap IVs
sort(IV_meanboot_PhyloFunc)/sum(IV_meanboot_PhyloFunc)
######hypothesis testing with bootstrap values for IVs#######
IV_boot_dataFrame<- data.frame(IV_values=c(IV_boot[,1:ncol(IV_boot)]),dimensions= rep(colnames(IV_boot), each=999), 
                               dimensions_general= rep(c("Phylo","Phylo","Phylo","Phylo","Func","Func","Func","tax"), each=999)) #data frame for linear model with IVs~dimensions
mod_aovIV<- aov(IV_values~dimensions, data= IV_boot_dataFrame) #linear model ANOVA type
mod_aovIV_2<- aov(IV_values~dimensions_general, data= IV_boot_dataFrame) #linear model ANOVA type merging dimensions
summary(mod_aovIV)
summary(mod_aovIV_2)
quartz()
plot(TukeyHSD(mod_aovIV), cex.axis=0.6, las= 1) #tukey test
quartz()
plot(TukeyHSD(mod_aovIV_2)) #tukey test for model 2
quartz()
par(mfrow= c(2, 2))
plot(mod_aovIV, which= c(1:4)) #model checking 

#####ploting IV profile for mammals#####

quartz()
par(las=2, mar=c(4,8,3,1))
barplot(height = sort(IV_meanboot),     
        ylim = c(0.0, 0.25),
        xlim= c(-0.2, 10),
        xaxs="i",
        yaxs="i",
        axis.lty=1,
        col= "white",
        axisnames= FALSE,
        axes= FALSE        
)
axis(side= 1, at= c(-0.2, 0.7, 1.9, 3.1, 4.3, 5.5, 6.7, 8, 9.1), 
     labels= (c("", "FDiv", "mntd", "FEve", "FD", "PEve", "PD", "Rich", "PSV" )))
par(las= 1)
axis(side= 2, at= axTicks(side= 2, axp= c(0.0,0.25,10)))
arrows(x0 = c(0.7, 1.9, 3.1, 4.3, 5.5, 6.7, 8, 9.1), y0 = IV_confint[1,names(sort(IV_meanboot))], 
       x1= c(0.7, 1.9, 3.1, 4.3, 5.5, 6.7, 8, 9.1), y1= IV_confint[2,names(sort(IV_meanboot))],
       code= 3, angle= 90, length = 0.1)

#Functional and taxonomic IVs
quartz()
par(las=2, mar=c(4,8,3,1))
barplot(height = sort(IV_meanboot_Func),     
        ylim = c(0.0, 0.50),
        xlim= c(-0.2, 5),
        xaxs="i",
        yaxs="i",
        axis.lty=1,
        col= "white",
        axisnames= FALSE,
        axes= FALSE        
)
axis(side= 1, at= c(-0.2, 0.7, 1.9, 3.1, 4.3), 
     labels= (c("", "FD", "FEve", "FDiv", "Rich" )))
par(las= 1)
axis(side= 2, at= axTicks(side= 2, axp= c(0.0,0.5,10)))
arrows(x0 = c(0.7, 1.9, 3.1, 4.3), y0 = IV_confint_Func[1,names(sort(IV_meanboot_Func))], 
       x1= c(0.7, 1.9, 3.1, 4.3), y1= IV_confint_Func[2,names(sort(IV_meanboot_Func))],
       code= 3, angle= 90, length = 0.1)

#Phylo and richness IVs
quartz()
par(las=2, mar=c(4,8,3,1))
barplot(height = sort(IV_meanboot_Phylo),     
        ylim = c(0.0, 0.40),
        xlim= c(-0.2, 6.5),
        xaxs="i",
        yaxs="i",
        axis.lty=1,
        col= "white",
        axisnames= FALSE,
        axes= FALSE        
)
axis(side= 1, at= c(-0.2, 0.7, 1.9, 3.1, 4.3,5.5), 
     labels= (c("", "mntd", "PEve", "PD", "Rich", "PSV" )))
par(las= 1)
axis(side= 2, at= axTicks(side= 2, axp= c(0.0,0.4,10)))
arrows(x0 = c(0.7, 1.9, 3.1, 4.3, 5.5), y0 = IV_confint_Phylo[1,names(sort(IV_meanboot_Phylo))], 
       x1= c(0.7, 1.9, 3.1, 4.3, 5.5), y1= IV_confint_Phylo[2,names(sort(IV_meanboot_Phylo))],
       code= 3, angle= 90, length = 0.1)

#Phylo and richness IVs
quartz()
par(las=2, mar=c(4,8,3,1))
barplot(height = sort(IV_meanboot_PhyloFunc),     
        ylim = c(0.0, 0.30),
        xlim= c(-0.2, 9),
        xaxs="i",
        yaxs="i",
        axis.lty=1,
        col= "white",
        axisnames= FALSE,
        axes= FALSE        
)
axis(side= 1, at= c(-0.2, 0.7, 1.9, 3.1, 4.3, 5.5, 6.7, 8), 
     labels= (c("", "FDiv", "mntd", "FEve", "FD", "PD", "Peve", "PSV" )))
par(las= 1)
axis(side= 2, at= axTicks(side= 2, axp= c(0.0,0.3,10)))
arrows(x0 = c(0.7, 1.9, 3.1, 4.3, 5.5, 6.7, 8), y0 = IV_confint_PhyloFunc[1,names(sort(IV_meanboot_PhyloFunc))], 
       x1= c(0.7, 1.9, 3.1, 4.3, 5.5, 6.7, 8), y1= IV_confint_PhyloFunc[2,names(sort(IV_meanboot_PhyloFunc))],
       code= 3, angle= 90, length = 0.1)
####calculating observed and bootstrap EE metric#####
EEobs_mammals<- dimensionality(matrix_Mobs, scale= TRUE, method= "standardize", evenness= "Camargo") #observed EE
EEobsAll_boot<- numeric(length= 999) # list with bootstrap IVs
sd(EEobsAll_boot)
EEobsFunc_boot<- numeric(length= 999) # list with bootstrap IVs
EEobsPhylo_boot<- numeric(length= 999) # list with bootstrap IVs
EEobsPhyloFunc_boot<- numeric(length= 999) # list with bootstrap IVs
for(i in 1:999){
  EEobsAll_boot[i]<- dimensionality(matrix.M = matrix_Mobs[sample(1:nrow(matrix_Mobs), size = nrow(matrix_Mobs), replace = TRUE), ], scale=T, method = "standardize", evenness = "Camargo")
}
EE_confint_All<- quantile(unlist(EEobsAll_boot), c(0.025,0.975)) #confidence interval for bootstrap IVs
EE_meanboot_All<- mean(unlist(EEobsAll_boot))

EEobs_matrixMFunc<- dimensionality(matrix_Mobs_Func, scale= T, method= "standardize", evenness= "Camargo") #EE for matrix M with functional metrics
for(i in 1:999){
  EEobsFunc_boot[i]<- dimensionality(matrix.M = matrix_Mobs_Func[sample(1:nrow(matrix_Mobs), size = nrow(matrix_Mobs), replace = TRUE), ], scale=T, method = "standardize", evenness = "Camargo")
}
EE_confint_Func<- quantile(unlist(EEobsFunc_boot), c(0.025,0.975)) #confidence interval for bootstrap IVs
EE_meanboot_Func<- mean(unlist(EEobsFunc_boot))

EEobs_matrixMPhylo<- dimensionality(matrix_Mobs_Phylo, scale= T, method= "standardize", evenness= "Camargo") #EE for matrix M with functional metrics
for(i in 1:999){
  EEobsPhylo_boot[i]<- dimensionality(matrix.M = matrix_Mobs_Phylo[sample(1:nrow(matrix_Mobs), size = nrow(matrix_Mobs), replace = TRUE), ], scale=T, method = "standardize", evenness = "Camargo")
}
EE_confint_Phylo<- quantile(unlist(EEobsPhylo_boot), c(0.025,0.975)) #confidence interval for bootstrap IVs
EE_meanboot_Phylo<- mean(unlist(EEobsPhylo_boot))

EEobs_matrixMPhyloFunc<- dimensionality(matrix_Mobs_PhyloFunc, scale= T, method= "standardize", evenness= "Camargo") #EE for matrix M with functional metrics
for(i in 1:999){
  EEobsPhyloFunc_boot[i]<- dimensionality(matrix.M = matrix_Mobs_PhyloFunc[sample(1:nrow(matrix_Mobs), size = nrow(matrix_Mobs_PhyloFunc), replace = TRUE), ], scale=T, method = "standardize", evenness = "Camargo")
}
EE_confint_PhyloFunc<- quantile(unlist(EEobsPhyloFunc_boot), c(0.025,0.975)) #confidence interval for bootstrap IVs
EE_meanboot_PhyloFunc<- mean(unlist(EEobsPhyloFunc_boot))

####ploting EE and IV equitability####

IV_booteveness_all<- apply(IV_boot_all, 1, function(i) camargo.eveness(n_spec = i, include_zeros = TRUE))
IV_booteveness_phylo<- apply(IV_boot_Phylo, 1, function(i) camargo.eveness(n_spec = i, include_zeros = TRUE))
IV_booteveness_func<- apply(IV_boot_Func, 1, function(i) camargo.eveness(n_spec = i, include_zeros = TRUE))
IV_booteveness_funcphylo<- apply(IV_boot_PhyloFunc, 1, function(i) camargo.eveness(n_spec = i,include_zeros = TRUE))
mean_bootEqIV_all<- mean(IV_booteveness_all)
sd(IV_booteveness_all)
mean_bootEqIV_phylo<- mean(IV_booteveness_phylo)
mean_bootEqIV_func<- mean(IV_booteveness_func)
mean_bootEqIV_phylofunc<- mean(IV_booteveness_funcphylo)
IV_allBoot<- matrix(c(mean_bootEqIV_all,quantile(IV_booteveness_all, c(0.025,0.975)),
         mean_bootEqIV_func,quantile(IV_booteveness_func, c(0.025,0.975)),
         mean_bootEqIV_phylo,quantile(IV_booteveness_phylo, c(0.025,0.975)),
         mean_bootEqIV_phylofunc,quantile(IV_booteveness_funcphylo, c(0.025,0.975))),nrow = 4, ncol = 3, byrow=TRUE, 
       dimnames= list(c("All","Func","Phylo","FuncPhylo"),c("mean","0.025","0.975")))
EE_allboot<- matrix(c(EE_meanboot_All,EE_confint_All,
                      EE_meanboot_Func,EE_confint_Func,
                      EE_meanboot_Phylo,EE_confint_Phylo,
                      EE_meanboot_PhyloFunc,EE_confint_PhyloFunc),nrow = 4, ncol = 3, byrow=TRUE, 
                    dimnames= list(c("All","Func","Phylo","FuncPhylo"),c("mean","0.025","0.975")))
quartz()
plot(EE_allboot[,1],IV_allBoot[,1],xlab = "EE", ylab= "IV equitability", xlim=c(0,1), ylim=c(0,1), pch=c(19,15,10,7), cex=2)
abline(v=0.5,h=0.5,lty=2)
arrows(x0 = EE_allboot[,2], x1 = EE_allboot[,3], y0 = IV_allBoot[,1],y1 = IV_allBoot[,1],code= 3, angle= 90, length = 0.1)
arrows(x0 = EE_allboot[,1], x1 = EE_allboot[,1], y0 = IV_allBoot[,2],y1 = IV_allBoot[,3],code= 3, angle= 90, length = 0.1)
legend(x = "topleft", legend= c("All metrics", "Func+rich", "Phylo+rich", "Phylo/Func"), pch = (c(19,15,10,7)))

###null EE metric model 1####
null_comm_model1<- replicate(999, sim3(t(comm))) #calculating null model 1 from Stevens and Tello (Ecography 2018) fixed richness
list_null_commMod1<- vector(mode="list", length= 999)

for(i in 1:length(list_null_commMod1)){
  comm_null<- t(null_comm_model1[,,i])
  colnames(comm_null)<- colnames(comm)
  list_null_commMod1[[i]]<- comm_null
} #null communities organized in list

list_matrixMobs_nullMod1_Phylo<- lapply(list_null_commMod1, function(i){ #null matrix M calculated according to null model 1
  matrixM_error_Phylo(i)
  } 
  ) #calculating null matrix M only with Phylo metrics ignoring error in dbFD function

list_matrixMobs_nullMod1_Func<- lapply(list_null_commMod1, function(i){ #null matrix M calculated according to null model 1
  matrixM_error_Func(i)
} 
) #calculating null matrix M only with Phylo metrics ignoring error in dbFD function

list_matrixMobs_nullMod1_all<- lapply(list_null_commMod1, function(i){ #null matrix M calculated according to null model 1
  matrixM_error_all(i)
} 
)


list_matrixMobs_nullMod1_all_noNA<- Filter(function(a) any(!is.na(a)), list_matrixMobs_nullMod1)
list_matrixMobs_nullMod1_Func_noNA<- Filter(function(a) any(!is.na(a)), list_matrixMobs_nullMod1_Func)
list_matrixMobs_nullMod1_Phylo_noNA<- Filter(function(a) any(!is.na(a)), list_matrixMobs_nullMod1_Phylo)
EEnullMod1_all<- lapply(list_matrixMobs_nullMod1_all_noNA, function(i) dimensionality(i, scale= T, method= "standardize", evenness= "Camargo")) #null values of EE metric according to null model 1
EEnullMod1_Func<- lapply(list_matrixMobs_nullMod1_Func_noNA, function(i) dimensionality(i, scale= T, method= "standardize", evenness= "Camargo")) #null values of EE metric according to null model 1
EEnullMod1_Phylo<- lapply(list_matrixMobs_nullMod1_Phylo_noNA, function(i) dimensionality(i, scale= T, method= "standardize", evenness= "Camargo")) #null values of EE metric according to null model 1

####null EE metric model 2####
list_null_commMod2<- vector(mode="list", length= 999)
for(i in 1:length(list_null_commMod2)){ #calculating null model 2 from Stevens and Tello (GEB 2014) fixed occurence and richness
  list_null_commMod2[[i]]<-sim9(t(comm), sim9, metric="c_score", nReps = 1000, saveSeed = FALSE,
                                burn_in = 0, suppressProg = TRUE)$Randomized.Data
}

list_null_commMod2_transp<- lapply(list_null_commMod2, function(i) t(i))
list_matrixMobs_nullMod2<- lapply(list_null_commMod2_transp, function(i){ #null matrix M calculated according to null model 2
  matrixM_error_all(i)
  } 
  ) #calculating null matrix M with all metrics 

list_matrixMobs_nullMod2_Phylo<- lapply(list_null_commMod2_transp, function(i){ #null matrix M calculated according to null model 2 only with Phylo metrics
  matrixM_error_Phylo(i)
} 
) #calculating null matrix M with Phylo metrics 

list_matrixMobs_nullMod2_Func<- lapply(list_null_commMod2_transp, function(i){ #null matrix M calculated according to null model 2 only with Phylo metrics
  matrixM_error_Func(i)
} 
) #calculating null matrix M with Func metrics 

EEnullMod2<- lapply(list_matrixMabs_nullMod2, function(i) dimensionality(i, scale= T, method= "standardize", evenness= "Camargo")) #null values of EE metric according to null model 2

####Graphics EE null models####
quartz()
hist(unlist(EEnullMod1), xlab= "EEs", ylab= "EEs frequency", cex= 1.3, main= "null model 1") #histogram of null EE according to null model 1
abline(v=quantile(sort(unlist(EEnullMod1)),c(0.025,0.975)),lty=3, lwd= 1.3) #observed 95% interval in histogram of null EE
abline(v=EEobs_mammals,lty=1, lwd= 1.3)
quartz()
hist(unlist(EEnullMod2), xlab= "EEs", ylab= "EEs frequency", cex= 1.3, main= "null model 2") #histogram of null EE according to null model 1
abline(v=quantile(sort(unlist(EEnullMod2)),c(0.025,0.975)),lty=3, lwd= 1.3) #observed 95% interval in histogram of null EE
abline(v=EEobs_mammals,lty=1, lwd= 1.3)



