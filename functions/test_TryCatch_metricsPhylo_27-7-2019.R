###TryCatch with Phylo metrics
matrixM_error_Phylo<- function(data){
  out<- tryCatch({
    comm_ord<- data[,match(rownames(as.matrix(cophenetic(phy))),colnames(data))]
    comm_ord_trait<- data[,match(rownames(as.matrix(distrait_orig)),colnames(data))]
    if(length(which(colSums(data)==0))>0){
      distrait<- distrait_orig
      matrix_phy<- as.matrix(cophenetic(phy))[-which(colSums(comm_ord)==0),-which(colSums(comm_ord)==0)]
      dist_phy<- dist(matrix_phy)
      distrait<- as.matrix(distrait)[-which(colSums(comm_ord_trait)==0),-which(colSums(comm_ord_trait)==0)]
      rownames(distrait)<- colnames(data)[-which(colSums(comm_ord_trait)==0)]
      colnames(distrait)<- colnames(data)[-which(colSums(comm_ord_trait)==0)]
      distrait<- as.dist(distrait)
      comm_ord<- comm_ord[,-which(colSums(comm_ord)==0)]
      comm_ord_trait<- comm_ord_trait[,-which(colSums(comm_ord_trait)==0)]
    }else{
      dist_phy<- cophenetic(phy)
      distrait<- as.matrix(distrait_orig)
      rownames(distrait)<- colnames(data)
      colnames(distrait)<- colnames(data)
      distrait<- as.dist(distrait)
    }
    PDfaith<-pd(data,phy)$PD #phylo diversity
    mntd<- mntd(samp = data, dis = cophenetic(phy))
    PSV<- psv(samp = data,tree = phy,compute.var=TRUE,scale.vcv=TRUE)$PSVs
    DBPhylo<- dbFD(x = dist_phy, a = comm_ord, calc.FRic = T, w.abun = FALSE, calc.FDiv = TRUE, calc.CWM = FALSE, calc.FGR = FALSE) #phylogenetic db measures (Vill?ger)
    #DBFunc<- dbFD(x = distrait, a = comm_ord_trait, calc.FRic = T, w.abun = FALSE, calc.FDiv = TRUE, calc.CWM = FALSE, calc.FGR = FALSE)
    Peve<- DBPhylo$FEve
    Peve[which(is.na(Peve))]<- 0 #Phylogenetic evenness
    #functional metrics
    #FDfaith<-pd(data,FD_tree)$PD #func diversity
    #FEve<- DBFunc$FEve #Functional evenness
    #FEve[which(is.na(FEve))]<- 0
    #FDiv<- DBFunc$FDiv #Functional divergence
    #FDiv[which(is.na(FDiv))]<- 0
    #taxonomic metric
    rich<- rowSums(data)
    matrix_Mobs_nullMod1<- as.matrix(data.frame(PDfaith, mntd, PSV, Peve, richness= rich))
    matrix_Mobs_nullMod1
  },
  error= function(err){
    return(NA)
  }
  )
  return(out)
}
