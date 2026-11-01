args<-commandArgs(trailingOnly=TRUE);root<-args[1]
suppressPackageStartupMessages(library(vegan))
t<-as.matrix(read.delim(gzfile(file.path(root,"TSS.tsv.gz")),row.names=1,check.names=FALSE))
y<-as.matrix(read.delim(gzfile(file.path(root,"CLR.tsv.gz")),row.names=1,check.names=FALSE))
m<-read.delim(file.path(root,"metadata.tsv"),row.names=1,check.names=FALSE)
N<-as.matrix(read.delim(file.path(root,"nuisance.tsv"),row.names=1,check.names=FALSE))
perms<-as.matrix(read.delim(gzfile(file.path(root,"permutations_1based.tsv.gz")),check.names=FALSE))
stopifnot(identical(rownames(t),rownames(m)),identical(rownames(y),rownames(m)),nrow(perms)==9999)
Alcohol<-as.numeric(m$analysis_group=="Alcohol");Z<-N[,-1,drop=FALSE]
results<-list()
for(metric in c("BrayCurtis","Aitchison")){
 d<-if(metric=="BrayCurtis") vegdist(t,"bray") else dist(y)
 ord<-if(ncol(Z)) dbrda(d~Alcohol+Condition(Z)) else dbrda(d~Alcohol)
 at<-anova(ord,permutations=perms[1:99,,drop=FALSE],by="margin",model="reduced")
 write.table(data.frame(F=as.numeric(attr(at,"F.perm"))),gzfile(file.path(root,paste0(metric,"_vegan_location_Fperm99.tsv.gz"))),sep="\t",quote=FALSE,row.names=FALSE)
 locF<-at["Alcohol","F"]
 locSS<-at["Alcohol",2];locR2<-locSS/ord$tot.chi
 bd<-betadisper(d,factor(m$analysis_group,levels=c("Health","Alcohol")),
               type="median",bias.adjust=TRUE,add=if(metric=="BrayCurtis")"lingoes" else FALSE)
 bt<-permutest(bd,permutations=perms)
 stopifnot(length(bt$perm)==9999)
 write.table(data.frame(F=as.numeric(bt$perm)),gzfile(file.path(root,paste0(metric,"_dispersion_Fperm.tsv.gz"))),sep="\t",quote=FALSE,row.names=FALSE)
 write.table(data.frame(sample_id=rownames(m),group=m$analysis_group,Study=m$Study,
                       distance_to_group_median=bd$distances),
             file.path(root,paste0(metric,"_dispersion_values.tsv")),sep="\t",quote=FALSE,row.names=FALSE)
 results[[metric]]<-data.frame(metric=metric,vegan_location_F=locF,vegan_location_R2=locR2,
              vegan_marginal_SS=locSS,vegan_total_inertia=ord$tot.chi,
              dispersion_F=bt$tab[1,"F"],
              dispersion_P=bt$tab[1,"Pr(>F)"],dispersion_B_actual=nrow(perms),
              dispersion_median=TRUE,dispersion_bias_adjust=TRUE,
              dispersion_distance_add=if(metric=="BrayCurtis")"lingoes" else "none",
              dispersion_method="actual_vegan_betadisper_spatial_median_restricted_residual_permutations")
}
write.table(do.call(rbind,results),file.path(root,"R_results.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
capture.output(sessionInfo(),file=file.path(root,"R_sessionInfo.txt"))
