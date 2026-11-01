args<-commandArgs(trailingOnly=TRUE);root<-args[1];seed<-as.integer(args[2])
suppressPackageStartupMessages(library(cluster))
x<-as.matrix(read.delim(file.path(root,"learning_CLR.tsv"),row.names=1,check.names=FALSE))
m<-read.delim(file.path(root,"learning_metadata.tsv"),stringsAsFactors=FALSE)
stopifnot(identical(rownames(x),m$sample_id))
set.seed(seed);n<-nrow(x);B<-200
groups<-split(seq_len(n),m$Study)
draws<-matrix(0L,B,n)
for(b in seq_len(B))draws[b,]<-unlist(lapply(groups,function(ix)sample(ix,length(ix),replace=TRUE)),use.names=FALSE)
write.table(draws,file.path(root,"bootstrap_draw_indices.tsv"),sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
states<-list();d<-dist(x)
for(k in 2:6){
 fit<-pam(d,k=k,diss=TRUE);med<-fit$id.med
 ord<-order(rownames(x)[med]);med<-med[ord];lab<-match(fit$clustering,ord)
 write.table(data.frame(sample_id=rownames(x),cluster=lab),file.path(root,paste0("original_K",k,".tsv")),sep="\t",row.names=FALSE,quote=FALSE)
 write.table(data.frame(cluster=seq_along(med),sample_id=rownames(x)[med]),file.path(root,paste0("medoids_K",k,".tsv")),sep="\t",row.names=FALSE,quote=FALSE)
 labels<-matrix(NA_integer_,B,n);medis<-matrix(NA_integer_,B,k)
 for(b in seq_len(B)){
  ix<-draws[b,]
  state<-tryCatch({
   if(length(unique(ix))<=k)stop("too_few_unique_draws")
   f<-pam(dist(x[ix,,drop=FALSE]),k=k,diss=TRUE);md<-ix[f$id.med]
   md<-md[order(rownames(x)[md])]
   if(length(unique(md))<k)stop("duplicate_medoid")
   dd<-sapply(md,function(j)rowSums(sweep(x,2,x[j,],"-")^2))
   la<-max.col(-dd,ties.method="first")
   if(length(unique(la[unique(ix)]))<k)stop("missing_cluster_on_shared_unique_IDs")
   labels[b,]<-la;medis[b,]<-md;"valid"
  },error=function(e)paste0("failed:",conditionMessage(e)))
  states[[length(states)+1]]<-data.frame(K=k,attempt=b,status=state)
 }
 write.table(labels,gzfile(file.path(root,paste0("bootstrap_labels_K",k,".tsv.gz"))),sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
 write.table(medis,file.path(root,paste0("bootstrap_medoids_K",k,".tsv")),sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
}
write.table(do.call(rbind,states),file.path(root,"bootstrap_status.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
capture.output(sessionInfo(),file=file.path(root,"R_sessionInfo.txt"))
