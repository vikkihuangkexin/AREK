#!/usr/bin/env Rscript
# Corrected entry: real vegan adonis2 and betadisper, with explicit design guards.
suppressPackageStartupMessages({library(vegan);library(permute)})
a<-commandArgs(TRUE);folder<-a[1];metapath<-a[2];groupcol<-a[3];subjectcol<-a[4];blockcol<-a[5];design<-a[6];covars<-if(a[7]=='')character() else strsplit(a[7],',',fixed=TRUE)[[1]];nperm<-as.integer(a[8]);set.seed(as.integer(a[9]))
readtab<-function(p)read.delim(p,row.names=1,check.names=FALSE,quote='',comment.char='',na.strings=c('NA',''))
alpha<-readtab(file.path(folder,'alpha.tsv'));meta<-if(metapath!='')readtab(metapath) else data.frame(row.names=rownames(alpha))
meta<-meta[match(rownames(alpha),rownames(meta)),,drop=FALSE];rownames(meta)<-rownames(alpha)
writeout<-function(x,name)write.table(x,file.path(folder,name),sep='\t',row.names=FALSE,quote=FALSE,na='NA')
reason<-if(design=='descriptive')'descriptive_design_no_group_inference' else '';ids<-rownames(alpha)
needed<-unique(c(groupcol,covars,if(subjectcol!='')subjectcol,if(blockcol!='')blockcol))
if(groupcol=='' || !groupcol%in%names(meta))reason<-'group_column_missing'
if(reason=='' && any(!needed%in%names(meta)))reason<-paste('metadata_columns_missing',paste(setdiff(needed,names(meta)),collapse=','))
if(reason==''){
 ok<-complete.cases(meta[,needed,drop=FALSE]);ids<-rownames(meta)[ok];meta<-meta[ok,,drop=FALSE];g<-factor(meta[[groupcol]])
 if(length(levels(g))<2)reason<-'single_arm_or_all_unknown'
 else if(min(table(g))<2)reason<-'fewer_than_two_samples_in_group'
 repeated<-subjectcol!='' && anyDuplicated(meta[[subjectcol]])>0
 if(repeated && design!='within_subject')reason<-'repeated_subjects: specify within_subject+block or use an explicit study-specific independent specimen selection'
 if(repeated && (blockcol=='' || !identical(as.character(meta[[blockcol]]),as.character(meta[[subjectcol]]))))reason<-'repeated_subjects_require_subject_permutation_block'
 if(repeated && any(vapply(split(as.character(g),meta[[subjectcol]]),function(x)length(unique(x))<2,logical(1))))reason<-'group_invariant_within_subject: within-subject permutation cannot estimate between-subject effect'
}
alpha_rows<-list()
for(metric in names(alpha)){
 p<-stat<-NA_real_;why<-reason;method<-'not_estimable';n<-length(ids)
 if(why=='')tryCatch({
  vals<-alpha[ids,metric];good<-is.finite(vals);gg<-droplevels(g[good]);vv<-vals[good];n<-length(vv)
  if(any(table(gg)<2)||nlevels(gg)<2)stop('insufficient_finite_values')
  if(length(covars))stop('alpha_covariate_adjusted_test_not_implemented; descriptive values retained')
  if(design=='within_subject'){
   if(nlevels(gg)!=2 || subjectcol=='')stop('paired_alpha_requires_two_groups_and_subject')
   sub<-as.character(meta[good,subjectcol]);tab<-table(sub,gg)
   if(any(tab!=1))stop('paired_alpha_requires_exactly_one_sample_per_subject_per_group')
   lev<-levels(gg);aa<-setNames(vv[gg==lev[1]],sub[gg==lev[1]]);bb<-setNames(vv[gg==lev[2]],sub[gg==lev[2]])
   test<-wilcox.test(aa[sort(names(aa))],bb[sort(names(aa))],paired=TRUE,exact=FALSE);method<-'paired_Wilcoxon_signed_rank'
  }else if(blockcol!='')stop('blocked_alpha_requires_explicit_paired_design')
  else if(nlevels(gg)==2){test<-wilcox.test(vv~gg,exact=FALSE);method<-'Mann_Whitney_asymptotic_R'}
  else{test<-kruskal.test(vv~gg);method<-'Kruskal_Wallis'}
  p<-test$p.value;stat<-unname(test$statistic)
 },error=function(e)why<<-conditionMessage(e))
 if(why=='' && (!is.finite(p)||!is.finite(stat)))why<-'nonfinite_statistic_or_p: constant_data_or_nonestimable'
 alpha_rows[[metric]]<-data.frame(metric=metric,method=method,statistic=stat,p_value=p,n_samples=n,reason=why)
}
aout<-do.call(rbind,alpha_rows);aout$FDR_p_value<-p.adjust(aout$p_value,'BH');writeout(aout,'alpha_tests.tsv')
perrows<-disprows<-list()
for(metric in c('BrayCurtis','Jaccard','Aitchison')){
 dd<-as.matrix(readtab(file.path(folder,paste0(metric,'.distance.tsv'))));why<-reason;p<-f<-r2<-NA_real_;dp<-df<-NA_real_;dwhy<-why
 if(why==''){
  d<-as.dist(dd[ids,ids,drop=FALSE]);dat<-data.frame(group=g)
  cvnames<-character()
  for(i in seq_along(covars)){nm<-paste0('cov',i);dat[[nm]]<-meta[[covars[i]]];if(is.character(dat[[nm]])){num<-suppressWarnings(as.numeric(dat[[nm]]));dat[[nm]]<-if(all(is.finite(num)))num else factor(dat[[nm]])};cvnames<-c(cvnames,nm)}
  ctrl<-how(nperm=nperm);if(blockcol!='')setBlocks(ctrl)<-factor(meta[[blockcol]])
  perms<-tryCatch(shuffleSet(nrow(dat),nset=nperm,control=ctrl),error=function(e)e)
  if(inherits(perms,'error'))why<-dwhy<-paste('permutations_unavailable',conditionMessage(perms))
  else if(!nrow(perms) || all(apply(perms,1,function(x)identical(as.character(g[x]),as.character(g)))))why<-dwhy<-'no_exchangeable_group_labels_within_blocks'
  else {
   tryCatch({
    form<-reformulate(c(cvnames,'group'),response='d');mm<-model.matrix(reformulate(c(cvnames,'group')),data=dat)
    if(qr(mm)$rank<ncol(mm)||nrow(mm)<=ncol(mm))stop('rank_deficient_or_no_residual_df')
    fit<-adonis2(form,data=dat,permutations=perms,by='margin');p<-fit['group','Pr(>F)'];f<-fit['group','F'];r2<-fit['group','R2']
   },error=function(e)why<<-conditionMessage(e))
   tryCatch({bd<-betadisper(d,g,type='median');bt<-permutest(bd,permutations=perms);df<-bt$tab[1,'F'];dp<-bt$tab[1,'Pr(>F)']},error=function(e)dwhy<<-conditionMessage(e))
  }
 }
 if(why=='' && (!is.finite(p)||!is.finite(f)))why<-'nonfinite_PERMANOVA: zero_distance_or_nonestimable'
 if(dwhy=='' && (!is.finite(dp)||!is.finite(df)))dwhy<-'nonfinite_PERMDISP: constant_dispersion_or_nonestimable'
 perrows[[metric]]<-data.frame(metric=metric,method='vegan_adonis2_marginal',F_statistic=f,R2=r2,p_value=p,n_samples=length(ids),reason=why)
 disprows[[metric]]<-data.frame(metric=metric,method='vegan_betadisper_spatial_median_permutest',F_statistic=df,p_value=dp,n_samples=length(ids),reason=dwhy)
}
po<-do.call(rbind,perrows);po$FDR_p_value<-p.adjust(po$p_value,'BH');writeout(po,'permanova.tsv')
bo<-do.call(rbind,disprows);bo$FDR_p_value<-p.adjust(bo$p_value,'BH');writeout(bo,'dispersion.tsv')
writeout(data.frame(sample_id=rownames(alpha),in_statistics=rownames(alpha)%in%ids,reason=ifelse(rownames(alpha)%in%ids,'included','metadata_complete_case_exclusion')),'statistics_sample_audit.tsv')
writeLines(c(capture.output(sessionInfo()),paste('design:',design,'group:',groupcol,'subject:',subjectcol,'block:',blockcol,'nperm:',nperm)),file.path(folder,'statistics_session.txt'))
