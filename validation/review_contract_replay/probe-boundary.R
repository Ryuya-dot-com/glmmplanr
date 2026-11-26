library(glmmplanr)
s <- rep(0:3,each=2L); i<-rep(0:1,4L); z<-rep(c(-1,1),4L)
a <- list(X=cbind(1,z),subj_idx=s,Z_subj=cbind(1,z),n_subj=4L,item_idx=i,Z_item=matrix(1,8,1),n_item=2L,testlet_idx=integer(),Z_testlet=matrix(numeric(),0,0),n_testlet=0L,chance=1/3,corr_mode=0L,regularized=TRUE,init_Sigma_item=matrix(.01),fix_cov=TRUE,max_em_iter=400L,return_ranef=TRUE)
rows<-list(); found<-list()
a$y<-as.integer(s==0L);a$max_pirls_iter<-44L
for(v in seq(.25,4,length.out=2001)) {
 a$init_Sigma_subj<-diag(c(.01,v));a$compute_joint_vcov<-FALSE
 f<-tryCatch(do.call(glmmplanr:::cpp_fit_crossed_glmm,a),error=function(e)e)
 a$compute_joint_vcov<-TRUE
 h<-tryCatch(do.call(glmmplanr:::cpp_fit_crossed_glmm,a),error=function(e)e)
 if(inherits(h,'error')) {found[[as.character(v)]]<-a;break}
 rows[[length(rows)+1L]]<-data.frame(v=v,plain=f$status,joint=h$status,step=h$inner_step)
}
saveRDS(list(rows=do.call(rbind,rows),found=found),commandArgs(TRUE)[1]);print(length(found));print(table(do.call(rbind,rows)$plain,do.call(rbind,rows)$joint))
