library(glmmplanr)
s <- rep(0:3,each=2L); i<-rep(0:1,4L); z<-rep(c(-1,1),4L)
a <- list(X=cbind(1,z),subj_idx=s,Z_subj=cbind(1,z),n_subj=4L,item_idx=i,Z_item=matrix(1,8,1),n_item=2L,testlet_idx=integer(),Z_testlet=matrix(numeric(),0,0),n_testlet=0L,chance=1/3,corr_mode=0L,regularized=TRUE,init_Sigma_item=matrix(.01),fix_cov=TRUE,max_em_iter=400L,return_ranef=TRUE)
rows<-list(); found<-list()
for(v in c(4,16,64)) for(code in 1:254) {
 a$y<-as.integer(intToBits(code))[1:8]; a$init_Sigma_subj<-diag(c(.01,v)); a$max_pirls_iter<-44L
 f<-tryCatch(do.call(glmmplanr:::cpp_fit_crossed_glmm,a),error=function(e)e)
 rows[[length(rows)+1L]]<-data.frame(code=code,v=v,status=if(inherits(f,'error')) conditionMessage(f) else f$status,step=if(inherits(f,'error')) NA else f$inner_step)
 if(inherits(f,'error')) found[[paste(code,v)]]<-a
}
saveRDS(list(rows=do.call(rbind,rows),found=found),commandArgs(TRUE)[1]); print(table(do.call(rbind,rows)$status))
