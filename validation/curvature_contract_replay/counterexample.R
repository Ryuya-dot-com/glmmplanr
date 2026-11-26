library(glmmplanr)
for (tau in 4) for(joint in c(FALSE,TRUE)) {
 ns<-4L; ni<-2L; s<-rep(0:(ns-1),each=ni); i<-rep(0:(ni-1),ns); z<-rep(c(-1,1),ns)
 a<-list(y=as.integer(s==0),X=cbind(1,z),subj_idx=as.integer(s),Z_subj=cbind(1,z),n_subj=ns,
  item_idx=as.integer(i),Z_item=matrix(1,length(i),1),n_item=ni,
  testlet_idx=integer(),Z_testlet=matrix(numeric(),0,0),n_testlet=0L,
  chance=1/3,corr_mode=0L,regularized=TRUE,init_Sigma_subj=diag(c(.01,tau^2)),init_Sigma_item=matrix(.01),
  fix_cov=TRUE,compute_joint_vcov=joint,return_ranef=TRUE,max_em_iter=400L,max_pirls_iter=44L)
 f<-do.call(glmmplanr:::cpp_fit_crossed_glmm,a)
 Ls<-t(chol(f$Sigma_subj)); Li<-t(chol(f$Sigma_item)); Z<-matrix(0,length(s),2*ns+ni)
 for(k in 1:2) Z[cbind(seq_along(s),2*s+k)]<-(a$Z_subj%*%Ls)[,k]
 Z[cbind(seq_along(s),2*ns+i+1)]<-(a$Z_item%*%Li)[,1]
 v<-c(as.vector(solve(Ls,t(f$u_subj))),as.vector(solve(Li,t(f$u_item))))
 eta<-drop(a$X%*%f$beta+Z%*%v); p<-plogis(eta); mu<-1/3+2/3*p
 dp<-2/3*p*(1-p); d2p<-dp*(1-2*p); y<-a$y
 w<-(dp^2-(y-mu)*d2p)/(mu*(1-mu))+(y-mu)*dp^2*(1-2*mu)/(mu^2*(1-mu)^2)
 H<-diag(length(v))+crossprod(Z,w*Z); ev<-eigen(H,symmetric=TRUE)
 phi<-function(vv) sum(dbinom(y,1,1/3+2/3*plogis(drop(a$X%*%f$beta+Z%*%vv)),log=TRUE))-.5*sum(vv^2)
 direction<-ev$vectors[,ncol(H)]; step<-.01
 cat('tau=',tau,'joint=',joint,'converged=',f$converged,'ll=',f$logLik,'beta=',f$beta,
  'minH=',min(ev$values),'inner=',f$inner_step,'phi_increase=',phi(v+step*direction)-phi(v),'\n')
 cat('finite_difference_curvature=', -(phi(v+step*direction)-2*phi(v)+phi(v-step*direction))/step^2, 'clamped_ll=',phi(v)-as.numeric(determinant(diag(length(v))+crossprod(Z,pmax(w,1e-6)*Z),logarithm=TRUE)$modulus)/2,'score_max=',max(abs(crossprod(Z,(y-mu)*dp/(mu*(1-mu)))-v)), '\n')
}
