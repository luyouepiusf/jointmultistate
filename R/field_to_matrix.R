field_to_matrix<-function(field){
  all_block_nrow<-apply(field,c(1,2),function(x)nrow(x[[1]]))
  all_block_ncol<-apply(field,c(1,2),function(x)ncol(x[[1]]))
  all_nonzero<-all_block_nrow>0&all_block_ncol>0
  block_nrow<-apply(all_block_nrow,1,max)
  block_ncol<-apply(all_block_ncol,2,max)
  result<-matrix(0,sum(block_nrow),sum(block_ncol))
  idx_row<-cumsum(block_nrow)
  idx_col<-cumsum(block_ncol)
  for(ii in 1:nrow(field)){
    for(jj in 1:ncol(field)){
      if(all_nonzero[ii,jj]){
        a_idx_row<-(idx_row[ii]-block_nrow[ii]+1):idx_row[ii]
        a_idx_col<-(idx_col[jj]-block_ncol[jj]+1):idx_col[jj]
        result[a_idx_row,a_idx_col]<-field[[ii,jj]]
      }
    }
  }
  return(list(
    mat=result,
    idx_row=rep(1:nrow(field),block_nrow),
    idx_col=rep(1:ncol(field),block_ncol)))
}
