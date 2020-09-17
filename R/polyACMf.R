PolychoricRM <- function (iRaw=NULL, IAdjust=0, NCore=2, estimate.acm=FALSE) {
   
   if (is.null(iRaw)) {
      return 
   }  else{
      nvar = ncol (iRaw)
      ncase = nrow (iRaw)
   }
   
   
   if (!is.integer(IAdjust)) {storage.mode(IAdjust) <- 'integer'}
   if (!is.integer(NCore)) {storage.mode(NCore) <- 'integer'}
   if (!is.integer(iRaw)) {storage.mode(iRaw) <- 'integer'}
   
   
   
   if (estimate.acm) {
      output =  .Call(c_polyACM_f, ncase, nvar, IAdjust, NCore, iRaw)
      names(output) = c('threshold','correlation','flag','ACM')
   } else {
      output = .Call(c_polyR_f, ncase, nvar, IAdjust, NCore, iRaw)
      names(output) = c('threshold','correlation','flag')
   }
   
   return(output)
}
