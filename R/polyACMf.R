
# EstimatepolyACM <- function (ncase, nvar, IAdjust, NCore, iRaw) {
  
#  if (!is.integer(ncase)) {storage.mode(ncase) <- 'integer'}
#  if (!is.integer(nvar)) {storage.mode(nvar) <- 'integer'}
#  if (!is.integer(IAdjust)) {storage.mode(IAdjust) <- 'integer'}
#  if (!is.integer(NCore)) {storage.mode(NCore) <- 'integer'}
#  if (!is.integer(iRaw)) {storage.mode(iRaw) <- 'integer'}

#  .Call(c_polyACM_f, ncase, nvar, IAdjust, NCore, iRaw)
  
# }
 

# EstimatepolyR <- function (ncase, nvar, IAdjust, NCore, iRaw) {
   
#   if (!is.integer(ncase)) {storage.mode(ncase) <- 'integer'}
#   if (!is.integer(nvar)) {storage.mode(nvar) <- 'integer'}
#   if (!is.integer(IAdjust)) {storage.mode(IAdjust) <- 'integer'}
#   if (!is.integer(NCore)) {storage.mode(NCore) <- 'integer'}
#   if (!is.integer(iRaw)) {storage.mode(iRaw) <- 'integer'}
   
#   .Call(c_polyR_f, ncase, nvar, IAdjust, NCore, iRaw)
   
# }
 
 
 PolychoricRM <- function (iRaw=NULL, IMissing=1, IAdjust=0, NCore=2, estimate.acm=FALSE) {
    
   if (is.null(iRaw)) {
     return 
     }  else{
       iRaw[is.na(iRaw)] <- - 999
       nvar = ncol (iRaw)
       ncase = nrow (iRaw)
     }
  
    if (!is.integer(IMissing)) {storage.mode(IMissing) <- 'integer'}
    if (!is.integer(IAdjust)) {storage.mode(IAdjust) <- 'integer'}
    if (!is.integer(NCore)) {storage.mode(NCore) <- 'integer'}
    if (!is.integer(iRaw)) {storage.mode(iRaw) <- 'integer'}

   
       
    if (estimate.acm) {
     output =  .Call(c_polyACM_f, ncase, nvar, IMissing, IAdjust, NCore, iRaw)
     names(output) = c('threshold','correlation','flag','ACM')
    } else {
      output = .Call(c_polyR_f, ncase, nvar, IMissing, IAdjust, NCore, iRaw)
      names(output) = c('threshold','correlation','flag')
    }
    
   return(output)
 }
 