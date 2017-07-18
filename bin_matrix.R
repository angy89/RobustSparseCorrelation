bin_matrix = function(R, th){
  R_bin = R
  R_bin[abs(R)<= th] = 0
  R_bin[abs(R)>th] = 1
  
  return(R_bin)
}
