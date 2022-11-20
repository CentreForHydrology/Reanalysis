halfdegree <-
function(x){
  # gets closest half degree to specified value
  y <- (x - 0.25) / 0.5
  y <- round(y)
  hd <- y * 0.5 + 0.25
  return(hd)
}
