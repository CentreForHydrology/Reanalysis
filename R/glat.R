glat <-
function(x){
  # gets integer location corresponding to half-degree of latitude
  # 1 = -89.75
  # 360 = 89.75
  y <- round((x + 89.75) / 0.5)
  y <- y + 1
  return(y)
}
