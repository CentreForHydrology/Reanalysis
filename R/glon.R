glon <-
function(x){
  # gets integer location corresponding to half-degree of longitude
  # 1 = -179.75
  # 720 = 179.75
  y <- round((x + 179.75) / 0.5)
  y <- y + 1
  return(y)
}
