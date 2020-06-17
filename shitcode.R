require (MASS)


sim ← function (n , sigma =2 .5 , pr = FALSE , prcor = FALSE ) {
  x1 ← rnorm (n )
  x2 ← x1 + 0 .5 * rnorm (n )
  x3 ← rnorm (n )
  x4 ← x3 + 1 .5 * rnorm (n )
  x5 ← x1 + rnorm (n)/1 .3
  x6 ← x2 + rnorm (n)/1 .3
  x7 ← x3 + x4 + rnorm (n)
  x8 ← x7 + 0 .5 * rnorm (n )
  if( prcor ) return ( round (cor( cbind ( x1 ,x2 ,x3 , x4 , x5 ,x6 ,x7 , x8 ) ) ,2) )
  lp ← x1 + x2 + .5*x3 + .4*x7
  y ← lp + sigma * rnorm (n )
  f ← lm( y ∼ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 )
  g ← stepAIC (f , trace =0)
  p ← g$ rank - 1
  xs ← if( p == 0) ’ none ’ else
    gsub ( ’[ \\+ x] ’ ,’’ , as.character ( formula ( g)) [3])
  if( pr ) print ( formula (g) , showEnv = FALSE )
  ssesw ← sum( resid ( g)∧ 2)
  s2s ← ssesw /g$ df.residual
  # S e t S S E s w / ( n - g d f - 1 ) = t r u e s i g m a ∧ 2
  gdf ← n - 1 - ssesw /( sigma∧ 2)
  # C o m p u t e r o o t m e a n s q u a r e d e r r o r a g a i n s t t r u e l i n e a r p r e d i c t o r
  rmse.full ← sqrt ( mean (( fitted ( f) - lp ) ∧ 2) )
  rmse.step ← sqrt ( mean (( fitted ( g) - lp ) ∧ 2) )
  list ( stats =c( n=n , vratio = s2s/( sigma∧ 2) ,
                   gdf = gdf , apparentdf =p , rmse.full = rmse.full , rmse.step = rmse.step ) ,
         xselected = xs )
}

rsim ← function (B , n) {
  xs ← character (B)
  r ← matrix (NA , nrow =B , ncol =6)
  for(i in 1: B) {
    w ← sim (n )
    r[i ,] ← w$ stats
    xs [ i] ← w$ xselected
  }
  colnames ( r) ← names (w$ stats )
  s ← apply (r , 2, median )
  p ← r[ , ’ apparentdf ’]
  s[ ’ apparentdf ’] ← mean (p )
  print ( round (s , 2) )
  print ( table ( p) )
  cat( ’ Prob [ correct model ]= ’ , round (sum( xs == ’ 1237 ’)/B , 2) , ’\n ’)
}

sim (50000 , prcor = TRUE )

set.seed (11)
rsim (100 , 20)  