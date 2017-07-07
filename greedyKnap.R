eval<-function(conf,w,v){
  Wm<-400
  vals<-v
  ws<-w
  cVal<-sum(conf*vals) #config value
  cW<-sum(conf*ws) #config weight
  lim<-Wm+mean(ws)
  #lim<-Wm+lightest-1 try the -1 so it may handle cases better where everything is the same weigh
  res<-numeric(3)
  res[1]<-ifelse(cW>lim,0,cVal) #if the weight is over the limit return 0, else return the value
  res[2]<-cW
  res[3]<-ifelse(cW>Wm && cW<lim,1,0) #if the weight is greater than Wm AND less than the limit return 1 else return 0
  return(res)
}
greedy<-function(w,v){
  vals<-v
  ws<-w
  rs<-vals/ws
  ind<-sort.int(rs,decreasing=T,index.return=T)$ix
  reOrderedWeights<-(ws[ind])
  curW<-0
  Wm<-400
  conf<-numeric(8)
  track<-1
  while(track<=8){
    if(curW+reOrderedWeights[track]>Wm){
      break
    }else{
    curW<-curW+reOrderedWeights[track]
    conf[ind[track]]<-1
    track<-track+1
    }
  }
  return(conf)
}
vals<-c(50,40,30,60,100,150,120,70)
ws<-c(50,30,60,50,50,50,50,200)
greedy(ws,vals)