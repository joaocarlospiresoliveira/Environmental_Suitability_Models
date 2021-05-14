
###################### Acknowledgments ##########################
### Dr. Matheus de Souza Lima-Ribeiro's team of Universidade Federal de Jataí.
### Dr. Diogo Souza Bezerra Rocha (Botanical Garden Research Institute / RJ).
### MSc. João Carlos Pires de Oliveira (PhD candidate)

## Install and Loading packages ####
# install.packages("rgdal")
# install.packages("sp")
# install.packages("raster")
# install.packages("maps")
# install.packages("mnormt") #mnormt / psych
# install.packages("psych")
# #install.packages("permut")#permut / vegan ... permut desatualizado
# install.packages("vegan")
# install.packages("dismo")
# install.packages("kernlab")
# install.packages("rJava")
# install.packages("randomForest")
# install.packages("earth")
# install.packages("mgcv")
# install.packages("nnet")
# install.packages("beepr")
# install.packages("doParallel")
# install.packages("biomod2")
# install.packages("sdmvspecies")
# install.packages("filesstrings")
# install.packages("dplyr")
library(sp)
library(rgdal)
library(tcltk2)
library(raster)
library(maps)
library(psych)
library(vegan)
library(mnormt)
library(dismo)
library(kernlab)
library(rJava)
library(randomForest)
library(earth)
library(mgcv)
library(nnet)
library(beepr)
library(doParallel)
library(biomod2)
library(sdmvspecies)
library(filesstrings)
library(dplyr)
# Criating output dir #
if (dir.exists("outputs") == F) {
  dir.create("outputs")
}

# Function to scale maps -- DON'T CHANGE###
rescMod = function(ras1, ras2, ras3){
  p1 <- rasterToPoints(stack(ras1))
  p2 <- rasterToPoints(stack(ras2))
  p3 <- rasterToPoints(stack(ras3))
  coord <- rasterToPoints(ras1)[, 1:2]
  id <- rep(c("cur", "fut.45", "fut.85"), each = nrow(p1))
  rd = rbind(p1, p2, p3)
  st = vegan::decostand(x = rd[,3], method = "range", MARGIN = 2)
  ens = data.frame(coord, Cur = st[id == "cur"], 
                   Fut.45 = st[id =="fut.45"],
                   Fut.85 = st[id =="fut.85"])
  ens.cur = raster::rasterize(ens[,c("x","y")],ras1[[1]], ens[,"Cur"])
  ens.fut.45 = raster::rasterize(ens[,c("x","y")],ras1[[1]], ens[,"Fut.45"])
  ens.fut.85 = raster::rasterize(ens[,c("x","y")],ras1[[1]], ens[,"Fut.85"])
  resul = stack(ens.cur, ens.fut.45, ens.fut.85)
  return(resul)
  rm(p1,p2,p3,coord,id,rd,st,ens,ens.fut.45, ens.fut.85)
  gc(reset = T, full = T)
}

# Function to scale maps -- DON'T CHANGE####
rescMod.One =  function(raster.layer) {
  if (!(class(raster.layer) %in% "RasterLayer")) {
    stop("raster.layer is not a RasterLayer objectect!")
  }
  min.value <- cellStats(raster.layer, min)
  if (min.value < 0) {
    raster.layer <- raster.layer + (0-min.value)
    max.value <- cellStats(raster.layer, max)
    raster.layer <- raster.layer/max.value
    return(raster.layer)
  } else {
    max.value <- cellStats(raster.layer, max)
    raster.layer <- raster.layer/max.value
    return(raster.layer)
  }
}

# Function to scale vars -- DON'T CHANGE###
rescVars = function(rast, nClimVars = NULL){
  p1 <- raster::rasterToPoints(stack(rast[[1:nClimVars]]))
  coord <- p1[, 1:2]
  id <- rep(paste0("PC", 1:ncol(p1[,-c(1:2)])), each = nrow(p1))
  rd = c(p1[,-c(1:2)][,1], p1[,-c(1:2)][,2], p1[,-c(1:2)][,3], p1[,-c(1:2)][,4], 
         p1[,-c(1:2)][,5], p1[,-c(1:2)][,6])
  st = vegan::decostand(x = rd, method = "range", MARGIN = 2)
  st = data.frame(lay = st,PC = id)
  ens = data.frame(coord,
                   PC1 = st[st[,"PC"] == "PC1",1],
                   PC2 = st[st[,"PC"] == "PC2",1],
                   PC3 = st[st[,"PC"] == "PC3",1],
                   PC4 = st[st[,"PC"] == "PC4",1],
                   PC5 = st[st[,"PC"] == "PC5",1],
                   PC6 = st[st[,"PC"] == "PC6",1])
  ens1.cur = raster::rasterize(ens[,c("x","y")], rast[[1]],ens[,(paste0('PC',1:6))])
  
  p2 <- raster::rasterToPoints(stack(rast[[(nClimVars+1):nlayers(rast)]]))
  coord <- p2[, 1:2]
  id <- rep(paste0("PC", 7:8), each = nrow(p2))
  rd2 = c(p2[,-c(1:2)][,1], p2[,-c(1:2)][,2])
  st2 = vegan::decostand(x = rd2, method = "range", MARGIN = 2)
  st2 = data.frame(lay = st2,PC = id)
  ens2 = data.frame(coord,
                    PC7 = st2[st2[,"PC"] == "PC7",1],
                    PC8 = st2[st2[,"PC"] == "PC8",1])
  ens2.cur = raster::rasterize(ens2[,c("x","y")], rast[[1]],ens2[,(paste0('PC',7:8))])
  
  resul = stack(ens1.cur,ens2.cur)
  return(resul)
  rm(p1,p2,coord,id,rd, rd2,st, st2,ens, ens2,ens1.cur, ens2.cur)
  gc(reset = T, full = T)
}

# Function to Evaluate All Models --- DON'T CHANGE  ###

eval.All.Model <- function(rast, dismoTestPrepared, tr ) {{
  p = raster::extract(rast, dismoTestPrepared[dismoTestPrepared[, 3] == 1, 1:2 ])
  a = raster::extract(rast, dismoTestPrepared[dismoTestPrepared[, 3] == 0, 1:2 ])}
  p <- stats::na.omit(p)
  a <- stats::na.omit(a)
  np <- length(p)
  na <- length(a)
  if (na == 0 | np == 0) {
    stop('cannot evaluate a model without absence and presence data that are not NA')}
  if (missing(tr)) {
    if (length(p) > 1000) {
      tr <- as.vector(quantile(p, 0:1000/1000))
    } else {
      tr <- p}
    if (length(a) > 1000) {
      tr <- c(tr, as.vector(quantile(a, 0:1000/1000)))
    } else {
      tr <- c(tr, a)}
    tr <- sort(unique( round(tr, 8)))
    tr <- c( tr - 0.0001, tr[length(tr)] + c(0, 0.0001))
  } else {
    tr <- sort(as.vector(tr))}
  N <- na + np
  xc <- new('ModelEvaluation')
  xc@presence = p
  xc@absence = a
  R <- sum(rank(c(p, a))[1:np]) - (np*(np+1)/2)
  xc@auc <- R / (as.numeric(na) * as.numeric(np))
  cr <- try( cor.test(c(p,a), c(rep(1, length(p)), rep(0, length(a))) ), silent=TRUE )
  if (class(cr) != 'try-error') {
    xc@cor <- cr$estimate
    xc@pcor <- cr$p.value}
  res <- matrix(ncol=4, nrow=length(tr))
  colnames(res) <- c('tp', 'fp', 'fn', 'tn')
  xc@t <- tr
  for (i in 1:length(tr)) {
    res[i,1] <- length(p[p>=tr[i]])  # a  true positives
    res[i,2] <- length(a[a>=tr[i]])  # b  false positives
    res[i,3] <- length(p[p<tr[i]])    # c  false negatives
    res[i,4] <- length(a[a<tr[i]])}    # d  true negatives
  xc@confusion = res
  a = res[,1]
  b = res[,2]
  c = res[,3]
  d = res[,4]
  # after Fielding and Bell	
  xc@np <- as.integer(np)
  xc@na <- as.integer(na)
  xc@prevalence = (a + c) / N
  xc@ODP = (b + d) / N
  xc@CCR = (a + d) / N
  xc@TPR = a / (a + c)
  xc@TNR = d / (b + d)
  xc@FPR = b / (b + d)
  xc@FNR = c/(a + c)
  xc@PPP = a/(a + b)
  xc@NPP = d/(c + d)
  xc@MCR = (b + c)/N
  xc@OR = (a*d)/(c*b)
  prA = (a+d)/N
  prY = (a+b)/N * (a+c)/N
  prN = (c+d)/N * (b+d)/N
  prE = prY + prN
  xc@kappa = (prA - prE) / (1-prE)
  return(xc)}

# Function to predict future predictions ###
preFut = function(rast, rast1 = NULL, model, GCM =  NULL) {
  if(missing(GCM)){
    GCM = 1
  }else{GCM = GCM}
  pre = foreach::foreach(
    gcm = 1:length(GCM),
    .combine = stack,
    .packages = c("raster","biomod2", 'sp', "kernlab","dismo",
                  "stats","randomForest", "nnet", "earth" )
  ) %dopar% {
    predict.enfa <-function (object.enfa,baseline.climate,new.climate,nf = 2,...) {
      m.baseline <- apply(slot(baseline.climate, "data"), 2, mean)
      sd.baseline <-apply(slot(baseline.climate, "data"), 2, sd)
      Zli <-object.enfa$li[, 1:nf]
      f1 <-function(x) {rep(x, object.enfa$pr)}
      Sli <- apply(Zli, 2, f1)
      m <- apply(Sli, 2, mean)
      cov <-t(as.matrix(Sli)) %*% as.matrix(Sli) / nrow(Sli)
      if (!missing("new.climate")) {
        new.climate.scale <- sweep(slot(new.climate, "data"), 2, m.baseline)
        new.climate.scale <-
          as.matrix(new.climate.scale) %*% diag(1 / sd.baseline)
        Zli <-
          new.climate.scale %*% as.matrix(object.enfa$co)}
      maha <- mahalanobis(Zli, center = m, cov = cov)
      map <-rasterize(data.frame(new.climate@coords),rast[[1]],maha) * -1
      return(invisible(map))
    }
    if (!"madifa" %in% class(model)) {
      raster::predict(rast[gcm][[1]], model)
    } else if ("madifa" %in% class(model) && "list" %in% class(rast1)) {
      climaPres.spdf = na.omit(data.frame(xyFromCell(rast, 1:ncell(rast)),
                                          raster::values(rast)))
      climaFut.spdf = na.omit(data.frame(xyFromCell(rast1[gcm][[1]], 
                                                    1:ncell(rast1[gcm][[1]])),
                                         raster::values(rast1[gcm][[1]])))
      suppressWarnings(gridded(climaPres.spdf) <-~ x + y)
      suppressWarnings(gridded(climaFut.spdf) <- ~ x + y)
      predict.enfa( object.enfa = model,baseline.climate =climaPres.spdf,
        new.climate = climaFut.spdf)}else{
          climaPres.spdf = na.omit(data.frame(xyFromCell(rast, 1:ncell(rast)),
                                              raster::values(rast)))
          climaFut.spdf = na.omit(data.frame(xyFromCell(rast1, 
                                                        1:ncell(rast1)),
                                             raster::values(rast1)))
          suppressWarnings(gridded(climaPres.spdf) <-~ x + y)
          suppressWarnings(gridded(climaFut.spdf) <- ~ x + y)
          predict.enfa( object.enfa = model,baseline.climate =climaPres.spdf,
                        new.climate = climaFut.spdf)}
  }
  pre = mean(pre)
  return(pre)
  rm(new.climate.scale,maha,map,climaPres.spdf,climaFut.spdf)
  gc(reset = T, full = T)
}

quiet <- function(x) { 
  sink(tempfile()) 
  on.exit(sink()) 
  invisible(force(x)) 
} 

### Whenever necessary:
# Parallel processing #
dc = detectCores()
getDoParWorkers()
cl <-
  parallel::makeCluster(3, outfile = paste0("./outputs/", "log_models.log"))
#cl <- parallel::makeCluster(10, type = "MPI", outfile = "./outputs/joao.log")
registerDoParallel(cl)
getDoParWorkers()

# Increased memory allocation
memory.limit(17592186044415) # or some other memory value (in kB)
# memory.limit(8062000000000) # or some other memory value (in kB)


## CHANGE raster TEMPORARY FILE DIRECTORY
## define the name of a temp directory where raster tmp files will be stored
raster_tmp_dir <- "raster_tmp"
if (dir.exists("raster_tmp") == F) {
  ## create the directory (only when starting modeling):
  dir.create(raster_tmp_dir,
             showWarnings = F,
             recursive = T)
}
#
# ## set raster options
rasterOptions(tmpdir = raster_tmp_dir)


##### Loading variaveis ####

bio.crop<-
  list.files(
    "./Present/PCA",
    full.names = TRUE,
    pattern = ".bil$"
  )
bio.crop
bio.crop <- stack(bio.crop)
# bio.crop <- disaggregate(bio.crop, fact = 5, fun = mean)
# bio.crop <- aggregate(bio.crop, fact = 5, fun = mean)
# bio.crop <- rescVars(bio.crop, nClimVars = 6)
names(bio.crop)<-
  c("PCA1", "PCA2", "PCA3","PCA4", "PCA5","PCA6", "PCA7", "PCA8")


#--------------------------------------------------#
#      LOADING FUTURE VARIABLES                ####
#------------------------------------------------#
###GCM 1: CanESM5

bio70_CA_45 <-
  list.files(
    "./Future/RCP45/CanESM5/ssp245/PCA",
    pattern = ".bil$",
    full.names = TRUE
  )
bio70_CA_45
bio70_CA_45 <- stack(bio70_CA_45)
# bio70_CA_45 <-disaggregate(bio70_CA_45, fact=5, fun=mean)
# bio70_CA_45 <-aggregate(bio70_CA_45, fact=5, fun=mean)
# bio70_CA_45<-rescVars(bio70_CA_45, nClimVars = 6)
bio70_CA_45
names(bio70_CA_45) <-
  c("PCA1", "PCA2", "PCA3","PCA4", "PCA5","PCA6", "PCA7", "PCA8")
names(bio70_CA_45)


###GCM 2: CNRM-CM6-1
bio70_CN_45 <-
  list.files(
    "./Future/RCP45/CNRM-CM6-1/ssp245/PCA",
    pattern = ".bil$",
    full.names = TRUE
  )
bio70_CN_45
bio70_CN_45 <- stack(bio70_CN_45)
# bio70_CN_45 <-disaggregate(bio70_CN_45, fact=5, fun=mean)
# bio70_CN_45 <-aggregate(bio70_CN_45, fact=5, fun=mean)
# bio70_CN_45<-rescVars(bio70_CN_45, nClimVars = 6)
bio70_CN_45
names(bio70_CN_45) <-
  c("PCA1", "PCA2", "PCA3","PCA4", "PCA5","PCA6", "PCA7", "PCA8")
names(bio70_CN_45)

###GCM 3: MIROC-ES2L

bio70_MI_45 <-
  list.files(
    "./Future/RCP45/MIROC-ES2L/ssp245/PCA",
    pattern = ".bil$",
    full.names = TRUE
  )
bio70_MI_45
bio70_MI_45 <- stack(bio70_MI_45)
# bio70_MI_45 <-disaggregate(bio70_MI_45, fact=5, fun=mean)
# bio70_MI_45 <-aggregate(bio70_MI_45, fact=5, fun=mean)
# bio70_MI_45<-rescVars(bio70_MI_45, nClimVars = 6)
bio70_MI_45
names(bio70_MI_45) <-
  c("PCA1", "PCA2", "PCA3","PCA4", "PCA5","PCA6", "PCA7", "PCA8")
names(bio70_MI_45)

# Loading raster layer with future projections 85 ###
#------------------------------------------------#
###GCM 1: CanESM5

bio70_CA_85 <-
  list.files(
    "./Future/RCP85/CanESM5/ssp585/PCA",
    pattern = ".bil$",
    full.names = TRUE
  )
bio70_CA_85
bio70_CA_85 <- stack(bio70_CA_85)
# bio70_CA_85 <-disaggregate(bio70_CA_85, fact=5, fun=mean)
# bio70_CA_85 <-aggregate(bio70_CA_85, fact=5, fun=mean)
# bio70_CA_85<-rescVars(bio70_CA_85, nClimVars = 6)
bio70_CA_85
names(bio70_CA_85) <-
  c("PCA1", "PCA2", "PCA3","PCA4", "PCA5","PCA6", "PCA7", "PCA8")
names(bio70_CA_85)


###GCM 2: CNRM-CM6-1
bio70_CN_85 <-
  list.files(
    "./Future/RCP85/CNRM-CM6-1/ssp585/PCA",
    pattern = ".bil$",
    full.names = TRUE
  )
bio70_CN_85
bio70_CN_85 <- stack(bio70_CN_85)
# bio70_CN_85 <-disaggregate(bio70_CN_85, fact=5, fun=mean)
# bio70_CN_85 <-aggregate(bio70_CN_85, fact=5, fun=mean)
# bio70_CN_85<-rescVars(bio70_CN_85, nClimVars = 6)
bio70_CN_85
names(bio70_CN_85) <-
  c("PCA1", "PCA2", "PCA3","PCA4", "PCA5","PCA6", "PCA7", "PCA8")
names(bio70_CN_85)

###GCM 3: MIROC-ES2L

bio70_MI_85 <-
  list.files(
    "./Future/RCP85/MIROC-ES2L/ssp585/PCA",
    pattern = ".bil$",
    full.names = TRUE
  )
bio70_MI_85
bio70_MI_85 <- stack(bio70_MI_85)
# bio70_MI_85 <-disaggregate(bio70_MI_85, fact=5, fun=mean)
# bio70_MI_85 <-aggregate(bio70_MI_85, fact=5, fun=mean)
# bio70_MI_85<-rescVars(bio70_MI_85, nClimVars = 6)
bio70_MI_85
names(bio70_MI_85) <-
  c("PCA1", "PCA2", "PCA3","PCA4", "PCA5","PCA6", "PCA7", "PCA8")
names(bio70_MI_85)

scenario.list.45 <- list(
  bio70_CA_45,
  bio70_CN_45,
  bio70_MI_45
)

scenario.list.85 <- list(
  bio70_CA_85,
  bio70_CN_85,
  bio70_MI_85
)

rm( bio70_CA_45,
    bio70_CN_45,
    bio70_MI_45,
    bio70_CA_85,
    bio70_CN_85,
    bio70_MI_85)

##### Occurrences ####
# Select you species data matrix # 
spp <- read.table(file.choose(), header = T, sep = ',')  

# plot all your occurence points #
# plot(bio.crop[[1]]) 
# points(spp[,c("lon","lat)])

dim(spp)
head(spp, 10)

table(spp$sp)

especies <- unique(spp$sp)
especies


# Creating VarImport Object ###
var = data.frame(
  varia = c(
    paste0("Clim", 1:6),
    paste0("Topo", 1:2),
    rep("BlockC1", 6),
    rep("BlockT2", 2)
  ),
  pos = c(1:8, 1:8),
  PC = c(1:8, 1:8)
)

vars <- unique(var$varia)

# Pseudo-absence Set

PAs <- 10
RUNs = 10

GCM45s <- c("bio70_CA_45", "bio70_CN_45", "bio70_MI_45")
GCM85s <- c("bio70_CA_85", "bio70_CN_85", "bio70_MI_85")

# If you wish to compute the importance of variables, put VarImport = TRUE,
# or set VarImport = FALSE not to do this.
VarImport = TRUE

# especie = especies[1]
#### Species Loop ####
# For sequential loop (One species) ###
# for (especie in especies[2:3]) {
                        
## For species in parallel ###
foreach(especie = especies, # For parallel looping (Multiple Species)
        .packages = c("raster", "biomod2",'sp',"sdmvspecies", "filesstrings",
                      "rgdal","maps","mnormt","kernlab","dismo","doParallel",
                      "stats","rJava","randomForest","nnet","psych", "earth"),
        .verbose = F,
        .errorhandling = "stop") %do% {

  if (dir.exists(paste0("./temp_output/",especie,"/")) == F) {
    dir.create(paste0("./temp_output/"))
    dir.create(paste0("./temp_output/",especie,"/"))
  }

  ini1 = Sys.time()
  print(paste0("Starting", " ", especie, " " ,"modeling"))
  
  # Creating empty objects to store results ###
  bioclim.e <-
    domain.e <-
    enfa.e <-
    glm.e <- gam.e <- mars.e <- maxent.e <- nnet.e <- svm.e <- rf.e <- NULL
  
  
  occs <- spp[spp$sp == especie, c("lon", "lat")]
  
  # Data checking and prepation
  ocor.val <- raster::extract(bio.crop, occs, cellnumbers = T)
  
  sum(is.na(ocor.val[, 1]))
  
  # plot(predictors[[1]], colNA = "red")
  
  # points(ocor.all[, -1], pch = 19, cex = 0.5)
  
  ocor.val <- cbind(occs, ocor.val)
  
  ocor.val <- na.omit(ocor.val)
  
  id <- duplicated(ocor.val[, "cells"]) # Checking dumplicate points
  
  sum(id == T)
  
  ocor <-
    ocor.val[id == F, c("lon", "lat")] # Removing duplicate points
  
  #------------------------------------------#
  #           SELECT PAs                  ###
  #----------------------------------------#
  
  try({
    coord1 = ocor
    sp::coordinates(coord1) <- ~ lon + lat
    raster::crs(coord1) <- raster::crs(bio.crop)
    
    dist.mean <- mean(sp::spDists(
      x = coord1,
      longlat = T,
      segments = FALSE
    ))
    dist.min = 5
    dist.min <-  min(sp::spDists(
      x = coord1,
      longlat = T,
      segments = F
    ))
    dist.min = 5
    
    write.table(
      c(dist.min, dist.mean),
      paste0('./outputs/',
             especie, "_","dist", ".csv"),
      row.names = F,
      sep = ";"
    )
  })
  
  PA.number <- nrow(ocor)
  PA.number #número de pontos de ocorrência espacialmente únicos
  
  diretorio = paste0("Occurrence.", especie)
  
  bioclim.var.All <- domain.var.All <- enfa.var.All <- glm.var.All <-
   mars.var.All <- maxent.var.All <- svm.var.All <- 
    nnet.var.All <- rf.var.All <- NULL
  
  # Loop PAs ####
  for (PA in seq(PAs)) {

   
    # Preparando
    invisible(capture.output(sel.PA <- biomod2::BIOMOD_FormatingData(
      resp.var = rep(1, times = nrow(ocor)),
      expl.var = raster::stack(bio.crop),
      resp.xy = ocor,
      resp.name = diretorio,
      PA.nb.rep = 1,
      #número de datasets de pseudoausências
      PA.nb.absences = PA.number,
      #= número de pseudoausências = número de pontos espacialmente únicos
      PA.strategy = "disk",
      PA.dist.min = dist.min * 1000,
      PA.dist.max = dist.mean * 1000,
      na.rm = TRUE
    )))
    
    li.p <- grep("pa", rownames(sel.PA@coord))
    
    
    pa <- sel.PA@coord[li.p,]
    
    invisible(capture.output(sel.back <- biomod2::BIOMOD_FormatingData(
      resp.var = rep(1, times = nrow(ocor)),
      expl.var = raster::stack(bio.crop),
      resp.xy = ocor,
      resp.name = diretorio,
      PA.nb.rep = 1,
      PA.nb.absences = 10000,
      PA.strategy = "disk",
      PA.dist.min = dist.min * 1000,
      PA.dist.max = dist.mean * 1000,
      na.rm = TRUE
    )))
    
    li.b <- grep("pa", rownames(sel.back@coord))
    
    
    back <- sel.back@coord[li.b,]
    
    
    rm(sel.back, sel.PA)
    # set.seed(0)
    areaToral <- nrow(rasterToPoints(bio.crop))
    
    # Loop RUN ####
    for (RUN in seq(RUNs)) {
      print(paste0(especie," ","PASet", PA," ","RUN", RUN))
      # Separating test/ training data
      
      
      id.training.pa <-
        sample(1:nrow(ocor), round(0.7 * nrow(ocor), 0)) # prepare data 70/30
      
      id.training.b <-
        sample(1:nrow(back), round(0.7 * nrow(back), 0)) # prepare data 70/30
      
      training.b <-
        na.omit(dismo::prepareData(bio.crop, p = ocor[id.training.pa, ], 
                                   b = back[id.training.b,], xy = T))
      # head(training.b)
      # tail(training.b)
      training.pa <-
        na.omit(dismo::prepareData(bio.crop, p = ocor[id.training.pa, ], 
                                   b = pa[id.training.pa, ], xy = T))
      test.pa <-
        na.omit(dismo::prepareData(bio.crop, p = ocor[-id.training.pa, ], 
                                   b = pa[-id.training.pa, ], xy = T))
      # test.back <-
      #   na.omit(dismo::prepareData(bio.crop, p = ocor[-id.training.pa, ], 
      #                              b = back[-id.training.b, ], xy = T))
      
      
      #### MODELING ###
      
      ##### Bioclim ####
      
      print(paste0(especie, " ","Bioclim"," ", "PA",  PA," ","RUN", RUN))
      
      bioclim_model <-
        dismo::bioclim(x = training.b[training.b[, "pb"] == 1,-1][,-c(1,2)])
      
      # dismo::response(bioclim_model)
      
      # Building Preojections ###
      
      # Current #
      bioclim_Cur <-
        raster::predict(bio.crop, bioclim_model)
      
      # Furure 45 #
      bioclim_Fut.45 <- preFut(rast = scenario.list.45,
                               model =  bioclim_model, GCM = GCM45s)
      
      # Furure 85 #
      bioclim_Fut.85 = preFut(rast = scenario.list.85,
                              model =  bioclim_model, GCM = GCM85s)
      
      
      bioclim_std = rescMod(bioclim_Cur, bioclim_Fut.45, bioclim_Fut.85)
      
      names(bioclim_std)<-c(paste0("Cur", ".",PA, ".", RUN), paste0("Fut.45", ".",PA, ".", RUN),
                            paste0("Fut.85", ".",PA, ".", RUN))
      
      bioclim_Cur <- subset(bioclim_std, grep(
        "Cur", names(bioclim_std)))
      writeRaster(
        bioclim_Cur, 
        paste0("./temp_output/",especie,"/", "bioclim_Cur","PA_",PA,"RUN_",RUN,".grd"),
        # formato = "GTiff",
        overwrite = T)
      
      
      bioclim_Fut.45 = subset(bioclim_std, grep(
        "Fut.45", names(bioclim_std)))
      writeRaster(
        bioclim_Fut.45,
        paste0("./temp_output/",especie,"/","bioclim_Fut.45","PA_",PA,"RUN_",RUN,".grd"),
        # formato = "GTiff",
        overwrite = T)
      
      bioclim_Fut.85= subset(bioclim_std, grep(
        "Fut.85", names(bioclim_std)))
      writeRaster(
        bioclim_Fut.85,
        paste0("./temp_output/",especie,"/","bioclim_Fut.85","PA_",PA,"RUN_",RUN,".grd"),
        # formato = "GTiff",
        overwrite = T)
      
      # Evaluating ###
      
      bioclim_eval <-
        eval.All.Model(subset(bioclim_std, grep("Cur", names(bioclim_std))), 
                       test.pa)
      
      bioclim_th.spec_sens <-
        dismo::threshold(bioclim_eval, "spec_sens")
      
      # bioclim_th.LPT5 <-
      #   quantile(raster::extract(bioclim_Cur, ocor[id.training.pa,]), 0.05, na.rm = T)
      # 
      # bioclim_th.VDl <-
      #   com.vdl(bioclim_Cur, test.pa, bioclim_eval)@maxVDl
      
      
      bioclim_eval.spec_sens <-
        eval.All.Model(subset(bioclim_std, grep("Cur", names(bioclim_std))), test.pa,
                       tr = bioclim_th.spec_sens)
      # bioclim_eval.LPT5 <-
      #   eval.All.Model(bioclim_Cur, test.pa,
      #                  tr = bioclim_th.LPT5)
      # bioclim_eval.VDl <-
      #   eval.All.Model(bioclim_Cur, test.pa, tr = bioclim_th.VDl)
      
      
      bioclim_e.spec_sens <- c(
        AUC.SpecSens = bioclim_eval.spec_sens@auc,
        TPR.SpecSens = bioclim_eval.spec_sens@TPR,
        TNR.SpecSens = bioclim_eval.spec_sens@TNR,
        thr.SpecSens = bioclim_th.spec_sens,
        # d.SpecSens = bioclim_eval.spec_sens@TPR * (1 - colSums(
        #   as.data.frame(
        #     bioclim_Cur >= bioclim_th.spec_sens,
        #     xy = F,
        #     na.rm = T
        #   )
        # ) / areaToral),
        TSS.SpecSens = (bioclim_eval.spec_sens@TPR + bioclim_eval.spec_sens@TNR)-1
        # TSSSpecSens = (
        #   bioclim_eval.spec_sens@TPR + bioclim_eval.spec_sens@TNR
        # ) - 1
      )
      
      # bioclim_e.LPT5 <- c(
      #   AUC.LPT5 = bioclim_eval.LPT5@auc,
      #   TPR.LPT5 = bioclim_eval.LPT5@TPR,
      #   TNR.LPT5 = bioclim_eval.LPT5@TNR,
      #   thr.LPT5 = bioclim_th.LPT5,
      #   d.LPT5 = bioclim_eval.LPT5@TPR * (1 - colSums(
      #     as.data.frame(
      #       bioclim_Cur >= bioclim_th.LPT5,
      #       xy = F,
      #       na.rm = T
      #     )
      #   ) / areaToral),
      #   TSS.LPT5 = bioclim_eval.LPT5@TSS
      #   # TSSLPT5 = (bioclim_eval.LPT5@TPR + bioclim_eval.LPT5@TNR) - 1
      # )
      # 
      # 
      # bioclim_e.VDl <- c(
      #   AUC.VDl = bioclim_eval.VDl@auc,
      #   TPR.VDl = bioclim_eval.VDl@TPR,
      #   TNR.VDl = bioclim_eval.VDl@TNR,
      #   thr.VDl = bioclim_th.VDl,
      #   d.VDl = bioclim_eval.VDl@TPR * (1 - colSums(
      #     as.data.frame(
      #       bioclim_Cur >= bioclim_th.VDl,
      #       xy = F,
      #       na.rm = T
      #     )
      #   ) / areaToral),
      #   TSS.VDl = bioclim_eval.VDl@TSS
      #   # TSSVDl = (bioclim_eval.VDl@TPR + bioclim_eval.VDl@TNR) - 1
      # )
      bioclim.e = rbind(bioclim.e, c(bioclim_e.spec_sens, PA = PA, RUN = RUN))#, bioclim_e.LPT5, bioclim_e.VDl))
      
      rownames(bioclim.e) = rep(paste0("bioclim"), 
                                nrow(bioclim.e))
      write.csv(bioclim.e,
                paste0("./temp_output/",especie,"/", "bioclim_eval.all.csv"), row.names = T)
      
      
      if(VarImport == T){
        
        bioclim_var.part.o  = foreach::foreach(
          vari = 1:length(vars),
          .combine = rbind,.packages = c("dismo"),.inorder = T) %dopar% {
            rmo = as.integer(subset(var, varia == vars[vari], select = pos)[[1]])
            
            bioclim_mod <-
              dismo::bioclim(x = bio.crop[[-(rmo)]] ,
                             p = training.pa[training.pa[, "pb"] == 1,c(1:2)])
            
            bioclim_VarImp1 <-(raster::predict(object = bio.crop[[-rmo]], 
                                               model = bioclim_mod))
            
            bioclim_VarImp1 = rescMod.One(bioclim_VarImp1)
            
            
            bioclim_obs <-
              as.data.frame(bioclim_Cur, xy = T)[,3]
            
            bioclim_sim <- 
              as.data.frame(bioclim_VarImp1, xy = T)[,3]
            
            er = na.omit(data.frame(sim = bioclim_sim, obs = bioclim_obs))
            
            bioclim_RMSE.VarImp1 <-
              rbind(MLmetrics::RMSE(er$sim, 
                                    er$obs))
            
            bioclim_NicheOverlap.VarImp1 <-
              rbind(1 - dismo::nicheOverlap(bioclim_Cur, 
                                            bioclim_VarImp1))
            
            cbind(bioclim_RMSE.VarImp1,bioclim_NicheOverlap.VarImp1
            )
            
          }
        bioclim_RMSE.name<-bioclim_NicheOverlap.name<-NULL
        
        for (vari in 1:length(vars)) {
          bioclim_RMSE.name  <-c(bioclim_RMSE.name,
                                 paste0(vars[vari],"_", "RMSE_", vari))
          
          bioclim_NicheOverlap.name  <-c(bioclim_NicheOverlap.name,
                                         paste0(vars[vari],"_", "NicheOvelap_", vari))
        }
        namesRow = c(bioclim_RMSE.name, bioclim_NicheOverlap.name)
        
        
        bioclim_var.part.o <-
          data.frame(metrics = namesRow,
                     imp = c(bioclim_var.part.o[,1], bioclim_var.part.o[,2]), 
                     PA = PA, RUN = RUN, alg = "bioclim")
        
        bioclim.var.All<-rbind(bioclim.var.All, 
                               bioclim_var.part.o)
      } else { next }
      
      
      remove(list = ls()[c(grep(
        "bioclim_", as.factor(ls())))])
      
      gc(reset = TRUE, full = T)
      
      # ##### Domain  ####
      # print(paste0(especie, " ","Domain"," ",  "PA",  PA," ","RUN", RUN))
      # 
      # domain_model <-
      #   domain(x = training.b[training.b[, "pb"] == 1, -c(1:3)], )
      # # dismo::response(domain_model)
      # 
      # # Building Preojections ###
      # 
      # # Current #
      # domain_Cur <-
      #   raster::predict(bio.crop, domain_model)
      # 
      # # Furure 45 #
      # domain_Fut.45 <- preFut(rast = scenario.list.45,
      #                        model =  domain_model, GCM = GCM45s)
      # 
      # # Furure 85 #
      # domain_Fut.85 = preFut(rast = scenario.list.85,
      #                       model =  domain_model, GCM = GCM85s)
      # 
      # domain_std = rescMod(domain_Cur, mean(domain_Fut.45), mean(domain_Fut.85))
      # 
      # names(domain_std)<-c(paste0("Cur", ".",PA, ".", RUN), paste0("Fut.45", ".",PA, ".", RUN),
      #                      paste0("Fut.85", ".",PA, ".", RUN))
      # 
      # domain_std = rescMod(domain_Cur, mean(domain_Fut.45), mean(domain_Fut.85))
      # 
      # domain_Cur <- subset(domain_std, grep(
      #   "Cur", names(domain_std)))
      # writeRaster(
      #   domain_Cur, 
      #   paste0("./temp_output/",especie,"/", "domain_Cur","PA_",PA,"RUN_",RUN,".grd"),
      #   # formato = "GTiff",
      #   overwrite = T)
      # 
      # 
      # domain_Fut.45 = subset(domain_std, grep(
      #   "Fut.45", names(domain_std)))
      # writeRaster(
      #   domain_Fut.45,
      #   paste0("./temp_output/",especie,"/","domain_Fut.45","PA_",PA,"RUN_",RUN,".grd"),
      #   # formato = "GTiff",
      #   overwrite = T)
      # 
      # domain_Fut.85= subset(domain_std, grep(
      #   "Fut.85", names(domain_std)))
      # writeRaster(
      #   domain_Fut.85,
      #   paste0("./temp_output/",especie,"/","domain_Fut.85","PA_",PA,"RUN_",RUN,".grd"),
      #   # formato = "GTiff",
      #   overwrite = T)
      #
      # # Evaluating ###
      # 
      # domain_eval <-
      #   eval.All.Model(subset(domain_std, grep("Cur", names(domain_std))),
      #                  test.pa)
      # 
      # domain_th.spec_sens <-
      #   dismo::threshold(domain_eval, "spec_sens")
      # 
      # # domain_th.LPT5 <-
      # #   quantile(raster::extract(domain_Cur, ocor[id.training.pa,]), 0.05, na.rm = T)
      # #
      # # domain_th.VDl <-
      # #   com.vdl(domain_Cur, test.pa, domain_eval)@maxVDl
      # 
      # 
      # domain_eval.spec_sens <-
      #   eval.All.Model(subset(domain_std, grep("Cur", names(domain_std))), test.pa,
      #                  tr = domain_th.spec_sens)
      # # domain_eval.LPT5 <-
      # #   eval.All.Model(domain_Cur, test.pa,
      # #                  tr = domain_th.LPT5)
      # # domain_eval.VDl <-
      # #   eval.All.Model(domain_Cur, test.pa, tr = domain_th.VDl)
      # 
      # 
      # domain_e.spec_sens <- c(
      #   AUC.SpecSens = domain_eval.spec_sens@auc,
      #   TPR.SpecSens = domain_eval.spec_sens@TPR,
      #   TNR.SpecSens = domain_eval.spec_sens@TNR,
      #   thr.SpecSens = domain_th.spec_sens,
      #   # d.SpecSens = domain_eval.spec_sens@TPR * (1 - colSums(
      #   #   as.data.frame(
      #   #     domain_Cur >= domain_th.spec_sens,
      #   #     xy = F,
      #   #     na.rm = T
      #   #   )
      #   # ) / areaToral),
      #   TSS.SpecSens = (domain_eval.spec_sens@TPR + domain_eval.spec_sens@TNR)-1
      #   # TSSSpecSens = (
      #   #   domain_eval.spec_sens@TPR + domain_eval.spec_sens@TNR
      #   # ) - 1
      # )
      # 
      # # domain_e.LPT5 <- c(
      # #   AUC.LPT5 = domain_eval.LPT5@auc,
      # #   TPR.LPT5 = domain_eval.LPT5@TPR,
      # #   TNR.LPT5 = domain_eval.LPT5@TNR,
      # #   thr.LPT5 = domain_th.LPT5,
      # #   d.LPT5 = domain_eval.LPT5@TPR * (1 - colSums(
      # #     as.data.frame(
      # #       domain_Cur >= domain_th.LPT5,
      # #       xy = F,
      # #       na.rm = T
      # #     )
      # #   ) / areaToral),
      # #   TSS.LPT5 = domain_eval.LPT5@TSS
      # #   # TSSLPT5 = (domain_eval.LPT5@TPR + domain_eval.LPT5@TNR) - 1
      # # )
      # #
      # #
      # # domain_e.VDl <- c(
      # #   AUC.VDl = domain_eval.VDl@auc,
      # #   TPR.VDl = domain_eval.VDl@TPR,
      # #   TNR.VDl = domain_eval.VDl@TNR,
      # #   thr.VDl = domain_th.VDl,
      # #   d.VDl = domain_eval.VDl@TPR * (1 - colSums(
      # #     as.data.frame(
      # #       domain_Cur >= domain_th.VDl,
      # #       xy = F,
      # #       na.rm = T
      # #     )
      # #   ) / areaToral),
      # #   TSS.VDl = domain_eval.VDl@TSS
      # #   # TSSVDl = (domain_eval.VDl@TPR + domain_eval.VDl@TNR) - 1
      # # )
      # domain.e = rbind(domain.e, c(domain_e.spec_sens, PA = PA, RUN = RUN))#, domain_e.LPT5, domain_e.VDl))
      # 
      # rownames(domain.e) = rep(paste0("domain"),
      #                          nrow(domain.e))
      # write.csv(domain.e,
      #           paste0("./temp_output/",especie,"/", "domain_eval.all.csv"), row.names = T)
      # 
      # domain.e = rbind(domain.e, c(domain_e.spec_sens, PA = PA, RUN = RUN))#, domain_e.LPT5, domain_e.VDl))
      # 
      # rownames(domain.e) = rep(paste0("domain"), 
      #                          nrow(domain.e))
      # write.csv(domain.e,
      #           paste0("./temp_output/",especie,"/", "domain_eval.all.csv"), row.names = T)
      #       if(VarImport == T){
      # domain_var.part.o  = foreach::foreach(
      #   vari = 1:length(vars),
      #   .combine = rbind,.packages = c("dismo"),.inorder = T) %dopar% {
      #     rmo = as.integer(subset(var, varia == vars[vari], select = pos)[[1]])
      #     
      #     domain_mod <-
      #       dismo::domain(x = bio.crop[[-(rmo)]] ,
      #                     p = training.pa[training.pa[, "pb"] == 1,c(1:2)])
      #     
      #     domain_VarImp1 <-(raster::predict(object = bio.crop[[-rmo]], 
      #                                       model = domain_mod))
      #     
      #     
      #     domain_obs <-
      #       as.data.frame(domain_Cur, xy = T)[,3]
      #     
      #     domain_sim <- 
      #       as.data.frame(domain_VarImp1, xy = T)[,3]
      #     
      #     er = na.omit(data.frame(sim = domain_sim, obs = domain_obs))
      #     
      #     domain_RMSE.VarImp1 <-
      #       rbind(MLmetrics::RMSE(er$sim, 
      #                             er$obs))
      #     
      # 
      # domain_NicheOverlap.VarImp1 <-
      #   rbind(1 - dismo::nicheOverlap(domain_Cur, 
      #                                 domain_VarImp1))
      #     
      #     
      #     cbind(domain_RMSE.VarImp1,domain_NicheOverlap.VarImp1
      #     )
      #     
      #   }
      # domain_RMSE.name<-domain_NicheOverlap.name<-NULL
      # 
      # for (vari in 1:length(vars)) {
      #   domain_RMSE.name  <-c(domain_RMSE.name,
      #                         paste0(vars[vari],"_", "RMSE_", vari))
      #   
      #   domain_NicheOverlap.name  <-c(domain_NicheOverlap.name,
      #                        paste0(vars[vari],"_", "Overlap_", vari))
      # }
      # namesRow = c(domain_RMSE.name, domain_NicheOverlap.name)
      # 
      # 
      # domain_var.part.o <-
      #   data.frame(metrics = namesRow,
      #              imp = c(domain_var.part.o[,1], domain_var.part.o[,2]), 
      #              PA = PA, RUN = RUN, alg = "domain")
      # 
      # domain.var.All<-rbind(domain.var.All, 
      #                       domain_var.part.o) } else { next }
      # 
      # remove(list = ls()[c(grep(
      #   "domain_", as.factor(ls())))])
      # 
      # gc(reset = TRUE, full = T)
      
      
      ##### ENFA  #### 
      print(paste0(especie, " ","ENFA"," ", "PA",  PA," ","RUN", RUN))
      
      climaPres <- raster::stack(bio.crop)
      # names(climaPres)# <- paste0("PC",1:6)
      
      climaPres.values <- raster::values(climaPres)
      climaPres.spdf <- na.omit(data.frame(xyFromCell(climaPres, 1:ncell(climaPres)), 
                                           climaPres.values))
      
      suppressWarnings(gridded(climaPres.spdf) <- ~x+y)
      climaPres <- raster::stack(climaPres.spdf)
      climaPres.values <- raster::values(climaPres)
      media.climaPres <- apply(slot(climaPres.spdf, "data"), 2, mean)
      sd.climaPres <- apply(slot(climaPres.spdf, "data"), 2, sd)
      climaPres.scale<- sweep(slot(climaPres.spdf, "data"),2, media.climaPres)
      climaPres.scale<- as.matrix(climaPres.scale) %*% diag(1/sd.climaPres)
      
      #adjustment of the ENFA model
      
      pr.cell <- raster::extract(climaPres, training.pa[training.pa[,"pb"]==1,1:2], cellnumber=T)
      pr <- data.frame(pr= rep(0, ncell(climaPres)), climaPres.values)
      pr[pr.cell[,"cells"], 1] <- 1
      pr <- na.omit(pr)
      pr <- pr[,1]
      enfa_model <- adehabitatHS::madifa(ade4::dudi.pca(climaPres.scale, 
                                                        center=F, scale=F, scannf=F), 
                                         pr, scannf=F)
      
      # Current #
      enfa_Cur <-  preFut(rast = bio.crop,rast1 = bio.crop,
                          model =  enfa_model)
      
      # Furure 45 #
      enfa_Fut.45 <- preFut(rast = bio.crop,rast1 = scenario.list.45,
                            model =  enfa_model, GCM = GCM45s)
      # Furure 85 #
      enfa_Fut.85 <- preFut(rast = bio.crop,rast1 = scenario.list.85,
                            model =  enfa_model, GCM = GCM85s)
      
      enfa_std = rescMod(enfa_Cur, enfa_Fut.45, enfa_Fut.85)
      
      names(enfa_std)<-c(paste0("Cur", ".",PA, ".", RUN), paste0("Fut.45", ".",PA, ".", RUN),
                         paste0("Fut.85", ".",PA, ".", RUN))
      
      enfa_Cur <- subset(enfa_std, grep(
        "Cur", names(enfa_std)))
      writeRaster(
        enfa_Cur, 
        paste0("./temp_output/",especie,"/", "enfa_Cur","PA_",PA,"RUN_",RUN,".grd"),
        # formato = "GTiff",
        overwrite = T)
      
      
      enfa_Fut.45 = subset(enfa_std, grep(
        "Fut.45", names(enfa_std)))
      writeRaster(
        enfa_Fut.45,
        paste0("./temp_output/",especie,"/","enfa_Fut.45","PA_",PA,"RUN_",RUN,".grd"),
        # formato = "GTiff",
        overwrite = T)
      
      enfa_Fut.85= subset(enfa_std, grep(
        "Fut.85", names(enfa_std)))
      writeRaster(
        enfa_Fut.85,
        paste0("./temp_output/",especie,"/","enfa_Fut.85","PA_",PA,"RUN_",RUN,".grd"),
        # formato = "GTiff",
        overwrite = T)
      
      # Evaluating ###
      
      enfa_eval <-
        eval.All.Model(subset(enfa_std, grep("Cur", names(enfa_std))), 
                       test.pa)
      
      enfa_th.spec_sens <-
        dismo::threshold(enfa_eval, "spec_sens")
      
      # enfa_th.LPT5 <-
      #   quantile(raster::extract(enfa_Cur, ocor[id.training.pa,]), 0.05, na.rm = T)
      # 
      # enfa_th.VDl <-
      #   com.vdl(enfa_Cur, test.pa, enfa_eval)@maxVDl
      
      
      enfa_eval.spec_sens <-
        eval.All.Model(subset(enfa_std, grep("Cur", names(enfa_std))), test.pa,
                       tr = enfa_th.spec_sens)
      # enfa_eval.LPT5 <-
      #   eval.All.Model(enfa_Cur, test.pa,
      #                  tr = enfa_th.LPT5)
      # enfa_eval.VDl <-
      #   eval.All.Model(enfa_Cur, test.pa, tr = enfa_th.VDl)
      
      
      enfa_e.spec_sens <- c(
        AUC.SpecSens = enfa_eval.spec_sens@auc,
        TPR.SpecSens = enfa_eval.spec_sens@TPR,
        TNR.SpecSens = enfa_eval.spec_sens@TNR,
        thr.SpecSens = enfa_th.spec_sens,
        # d.SpecSens = enfa_eval.spec_sens@TPR * (1 - colSums(
        #   as.data.frame(
        #     enfa_Cur >= enfa_th.spec_sens,
        #     xy = F,
        #     na.rm = T
        #   )
        # ) / areaToral),
        TSS.SpecSens = (enfa_eval.spec_sens@TPR + enfa_eval.spec_sens@TNR)-1
        # TSSSpecSens = (
        #   enfa_eval.spec_sens@TPR + enfa_eval.spec_sens@TNR
        # ) - 1
      )
      
      # enfa_e.LPT5 <- c(
      #   AUC.LPT5 = enfa_eval.LPT5@auc,
      #   TPR.LPT5 = enfa_eval.LPT5@TPR,
      #   TNR.LPT5 = enfa_eval.LPT5@TNR,
      #   thr.LPT5 = enfa_th.LPT5,
      #   d.LPT5 = enfa_eval.LPT5@TPR * (1 - colSums(
      #     as.data.frame(
      #       enfa_Cur >= enfa_th.LPT5,
      #       xy = F,
      #       na.rm = T
      #     )
      #   ) / areaToral),
      #   TSS.LPT5 = enfa_eval.LPT5@TSS
      #   # TSSLPT5 = (enfa_eval.LPT5@TPR + enfa_eval.LPT5@TNR) - 1
      # )
      # 
      # 
      # enfa_e.VDl <- c(
      #   AUC.VDl = enfa_eval.VDl@auc,
      #   TPR.VDl = enfa_eval.VDl@TPR,
      #   TNR.VDl = enfa_eval.VDl@TNR,
      #   thr.VDl = enfa_th.VDl,
      #   d.VDl = enfa_eval.VDl@TPR * (1 - colSums(
      #     as.data.frame(
      #       enfa_Cur >= enfa_th.VDl,
      #       xy = F,
      #       na.rm = T
      #     )
      #   ) / areaToral),
      #   TSS.VDl = enfa_eval.VDl@TSS
      #   # TSSVDl = (enfa_eval.VDl@TPR + enfa_eval.VDl@TNR) - 1
      # )
      enfa.e = rbind(enfa.e, c(enfa_e.spec_sens, PA = PA, RUN = RUN))#, enfa_e.LPT5, enfa_e.VDl))
      
      rownames(enfa.e) = rep(paste0("enfa"), 
                             nrow(enfa.e))
      
      write.csv(enfa.e,
                paste0("./temp_output/",especie,"/", "enfa_eval.all.csv"), row.names = T)
      
      if(VarImport == T){
        
        enfa_var.part.o  = foreach::foreach(
          vari = 1:length(vars),
          .combine = rbind,.packages = c("dismo", "foreach", "adehabitatHS"),
          .inorder = T) %dopar% {
            rmo = as.integer(subset(var, varia == vars[vari], select = pos)[[1]])
            
            
            enfa_climaPr <- raster::stack(bio.crop[[-rmo]])
            names(enfa_climaPr)# <- paste0("PC",1:6)
            
            enfa_climaPr.values <- raster::values(enfa_climaPr)
            enfa_climaPr.spdf <-
              na.omit(data.frame(xyFromCell(enfa_climaPr, 1:ncell(enfa_climaPr)), 
                                 enfa_climaPr.values))
            
            enfa_climaPr<-rasterize(enfa_climaPr.spdf[,1:2],enfa_Cur[[1]],
                                    enfa_climaPr.spdf[,-c(1:2)])
            
            suppressWarnings(gridded(enfa_climaPr.spdf)<-~x+y)
            
            enfa_climaPr <- raster::stack(enfa_climaPr.spdf)
            enfa_climaPr.values <- raster::values(enfa_climaPr)
            media.enfa_climaPr <- apply(slot(enfa_climaPr.spdf, "data"), 2, mean)
            sd.enfa_climaPr <- apply(slot(enfa_climaPr.spdf, "data"), 2, sd)
            enfa_climaPr.scale<- sweep(slot(enfa_climaPr.spdf, "data"),2, 
                                       media.enfa_climaPr)
            enfa_climaPr.scale <-
              as.matrix(enfa_climaPr.scale) %*% diag(1 / sd.enfa_climaPr)
            
            enfa_pr.cel <- raster::extract(enfa_climaPr, ocor, cellnumber=T)
            enfa_pr1 <-
              data.frame(enfa_pr1 = rep(0, ncell(enfa_climaPr)), enfa_climaPr.values)
            enfa_pr1[enfa_pr.cel[,"cells"], 1] <- 1
            enfa_pr1 <- na.omit(enfa_pr1)
            enfa_pr1 <- enfa_pr1[,1]
            
            enfa_mod <- adehabitatHS::madifa(ade4::dudi.pca(enfa_climaPr.scale, 
                                                            center=F, scale=F, scannf=F), 
                                             enfa_pr1, scannf=F)

            
            enfa_VarImp1 <-
              extend(rescMod.One(preFut(model = enfa_mod, 
                                        enfa_climaPr, enfa_climaPr)),
                     enfa_Cur@extent)
            
            
            enfa_obs <-
              as.data.frame(enfa_Cur, xy = T)[,3]
            
            enfa_sim <- 
              as.data.frame(enfa_VarImp1, xy = T)[,3]
            
            er = na.omit(data.frame(sim = enfa_sim, obs = enfa_obs))
            
            enfa_RMSE.VarImp1 <-
              rbind(MLmetrics::RMSE(er$sim, 
                                    er$obs))
            
            
            enfa_NicheOverlap.VarImp1 <-
              rbind(1 - dismo::nicheOverlap(enfa_Cur, 
                                            enfa_VarImp1))
            
            cbind(enfa_RMSE.VarImp1,enfa_NicheOverlap.VarImp1
            )
            
          }
        enfa_RMSE.name<-enfa_NicheOverlap.name<-NULL
        
        for (vari in 1:length(vars)) {
          enfa_RMSE.name  <-c(enfa_RMSE.name,
                              paste0(vars[vari],"_", "RMSE_", vari))
          
          enfa_NicheOverlap.name  <-c(enfa_NicheOverlap.name,
                                      paste0(vars[vari],"_", "Overlap_", vari))
        }
        namesRow = c(enfa_RMSE.name, enfa_NicheOverlap.name)
        
        
        enfa_var.part.o <-
          data.frame(metrics = namesRow,
                     imp = c(enfa_var.part.o[,1], enfa_var.part.o[,2]), 
                     PA = PA, RUN = RUN, alg = "enfa")
        
        enfa.var.All<-rbind(enfa.var.All, 
                            enfa_var.part.o)
      } else{ next }
      
      rm(climaPres.values,climaPres.spdf, media.climaPres,sd.climaPres,
         climaPres.scale, pr.cell,pr)
      remove(list = ls()[c(grep(
        "enfa_", as.factor(ls())))])
      gc(reset = TRUE, full = T)
      
      
      ##### GLM  ####
      print(paste0(especie, " ","GLM"," ",  "PA",  PA," ","RUN", RUN))
      
      pb <- training.pa$pb
      
      glm_model <- glm(
        formula = pb ~ .,
        data = training.pa[, -c(1:3)],
        type = 'quadratic',
        interaction.level = 0,
        myFormula = NULL,
        test = 'AIC',
        family = binomial(link = 'logit'),
        control = glm.control(
          epsilon = 1e-08,
          maxit = 50,
          trace = FALSE
        )
      )
      
      # Building Preojections ###
      
      # Current #
      glm_Cur <-
        raster::predict(bio.crop, glm_model)
      
      # Furure 45 #
      glm_Fut.45 <- preFut(rast = scenario.list.45,
                           model =  glm_model, GCM = GCM45s)
      
      # Furure 85 #
      glm_Fut.85 = preFut(rast = scenario.list.85,
                          model =  glm_model, GCM = GCM85s)
      
      glm_std = rescMod(glm_Cur, glm_Fut.45, glm_Fut.85)
      
      names(glm_std)<-c(paste0("Cur", ".",PA, ".", RUN), paste0("Fut.45", ".",PA, ".", RUN),
                        paste0("Fut.85", ".",PA, ".", RUN))
      
      glm_Cur <- subset(glm_std, grep(
        "Cur", names(glm_std)))
      writeRaster(
        glm_Cur, 
        paste0("./temp_output/",especie,"/", "glm_Cur","PA_",PA,"RUN_",RUN,".grd"),
        # formato = "GTiff",
        overwrite = T)
      
      
      glm_Fut.45 = subset(glm_std, grep(
        "Fut.45", names(glm_std)))
      writeRaster(
        glm_Fut.45,
        paste0("./temp_output/",especie,"/","glm_Fut.45","PA_",PA,"RUN_",RUN,".grd"),
        # formato = "GTiff",
        overwrite = T)
      
      glm_Fut.85= subset(glm_std, grep(
        "Fut.85", names(glm_std)))
      writeRaster(
        glm_Fut.85,
        paste0("./temp_output/",especie,"/","glm_Fut.85","PA_",PA,"RUN_",RUN,".grd"),
        # formato = "GTiff",
        overwrite = T)
      
      # Evaluating ###
      
      glm_eval <-
        eval.All.Model(subset(glm_std, grep("Cur", names(glm_std))), 
                       test.pa)
      
      glm_th.spec_sens <-
        dismo::threshold(glm_eval, "spec_sens")
      
      # glm_th.LPT5 <-
      #   quantile(raster::extract(glm_Cur, ocor[id.training.pa,]), 0.05, na.rm = T)
      # 
      # glm_th.VDl <-
      #   com.vdl(glm_Cur, test.pa, glm_eval)@maxVDl
      
      
      glm_eval.spec_sens <-
        eval.All.Model(subset(glm_std, grep("Cur", names(glm_std))), test.pa,
                       tr = glm_th.spec_sens)
      # glm_eval.LPT5 <-
      #   eval.All.Model(glm_Cur, test.pa,
      #                  tr = glm_th.LPT5)
      # glm_eval.VDl <-
      #   eval.All.Model(glm_Cur, test.pa, tr = glm_th.VDl)
      
      
      glm_e.spec_sens <- c(
        AUC.SpecSens = glm_eval.spec_sens@auc,
        TPR.SpecSens = glm_eval.spec_sens@TPR,
        TNR.SpecSens = glm_eval.spec_sens@TNR,
        thr.SpecSens = glm_th.spec_sens,
        # d.SpecSens = glm_eval.spec_sens@TPR * (1 - colSums(
        #   as.data.frame(
        #     glm_Cur >= glm_th.spec_sens,
        #     xy = F,
        #     na.rm = T
        #   )
        # ) / areaToral),
        TSS.SpecSens = (glm_eval.spec_sens@TPR + glm_eval.spec_sens@TNR)-1
        # TSSSpecSens = (
        #   glm_eval.spec_sens@TPR + glm_eval.spec_sens@TNR
        # ) - 1
      )
      
      # glm_e.LPT5 <- c(
      #   AUC.LPT5 = glm_eval.LPT5@auc,
      #   TPR.LPT5 = glm_eval.LPT5@TPR,
      #   TNR.LPT5 = glm_eval.LPT5@TNR,
      #   thr.LPT5 = glm_th.LPT5,
      #   d.LPT5 = glm_eval.LPT5@TPR * (1 - colSums(
      #     as.data.frame(
      #       glm_Cur >= glm_th.LPT5,
      #       xy = F,
      #       na.rm = T
      #     )
      #   ) / areaToral),
      #   TSS.LPT5 = glm_eval.LPT5@TSS
      #   # TSSLPT5 = (glm_eval.LPT5@TPR + glm_eval.LPT5@TNR) - 1
      # )
      # 
      # 
      # glm_e.VDl <- c(
      #   AUC.VDl = glm_eval.VDl@auc,
      #   TPR.VDl = glm_eval.VDl@TPR,
      #   TNR.VDl = glm_eval.VDl@TNR,
      #   thr.VDl = glm_th.VDl,
      #   d.VDl = glm_eval.VDl@TPR * (1 - colSums(
      #     as.data.frame(
      #       glm_Cur >= glm_th.VDl,
      #       xy = F,
      #       na.rm = T
      #     )
      #   ) / areaToral),
      #   TSS.VDl = glm_eval.VDl@TSS
      #   # TSSVDl = (glm_eval.VDl@TPR + glm_eval.VDl@TNR) - 1
      # )
      glm.e = rbind(glm.e, c(glm_e.spec_sens, PA = PA, RUN = RUN))#, glm_e.LPT5, glm_e.VDl))
      
      rownames(glm.e) = rep(paste0("glm"), 
                            nrow(glm.e))
      
      write.csv(glm.e,
                paste0("./temp_output/",especie,"/", "glm_eval.all.csv"), row.names = T)
      
      if(VarImport == T){
        
        glm_var.part.o  = foreach::foreach(
          vari = 1:length(vars),
          .combine = rbind,.packages = c("dismo", "foreach", "adehabitatHS"),
          .inorder = T) %dopar% {
            rmo = as.integer(subset(var, varia == vars[vari], select = pos)[[1]])
            
            glm_Data = training.pa[, -c(1:3)][, -rmo]  
            
            glm_mod <- glm(
              formula = pb ~ .,
              data = glm_Data,
              type = 'quadratic',
              interaction.level = 0,
              myFormula = NULL,
              test = 'AIC',
              family = binomial(link = 'logit'),
              control = glm.control(
                epsilon = 1e-08,
                maxit = 50,
                trace = FALSE
              )
            )
            glm_PCA <- raster::stack(bio.crop[[-rmo]])
            
            glm_VarImp1 <-rescMod.One(raster::predict(glm_PCA, glm_mod))
            
            
            glm_obs <-
              as.data.frame(glm_Cur, xy = T)[,3]
            
            glm_sim <- 
              as.data.frame(glm_VarImp1, xy = T)[,3]
            
            er = na.omit(data.frame(sim = glm_sim, obs = glm_obs))
            
            glm_RMSE.VarImp1 <-
              rbind(MLmetrics::RMSE(er$sim, 
                                    er$obs))
            
            
            glm_NicheOverlap.VarImp1 <-
              rbind(1 - dismo::nicheOverlap(glm_Cur, 
                                            glm_VarImp1))
            
            cbind(glm_RMSE.VarImp1,glm_NicheOverlap.VarImp1
            )
            
          }
        glm_RMSE.name<-glm_NicheOverlap.name<-NULL
        
        for (vari in 1:length(vars)) {
          glm_RMSE.name  <-c(glm_RMSE.name,
                             paste0(vars[vari],"_", "RMSE_", vari))
          
          glm_NicheOverlap.name  <-c(glm_NicheOverlap.name,
                                     paste0(vars[vari],"_", "Overlap_", vari))
        }
        namesRow = c(glm_RMSE.name, glm_NicheOverlap.name)
        
        
        glm_var.part.o <-
          data.frame(metrics = namesRow,
                     imp = c(glm_var.part.o[,1], glm_var.part.o[,2]), 
                     PA = PA, RUN = RUN, alg = "glm")
        
        glm.var.All<-rbind(glm.var.All, 
                           glm_var.part.o)
      } else {next}
      
      remove(list = ls()[c(grep(
        "glm_", as.factor(ls())))])
      gc(reset = TRUE, full = T)
      
      ##### MARS ####
      print(paste0(especie, " ","MARS"," ", "PA",  PA," ","RUN", RUN))
      
      mars_model <- earth::earth(
        pb ~ .,
        data = training.pa[,-c(1:2)],
        # type = 'simple',
        # interaction.level = 0,
        penalty = 1,
        nprune = NULL,
        pmethod = 'backward'
      )
      
      # Building Preojections ###
      
      # Current #
      mars_Cur <-
        raster::predict(bio.crop, mars_model)
      
      # Furure 45 #
      mars_Fut.45 <- preFut(rast = scenario.list.45,
                            model =  mars_model, GCM = GCM45s)
      
      # Furure 85 #
      mars_Fut.85 = preFut(rast = scenario.list.85,
                           model =  mars_model, GCM = GCM85s)
      
      mars_std = rescMod(mars_Cur, mars_Fut.45, mars_Fut.85)
      
      names(mars_std)<-c(paste0("Cur", ".",PA, ".", RUN), paste0("Fut.45", ".",PA, ".", RUN),
                         paste0("Fut.85", ".",PA, ".", RUN))
      
      mars_Cur <- subset(mars_std, grep(
        "Cur", names(mars_std)))
      writeRaster(
        mars_Cur, 
        paste0("./temp_output/",especie,"/", "mars_Cur","PA_",PA,"RUN_",RUN,".grd"),
        # formato = "GTiff",
        overwrite = T)
      
      
      mars_Fut.45 = subset(mars_std, grep(
        "Fut.45", names(mars_std)))
      writeRaster(
        mars_Fut.45,
        paste0("./temp_output/",especie,"/","mars_Fut.45","PA_",PA,"RUN_",RUN,".grd"),
        # formato = "GTiff",
        overwrite = T)
      
      mars_Fut.85 = subset(mars_std, grep(
        "Fut.85", names(mars_std)))
      writeRaster(
        mars_Fut.85,
        paste0("./temp_output/",especie,"/","mars_Fut.85","PA_",PA,"RUN_",RUN,".grd"),
        # formato = "GTiff",
        overwrite = T)
      
      # Evaluating ###
      
      mars_eval <-
        eval.All.Model(subset(mars_std, grep("Cur", names(mars_std))), 
                       test.pa)
      
      mars_th.spec_sens <-
        dismo::threshold(mars_eval, "spec_sens")
      
      # mars_th.LPT5 <-
      #   quantile(raster::extract(mars_Cur, ocor[id.training.pa,]), 0.05, na.rm = T)
      # 
      # mars_th.VDl <-
      #   com.vdl(mars_Cur, test.pa, mars_eval)@maxVDl
      
      
      mars_eval.spec_sens <-
        eval.All.Model(subset(mars_std, grep("Cur", names(mars_std))), test.pa,
                       tr = mars_th.spec_sens)
      # mars_eval.LPT5 <-
      #   eval.All.Model(mars_Cur, test.pa,
      #                  tr = mars_th.LPT5)
      # mars_eval.VDl <-
      #   eval.All.Model(mars_Cur, test.pa, tr = mars_th.VDl)
      
      
      mars_e.spec_sens <- c(
        AUC.SpecSens = mars_eval.spec_sens@auc,
        TPR.SpecSens = mars_eval.spec_sens@TPR,
        TNR.SpecSens = mars_eval.spec_sens@TNR,
        thr.SpecSens = mars_th.spec_sens,
        # d.SpecSens = mars_eval.spec_sens@TPR * (1 - colSums(
        #   as.data.frame(
        #     mars_Cur >= mars_th.spec_sens,
        #     xy = F,
        #     na.rm = T
        #   )
        # ) / areaToral),
        TSS.SpecSens = (mars_eval.spec_sens@TPR + mars_eval.spec_sens@TNR)-1
        # TSSSpecSens = (
        #   mars_eval.spec_sens@TPR + mars_eval.spec_sens@TNR
        # ) - 1
      )
      
      # mars_e.LPT5 <- c(
      #   AUC.LPT5 = mars_eval.LPT5@auc,
      #   TPR.LPT5 = mars_eval.LPT5@TPR,
      #   TNR.LPT5 = mars_eval.LPT5@TNR,
      #   thr.LPT5 = mars_th.LPT5,
      #   d.LPT5 = mars_eval.LPT5@TPR * (1 - colSums(
      #     as.data.frame(
      #       mars_Cur >= mars_th.LPT5,
      #       xy = F,
      #       na.rm = T
      #     )
      #   ) / areaToral),
      #   TSS.LPT5 = mars_eval.LPT5@TSS
      #   # TSSLPT5 = (mars_eval.LPT5@TPR + mars_eval.LPT5@TNR) - 1
      # )
      # 
      # 
      # mars_e.VDl <- c(
      #   AUC.VDl = mars_eval.VDl@auc,
      #   TPR.VDl = mars_eval.VDl@TPR,
      #   TNR.VDl = mars_eval.VDl@TNR,
      #   thr.VDl = mars_th.VDl,
      #   d.VDl = mars_eval.VDl@TPR * (1 - colSums(
      #     as.data.frame(
      #       mars_Cur >= mars_th.VDl,
      #       xy = F,
      #       na.rm = T
      #     )
      #   ) / areaToral),
      #   TSS.VDl = mars_eval.VDl@TSS
      #   # TSSVDl = (mars_eval.VDl@TPR + mars_eval.VDl@TNR) - 1
      # )
      mars.e = rbind(mars.e, c(mars_e.spec_sens, PA = PA, RUN = RUN))#, mars_e.LPT5, mars_e.VDl))
      
      rownames(mars.e) = rep(paste0("mars"), 
                             nrow(mars.e))
      
      write.csv(mars.e,
                paste0("./temp_output/",especie,"/", "mars_eval.all.csv"), row.names = T)
      
      if(VarImport == T){
        
        mars_var.part.o  = foreach::foreach(
          vari = 1:length(vars),
          .combine = rbind,.packages = c("dismo", "foreach", "adehabitatHS"),
          .inorder = T) %dopar% {
            rmo = as.integer(subset(var, varia == vars[vari], select = pos)[[1]])
            
            mars_Data = training.pa[, -c(1:3)][, -rmo]  
            
            mars_mod <- earth::earth(
              formula = pb ~ .,
              data = mars_Data,
              # type = 'simple',
              # interaction.level = 0,
              penalty = 1,
              nprune = NULL,
              pmethod = 'backward'
            )
            mars_PCA <- raster::stack(bio.crop[[-rmo]])
            
            mars_VarImp1 <-rescMod.One(raster::predict(mars_PCA, mars_mod))
            
            
            mars_obs <-
              as.data.frame(mars_Cur, xy = T)[,3]
            
            mars_sim <- 
              as.data.frame(mars_VarImp1, xy = T)[,3]
            
            er = na.omit(data.frame(sim = mars_sim, obs = mars_obs))
            
            mars_RMSE.VarImp1 <-
              rbind(MLmetrics::RMSE(er$sim, 
                                    er$obs))
            
            mars_NicheOverlap.VarImp1 <-
              rbind(1 - dismo::nicheOverlap(mars_Cur, 
                                            mars_VarImp1))
            
            cbind(mars_RMSE.VarImp1,mars_NicheOverlap.VarImp1
            )
            
          }
        mars_RMSE.name<-mars_NicheOverlap.name<-NULL
        
        for (vari in 1:length(vars)) {
          mars_RMSE.name  <-c(mars_RMSE.name,
                              paste0(vars[vari],"_", "RMSE_", vari))
          
          mars_NicheOverlap.name  <-c(mars_NicheOverlap.name,
                                      paste0(vars[vari],"_", "Overlap_", vari))
        }
        namesRow = c(mars_RMSE.name, mars_NicheOverlap.name)
        
        
        mars_var.part.o <-
          data.frame(metrics = namesRow,
                     imp = c(mars_var.part.o[,1], mars_var.part.o[,2]), 
                     PA = PA, RUN = RUN, alg = "mars")
        
        mars.var.All<-rbind(mars.var.All, 
                            mars_var.part.o)
      } else { next }
      
      remove(list = ls()[c(grep(
        "mars_", as.factor(ls())))])
      gc(reset = TRUE, full = T)
      
      
      ##### Maxent #### -> ku.enm ####
      print(paste0(especie, " ","MAXENT"," ",  "PA",  PA," ","RUN", RUN))
      
      maxent_model <- quiet(maxent(
        x = training.b[, -c(1:3)],
        p = training.b[, 'pb'],
        memory_allocated = 512,
        background_data_dir = 'default',
        maximumbackground = 'default',
        maximumiterations = 200,
        visible = FALSE,
        linear = TRUE,
        quadratic = TRUE,
        product = TRUE,
        threshold = TRUE,
        hinge = TRUE,
        lq2lqptthreshold = 80,
        l2lqthreshold = 10,
        hingethreshold = 15,
        beta_threshold = -1,
        beta_categorical = -1,
        beta_lqp = -1,
        beta_hinge = -1,
        betamultiplier = 1,
        defaultprevalence = 0.5
      ))
      # Building Preojections ###
      
      # Current #
      maxent_Cur <-
        quiet(raster::predict(bio.crop, maxent_model))
      
      
      # Furure 45 #
      maxent_Fut.45 <- quiet(preFut(rast = scenario.list.45,
                                    model =  maxent_model, GCM = GCM45s))
      
      # Furure 85 #
      maxent_Fut.85 = quiet(preFut(rast = scenario.list.85,
                                   model =  maxent_model, GCM = GCM85s))
      
      maxent_std = rescMod(maxent_Cur, maxent_Fut.45, maxent_Fut.85)
      
      names(maxent_std)<-c(paste0("Cur", ".",PA, ".", RUN), paste0("Fut.45", ".",PA, ".", RUN),
                           paste0("Fut.85", ".",PA, ".", RUN))
      
      maxent_Cur <- subset(maxent_std, grep(
        "Cur", names(maxent_std)))
      writeRaster(
        maxent_Cur, 
        paste0("./temp_output/",especie,"/", "maxent_Cur","PA_",PA,"RUN_",RUN,".grd"),
        # formato = "GTiff",
        overwrite = T)
      
      
      maxent_Fut.45 = subset(maxent_std, grep(
        "Fut.45", names(maxent_std)))
      writeRaster(
        maxent_Fut.45,
        paste0("./temp_output/",especie,"/","maxent_Fut.45","PA_",PA,"RUN_",RUN,".grd"),
        # formato = "GTiff",
        overwrite = T)
      
      maxent_Fut.85= subset(maxent_std, grep(
        "Fut.85", names(maxent_std)))
      writeRaster(
        maxent_Fut.85,
        paste0("./temp_output/",especie,"/","maxent_Fut.85","PA_",PA,"RUN_",RUN,".grd"),
        # formato = "GTiff",
        overwrite = T)
      
      # Evaluating ###
      
      maxent_eval <-
        eval.All.Model(subset(maxent_std, grep("Cur", names(maxent_std))), 
                       test.pa)
      
      maxent_th.spec_sens <-
        dismo::threshold(maxent_eval, "spec_sens")
      
      # maxent_th.LPT5 <-
      #   quantile(raster::extract(maxent_Cur, ocor[id.training.pa,]), 0.05, na.rm = T)
      # 
      # maxent_th.VDl <-
      #   com.vdl(maxent_Cur, test.pa, maxent_eval)@maxVDl
      
      
      maxent_eval.spec_sens <-
        eval.All.Model(subset(maxent_std, grep("Cur", names(maxent_std))), test.pa,
                       tr = maxent_th.spec_sens)
      # maxent_eval.LPT5 <-
      #   eval.All.Model(maxent_Cur, test.pa,
      #                  tr = maxent_th.LPT5)
      # maxent_eval.VDl <-
      #   eval.All.Model(maxent_Cur, test.pa, tr = maxent_th.VDl)
      
      
      maxent_e.spec_sens <- c(
        AUC.SpecSens = maxent_eval.spec_sens@auc,
        TPR.SpecSens = maxent_eval.spec_sens@TPR,
        TNR.SpecSens = maxent_eval.spec_sens@TNR,
        thr.SpecSens = maxent_th.spec_sens,
        # d.SpecSens = maxent_eval.spec_sens@TPR * (1 - colSums(
        #   as.data.frame(
        #     maxent_Cur >= maxent_th.spec_sens,
        #     xy = F,
        #     na.rm = T
        #   )
        # ) / areaToral),
        TSS.SpecSens = (maxent_eval.spec_sens@TPR + maxent_eval.spec_sens@TNR)-1
        # TSSSpecSens = (
        #   maxent_eval.spec_sens@TPR + maxent_eval.spec_sens@TNR
        # ) - 1
      )
      
      # maxent_e.LPT5 <- c(
      #   AUC.LPT5 = maxent_eval.LPT5@auc,
      #   TPR.LPT5 = maxent_eval.LPT5@TPR,
      #   TNR.LPT5 = maxent_eval.LPT5@TNR,
      #   thr.LPT5 = maxent_th.LPT5,
      #   d.LPT5 = maxent_eval.LPT5@TPR * (1 - colSums(
      #     as.data.frame(
      #       maxent_Cur >= maxent_th.LPT5,
      #       xy = F,
      #       na.rm = T
      #     )
      #   ) / areaToral),
      #   TSS.LPT5 = maxent_eval.LPT5@TSS
      #   # TSSLPT5 = (maxent_eval.LPT5@TPR + maxent_eval.LPT5@TNR) - 1
      # )
      # 
      # 
      # maxent_e.VDl <- c(
      #   AUC.VDl = maxent_eval.VDl@auc,
      #   TPR.VDl = maxent_eval.VDl@TPR,
      #   TNR.VDl = maxent_eval.VDl@TNR,
      #   thr.VDl = maxent_th.VDl,
      #   d.VDl = maxent_eval.VDl@TPR * (1 - colSums(
      #     as.data.frame(
      #       maxent_Cur >= maxent_th.VDl,
      #       xy = F,
      #       na.rm = T
      #     )
      #   ) / areaToral),
      #   TSS.VDl = maxent_eval.VDl@TSS
      #   # TSSVDl = (maxent_eval.VDl@TPR + maxent_eval.VDl@TNR) - 1
      # )
      maxent.e = rbind(maxent.e, c(maxent_e.spec_sens, PA = PA, RUN = RUN))#, maxent_e.LPT5, maxent_e.VDl))
      
      rownames(maxent.e) = rep(paste0("maxent"), 
                               nrow(maxent.e))
      
      write.csv(maxent.e,
                paste0("./temp_output/",especie,"/", "maxent_eval.all.csv"), row.names = T)
      
      if(VarImport == T){
        
        maxent_var.part.o  = foreach::foreach(
          vari = 1:length(vars),
          .combine = rbind,.packages = c("dismo", "foreach", "adehabitatHS"),
          .inorder = T) %dopar% {
            rmo = as.integer(subset(var, varia == vars[vari], select = pos)[[1]])
            
            maxent_mod <- quiet(maxent(
              x = training.b[,-c(1:3)][, -rmo],
              p = training.b[, 'pb'],
              memory_allocated = 512,
              background_data_dir = 'default',
              maximumbackground = 'default',
              maximumiterations = 200,
              visible = FALSE,
              linear = TRUE,
              quadratic = TRUE,
              product = TRUE,
              threshold = TRUE,
              hinge = TRUE,
              lq2lqptthreshold = 80,
              l2lqthreshold = 10,
              hingethreshold = 15,
              beta_threshold = -1,
              beta_categorical = -1,
              beta_lqp = -1,
              beta_hinge = -1,
              betamultiplier = 1,
              defaultprevalence = 0.5
            ))
            maxent_PCA <- raster::stack(bio.crop[[-rmo]])
            
            maxent_VarImp1 <-quiet(rescMod.One(raster::predict(maxent_PCA, 
                                                                        maxent_mod)))
            
            
            maxent_obs <-
              as.data.frame(maxent_Cur, xy = T)[,3]
            
            maxent_sim <- 
              as.data.frame(maxent_VarImp1, xy = T)[,3]
            
            er = na.omit(data.frame(sim = maxent_sim, obs = maxent_obs))
            
            maxent_RMSE.VarImp1 <-
              rbind(MLmetrics::RMSE(er$sim, 
                                    er$obs))
            
            
            maxent_NicheOverlap.VarImp1 <-
              rbind(1 - dismo::nicheOverlap(maxent_Cur, 
                                            maxent_VarImp1))
            
            cbind(maxent_RMSE.VarImp1,maxent_NicheOverlap.VarImp1
            )
            
          }
        maxent_RMSE.name<-maxent_NicheOverlap.name<-NULL
        
        for (vari in 1:length(vars)) {
          maxent_RMSE.name  <-c(maxent_RMSE.name,
                                paste0(vars[vari],"_", "RMSE_", vari))
          
          maxent_NicheOverlap.name  <-c(maxent_NicheOverlap.name,
                                        paste0(vars[vari],"_", "Overlap_", vari))
        }
        namesRow = c(maxent_RMSE.name, maxent_NicheOverlap.name)
        
        
        maxent_var.part.o <-
          data.frame(metrics = namesRow,
                     imp = c(maxent_var.part.o[,1], maxent_var.part.o[,2]), 
                     PA = PA, RUN = RUN, alg = "maxent")
        
        maxent.var.All<-rbind(maxent.var.All, 
                              maxent_var.part.o)
      } else { next }
      
      remove(list = ls()[c(grep(
        "maxent_", as.factor(ls())))])
      gc(reset = TRUE, full = T)
      
      
      # ##### SVM ####
      print(paste0(especie, " ","SVM"," ",  "PA",  PA," ","RUN", RUN))
      
      svm_model <- ksvm(pb ~ ., data = training.b[,-c(1:2)])
      
      # Building Preojections ###
      
      # Current #
      svm_Cur <-
        raster::predict(bio.crop, svm_model)
      
      
      # Furure 45 #
      svm_Fut.45 <- preFut(rast = scenario.list.45,
                           model =  svm_model, GCM = GCM45s)
      
      # Furure 85 #
      svm_Fut.85 = preFut(rast = scenario.list.85,
                          model =  svm_model, GCM = GCM85s)
      
      svm_std = rescMod(svm_Cur, svm_Fut.45, svm_Fut.85)
      
      names(svm_std)<-c(paste0("Cur", ".",PA, ".", RUN), paste0("Fut.45", ".",PA, ".", RUN),
                        paste0("Fut.85", ".",PA, ".", RUN))
      
      svm_Cur <- subset(svm_std, grep(
        "Cur", names(svm_std)))
      writeRaster(
        svm_Cur, 
        paste0("./temp_output/",especie,"/", "svm_Cur","PA_",PA,"RUN_",RUN,".grd"),
        # formato = "GTiff",
        overwrite = T)
      
      
      svm_Fut.45 = subset(svm_std, grep(
        "Fut.45", names(svm_std)))
      writeRaster(
        svm_Fut.45,
        paste0("./temp_output/",especie,"/","svm_Fut.45","PA_",PA,"RUN_",RUN,".grd"),
        # formato = "GTiff",
        overwrite = T)
      
      svm_Fut.85= subset(svm_std, grep(
        "Fut.85", names(svm_std)))
      writeRaster(
        svm_Fut.85,
        paste0("./temp_output/",especie,"/","svm_Fut.85","PA_",PA,"RUN_",RUN,".grd"),
        # formato = "GTiff",
        overwrite = T)
      
      # Evaluating ###
      
      svm_eval <-
        eval.All.Model(subset(svm_std, grep("Cur", names(svm_std))), 
                       test.pa)
      
      svm_th.spec_sens <-
        dismo::threshold(svm_eval, "spec_sens")
      
      # svm_th.LPT5 <-
      #   quantile(raster::extract(svm_Cur, ocor[id.training.b,]), 0.05, na.rm = T)
      # 
      # svm_th.VDl <-
      #   com.vdl(svm_Cur, test.pa, svm_eval)@maxVDl
      
      
      svm_eval.spec_sens <-
        eval.All.Model(subset(svm_std, grep("Cur", names(svm_std))), test.pa,
                       tr = svm_th.spec_sens)
      # svm_eval.LPT5 <-
      #   eval.All.Model(svm_Cur, test.pa,
      #                  tr = svm_th.LPT5)
      # svm_eval.VDl <-
      #   eval.All.Model(svm_Cur, test.pa, tr = svm_th.VDl)
      
      
      svm_e.spec_sens <- c(
        AUC.SpecSens = svm_eval.spec_sens@auc,
        TPR.SpecSens = svm_eval.spec_sens@TPR,
        TNR.SpecSens = svm_eval.spec_sens@TNR,
        thr.SpecSens = svm_th.spec_sens,
        # d.SpecSens = svm_eval.spec_sens@TPR * (1 - colSums(
        #   as.data.frame(
        #     svm_Cur >= svm_th.spec_sens,
        #     xy = F,
        #     na.rm = T
        #   )
        # ) / areaToral),
        TSS.SpecSens = (svm_eval.spec_sens@TPR + svm_eval.spec_sens@TNR)-1
        # TSSSpecSens = (
        #   svm_eval.spec_sens@TPR + svm_eval.spec_sens@TNR
        # ) - 1
      )
      
      # svm_e.LPT5 <- c(
      #   AUC.LPT5 = svm_eval.LPT5@auc,
      #   TPR.LPT5 = svm_eval.LPT5@TPR,
      #   TNR.LPT5 = svm_eval.LPT5@TNR,
      #   thr.LPT5 = svm_th.LPT5,
      #   d.LPT5 = svm_eval.LPT5@TPR * (1 - colSums(
      #     as.data.frame(
      #       svm_Cur >= svm_th.LPT5,
      #       xy = F,
      #       na.rm = T
      #     )
      #   ) / areaToral),
      #   TSS.LPT5 = svm_eval.LPT5@TSS
      #   # TSSLPT5 = (svm_eval.LPT5@TPR + svm_eval.LPT5@TNR) - 1
      # )
      # 
      # 
      # svm_e.VDl <- c(
      #   AUC.VDl = svm_eval.VDl@auc,
      #   TPR.VDl = svm_eval.VDl@TPR,
      #   TNR.VDl = svm_eval.VDl@TNR,
      #   thr.VDl = svm_th.VDl,
      #   d.VDl = svm_eval.VDl@TPR * (1 - colSums(
      #     as.data.frame(
      #       svm_Cur >= svm_th.VDl,
      #       xy = F,
      #       na.rm = T
      #     )
      #   ) / areaToral),
      #   TSS.VDl = svm_eval.VDl@TSS
      #   # TSSVDl = (svm_eval.VDl@TPR + svm_eval.VDl@TNR) - 1
      # )
      svm.e = rbind(svm.e, c(svm_e.spec_sens, PA = PA, RUN = RUN))#, svm_e.LPT5, svm_e.VDl))
      
      rownames(svm.e) = rep(paste0("svm"), 
                            nrow(svm.e))
      
      write.csv(svm.e,
                paste0("./temp_output/",especie,"/", "svm_eval.all.csv"), row.names = T)
      
      if(VarImport == T){
        
        svm_var.part.o  = foreach::foreach(
          vari = 1:length(vars),
          .combine = rbind,.packages = c("dismo", "foreach", "adehabitatHS"),
          .inorder = T) %dopar% {
            rmo = as.integer(subset(var, varia == vars[vari], select = pos)[[1]])
            
            rp <- 1 + rmo
            
            svmTraining = training.b[,-c(1:2)]
            
            svm_mod <-  kernlab::ksvm( pb ~ ., data = svmTraining[, -rp])
            
            svm_PCA <- raster::stack(bio.crop[[-rmo]])
            
            svm_VarImp1 <-(rescMod.One(raster::predict(svm_PCA,
                                                                svm_mod)))
            
            
            svm_obs <-
              as.data.frame(svm_Cur, xy = T)[,3]
            
            svm_sim <- 
              as.data.frame(svm_VarImp1, xy = T)[,3]
            
            er = na.omit(data.frame(sim = svm_sim, obs = svm_obs))
            
            svm_RMSE.VarImp1 <-
              rbind(MLmetrics::RMSE(er$sim, 
                                    er$obs))
            
            
            svm_NicheOverlap.VarImp1 <-
              rbind(1 - dismo::nicheOverlap(svm_Cur, 
                                            svm_VarImp1))
            
            cbind(svm_RMSE.VarImp1,svm_NicheOverlap.VarImp1
            )
            
          }
        svm_RMSE.name<-svm_NicheOverlap.name<-NULL
        
        for (vari in 1:length(vars)) {
          svm_RMSE.name  <-c(svm_RMSE.name,
                             paste0(vars[vari],"_", "RMSE_", vari))
          
          svm_NicheOverlap.name  <-c(svm_NicheOverlap.name,
                                     paste0(vars[vari],"_", "Overlap_", vari))
        }
        namesRow = c(svm_RMSE.name, svm_NicheOverlap.name)
        
        
        svm_var.part.o <-
          data.frame(metrics = namesRow,
                     imp = c(svm_var.part.o[,1], svm_var.part.o[,2]), 
                     PA = PA, RUN = RUN, alg = "svm")
        
        svm.var.All<-rbind(svm.var.All, 
                           svm_var.part.o)
      } else { next }
      
      remove(list = ls()[c(grep(
        "svm_", as.factor(ls())))])
      gc(reset = TRUE, full = T)
      
      
      ##### ANN ####
      print(paste0(especie, " ","ANN"," ", "PA",  PA," ","RUN", RUN))
      
      nnet_model <- nnet(
        pb ~ . ,
        data = training.pa[,-c(1:2)],
        NbCV = 5,
        size = 1,
        decay = 0.05,
        rang = 0.1,
        maxit = 200, 
        trace = F
      )
      
      # Building Preojections ###
      
      # Current #
      nnet_Cur <-
        raster::predict(bio.crop, nnet_model)
      
      
      # Furure 45 #
      nnet_Fut.45 <- preFut(rast = scenario.list.45,
                            model =  nnet_model, GCM = GCM45s)
      
      # Furure 85 #
      nnet_Fut.85 = preFut(rast = scenario.list.85,
                           model =  nnet_model, GCM = GCM85s)
      
      nnet_std = rescMod(nnet_Cur, nnet_Fut.45, nnet_Fut.85)
      
      names(nnet_std)<-c(paste0("Cur", ".",PA, ".", RUN), paste0("Fut.45", ".",PA, ".", RUN),
                         paste0("Fut.85", ".",PA, ".", RUN))
      
      nnet_Cur <- subset(nnet_std, grep(
        "Cur", names(nnet_std)))
      writeRaster(
        nnet_Cur, 
        paste0("./temp_output/",especie,"/", "nnet_Cur","PA_",PA,"RUN_",RUN,".grd"),
        # formato = "GTiff",
        overwrite = T)
      
      
      nnet_Fut.45 = subset(nnet_std, grep(
        "Fut.45", names(nnet_std)))
      writeRaster(
        nnet_Fut.45,
        paste0("./temp_output/",especie,"/","nnet_Fut.45","PA_",PA,"RUN_",RUN,".grd"),
        # formato = "GTiff",
        overwrite = T)
      
      nnet_Fut.85= subset(nnet_std, grep(
        "Fut.85", names(nnet_std)))
      writeRaster(
        nnet_Fut.85,
        paste0("./temp_output/",especie,"/","nnet_Fut.85","PA_",PA,"RUN_",RUN,".grd"),
        # formato = "GTiff",
        overwrite = T)
      
      # Evaluating ###
      
      nnet_eval <-
        eval.All.Model(subset(nnet_std, grep("Cur", names(nnet_std))), 
                       test.pa)
      
      nnet_th.spec_sens <-
        dismo::threshold(nnet_eval, "spec_sens")
      
      # nnet_th.LPT5 <-
      #   quantile(raster::extract(nnet_Cur, ocor[id.training.pa,]), 0.05, na.rm = T)
      # 
      # nnet_th.VDl <-
      #   com.vdl(nnet_Cur, test.pa, nnet_eval)@maxVDl
      
      
      nnet_eval.spec_sens <-
        eval.All.Model(subset(nnet_std, grep("Cur", names(nnet_std))), test.pa,
                       tr = nnet_th.spec_sens)
      # nnet_eval.LPT5 <-
      #   eval.All.Model(nnet_Cur, test.pa,
      #                  tr = nnet_th.LPT5)
      # nnet_eval.VDl <-
      #   eval.All.Model(nnet_Cur, test.pa, tr = nnet_th.VDl)
      
      
      nnet_e.spec_sens <- c(
        AUC.SpecSens = nnet_eval.spec_sens@auc,
        TPR.SpecSens = nnet_eval.spec_sens@TPR,
        TNR.SpecSens = nnet_eval.spec_sens@TNR,
        thr.SpecSens = nnet_th.spec_sens,
        # d.SpecSens = nnet_eval.spec_sens@TPR * (1 - colSums(
        #   as.data.frame(
        #     nnet_Cur >= nnet_th.spec_sens,
        #     xy = F,
        #     na.rm = T
        #   )
        # ) / areaToral),
        TSS.SpecSens = (nnet_eval.spec_sens@TPR + nnet_eval.spec_sens@TNR)-1
        # TSSSpecSens = (
        #   nnet_eval.spec_sens@TPR + nnet_eval.spec_sens@TNR
        # ) - 1
      )
      
      # nnet_e.LPT5 <- c(
      #   AUC.LPT5 = nnet_eval.LPT5@auc,
      #   TPR.LPT5 = nnet_eval.LPT5@TPR,
      #   TNR.LPT5 = nnet_eval.LPT5@TNR,
      #   thr.LPT5 = nnet_th.LPT5,
      #   d.LPT5 = nnet_eval.LPT5@TPR * (1 - colSums(
      #     as.data.frame(
      #       nnet_Cur >= nnet_th.LPT5,
      #       xy = F,
      #       na.rm = T
      #     )
      #   ) / areaToral),
      #   TSS.LPT5 = nnet_eval.LPT5@TSS
      #   # TSSLPT5 = (nnet_eval.LPT5@TPR + nnet_eval.LPT5@TNR) - 1
      # )
      # 
      # 
      # nnet_e.VDl <- c(
      #   AUC.VDl = nnet_eval.VDl@auc,
      #   TPR.VDl = nnet_eval.VDl@TPR,
      #   TNR.VDl = nnet_eval.VDl@TNR,
      #   thr.VDl = nnet_th.VDl,
      #   d.VDl = nnet_eval.VDl@TPR * (1 - colSums(
      #     as.data.frame(
      #       nnet_Cur >= nnet_th.VDl,
      #       xy = F,
      #       na.rm = T
      #     )
      #   ) / areaToral),
      #   TSS.VDl = nnet_eval.VDl@TSS
      #   # TSSVDl = (nnet_eval.VDl@TPR + nnet_eval.VDl@TNR) - 1
      # )
      nnet.e = rbind(nnet.e, c(nnet_e.spec_sens, PA = PA, RUN = RUN))#, nnet_e.LPT5, nnet_e.VDl))
      
      rownames(nnet.e) = rep(paste0("nnet"), 
                             nrow(nnet.e))
      
      write.csv(nnet.e,
                paste0("./temp_output/",especie,"/", "nnet_eval.all.csv"), row.names = T)
      
      if(VarImport == T){
        
        nnet_var.part.o  = foreach::foreach(
          vari = 1:length(vars),
          .combine = rbind,.packages = c("dismo", "foreach", "adehabitatHS"),
          .inorder = T) %dopar% {
            rmo = as.integer(subset(var, varia == vars[vari], select = pos)[[1]])
            
            nnet_mod <-
              nnet::nnet(
                pb ~ . ,
                data = training.pa[, -c(1:3)][, -rmo],
                NbCV = 5,
                size = 1,
                decay = 0.05,
                rang = 0.1,
                maxit = 200, 
                trace = F
              )
            
            nnet_PCA <- raster::stack(bio.crop[[-rmo]])
            
            nnet_VarImp1 <-(rescMod.One(raster::predict(nnet_PCA, 
                                                                 nnet_mod)))
            
            
            nnet_obs <-
              as.data.frame(nnet_Cur, xy = T)[,3]
            
            nnet_sim <- 
              as.data.frame(nnet_VarImp1, xy = T)[,3]
            
            er = na.omit(data.frame(sim = nnet_sim, obs = nnet_obs))
            
            nnet_RMSE.VarImp1 <-
              rbind(MLmetrics::RMSE(er$sim, 
                                    er$obs))
            
            
            nnet_NicheOverlap.VarImp1 <-
              rbind(1 - dismo::nicheOverlap(nnet_Cur, 
                                            nnet_VarImp1))
            
            cbind(nnet_RMSE.VarImp1,nnet_NicheOverlap.VarImp1
            )
            
          }
        nnet_RMSE.name<-nnet_NicheOverlap.name<-NULL
        
        for (vari in 1:length(vars)) {
          nnet_RMSE.name  <-c(nnet_RMSE.name,
                              paste0(vars[vari],"_", "RMSE_", vari))
          
          nnet_NicheOverlap.name  <-c(nnet_NicheOverlap.name,
                                      paste0(vars[vari],"_", "Overlap_", vari))
        }
        namesRow = c(nnet_RMSE.name, nnet_NicheOverlap.name)
        
        
        nnet_var.part.o <-
          data.frame(metrics = namesRow,
                     imp = c(nnet_var.part.o[,1], nnet_var.part.o[,2]), 
                     PA = PA, RUN = RUN, alg = "nnet")
        
        nnet.var.All<-rbind(nnet.var.All, 
                            nnet_var.part.o)
      } else { next }
      
      remove(list = ls()[c(grep(
        "nnet_", as.factor(ls())))])
      gc(reset = TRUE, full = T)
      
      
      # ##### Random Forest ####
      print(paste0(especie, " ","Random Forest"," ", "PA",  PA," ","RUN", RUN))
      
      suppressWarnings(rf_model <-
                         randomForest::randomForest(
                           pb ~ .,
                           data = training.pa[,-c(1:2)],
                           ntree = 500,
                           nodesize = 5, importance = T
                         ))
      # Building Preojections ###
      
      # Current #
      rf_Cur <-
        raster::predict(bio.crop, rf_model)
      
      
      # Furure 45 #
      rf_Fut.45 <- preFut(rast = scenario.list.45,
                          model =  rf_model, GCM = GCM45s)
      
      # Furure 85 #
      rf_Fut.85 = preFut(rast = scenario.list.85,
                         model =  rf_model, GCM = GCM85s)
      
      rf_std = rescMod(rf_Cur, rf_Fut.45, rf_Fut.85)
      
      names(rf_std)<-c(paste0("Cur", ".",PA, ".", RUN), paste0("Fut.45", ".",PA, ".", RUN),
                       paste0("Fut.85", ".",PA, ".", RUN))
      
      rf_Cur <- subset(rf_std, grep(
        "Cur", names(rf_std)))
      writeRaster(
        rf_Cur, 
        paste0("./temp_output/",especie,"/", "rf_Cur","PA_",PA,"RUN_",RUN,".grd"),
        # formato = "GTiff",
        overwrite = T)
      
      
      rf_Fut.45 = subset(rf_std, grep(
        "Fut.45", names(rf_std)))
      writeRaster(
        rf_Fut.45,
        paste0("./temp_output/",especie,"/","rf_Fut.45","PA_",PA,"RUN_",RUN,".grd"),
        # formato = "GTiff",
        overwrite = T)
      
      rf_Fut.85= subset(rf_std, grep(
        "Fut.85", names(rf_std)))
      writeRaster(
        rf_Fut.85,
        paste0("./temp_output/",especie,"/","rf_Fut.85","PA_",PA,"RUN_",RUN,".grd"),
        # formato = "GTiff",
        overwrite = T)
      
      # Evaluating ###
      
      rf_eval <-
        eval.All.Model(subset(rf_std, grep("Cur", names(rf_std))), 
                       test.pa)
      
      rf_th.spec_sens <-
        dismo::threshold(rf_eval, "spec_sens")
      
      # rf_th.LPT5 <-
      #   quantile(raster::extract(rf_Cur, ocor[id.training.pa,]), 0.05, na.rm = T)
      # 
      # rf_th.VDl <-
      #   com.vdl(rf_Cur, test.pa, rf_eval)@maxVDl
      
      
      rf_eval.spec_sens <-
        eval.All.Model(subset(rf_std, grep("Cur", names(rf_std))), test.pa,
                       tr = rf_th.spec_sens)
      # rf_eval.LPT5 <-
      #   eval.All.Model(rf_Cur, test.pa,
      #                  tr = rf_th.LPT5)
      # rf_eval.VDl <-
      #   eval.All.Model(rf_Cur, test.pa, tr = rf_th.VDl)
      
      
      rf_e.spec_sens <- c(
        AUC.SpecSens = rf_eval.spec_sens@auc,
        TPR.SpecSens = rf_eval.spec_sens@TPR,
        TNR.SpecSens = rf_eval.spec_sens@TNR,
        thr.SpecSens = rf_th.spec_sens,
        # d.SpecSens = rf_eval.spec_sens@TPR * (1 - colSums(
        #   as.data.frame(
        #     rf_Cur >= rf_th.spec_sens,
        #     xy = F,
        #     na.rm = T
        #   )
        # ) / areaToral),
        TSS.SpecSens = (rf_eval.spec_sens@TPR + rf_eval.spec_sens@TNR)-1
        # TSSSpecSens = (
        #   rf_eval.spec_sens@TPR + rf_eval.spec_sens@TNR
        # ) - 1
      )
      
      # rf_e.LPT5 <- c(
      #   AUC.LPT5 = rf_eval.LPT5@auc,
      #   TPR.LPT5 = rf_eval.LPT5@TPR,
      #   TNR.LPT5 = rf_eval.LPT5@TNR,
      #   thr.LPT5 = rf_th.LPT5,
      #   d.LPT5 = rf_eval.LPT5@TPR * (1 - colSums(
      #     as.data.frame(
      #       rf_Cur >= rf_th.LPT5,
      #       xy = F,
      #       na.rm = T
      #     )
      #   ) / areaToral),
      #   TSS.LPT5 = rf_eval.LPT5@TSS
      #   # TSSLPT5 = (rf_eval.LPT5@TPR + rf_eval.LPT5@TNR) - 1
      # )
      # 
      # 
      # rf_e.VDl <- c(
      #   AUC.VDl = rf_eval.VDl@auc,
      #   TPR.VDl = rf_eval.VDl@TPR,
      #   TNR.VDl = rf_eval.VDl@TNR,
      #   thr.VDl = rf_th.VDl,
      #   d.VDl = rf_eval.VDl@TPR * (1 - colSums(
      #     as.data.frame(
      #       rf_Cur >= rf_th.VDl,
      #       xy = F,
      #       na.rm = T
      #     )
      #   ) / areaToral),
      #   TSS.VDl = rf_eval.VDl@TSS
      #   # TSSVDl = (rf_eval.VDl@TPR + rf_eval.VDl@TNR) - 1
      # )
      rf.e = rbind(rf.e, c(rf_e.spec_sens, PA = PA, RUN = RUN))#, rf_e.LPT5, rf_e.VDl))
      
      rownames(rf.e) = rep(paste0("rf"), 
                           nrow(rf.e))
      
      write.csv(rf.e,
                paste0("./temp_output/",especie,"/", "rf_eval.all.csv"), row.names = T)
      
      if(VarImport == T){
        
        rf_var.part.o  = foreach::foreach(
          vari = 1:length(vars),
          .combine = rbind,.packages = c("dismo", "foreach", "adehabitatHS"),
          .inorder = T) %dopar% {
            rmo = as.integer(subset(var, varia == vars[vari], select = pos)[[1]])
            
            suppressWarnings(rf_mod <-
                               randomForest::randomForest(
                                 pb ~ .,
                                 data =training.pa[,-c(1:3)][,-rmo],
                                 ntree = 500,
                                 nodesize = 5, importance = T
                               ))
            
            rf_PCA <- raster::stack(bio.crop[[-rmo]])
            
            rf_VarImp1 <-(rescMod.One(raster::predict(rf_PCA, 
                                                               rf_mod)))
            
            
            rf_obs <-
              as.data.frame(rf_Cur, xy = T)[,3]
            
            rf_sim <- 
              as.data.frame(rf_VarImp1, xy = T)[,3]
            
            er = na.omit(data.frame(sim = rf_sim, obs = rf_obs))
            
            rf_RMSE.VarImp1 <-
              rbind(MLmetrics::RMSE(er$sim, 
                                    er$obs))
            
            
            rf_NicheOverlap.VarImp1 <-
              rbind(1 - dismo::nicheOverlap(rf_Cur, 
                                            rf_VarImp1))
            
            cbind(rf_RMSE.VarImp1,rf_NicheOverlap.VarImp1
            )
            
          }
        rf_RMSE.name<-rf_NicheOverlap.name<-NULL
        
        for (vari in 1:length(vars)) {
          rf_RMSE.name  <-c(rf_RMSE.name,
                            paste0(vars[vari],"_", "RMSE_", vari))
          
          rf_NicheOverlap.name  <-c(rf_NicheOverlap.name,
                                    paste0(vars[vari],"_", "Overlap_", vari))
        }
        namesRow = c(rf_RMSE.name, rf_NicheOverlap.name)
        
        
        rf_var.part.o <-
          data.frame(metrics = namesRow,
                     imp = c(rf_var.part.o[,1], rf_var.part.o[,2]), 
                     PA = PA, RUN = RUN, alg = "rf")
        
        rf.var.All<-rbind(rf.var.All, 
                          rf_var.part.o)
      } else { next }
      
      remove(list = ls()[c(grep(
        "rf_", as.factor(ls())))])
      gc(reset = TRUE, full = T)
      
    }##### End RUN loop ###
    gc(reset = TRUE, full = TRUE)
    # print(paste0(especie," ","PA set"," ", PA))
  } ##### End PAs loop ###
  gc(reset = TRUE, full = TRUE)
  
  
  #### Writing All Reults ####
  
  ##### save VarImpor ####
  if(VarImport == T){
    
  varImport = rbind(
    bioclim.var.All,
    enfa.var.All,
    glm.var.All,
    mars.var.All,
    maxent.var.All,
    svm.var.All,
    nnet.var.All,
    rf.var.All
  )
  
  write.csv(
    varImport,
    paste0("./outputs/",especie,"_","varImp.All.csv"),
    row.names = F)
  
  rm(varImport,
     bioclim.var.All,
     enfa.var.All,
     glm.var.All,
     mars.var.All,
     maxent.var.All,
     svm.var.All,
     nnet.var.All,
     rf.var.All)
  
  } else { next }
  

  
  # Writing Evaluation Resultes #
  {
  # bioclim #
  bioclim_e = read.table(paste0("./temp_output/",especie,"/", "bioclim_eval.all.csv"), 
             header = T, sep = ',')
  
  # domain #
  # domain.e = read.table(paste0("./temp_output/",especie,"/", "domain_eval.all.csv"), 
  #                       header = T, sep = ',')
  
  # enfa #
  enfa_e = read.table(paste0("./temp_output/",especie,"/", "enfa_eval.all.csv"), 
                      header = T, sep = ',')
  
  # glm #
  glm_e = read.table(paste0("./temp_output/",especie,"/", "glm_eval.all.csv"), 
                     header = T, sep = ',')
  
  # mars #
  mars_e = read.table(paste0("./temp_output/",especie,"/", "mars_eval.all.csv"), 
                      header = T, sep = ',')
  
  # maxent #
  maxent_e = read.table(paste0("./temp_output/",especie,"/", "maxent_eval.all.csv"), 
                        header = T, sep = ',')
  
  # svm #
  svm_e = read.table(paste0("./temp_output/",especie,"/", "svm_eval.all.csv"), 
                     header = T, sep = ',')
  
  # nnet #
  nnet_e = read.table(paste0("./temp_output/",especie,"/", "nnet_eval.all.csv"), 
                      header = T, sep = ',')
  
  # rf #
  rf_e = read.table(paste0("./temp_output/",especie,"/", "rf_eval.all.csv"), 
                    header = T, sep = ',')

  Evaluation.all <- rbind(bioclim_e, enfa_e,  glm_e, mars_e, maxent_e, svm_e, nnet_e, rf_e)

  rm(bioclim_e, enfa_e,  glm_e, mars_e, maxent_e, svm_e, nnet_e, rf_e)  
  
  Evaluation.all <- cbind(Evaluation.all, ID = 1:nrow(Evaluation.all))

  colnames(Evaluation.all)<-c("Model","AUC","TPR","TNR","threshold","TSS","PA","RUN", "ID")
  
  ## Write Evaluation ##
  write.csv(Evaluation.all,
            paste0("./outputs/", especie, "_", "Evaluation.all.csv"), row.names = T)
  
  ## Write Thresholds ##
  th.all <- cbind(th = Evaluation.all[,"threshold"])
  write.csv(th.all,
            paste0("./outputs/", especie, "_", "th.all.csv"), row.names = T)
  
  th.mean <- mean(th.all)
  write.csv(th.mean,
            paste0("./outputs/", especie, "_", "th.mean.csv"), row.names = T)
  
  # Selectning "good models" #
  
  sel2 = Evaluation.all[Evaluation.all[, "TSS"] > 0.400, ]
  sel2 <- na.omit(sel2)
  write.csv(sel2,
            paste0("./outputs/", especie, "_", "Selected.Models.csv"), row.names = T)
  rm(Evaluation.all)
  gc(reset = TRUE, full = TRUE)
  }
  ## Write Rasters ##
  {
  # Current ####
  # bioclim ###
  bioclim.cur = stack()
  for(PA in 1:PAs){
      for (RUN in 1:RUNs) {
        bioclim.cur = stack(bioclim.cur, stack(list.files(paste0("./temp_output/",especie,"/"), 
                                                    paste0('bioclim_Cur',"PA_",PA,"RUN_",RUN,".grd"),
                                       full.names = T)))}}
  names(bioclim.cur) <- paste0("bioclim.cur.", 1:nlayers(bioclim.cur))

  # # domain ###
  # domain.cur = stack()
  # for(PA in 1:PAs){
  #   for (RUN in 1:RUNs) {
  #     domain.cur = stack(domain.cur, stack(list.files(paste0("./temp_output/",especie,"/"), 
  #                                                     paste0('domain_Cur',"PA_",PA,"RUN_",RUN,".grd"),
  #                                                     full.names = T)))}}
  # names(domain.cur) <- paste0("domain.cur.", 1:nlayers(domain.cur))
 
  # enfa ###
  enfa.cur = stack()
  for(PA in 1:PAs){
    for (RUN in 1:RUNs) {
      enfa.cur = stack(enfa.cur, stack(list.files(paste0("./temp_output/",especie,"/"), 
                                                  paste0('enfa_Cur',"PA_",PA,"RUN_",RUN,".grd"),
                                                  full.names = T)))}}
  names(enfa.cur) <- paste0("enfa.cur.", 1:nlayers(enfa.cur))
  
  # glm ###
  glm.cur = stack()
  for(PA in 1:PAs){
    for (RUN in 1:RUNs) {
      glm.cur = stack(glm.cur, stack(list.files(paste0("./temp_output/",especie,"/"), 
                                                paste0('glm_Cur',"PA_",PA,"RUN_",RUN,".grd"),
                                                full.names = T)))}}
  names(glm.cur) <- paste0("glm.cur.", 1:nlayers(glm.cur))
 
  # mars ###
  mars.cur = stack()
  for(PA in 1:PAs){
    for (RUN in 1:RUNs) {
      mars.cur = stack(mars.cur, stack(list.files(paste0("./temp_output/",especie,"/"), 
                                                  paste0('mars_Cur',"PA_",PA,"RUN_",RUN,".grd"),
                                                  full.names = T)))}}
  names(mars.cur) <- paste0("mars.cur.", 1:nlayers(mars.cur))
  
  # maxent ###
  maxent.cur = stack()
  for(PA in 1:PAs){
    for (RUN in 1:RUNs) {
      maxent.cur = stack(maxent.cur, stack(list.files(paste0("./temp_output/",especie,"/"), 
                                                      paste0('maxent_Cur',"PA_",PA,"RUN_",RUN,".grd"),
                                                      full.names = T)))}}
  names(maxent.cur) <- paste0("maxent.cur.", 1:nlayers(maxent.cur))
 
  # svm ###
  svm.cur = stack()
  for(PA in 1:PAs){
    for (RUN in 1:RUNs) {
      svm.cur = stack(svm.cur, stack(list.files(paste0("./temp_output/",especie,"/"), 
                                                paste0('svm_Cur',"PA_",PA,"RUN_",RUN,".grd"),
                                                full.names = T)))}}
  names(svm.cur) <- paste0("svm.cur.", 1:nlayers(svm.cur))

  # nnet ###
  nnet.cur = stack()
  for(PA in 1:PAs){
    for (RUN in 1:RUNs) {
      nnet.cur = stack(nnet.cur, stack(list.files(paste0("./temp_output/",especie,"/"), 
                                                  paste0('nnet_Cur',"PA_",PA,"RUN_",RUN,".grd"),
                                                  full.names = T)))}}
  names(nnet.cur) <- paste0("nnet.cur.", 1:nlayers(nnet.cur))
  
  # rf ###
  rf.cur = stack()
  for(PA in 1:PAs){
    for (RUN in 1:RUNs) {
      rf.cur = stack(rf.cur, stack(list.files(paste0("./temp_output/",especie,"/"), 
                                              paste0('rf_Cur',"PA_",PA,"RUN_",RUN,".grd"),
                                              full.names = T)))}}
  names(rf.cur) <- paste0("rf.cur.", 1:nlayers(rf.cur))
 
  # Current Ensemble ####
  gc(reset = T, full = T)
  Current.mean <- mean(subset(stack(bioclim.cur, enfa.cur, glm.cur, mars.cur, 
                                    maxent.cur, svm.cur, nnet.cur, rf.cur),
                              sel2[, "ID"]))
  writeRaster(
    Current.mean,
    paste0("./outputs/", especie, "_", "Current.mean.tif"), formtat = "GTiff", overwrite = T)
  
  # Binary Transformation #
  Current.bin <- biomod2::BinaryTransformation(Current.mean, th.mean)
  rm(Current.mean)
  gc(reset = T, full = T)
  writeRaster(
    Current.bin,
    paste0("./outputs/", especie, "_", "Current.bin.tif"), formtat = "GTiff", overwrite = T)
  rm(Current.bin,bioclim.cur, enfa.cur, glm.cur, mars.cur, 
     maxent.cur, svm.cur, nnet.cur, rf.cur)
  gc(reset = TRUE, full = TRUE)}
  # Future 45 ####
  {
    # bioclim 45 ###
    
    bioclim.fut.45 = stack()
    for(PA in 1:PAs){
      for (RUN in 1:RUNs) {
        bioclim.fut.45 = stack(bioclim.fut.45, 
                               stack(list.files(paste0("./temp_output/",especie,"/"), 
                                                paste0('bioclim_Fut.45',"PA_",PA,"RUN_",RUN,".grd"),
                                                full.names = T)))}}
    names(bioclim.fut.45) <- paste0("bioclim.fut.45.", 1:nlayers(bioclim.fut.45))
    
    
    # domain 45 ###
    # domain.fut.45 = stack()
    # for(PA in 1:PAs){
    #   for (RUN in 1:RUNs) {
    #     domain.fut.45 = stack(domain.fut.45, 
    #                           stack(list.files(paste0("./temp_output/",especie,"/"), 
    #                                            paste0('domain_Fut.45',"PA_",PA,"RUN_",RUN,".grd"),
    #                                            full.names = T)))}}
    # names(domain.fut.45) <- paste0("domain.fut.45.", 1:nlayers(domain.fut.45))
    
    # enfa 45 ###
    
    enfa.fut.45 = stack()
    for(PA in 1:PAs){
      for (RUN in 1:RUNs) {
        enfa.fut.45 = stack(enfa.fut.45, 
                            stack(list.files(paste0("./temp_output/",especie,"/"), 
                                             paste0('enfa_Fut.45',"PA_",PA,"RUN_",RUN,".grd"),
                                             full.names = T)))}}
    names(enfa.fut.45) <- paste0("enfa.fut.45.", 1:nlayers(enfa.fut.45))
    
    # glm 45 ###
    
    glm.fut.45 = stack()
    for(PA in 1:PAs){
      for (RUN in 1:RUNs) {
        glm.fut.45 = stack(glm.fut.45, 
                           stack(list.files(paste0("./temp_output/",especie,"/"), 
                                            paste0('glm_Fut.45',"PA_",PA,"RUN_",RUN,".grd"),
                                            full.names = T)))}}
    names(glm.fut.45) <- paste0("glm.fut.45.", 1:nlayers(glm.fut.45))
    
    # mars 45 ###
    
    mars.fut.45 = stack()
    for(PA in 1:PAs){
      for (RUN in 1:RUNs) {
        mars.fut.45 = stack(mars.fut.45, 
                            stack(list.files(paste0("./temp_output/",especie,"/"), 
                                             paste0('mars_Fut.45',"PA_",PA,"RUN_",RUN,".grd"),
                                             full.names = T)))}}
    names(mars.fut.45) <- paste0("mars.fut.45.", 1:nlayers(mars.fut.45))
    
    
    # maxent 45 ###
    
    maxent.fut.45 = stack()
    for(PA in 1:PAs){
      for (RUN in 1:RUNs) {
        maxent.fut.45 = stack(maxent.fut.45, 
                              stack(list.files(paste0("./temp_output/",especie,"/"), 
                                               paste0('maxent_Fut.45',"PA_",PA,"RUN_",RUN,".grd"),
                                               full.names = T)))}}
    names(maxent.fut.45) <- paste0("maxent.fut.45.", 1:nlayers(maxent.fut.45))
    
    # svm 45 ###
    
    svm.fut.45 = stack()
    for(PA in 1:PAs){
      for (RUN in 1:RUNs) {
        svm.fut.45 = stack(svm.fut.45, 
                           stack(list.files(paste0("./temp_output/",especie,"/"), 
                                            paste0('svm_Fut.45',"PA_",PA,"RUN_",RUN,".grd"),
                                            full.names = T)))}}
    names(svm.fut.45) <- paste0("svm.fut.45.", 1:nlayers(svm.fut.45))
    
    # nnet 45 ###
    
    nnet.fut.45 = stack()
    for(PA in 1:PAs){
      for (RUN in 1:RUNs) {
        nnet.fut.45 = stack(nnet.fut.45, 
                            stack(list.files(paste0("./temp_output/",especie,"/"), 
                                             paste0('nnet_Fut.45',"PA_",PA,"RUN_",RUN,".grd"),
                                             full.names = T)))}}
    
    names(nnet.fut.45) <- paste0("nnet.fut.45.", 1:nlayers(nnet.fut.45))
    
    # rf 45 ###
    rf.fut.45 = stack()
    for(PA in 1:PAs){
      for (RUN in 1:RUNs) {
        rf.fut.45 = stack(rf.fut.45, 
                          stack(list.files(paste0("./temp_output/",especie,"/"), 
                                           paste0('rf_Fut.45',"PA_",PA,"RUN_",RUN,".grd"),
                                           full.names = T)))}}
    
    names(rf.fut.45) <- paste0("rf.fut.85.", 1:nlayers(rf.fut.45))
    
    
    ## Future Ensemble 45 ####
    
    Future.45.mean <- mean(subset(stack(bioclim.fut.45,enfa.fut.45, glm.fut.45, mars.fut.45,
                                        maxent.fut.45, svm.fut.45, nnet.fut.45, rf.fut.45)
                                  , sel2[, "ID"]))
    rm(bioclim.fut.45,enfa.fut.45, glm.fut.45, mars.fut.45,
       maxent.fut.45, svm.fut.45, nnet.fut.45, rf.fut.45)
    
    writeRaster(
      Future.45.mean,
      paste0("./outputs/", especie, "_", "Future.45.mean.tif"), formtat = "GTiff", overwrite = T)
    
    # Binary Transformation #
    Future.45.bin <- biomod2::BinaryTransformation(Future.45.mean, th.mean)
    rm(Future.45.mean)
    gc(reset = T)
    writeRaster(
      Future.45.bin,
      paste0("./outputs/", especie, "_", "Future.45.bin.tif"), formtat = "GTiff", overwrite = T)
    rm(Future.45.bin)
    gc(reset = TRUE)}
  
  # Future 85 ####
  {
    # bioclim 85 ###
    
    bioclim.fut.85 = stack()
    for(PA in 1:PAs){
      for (RUN in 1:RUNs) {
        bioclim.fut.85 = stack(bioclim.fut.85, 
                               stack(list.files(paste0("./temp_output/",especie,"/"), 
                                                paste0('bioclim_Fut.85',"PA_",PA,"RUN_",RUN,".grd"),
                                                full.names = T)))}}
    names(bioclim.fut.85) <- paste0("bioclim.fut.85.", 1:nlayers(bioclim.fut.85))
    
    # domain 85 ###
    # domain.fut.85 = stack()
    # for(PA in 1:PAs){
    #   for (RUN in 1:RUNs) {
    #     domain.fut.85 = stack(domain.fut.85, 
    #                           stack(list.files(paste0("./temp_output/",especie,"/"), 
    #                                            paste0('domain_Fut.85',"PA_",PA,"RUN_",RUN,".grd"),
    #                                            full.names = T)))}}
    # names(domain.fut.85) <- paste0("domain.fut.85.", 1:nlayers(domain.fut.85))
    
    # enfa 85 ###
    
    enfa.fut.85 = stack()
    for(PA in 1:PAs){
      for (RUN in 1:RUNs) {
        enfa.fut.85 = stack(enfa.fut.85, 
                            stack(list.files(paste0("./temp_output/",especie,"/"), 
                                             paste0('enfa_Fut.85',"PA_",PA,"RUN_",RUN,".grd"),
                                             full.names = T)))}}
    names(enfa.fut.85) <- paste0("enfa.fut.85.", 1:nlayers(enfa.fut.85))
    
    # glm 85 ###
    
    glm.fut.85 = stack()
    for(PA in 1:PAs){
      for (RUN in 1:RUNs) {
        glm.fut.85 = stack(glm.fut.85, 
                           stack(list.files(paste0("./temp_output/",especie,"/"), 
                                            paste0('glm_Fut.85',"PA_",PA,"RUN_",RUN,".grd"),
                                            full.names = T)))}}
    names(glm.fut.85) <- paste0("glm.fut.85.", 1:nlayers(glm.fut.85))
    
    # mars 85 ###
    
    mars.fut.85 = stack()
    for(PA in 1:PAs){
      for (RUN in 1:RUNs) {
        mars.fut.85 = stack(mars.fut.85, 
                            stack(list.files(paste0("./temp_output/",especie,"/"), 
                                             paste0('mars_Fut.85',"PA_",PA,"RUN_",RUN,".grd"),
                                             full.names = T)))}}
    names(mars.fut.85) <- paste0("mars.fut.85.", 1:nlayers(mars.fut.85))
    
    
    # maxent 85 ###
    
    maxent.fut.85 = stack()
    for(PA in 1:PAs){
      for (RUN in 1:RUNs) {
        maxent.fut.85 = stack(maxent.fut.85, 
                              stack(list.files(paste0("./temp_output/",especie,"/"), 
                                               paste0('maxent_Fut.85',"PA_",PA,"RUN_",RUN,".grd"),
                                               full.names = T)))}}
    names(maxent.fut.85) <- paste0("maxent.fut.85.", 1:nlayers(maxent.fut.85))
    
    # svm 85 ###
    
    svm.fut.85 = stack()
    for(PA in 1:PAs){
      for (RUN in 1:RUNs) {
        svm.fut.85 = stack(svm.fut.85, 
                           stack(list.files(paste0("./temp_output/",especie,"/"), 
                                            paste0('svm_Fut.85',"PA_",PA,"RUN_",RUN,".grd"),
                                            full.names = T)))}}
    names(svm.fut.85) <- paste0("svm.fut.85.", 1:nlayers(svm.fut.85))
    
    # nnet 85 ###
    
    nnet.fut.85 = stack()
    for(PA in 1:PAs){
      for (RUN in 1:RUNs) {
        nnet.fut.85 = stack(nnet.fut.85, 
                            stack(list.files(paste0("./temp_output/",especie,"/"), 
                                             paste0('nnet_Fut.85',"PA_",PA,"RUN_",RUN,".grd"),
                                             full.names = T)))}}
    names(nnet.fut.85) <- paste0("nnet.fut.85.", 1:nlayers(nnet.fut.85))
    
    # rf 85 ###
    
    rf.fut.85 = stack()
    for(PA in 1:PAs){
      for (RUN in 1:RUNs) {
        rf.fut.85 = stack(rf.fut.85, 
                          stack(list.files(paste0("./temp_output/",especie,"/"), 
                                           paste0('rf_Fut.85',"PA_",PA,"RUN_",RUN,".grd"),
                                           full.names = T)))}}
    
    names(rf.fut.85) <- paste0("rf.fut.85.", 1:nlayers(rf.fut.85))
    
    
    # Future Ensemble 85 ####
    gc(reset = T)
    
    Future.85.mean <- mean(subset(stack(bioclim.fut.85,enfa.fut.85, glm.fut.85, mars.fut.85,
                                        maxent.fut.85, svm.fut.85, nnet.fut.85, rf.fut.85)
                                  , sel2[, "ID"]))
    rm(bioclim.fut.85, enfa.fut.85, glm.fut.85, mars.fut.85,
       maxent.fut.85, svm.fut.85, nnet.fut.85, rf.fut.85)
    gc(reset = T)
    
    writeRaster(
      Future.85.mean,
      paste0("./outputs/", especie, "_", "Future.85.mean.tif"), formtat = "GTiff", overwrite = T)
    # Binary Transformation #
    Future.85.bin <- biomod2::BinaryTransformation(Future.85.mean, th.mean)
    rm(Future.85.mean)
    gc(reset = T)
    writeRaster(
      Future.85.bin,
      paste0("./outputs/", especie, "_", "Future.85.bin.tif"), formtat = "GTiff", overwrite = T)
    }
  
 
  rm(Future.85.bin, sel2)
  gc(reset = TRUE, full = TRUE)
  # unlink(especie,recursive = T, force = T)
  #--------------------#
  # Move the files  ###
  #------------------#
  file.move((list.files("./outputs/", paste0(especie, "_", "Current"),
                        full.names = TRUE)), (paste0("./outputs/", especie, '.', 'Presente')), 
                                                                                    overwrite = TRUE)
  
  filesstrings::file.move((list.files(
    "./outputs/", paste0(especie, "_", "Future.45"),
    full.names = TRUE
  )), (paste0("./outputs/", especie, '.', 'RCP4.5')), overwrite = TRUE)
  
  
  filesstrings::file.move((list.files("./outputs/",  paste0(especie, "_", "Future.85"),
                        full.names = TRUE)), (paste0("./outputs/", especie, '.', 'RCP8.5')), 
            overwrite = TRUE)
  
  filesstrings::file.move((list.files(
    "./outputs/", paste0(especie, "_"),
    full.names = TRUE
  )), (paste0("./outputs/", especie)), overwrite = TRUE)
  
  
  # Time Compute ###
  sink("./outputs/tempo.txt", append = T)
  print(Sys.time() - ini1)
  sink()
  print(paste0("Finishing", " ",especie," " ,"modeling"))
}##### End species loop ####

beep(sound = 2)
beep(sound = 2)
base::quit(save = "yes")
