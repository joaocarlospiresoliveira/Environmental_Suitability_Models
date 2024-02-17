###################### Acknowledgments ##########################
### Dr. Matheus de Souza Lima-Ribeiro's team of Universidade Federal de Jataí.
### Dr. Diogo Souza Bezerra Rocha (International Institute for Sustainability).
### Dr. João Carlos Pires de Oliveira

############################## ATENÇÃO #########################################
# Antes de modelar coloque dentro da pasta de trabalho (work directory - WD) uma 
# pasta com nome 'vars'. a pasta 'vars', deve ter as pastas 'Presente' e 'Future'.
###--------------------------------------------------------------------------###

############################# CHECAGEM ########################################
## 1 - Verificar se a pastar 'vars' e suas subpastas estão na Work Directory - WD
## 2 - Planilha de ocorrências deve conter apenas 3 colunas (sp, lon, lat)
## 3 - Coloque a palilha com as ocorrências nomaeda como 'spp.csv' na WD


## Install and Library packages ####
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
# install.packages("sdmvspecies")
# install.packages("filesstrings")
# install.packages("dplyr")
# install.packages("mda")
library(terra)
library(sp)
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
library(sdmvspecies)
library(filesstrings)
library(dplyr)
library(mda)
# Criating output dir #
if (dir.exists("outputs") == F) {
  dir.create("outputs")
}
# Carregando funcoes necessarias #
# Function to predict future predictions ###
preFut = function(rast, rast1 = NULL, model, GCM =  NULL) {
  if(missing(GCM)){
    GCM = 1
  }else{GCM = GCM}
  pre = foreach::foreach(
    gcm = 1:length(GCM),
    .combine = c,
    .packages = c("terra", 'sp', "kernlab","dismo",
                  "stats","randomForest", "nnet", "earth" )
  ) %do% {
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
      map <-rasterize(as.matrix(new.climate@coords),rast[[1]],maha) * -1
      return(invisible(map))
    }
    if (!"madifa" %in% class(model)) {
      terra::predict(rast[gcm][[1]], model, na.rm = T)
    } else if ("madifa" %in% class(model) && "list" %in% class(rast1)) {
      climaPres.spdf = na.omit(data.frame(xyFromCell(rast, 1:ncell(rast)),
                                          terra::values(rast)))
      climaFut.spdf = na.omit(data.frame(xyFromCell(rast1[gcm][[1]], 
                                                    1:ncell(rast1[gcm][[1]])),
                                         terra::values(rast1[gcm][[1]])))
      suppressWarnings(gridded(climaPres.spdf) <-~ x + y)
      suppressWarnings(gridded(climaFut.spdf) <- ~ x + y)
      predict.enfa( object.enfa = model,baseline.climate =climaPres.spdf,
                    new.climate = climaFut.spdf)}else{
                      climaPres.spdf = na.omit(data.frame(xyFromCell(rast, 1:ncell(rast)),
                                                          terra::values(rast)))
                      climaFut.spdf = na.omit(data.frame(xyFromCell(rast1, 
                                                                    1:ncell(rast1)),
                                                         terra::values(rast1)))
                      suppressWarnings(gridded(climaPres.spdf) <-~ x + y)
                      suppressWarnings(gridded(climaFut.spdf) <- ~ x + y)
                      predict.enfa( object.enfa = model,baseline.climate =climaPres.spdf,
                                    new.climate = climaFut.spdf)}
  }
  pre = mean(pre)
  names(pre) = c("lyr1")
  return(pre)
  rm(new.climate.scale,maha,map,climaPres.spdf,climaFut.spdf)
  gc(reset = T, full = T)
}

# Function to scale maps -- DON'T CHANGE####
eval.All.Model <- function(rast, dismoTestPrepared, tr ) {
  p = terra::extract(rast, dismoTestPrepared[dismoTestPrepared[, 3] == 1, 1:2 ], ID = F)
  a = terra::extract(rast, dismoTestPrepared[dismoTestPrepared[, 3] == 0, 1:2 ], ID = F)
  p <- stats::na.omit(p[[1]])
  a <- stats::na.omit(a[[1]])
  np <- length(p)
  na <- length(a)
  if (na == 0 | np == 0) {
    stop('cannot evaluate a model without absence and presence data that are not NA')
  }
  if (missing(tr)) {
    if (length(p) > 1000) {
      tr <- as.vector(quantile(p, 0:1000/1000))
    }
    else {
      tr <- p
    }
    if (length(a) > 1000) {
      tr <- c(tr, as.vector(quantile(a, 0:1000/1000)))
    }
    else {
      tr <- c(tr, a)
    }
    tr <- sort(unique(round(tr, 8)))
    tr <- c(tr - 1e-04, tr[length(tr)] + c(0, 1e-04))
  }
  else {
    tr <- sort(as.vector(tr))
  }
  N <- na + np
  xc <- new('ModelEvaluation')
  xc@presence = p
  xc@absence = a
  R <- sum(rank(c(p, a))[1:np]) - (np * (np + 1)/2)
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


# Function to scale maps -- DON'T CHANGE###
rescMod = function(ras1, ras2, ras3){
  names(ras1) = names(ras2) = names(ras3) = c("lyr1")
  p1 <- terra::as.data.frame(ras1)
  p2 <- terra::as.data.frame(ras2)
  p3 <- terra::as.data.frame(ras3)
  coord <- terra::as.data.frame(ras1, xy = T)[, 1:2]
  id <- rep(c("cur", "fut.45", "fut.85"), each = nrow(p1))
  rd = rbind(p1, p2, p3)
  st = vegan::decostand(x = rd[,1], method = "range", MARGIN = 2)
  ens = as.matrix(data.frame(coord, Cur = st[id == "cur"], 
                             Fut.45 = st[id =="fut.45"],
                             Fut.85 = st[id =="fut.85"]))
  ens.cur = terra::rasterize(ens[,c("x","y")],ras1[[1]], ens[,"Cur"])
  ens.fut.45 = terra::rasterize(ens[,c("x","y")],ras1[[1]], ens[,"Fut.45"])
  ens.fut.85 = terra::rasterize(ens[,c("x","y")],ras1[[1]], ens[,"Fut.85"])
  resul = c(ens.cur, ens.fut.45, ens.fut.85)
  names(resul)<- c("Cur", "Fut.45", "Fut.85")
  return(resul)
  rm(p1,p2,p3,coord,id,rd,st,ens,ens.fut.45, ens.fut.85)
  gc(reset = T)
}

# Function to Evaluate All Models --- DON'T CHANGE  ####

eval.All.Model <- function(rast, dismoTestPrepared, tr ) {
  p = terra::extract(rast, dismoTestPrepared[dismoTestPrepared[, 3] == 1, 1:2 ], ID = F)
  a = terra::extract(rast, dismoTestPrepared[dismoTestPrepared[, 3] == 0, 1:2 ], ID = F)
  p <- stats::na.omit(p[[1]])
  a <- stats::na.omit(a[[1]])
  np <- length(p)
  na <- length(a)
  if (na == 0 | np == 0) {
    stop('cannot evaluate a model without absence and presence data that are not NA')
  }
  if (missing(tr)) {
    if (length(p) > 1000) {
      tr <- as.vector(quantile(p, 0:1000/1000))
    }
    else {
      tr <- p
    }
    if (length(a) > 1000) {
      tr <- c(tr, as.vector(quantile(a, 0:1000/1000)))
    }
    else {
      tr <- c(tr, a)
    }
    tr <- sort(unique(round(tr, 8)))
    tr <- c(tr - 1e-04, tr[length(tr)] + c(0, 1e-04))
  }
  else {
    tr <- sort(as.vector(tr))
  }
  N <- na + np
  xc <- new('ModelEvaluation')
  xc@presence = p
  xc@absence = a
  R <- sum(rank(c(p, a))[1:np]) - (np * (np + 1)/2)
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

quiet <- function(x) { 
  sink(tempfile()) 
  on.exit(sink()) 
  invisible(force(x)) 
} 

### Whenever necessary:
# Parallel processing #
detectCores() # Identifica o numero de nucleos do processador
getDoParWorkers() # mostra quantos nucleos o R esta utilizando
cl <- # seleciona X nucleos do processador
  parallel::makeCluster(5, outfile = paste0("./outputs/", "log_models.log"))
#cl <- parallel::makeCluster(10, type = "MPI", outfile = "./outputs/joao.log")
registerDoParallel(cl) # registra o paralelismo
getDoParWorkers()

# # Increased memory allocation
# memory.limit(17592186044415) # or some other memory value (in kB)
# # memory.limit(8062000000000) # or some other memory value (in kB)


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


## Occurrences ####
# Select you species data matrix # 
spp <- read.csv("./Araucaria-ok.csv", header = T)  

# plot all your occurence points #
# plot(bio.crop[[1]]) 
 # points(spp[,c("lon","lat")], pch = 19, col = "red") # adiciona pontos de ocorrencia no mapa


dim(spp) # mostra numero de linha e colunas da matriz
head(spp, 10) # olhar o cabecario da matriz

table(spp$sp) # checando o numero de especies

especies <- unique(spp$sp) # cria objeto com os nomes das especies
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

PAs <- 5 # numero de conjuntos de pseudo-ausencias
RUNs = 5 # numero de repeticoes

GCM45s <- c("bio70_CA_45", "bio70_CN_45", "bio70_MI_45")
GCM85s <- c("bio70_CA_85", "bio70_CN_85", "bio70_MI_85")

# If you wish to compute the importance of variables, put VarImport = TRUE,
# or set VarImport = FALSE not to do this.
VarImport = FALSE

# especie = especies[1]
## Species Loop ####
# For sequential loop (One species) ###
# for (especie in especies[2:3]) {

## For species in parallel ### 
foreach(especie = especies, # For parallel looping (Multiple Species)
        .packages = c("terra",'sp',"sdmvspecies", "filesstrings",
                      "rgdal","maps","mnormt","kernlab","dismo","doParallel",
                      "stats","rJava","randomForest","nnet","psych", "earth"),
        .verbose = F,
        .errorhandling = "stop") %dopar% {
          # criando pasta de outputs temporarios #
          if (dir.exists(paste0("./temp_output/",especie,"/")) == F) {
            dir.create(paste0("./temp_output/"))
            dir.create(paste0("./temp_output/",especie,"/"))
          }
          
          ini1 = Sys.time()
          print(paste0("Starting", " ", especie, " " ,"modeling"))
          
          ##### Loading variaveis ####
          
          bio.crop<-
            list.files(
              "./vars/Present/PCA",  pattern = ".grd$",
              full.names = TRUE
            )
          bio.crop <- terra::rast(bio.crop) # para rodar somente com dados de clima substituir "bio.crop" por "bio.crop[-c(7,8)]"
          names(bio.crop)<- c("PCA1", "PCA2", "PCA3","PCA4", "PCA5","PCA6", "PCA7", "PCA8") # quando rodar apenas com clima remover o trecho (, "PCA7", "PCA8") dessa linha 
          names(bio.crop)
          
          #--------------------------------------------------#
          ####      LOADING FUTURE VARIABLES             ####
          #------------------------------------------------#
          ###GCM 1: CanESM5
          
          bio70_CA_45 <-
            list.files(
              "./vars/Future/CanESM5/ssp245/PCA",
              pattern = ".grd$",
              full.names = TRUE
            )
          bio70_CA_45 <- terra::rast(bio70_CA_45)
          names(bio70_CA_45) <-
            c("PCA1", "PCA2", "PCA3","PCA4", "PCA5","PCA6", "PCA7", "PCA8")
          names(bio70_CA_45)
          
          
          ###GCM 2: CNRM-CM6-1
          bio70_CN_45 <-
            list.files(
              "./vars/Future/CNRM-CM6-1/ssp245/PCA",
              pattern = ".grd$",
              full.names = TRUE
            )
          bio70_CN_45 <- terra::rast(bio70_CN_45)
          names(bio70_CN_45) <-
            c("PCA1", "PCA2", "PCA3","PCA4", "PCA5","PCA6", "PCA7", "PCA8")
          names(bio70_CN_45)
          
          ###GCM 3: MIROC-ES2L
          
          bio70_MI_45 <-
            list.files(
              "./vars/Future/MIROC-ES2L/ssp245/PCA",
              pattern = ".grd$",
              full.names = TRUE
            )
          bio70_MI_45 <- terra::rast(bio70_MI_45)
          names(bio70_MI_45) <-
            c("PCA1", "PCA2", "PCA3","PCA4", "PCA5","PCA6", "PCA7", "PCA8")
          names(bio70_MI_45)
          
          # Loading raster layer with future projections 85 ###
          #------------------------------------------------#
          
          ###GCM 1: CanESM5
          bio70_CA_85 <-
            list.files(
              "./vars/Future/CanESM5/ssp585/PCA",
              pattern = ".grd$",
              full.names = TRUE
            )
          bio70_CA_85 <- terra::rast(bio70_CA_85)
          names(bio70_CA_85) <-
            c("PCA1", "PCA2", "PCA3","PCA4", "PCA5","PCA6", "PCA7", "PCA8")
          names(bio70_CA_85)
          
          
          ###GCM 2: CNRM-CM6-1
          bio70_CN_85 <-
            list.files(
              "./vars/Future/CNRM-CM6-1/ssp585/PCA",
              pattern = ".grd$",
              full.names = TRUE
            )
          bio70_CN_85 <- terra::rast(bio70_CN_85)
          names(bio70_CN_85) <-
            c("PCA1", "PCA2", "PCA3","PCA4", "PCA5","PCA6", "PCA7", "PCA8")
          names(bio70_CN_85)
          
          ###GCM 3: MIROC-ES2L
          bio70_MI_85 <-
            list.files(
              "./vars/Future/MIROC-ES2L/ssp585/PCA",
              pattern = ".grd$",
              full.names = TRUE
            )
          bio70_MI_85 <- terra::rast(bio70_MI_85)
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
          # remover objetos desnecessarios #
          rm( bio70_CA_45,
              bio70_CN_45,
              bio70_MI_45,
              bio70_CA_85,
              bio70_CN_85,
              bio70_MI_85)
          
          
          # Creating empty objects to store evaluation models ###
          bioclim.e <-
            # domain.e <-
            enfa.e <-
            glm.e <- mars.e <- maxent.e <- nnet.e <- svm.e <- rf.e <- NULL
          
          ### Checagem das coordenadas ###
          occs <- spp[spp$sp == especie, c("lon", "lat")] # selecionar apenas lon e lat
          
          # Data checking and preparation
          ocor.val <- terra::extract(bio.crop, occs, cells = T)
          
          sum(is.na(ocor.val[, 1])) # retorna o numero de linhas que contem NA
          
          # plot(bio.crop[[1]], colNA = "red")
          
          # points(occs, pch = 19, cex = 0.5)
          
          ocor.val <- cbind(occs, ocor.val) # unir coordenadas com a matriz de dados
          
          ocor.val <- na.omit(ocor.val) # exclui linhas com NA
          
          id <- duplicated(ocor.val[, "cell"]) # Checking dumplicate points
          
          sum(id == T) # somar o numero de duplicatas
          
          ocor <-
            ocor.val[id == F, c("lon", "lat")] # Removing duplicate points
          
          # Salvando matriz com pontos espacialmente unicos #
          
          unique_points <- data.frame(sp = especie, ocor)
          write.csv(unique_points, paste0(especie, "_spatial_uniqe.csv"),
                    row.names = F)
          #------------------------------------------#
          #           SELECT PAs                  ###
          #----------------------------------------#
          # computando distancias para limitar o buffer de amostragem das pseudo-ausencias #
          try({
            coord1 = ocor
            sp::coordinates(coord1) <- ~ lon + lat
            terra::crs(coord1) <- terra::crs(bio.crop)
            
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
          PA.number #numero de pontos de ocorrencia espacialmente unicos
          
          diretorio = paste0("Occurrence.", especie)
          
          bioclim.var.All <- domain.var.All <- enfa.var.All <- glm.var.All <-
            mars.var.All <- maxent.var.All <- svm.var.All <- 
            nnet.var.All <- rf.var.All <- NULL
          
          #### Loop PAs ####
          for (PA in seq(PAs)) {
            
            p = rep(1, times = nrow(ocor))
            # gerar pseudo-ausencias 
            invisible(capture.output(sel.PA <- biomod2::BIOMOD_FormatingData(
              resp.var = p,
              expl.var = bio.crop,
              resp.xy = ocor,
              resp.name = diretorio,
              PA.nb.rep = 1,
              #numero de datasets de pseudo-ausencias
              PA.nb.absences = PA.number,
              #= numero de pseudo-ausencias = numero de pontos espacialmente unicos
              PA.strategy = "disk",
              PA.dist.min = dist.min * 1000,
              PA.dist.max = dist.mean * 1000,
              na.rm = TRUE
            )))
            
            li.p <- grep("pa", rownames(sel.PA@coord))
            
            
            pa <- sel.PA@coord[li.p,]
            # gerar o background #
            invisible(capture.output(sel.back <- biomod2::BIOMOD_FormatingData(
              resp.var = p,
              expl.var = bio.crop,
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
            areaToral <- nrow(terra::as.data.frame(bio.crop))
            
            #### Loop RUN ####
            for (RUN in seq(RUNs)) {
              print(paste0(especie," ","PASet", PA," ","RUN", RUN))
              # Separating test/ training data
              
              
              id.training.pa <-
                sample(1:nrow(ocor), round(0.7 * nrow(ocor), 0)) # prepare data 70/30
              
              
              
              # preparando dados de treino #
              training.b <- # preparando dados para modelar com background
                na.omit(dismo::prepareData(bio.crop, p = ocor[id.training.pa, ], 
                                           b = back, xy = T)[,-4])
              # head(training.b)
              # tail(training.b)
              training.pa <- # preparando dados para modelar com pseudo-ausencias
                na.omit(dismo::prepareData(bio.crop, p = ocor[id.training.pa, ], 
                                           b = pa[id.training.pa, ], xy = T)[,-4])
              # preparando dados de teste # 
              test.pa <-
                na.omit(dismo::prepareData(bio.crop, p = ocor[-id.training.pa, ], 
                                           b = pa[-id.training.pa, ], xy = T)[,-4])
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
                terra::predict(bio.crop, 
                               bioclim_model)
              
              # Future 45 #
              bioclim_Fut.45 <- preFut(rast = scenario.list.45,
                                       model =  bioclim_model, GCM = GCM45s)
              
              # Future 85 #
              bioclim_Fut.85 = preFut(rast = scenario.list.85,
                                      model =  bioclim_model, GCM = GCM85s)
              
              # padronizando projecoes # 
              bioclim_std = rescMod(bioclim_Cur, bioclim_Fut.45, bioclim_Fut.85)
              # renomeando os mapas #
              names(bioclim_std)<-c(paste0("Cur", ".",PA, ".", RUN), paste0("Fut.45", ".",PA, ".", RUN),
                                    paste0("Fut.85", ".",PA, ".", RUN))
              # salvando resultados temporarios #
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
              # avaliou o modelo utilizando multiplos limiares/thresholds
              bioclim_eval <-
                eval.All.Model(subset(bioclim_std, grep("Cur", names(bioclim_std))), 
                               test.pa)
              # escolhendo o limiar/threshold que maximiza acertos especificidade e sensibilidade
              bioclim_th.spec_sens <-
                dismo::threshold(bioclim_eval, "spec_sens")
              
              # bioclim_th.LPT5 <-
              #   quantile(terra::extract(bioclim_Cur, ocor[id.training.pa,]), 0.05, na.rm = T)
              # 
              # bioclim_th.VDl <-
              #   com.vdl(bioclim_Cur, test.pa, bioclim_eval)@maxVDl
              
              # cumputa valores das metricas para o limiar/threshold selecionado
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
                    
                    bioclim_VarImp1 <-(terra::predict(object = bio.crop[[-rmo]], 
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
              }
              
              
              remove(list = ls()[c(grep(
                "bioclim_", as.factor(ls())))])
              
              gc(reset = TRUE, full = T)
              
               ##### Domain  ####
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
              #   terra::predict(bio.crop, domain_model)
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
              # #   quantile(terra::extract(domain_Cur, ocor[id.training.pa,]), 0.05, na.rm = T)
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
              #     domain_VarImp1 <-(terra::predict(object = bio.crop[[-rmo]], 
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
              #                       domain_var.part.o) }
              # 
              # remove(list = ls()[c(grep(
              #   "domain_", as.factor(ls())))])
              # 
              # gc(reset = TRUE, full = T)
              
              
              ##### ENFA  #### 
              print(paste0(especie, " ","ENFA"," ", "PA",  PA," ","RUN", RUN))
              
              climaPres <- (bio.crop)
              # names(climaPres)# <- paste0("PC",1:6)
              
              climaPres.values <- terra::values(climaPres) # obtendo os valores das PCAs
              climaPres.spdf <- na.omit(data.frame(terra::xyFromCell(climaPres, 1:ncell(climaPres)), 
                                                   climaPres.values))
              
              suppressWarnings(gridded(climaPres.spdf) <- ~x+y)
              climaPres <- rast(climaPres.spdf)
              climaPres.values <- terra::values(climaPres)
              media.climaPres <- apply(slot(climaPres.spdf, "data"), 2, mean)
              sd.climaPres <- apply(slot(climaPres.spdf, "data"), 2, sd)
              climaPres.scale<- sweep(slot(climaPres.spdf, "data"),2, media.climaPres)
              climaPres.scale<- as.matrix(climaPres.scale) %*% diag(1/sd.climaPres)
              
              #adjustment of the ENFA model
              
              pr.cell <- terra::extract(climaPres, training.pa[training.pa[,"pb"]==1,1:2], cells=T, ID = F)
              pr <- data.frame(pr= rep(0, ncell(climaPres)), climaPres.values)
              pr[pr.cell[,"cell"], 1] <- 1
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
              #   quantile(terra::extract(enfa_Cur, ocor[id.training.pa,]), 0.05, na.rm = T)
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
                    
                    
                    enfa_climaPr <- terra::c(bio.crop[[-rmo]])
                    names(enfa_climaPr)# <- paste0("PC",1:6)
                    
                    enfa_climaPr.values <- terra::values(enfa_climaPr)
                    enfa_climaPr.spdf <-
                      na.omit(data.frame(xyFromCell(enfa_climaPr, 1:ncell(enfa_climaPr)), 
                                         enfa_climaPr.values))
                    
                    enfa_climaPr<-rasterize(enfa_climaPr.spdf[,1:2],enfa_Cur[[1]],
                                            enfa_climaPr.spdf[,-c(1:2)])
                    
                    suppressWarnings(gridded(enfa_climaPr.spdf)<-~x+y)
                    
                    enfa_climaPr <- terra::c(enfa_climaPr.spdf)
                    enfa_climaPr.values <- terra::values(enfa_climaPr)
                    media.enfa_climaPr <- apply(slot(enfa_climaPr.spdf, "data"), 2, mean)
                    sd.enfa_climaPr <- apply(slot(enfa_climaPr.spdf, "data"), 2, sd)
                    enfa_climaPr.scale<- sweep(slot(enfa_climaPr.spdf, "data"),2, 
                                               media.enfa_climaPr)
                    enfa_climaPr.scale <-
                      as.matrix(enfa_climaPr.scale) %*% diag(1 / sd.enfa_climaPr)
                    
                    enfa_pr.cel <- terra::extract(enfa_climaPr, ocor, cells=T)
                    enfa_pr1 <-
                      data.frame(enfa_pr1 = rep(0, ncell(enfa_climaPr)), enfa_climaPr.values)
                    enfa_pr1[enfa_pr.cel[,"cell"], 1] <- 1
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
              }
              
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
                terra::predict(bio.crop, glm_model)
              
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
              #   quantile(terra::extract(glm_Cur, ocor[id.training.pa,]), 0.05, na.rm = T)
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
                    glm_PCA <- terra::c(bio.crop[[-rmo]])
                    
                    glm_VarImp1 <-rescMod.One(terra::predict(glm_PCA, glm_mod))
                    
                    
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
              }
              
              remove(list = ls()[c(grep(
                "glm_", as.factor(ls())))])
              gc(reset = TRUE, full = T)
              
              ##### MARS ####
              print(paste0(especie, " ","MARS"," ", "PA",  PA," ","RUN", RUN))
              
              pb <- training.pa$pb
              
              mars_model <- earth::earth(
                pb ~ .,
                data = training.pa[,-c(1:3)],
                # type = 'simple',
                # interaction.level = 0,
                penalty = 1,
                nprune = NULL,
                pmethod = 'backward'
              )
              
              # Building Preojections ###
              
              # Current #
              mars_Cur <-
                terra::predict(bio.crop, mars_model)
              
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
              #   quantile(terra::extract(mars_Cur, ocor[id.training.pa,]), 0.05, na.rm = T)
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
                    mars_PCA <- terra::c(bio.crop[[-rmo]])
                    
                    mars_VarImp1 <-rescMod.One(terra::predict(mars_PCA, mars_mod))
                    
                    
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
              }
              
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
                defaultprevalence = 0.5,
                na.rm = TRUE
              ))
              # Building Preojections ###
              
              # Current #
              maxent_Cur <-
                quiet(terra::predict(bio.crop, maxent_model, na.rm = T))
              
              
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
              #   quantile(terra::extract(maxent_Cur, ocor[id.training.pa,]), 0.05, na.rm = T)
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
                    maxent_PCA <- terra::c(bio.crop[[-rmo]])
                    
                    maxent_VarImp1 <-quiet(rescMod.One(terra::predict(maxent_PCA, 
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
              }
              
              remove(list = ls()[c(grep(
                "maxent_", as.factor(ls())))])
              gc(reset = TRUE, full = T)
              
              
              ##### SVM ####
              print(paste0(especie, " ","SVM"," ",  "PA",  PA," ","RUN", RUN))
              
              svm_model <- ksvm(pb ~ ., data = training.b[,-c(1:2)])
              
              # Building Preojections ###
              
              # Current #
              svm_Cur <-
                terra::predict(bio.crop, svm_model, na.rm = T)
              
              
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
              #   quantile(terra::extract(svm_Cur, ocor[id.training.b,]), 0.05, na.rm = T)
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
                    
                    svm_PCA <- terra::c(bio.crop[[-rmo]])
                    
                    svm_VarImp1 <-(rescMod.One(terra::predict(svm_PCA,
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
              }
              
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
                terra::predict(bio.crop, nnet_model, na.rm = T)
              
              
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
              #   quantile(terra::extract(nnet_Cur, ocor[id.training.pa,]), 0.05, na.rm = T)
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
                    
                    nnet_PCA <- terra::c(bio.crop[[-rmo]])
                    
                    nnet_VarImp1 <-(rescMod.One(terra::predict(nnet_PCA, 
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
              }
              
              remove(list = ls()[c(grep(
                "nnet_", as.factor(ls())))])
              gc(reset = TRUE, full = T)
              
              
              ##### Random Forest ####
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
                terra::predict(bio.crop, rf_model, na.rm = T)
              
              
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
              #   quantile(terra::extract(rf_Cur, ocor[id.training.pa,]), 0.05, na.rm = T)
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
                    
                    rf_PCA <- terra::c(bio.crop[[-rmo]])
                    
                    rf_VarImp1 <-(rescMod.One(terra::predict(rf_PCA, 
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
              }
              
              remove(list = ls()[c(grep(
                "rf_", as.factor(ls())))])
              gc(reset = TRUE, full = T)
              
            }##### End RUN loop ###
            gc(reset = TRUE, full = TRUE)
            # print(paste0(especie," ","PA set"," ", PA))
          } ##### End PAs loop ###
          gc(reset = TRUE, full = TRUE)
          
         
          ## Writing All Reults ####
          
          ## save VarImpor ####
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
            
          }
          
          
          
          # Writing Evaluation Results #
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
            # unindo matrizes de avaliacao #
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
            
     
            
            # Selectning "good models" #
            
            sel2 = Evaluation.all[Evaluation.all[, "TSS"] > 0.400, ]
            sel2 <- na.omit(sel2)
            write.csv(sel2,
                      paste0("./outputs/", especie, "_", "Selected.Models.csv"), row.names = T)
            
            
            rm(Evaluation.all)
            gc(reset = TRUE, full = TRUE)
          }
          ## Write Rasters ##
          # {
            # Current ####
            # bioclim ###
            bioclim.cur = c()
            for(PA in 1:PAs){
              for (RUN in 1:RUNs) {
                bioclim.cur = c(bioclim.cur,c(list.files(paste0("./temp_output/",especie,"/"), 
                                                         paste0('bioclim_Cur',"PA_",PA,"RUN_",RUN,".grd$"),
                                                         full.names = T)))
              }}
            bioclim.cur = rast(bioclim.cur)
            names(bioclim.cur) <- paste0("bioclim.cur.", 1:nlyr(bioclim.cur))
            
            # # domain ###
            # domain.cur = c()
            # for(PA in 1:PAs){
            #   for (RUN in 1:RUNs) {
            #     domain.cur = c(domain.cur,c(list.files(paste0("./temp_output/",especie,"/"), 
            #                                            paste0('domain_Cur',"PA_",PA,"RUN_",RUN,".grd$"),
            #                                            full.names = T)))
            #   }}
            # domain.cur = rast(domain.cur)
            # names(domain.cur) <- paste0("domain.cur.", 1:nlyr(domain.cur))
            
            # enfa ###
            enfa.cur = c()
            for(PA in 1:PAs){
              for (RUN in 1:RUNs) {
                enfa.cur = c(enfa.cur,c(list.files(paste0("./temp_output/",especie,"/"), 
                                                   paste0('enfa_Cur',"PA_",PA,"RUN_",RUN,".grd$"),
                                                   full.names = T)))
              }}
            enfa.cur = rast(enfa.cur)
            names(enfa.cur) <- paste0("enfa.cur.", 1:nlyr(enfa.cur))
            
            # glm ###
            glm.cur = c()
            for(PA in 1:PAs){
              for (RUN in 1:RUNs) {
                glm.cur = c(glm.cur,c(list.files(paste0("./temp_output/",especie,"/"), 
                                                 paste0('glm_Cur',"PA_",PA,"RUN_",RUN,".grd$"),
                                                 full.names = T)))
              }}
            glm.cur = rast(glm.cur)
            names(glm.cur) <- paste0("glm.cur.", 1:nlyr(glm.cur))
            
            # mars ###
            mars.cur = c()
            for(PA in 1:PAs){
              for (RUN in 1:RUNs) {
                mars.cur = c(mars.cur,c(list.files(paste0("./temp_output/",especie,"/"), 
                                                   paste0('mars_Cur',"PA_",PA,"RUN_",RUN,".grd$"),
                                                   full.names = T)))
              }}
            mars.cur = rast(mars.cur)
            names(mars.cur) <- paste0("mars.cur.", 1:nlyr(mars.cur))
            
            # maxent ###
            maxent.cur = c()
            for(PA in 1:PAs){
              for (RUN in 1:RUNs) {
                maxent.cur = c(maxent.cur,c(list.files(paste0("./temp_output/",especie,"/"), 
                                                       paste0('maxent_Cur',"PA_",PA,"RUN_",RUN,".grd$"),
                                                       full.names = T)))
              }}
            maxent.cur = rast(maxent.cur)
            names(maxent.cur) <- paste0("maxent.cur.", 1:nlyr(maxent.cur))
            
            # svm ###
            svm.cur = c()
            for(PA in 1:PAs){
              for (RUN in 1:RUNs) {
                svm.cur = c(svm.cur,c(list.files(paste0("./temp_output/",especie,"/"), 
                                                 paste0('svm_Cur',"PA_",PA,"RUN_",RUN,".grd$"),
                                                 full.names = T)))
              }}
            svm.cur = rast(svm.cur)
            names(svm.cur) <- paste0("svm.cur.", 1:nlyr(svm.cur))
            
            # nnet ###
            nnet.cur = c()
            for(PA in 1:PAs){
              for (RUN in 1:RUNs) {
                nnet.cur = c(nnet.cur,c(list.files(paste0("./temp_output/",especie,"/"), 
                                                   paste0('nnet_Cur',"PA_",PA,"RUN_",RUN,".grd$"),
                                                   full.names = T)))
              }}
            nnet.cur = rast(nnet.cur)
            names(nnet.cur) <- paste0("nnet.cur.", 1:nlyr(nnet.cur))
            
            # rf ###
            rf.cur = c()
            for(PA in 1:PAs){
              for (RUN in 1:RUNs) {
                rf.cur = c(rf.cur,c(list.files(paste0("./temp_output/",especie,"/"), 
                                               paste0('rf_Cur',"PA_",PA,"RUN_",RUN,".grd$"),
                                               full.names = T)))
              }}
            rf.cur = rast(rf.cur)
            names(rf.cur) <- paste0("rf.cur.", 1:nlyr(rf.cur))
            
            # Current Ensemble ####
            Current.all = c(bioclim.cur, enfa.cur, glm.cur, mars.cur, 
                                maxent.cur, svm.cur, nnet.cur, rf.cur)
            # fazendo a media ponderada pelo TSS #
            Current.mean <- terra::weighted.mean(subset(Current.all,
                                        sel2[, "ID"]), sel2$TSS)
            writeRaster(
              Current.mean,
              paste0("./outputs/", especie, "_", "Current.mean.tif"), overwrite = T)
            
            # Binary Transformation #
            Current.bin <- mean(subset(Current.all,
                                                                              sel2[, "ID"]) >= 
                                                                       sel2$threshold
                                )
            # Current.bin = c(Current.bin == terra::maxValue(Current.bin))
            
            rm(Current.mean, Current.all)
            gc(reset = T, full = T)
            writeRaster(
              Current.bin,
              paste0("./outputs/", especie, "_", "Current.bin.tif"), overwrite = T)
            rm(Current.bin,bioclim.cur, enfa.cur, glm.cur, mars.cur, 
               maxent.cur, svm.cur, nnet.cur, rf.cur)
            gc(reset = TRUE, full = TRUE)
            # }
          # Future 45 ####
          # {
          # bioclim 45 ###
          
          bioclim.fut.45 = c()
          for(PA in 1:PAs){
            for (RUN in 1:RUNs) {
              bioclim.fut.45 = c(bioclim.fut.45, 
                                 c(list.files(paste0("./temp_output/",especie,"/"), 
                                              paste0('bioclim_Fut.45',"PA_",PA,"RUN_",RUN,".grd$"),
                                              full.names = T)))}}
          bioclim.fut.45 = rast(bioclim.fut.45)
          names(bioclim.fut.45) <- paste0("bioclim.fut.45.", 1:nlyr(bioclim.fut.45))
          
          
          # # domain 45 ###
          
          # domain.fut.45 = c()
          # for(PA in 1:PAs){
          #   for (RUN in 1:RUNs) {
          #     domain.fut.45 = c(domain.fut.45, 
          #                       c(list.files(paste0("./temp_output/",especie,"/"), 
          #                                    paste0('domain_Fut.45',"PA_",PA,"RUN_",RUN,".grd$"),
          #                                    full.names = T)))}}
          # domain.fut.45 = rast(domain.fut.45)
          # names(domain.fut.45) <- paste0("domain.fut.45.", 1:nlyr(domain.fut.45))
          
          # enfa 45 ###
          
          enfa.fut.45 = c()
          for(PA in 1:PAs){
            for (RUN in 1:RUNs) {
              enfa.fut.45 = c(enfa.fut.45, 
                              c(list.files(paste0("./temp_output/",especie,"/"), 
                                           paste0('enfa_Fut.45',"PA_",PA,"RUN_",RUN,".grd$"),
                                           full.names = T)))}}
          enfa.fut.45 = rast(enfa.fut.45)
          names(enfa.fut.45) <- paste0("enfa.fut.45.", 1:nlyr(enfa.fut.45))
          
          # glm 45 ###
          
          glm.fut.45 = c()
          for(PA in 1:PAs){
            for (RUN in 1:RUNs) {
              glm.fut.45 = c(glm.fut.45, 
                             c(list.files(paste0("./temp_output/",especie,"/"), 
                                          paste0('glm_Fut.45',"PA_",PA,"RUN_",RUN,".grd$"),
                                          full.names = T)))}}
          glm.fut.45 = rast(glm.fut.45)
          names(glm.fut.45) <- paste0("glm.fut.45.", 1:nlyr(glm.fut.45))
          
          # mars 45 ###
          
          mars.fut.45 = c()
          for(PA in 1:PAs){
            for (RUN in 1:RUNs) {
              mars.fut.45 = c(mars.fut.45, 
                              c(list.files(paste0("./temp_output/",especie,"/"), 
                                           paste0('mars_Fut.45',"PA_",PA,"RUN_",RUN,".grd$"),
                                           full.names = T)))}}
          mars.fut.45 = rast(mars.fut.45)
          names(mars.fut.45) <- paste0("mars.fut.45.", 1:nlyr(mars.fut.45))
          
          # maxent 45 ###
          
          maxent.fut.45 = c()
          for(PA in 1:PAs){
            for (RUN in 1:RUNs) {
              maxent.fut.45 = c(maxent.fut.45, 
                                c(list.files(paste0("./temp_output/",especie,"/"), 
                                             paste0('maxent_Fut.45',"PA_",PA,"RUN_",RUN,".grd$"),
                                             full.names = T)))}}
          maxent.fut.45 = rast(maxent.fut.45)
          names(maxent.fut.45) <- paste0("maxent.fut.45.", 1:nlyr(maxent.fut.45))
          
          # svm 45 ###
          
          svm.fut.45 = c()
          for(PA in 1:PAs){
            for (RUN in 1:RUNs) {
              svm.fut.45 = c(svm.fut.45, 
                             c(list.files(paste0("./temp_output/",especie,"/"), 
                                          paste0('svm_Fut.45',"PA_",PA,"RUN_",RUN,".grd$"),
                                          full.names = T)))}}
          svm.fut.45 = rast(svm.fut.45)
          names(svm.fut.45) <- paste0("svm.fut.45.", 1:nlyr(svm.fut.45))
          
          # nnet 45 ###
          
          nnet.fut.45 = c()
          for(PA in 1:PAs){
            for (RUN in 1:RUNs) {
              nnet.fut.45 = c(nnet.fut.45, 
                              c(list.files(paste0("./temp_output/",especie,"/"), 
                                           paste0('nnet_Fut.45',"PA_",PA,"RUN_",RUN,".grd$"),
                                           full.names = T)))}}
          nnet.fut.45 = rast(nnet.fut.45)
          names(nnet.fut.45) <- paste0("nnet.fut.45.", 1:nlyr(nnet.fut.45))
          
          # rf 45 ###
          
          rf.fut.45 = c()
          for(PA in 1:PAs){
            for (RUN in 1:RUNs) {
              rf.fut.45 = c(rf.fut.45, 
                            c(list.files(paste0("./temp_output/",especie,"/"), 
                                         paste0('rf_Fut.45',"PA_",PA,"RUN_",RUN,".grd$"),
                                         full.names = T)))}}
          rf.fut.45 = rast(rf.fut.45)
          names(rf.fut.45) <- paste0("rf.fut.45.", 1:nlyr(rf.fut.45))
            
            ## Future Ensemble 45 ####
            
            Future.45.all = c(bioclim.fut.45,enfa.fut.45, glm.fut.45, mars.fut.45,
                                  maxent.fut.45, svm.fut.45, nnet.fut.45, rf.fut.45)
            
            Future.45.mean <- terra::weighted.mean(subset(Future.45.all, 
                                                           sel2[, "ID"]), sel2$TSS)
            
            rm(bioclim.fut.45,enfa.fut.45, glm.fut.45, mars.fut.45,
               maxent.fut.45, svm.fut.45, nnet.fut.45, rf.fut.45)
            
            writeRaster(
              Future.45.mean,
              paste0("./outputs/", especie, "_", "Future.45.mean.tif"), overwrite = T)
            
            # Binary Transformation #
            Future.45.bin <- terra::mean(subset(Future.45.all,
                                                  sel2[, "ID"]) >= sel2$threshold
                                          )
            
            # Future.45.bin = c(Future.45.bin == terra::maxValue(Future.45.bin))
            

            writeRaster(
              Future.45.bin,
              paste0("./outputs/", especie, "_", "Future.45.bin.tif"), overwrite = T)
            rm(Future.45.mean, Future.45.bin, Future.45.all)
            gc(reset = TRUE)
            # }
          
            # Future 85 ####
            # {
            # bioclim 85 ###
            
            bioclim.fut.85 = c()
            for(PA in 1:PAs){
              for (RUN in 1:RUNs) {
                bioclim.fut.85 = c(bioclim.fut.85, 
                                   c(list.files(paste0("./temp_output/",especie,"/"), 
                                                paste0('bioclim_Fut.85',"PA_",PA,"RUN_",RUN,".grd$"),
                                                full.names = T)))}}
            bioclim.fut.85 = rast(bioclim.fut.85)
            names(bioclim.fut.85) <- paste0("bioclim.fut.85.", 1:nlyr(bioclim.fut.85))
            
            
            # # domain 85 ###
            
            # domain.fut.85 = c()
            # for(PA in 1:PAs){
            #   for (RUN in 1:RUNs) {
            #     domain.fut.85 = c(domain.fut.85, 
            #                       c(list.files(paste0("./temp_output/",especie,"/"), 
            #                                    paste0('domain_Fut.85',"PA_",PA,"RUN_",RUN,".grd$"),
            #                                    full.names = T)))}}
            # domain.fut.85 = rast(domain.fut.85)
            # names(domain.fut.85) <- paste0("domain.fut.85.", 1:nlyr(domain.fut.85))
            
            # enfa 85 ###
            
            enfa.fut.85 = c()
            for(PA in 1:PAs){
              for (RUN in 1:RUNs) {
                enfa.fut.85 = c(enfa.fut.85, 
                                c(list.files(paste0("./temp_output/",especie,"/"), 
                                             paste0('enfa_Fut.85',"PA_",PA,"RUN_",RUN,".grd$"),
                                             full.names = T)))}}
            enfa.fut.85 = rast(enfa.fut.85)
            names(enfa.fut.85) <- paste0("enfa.fut.85.", 1:nlyr(enfa.fut.85))
            
            # glm 85 ###
            
            glm.fut.85 = c()
            for(PA in 1:PAs){
              for (RUN in 1:RUNs) {
                glm.fut.85 = c(glm.fut.85, 
                               c(list.files(paste0("./temp_output/",especie,"/"), 
                                            paste0('glm_Fut.85',"PA_",PA,"RUN_",RUN,".grd$"),
                                            full.names = T)))}}
            glm.fut.85 = rast(glm.fut.85)
            names(glm.fut.85) <- paste0("glm.fut.85.", 1:nlyr(glm.fut.85))
            
            # mars 85 ###
            
            mars.fut.85 = c()
            for(PA in 1:PAs){
              for (RUN in 1:RUNs) {
                mars.fut.85 = c(mars.fut.85, 
                                c(list.files(paste0("./temp_output/",especie,"/"), 
                                             paste0('mars_Fut.85',"PA_",PA,"RUN_",RUN,".grd$"),
                                             full.names = T)))}}
            mars.fut.85 = rast(mars.fut.85)
            names(mars.fut.85) <- paste0("mars.fut.85.", 1:nlyr(mars.fut.85))
            
            # maxent 85 ###
            
            maxent.fut.85 = c()
            for(PA in 1:PAs){
              for (RUN in 1:RUNs) {
                maxent.fut.85 = c(maxent.fut.85, 
                                  c(list.files(paste0("./temp_output/",especie,"/"), 
                                               paste0('maxent_Fut.85',"PA_",PA,"RUN_",RUN,".grd$"),
                                               full.names = T)))}}
            maxent.fut.85 = rast(maxent.fut.85)
            names(maxent.fut.85) <- paste0("maxent.fut.85.", 1:nlyr(maxent.fut.85))
            
            # svm 85 ###
            
            svm.fut.85 = c()
            for(PA in 1:PAs){
              for (RUN in 1:RUNs) {
                svm.fut.85 = c(svm.fut.85, 
                               c(list.files(paste0("./temp_output/",especie,"/"), 
                                            paste0('svm_Fut.85',"PA_",PA,"RUN_",RUN,".grd$"),
                                            full.names = T)))}}
            svm.fut.85 = rast(svm.fut.85)
            names(svm.fut.85) <- paste0("svm.fut.85.", 1:nlyr(svm.fut.85))
            
            # nnet 85 ###
            
            nnet.fut.85 = c()
            for(PA in 1:PAs){
              for (RUN in 1:RUNs) {
                nnet.fut.85 = c(nnet.fut.85, 
                                c(list.files(paste0("./temp_output/",especie,"/"), 
                                             paste0('nnet_Fut.85',"PA_",PA,"RUN_",RUN,".grd$"),
                                             full.names = T)))}}
            nnet.fut.85 = rast(nnet.fut.85)
            names(nnet.fut.85) <- paste0("nnet.fut.85.", 1:nlyr(nnet.fut.85))
            
            # rf 85 ###
            
            rf.fut.85 = c()
            for(PA in 1:PAs){
              for (RUN in 1:RUNs) {
                rf.fut.85 = c(rf.fut.85, 
                              c(list.files(paste0("./temp_output/",especie,"/"), 
                                           paste0('rf_Fut.85',"PA_",PA,"RUN_",RUN,".grd$"),
                                           full.names = T)))}}
            rf.fut.85 = rast(rf.fut.85)
            names(rf.fut.85) <- paste0("rf.fut.85.", 1:nlyr(rf.fut.85))
          
            ## Future Ensemble 85 ####
            
            Future.85.all = c(bioclim.fut.85,enfa.fut.85, glm.fut.85, mars.fut.85,
                                  maxent.fut.85, svm.fut.85, nnet.fut.85, rf.fut.85)
            
            Future.85.mean <- terra::weighted.mean(subset(Future.85.all, 
                                                           sel2[, "ID"]), sel2$TSS)
            
            rm(bioclim.fut.85,enfa.fut.85, glm.fut.85, mars.fut.85,
               maxent.fut.85, svm.fut.85, nnet.fut.85, rf.fut.85)
            
            writeRaster(
              Future.85.mean,
              paste0("./outputs/", especie, "_", "Future.85.mean.tif"), 
              overwrite = T)
            
            # Binary Transformation #
            Future.85.bin <- terra::mean(subset(Future.85.all,
                                                 sel2[, "ID"]) >= sel2$threshold
            )
            
            # Future.85.bin = c(Future.85.bin == terra::maxValue(Future.85.bin))
  
            writeRaster(
              Future.85.bin,
              paste0("./outputs/", especie, "_", "Future.85.bin.tif"), 
              overwrite = T)
            
            rm(Future.85.mean, Future.85.bin, Future.85.all)
            gc(reset = T)
          # }
          
          
          rm(sel2)
          gc(reset = TRUE, full = TRUE)
          # unlink(paste0("./temp_output/" ,especie),recursive = T, force = T)
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
# base::quit(save = "yes")
