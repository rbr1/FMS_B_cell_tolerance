Get_polygons<-function(in_dir){
  p_file<-concat(c(in_dir,"MFI_gates.txt"))
  pp <- as.matrix(read.csv(p_file, head=TRUE, sep="\t"))
  polygon_header<-as.character(pp[,1])
  x_poly<-as.numeric(pp[,2])
  y_poly<-as.numeric(pp[,3])
  polygon<-as.matrix(cbind(x_poly, y_poly))
  dimnames(polygon)[[2]]<-c("x","y")
  return(list(polygon_header, polygon))
}
###############################################################
concat = function(v) {
  res = ""
  for (i in 1:length(v)){
    res = paste(res,v[i],sep="")
  }
  res
}
###############################################################
add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha)) }
###############################################################
Get_files_data<-function(in_dir, file_table){
  p <- as.matrix(read.csv(file_table, head=TRUE, sep="\t"))
  p=p[which(p[,"Experiment.type"] %in% c("FMS data")),]
  locations = gsub(".//",in_dir,p[,"Location"], fixed = T)
  samples= p[,"Day_sample"]
  names(locations) = samples
  patient = p[,"Patient"]
  return(cbind(locations, patient))
}
Get_files_data_all<-function(in_dir, file_table){
  p <- as.matrix(read.csv(file_table, head=TRUE, sep="\t"))
  p=p[which(p[,"Experiment.type"] %in% c("FMS data", "Proliferation check")),]
  locations = gsub(".//",in_dir,p[,"Location"], fixed = T)
  samples= p[,"Day_sample"]
  names(locations) = samples
  patient = p[,"Patient"]
  return(cbind(locations, patient))
}
in_dir = "~/Google_Drive/Projects/Fibromyalgia/Cell_culture/FMS/"
file_table = concat(c(in_dir, "Files_information.txt"))

locations_patient = Get_files_data_all(in_dir, file_table)
locations = locations_patient[,1]
patient = locations_patient[,2]
patients = sort(unique(patient))

###############################################################
library(flowCore)

list<-Get_polygons(in_dir)
polygon_header<-list[[1]]
polygon<-list[[2]]

cex = 0.9
cols<-NULL
alpha<-0.5
cols<-c(cols, rgb(rbind(c(1,0,0)), alpha=alpha))
cols<-c(cols, rgb(rbind(c(0.4,0.9,0.4)), alpha=0.8))
cols<-c(cols, rgb(rbind(c(0,0.3,1)), alpha=alpha/2))
cols<-c(cols, rgb(rbind(c(0.5,0.5,0.5)), alpha=alpha))
cols<-c(cols, rgb(rbind(c(0,1,1)), alpha=alpha))

library(mgcv)

###### Gating
Get_lymphocytes<-function(in_dir, expr_matrix, sample){
  dat<-expr_matrix
  main = sample
  out_file_tmp = concat(c(in_dir,"FCS_matrices/tmp.txt"))
  
  cex = 0.9
  plot(dat[,"FSC-A"], dat[,"SSC-A"], pch=20, col=cols[4],xlab="FSC-H", ylab="SSC-A", main=main,log="",cex.axis = cex, cex.lab = cex, cex = 0.3, cex.main = cex)
  # coords <- locator() ### click points on graph then click escape 
  # write.table(cbind(coords$x, coords$y), file = out_file_tmp, append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = FALSE,col.names = TRUE, qmethod = c("escape", "double"),fileEncoding = "")
  
  poly1=polygon[which(polygon_header== "Lymph1"),]
  ioall <- in.out(poly1,dat[,c("FSC-A", "SSC-A")])
  points(dat[ioall,c("FSC-A", "SSC-A")], pch=16, col="white") 
  points(dat[ioall,c("FSC-A", "SSC-A")], pch=16, col=cols[1], cex = 0.3) 
  polygon(poly1,lwd = 2) 
  perc_captured_1<-signif(length(which(ioall==TRUE))*100/length(ioall), digits=3)
  text(mean(poly1[,1]), mean(poly1[,2]), labels = concat(c(perc_captured_1 ," %\nlymphocytes" )),offset =0.5, pos=1,cex = cex, font = 4)
  
  cells<-dat[ioall,]
  return(cells)
}

### filter lymphocytes
fileout=concat(c(in_dir, "FCS_matrices/Filtering_step1_lymphocytes.jpeg"))
w = 1800
jpeg(file=fileout, height=w*14, width=w*12 , res=600)
par(mfrow= c(17,12), mar = c(4, 4, 3, 0.1))

lymphocytes = NULL
all_cells_uncompensated = NULL
all_cells_unnorm = NULL

### new: 
#FMHP020
#FMHP028
#FMSP066
#FMSP067

which(names(locations)=="Day_5_FMHP020_A")

for(ind in c(1:length(locations))){
  sample = names(locations)[ind]
  fcs_data <- read.FCS(locations[ind])

  comp_matrix <- keyword(fcs_data)$SPILL
  expr_matrix_uncompensated <- exprs(fcs_data)

  params <- pData(parameters(fcs_data))  # Get the parameter metadata
  
  marker_names <- ifelse(!is.na(params$desc), params$desc, params$name)
  cbind(marker_names, colnames(expr_matrix_uncompensated))
  colnames(expr_matrix_uncompensated) <- marker_names
  rownames(expr_matrix_uncompensated) = paste("c",c(1:length(expr_matrix_uncompensated[,1])), sep = "")
  
  all_cells_uncompensated = c(all_cells_uncompensated, list(expr_matrix_uncompensated))

  compensated_data <- compensate(fcs_data, comp_matrix)
  
  # Define the logicle transformation
  lgcl <- logicleTransform(w = 0.5, t = 262144, m = 4.5, a = 0)
  
  expr_matrix <- exprs(compensated_data)
  params <- pData(parameters(compensated_data))  # Get the parameter metadata
  
  marker_names <- ifelse(!is.na(params$desc), params$desc, params$name)
  cbind(marker_names, colnames(expr_matrix))
  colnames(expr_matrix) <- marker_names
  rownames(expr_matrix) = paste("c",c(1:length(expr_matrix[,1])), sep = "")
  
  all_cells_unnorm = c(all_cells_unnorm, list(expr_matrix))
 
  print (concat(c(ind,"/", length(locations)," done")))

  #transform_columns = c("CD21","IgD","CD3/CD10","CD38","CD19","IgM","CD24","Proliferation Dye" ,"CD27","CD22","CD14")
  transform_columns = c("CD21","IgD","CD3/CD10","CD38","CD19","IgM","CD24","CD27","CD22","CD14")
  for(i in c(1:length(transform_columns))){
    expr_matrix[,transform_columns[i]] = lgcl(expr_matrix[,transform_columns[i]])
  }
  ## get lymphocytes and singlets
  cells = Get_lymphocytes(in_dir, expr_matrix, sample)
  lymphocytes = c(lymphocytes, list(cells))
  print (concat(c(ind,"/", length(locations)," done")))
}
dev.off()   
names(lymphocytes) = names(locations)
names(all_cells_uncompensated) = names(locations)
names(all_cells_unnorm)= names(locations)
saveRDS(file = concat(c(in_dir, "FCS_matrices/Filtering_step1_lymphocytes.rds")), lymphocytes)
saveRDS(file = concat(c(in_dir, "FCS_matrices/Filtering_step1_all_cells_uncompensated.rds")), all_cells_uncompensated)
saveRDS(file = concat(c(in_dir, "FCS_matrices/Filtering_step1_all_cells_untransformed.rds")), all_cells_unnorm)

### filter single cells
lymphocytes = readRDS(file = concat(c(in_dir, "FCS_matrices/Filtering_step1_lymphocytes.rds")))
names(lymphocytes) = names(locations)

names(lymphocytes)[which(names(lymphocytes) %in% rownames(thresholds)==F)]

Get_singlets<-function(in_dir, mat, sample){
  main = sample
  out_file_tmp = concat(c(in_dir,"FCS_matrices/tmp.txt"))

  cells = mat
  x = "SSC-A"
  y = "SSC-H"
  plot(cells[,x], cells[,y], pch=20, col=cols[4],xlab="SSC-A", ylab="SSC-H", main=main,log="",cex.axis = cex, cex.lab = cex, cex = 0.3, cex.main = cex)
  # coords <- locator() ### click points on graph then click escape 
  # write.table(cbind(coords$x, coords$y), file = out_file_tmp, append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = FALSE,col.names = TRUE, qmethod = c("escape", "double"),fileEncoding = "")
  poly1=polygon[which(polygon_header== "SSC-A_SSC-H"),]
  ioall <- in.out(poly1,cells[,c(x,y)])
  points(cells[ioall,c(x,y)], pch=16, col="white") 
  points(cells[ioall,c(x,y)], pch=16, col=cols[2], cex = 0.3) 
  polygon(poly1,lwd = 2) 
  perc_captured_1<-signif(length(which(ioall==TRUE))*100/length(ioall), digits=3)
  text(mean(poly1[,1])/2, mean(poly1[,2]), labels = concat(c(perc_captured_1 ," %\nsinglets")),offset =0.5, pos=1,cex = cex, font = 4)
  cells1<-cells[ioall,]
  
  return(cells1)
  
}

fileout=concat(c(in_dir, "FCS_matrices/Filtering_step1_singlets.jpeg"))
w = 1800
jpeg(file=fileout, height=w*14, width=w*12 , res=600)
par(mfrow= c(17,12), mar = c(4, 4, 3, 0.1))

singlets = NULL
for(ind in c(1:length(lymphocytes))){
  sample = concat(c(ind,"  ", names(locations)[ind]))
  mat = lymphocytes[[ind]]
  cells1 = Get_singlets(in_dir, mat, sample)
  singlets = c(singlets, list(cells1))
  print (concat(c(ind,"/", length(locations)," done")))
}
dev.off()
names(singlets) = names(lymphocytes)

saveRDS(file = concat(c(in_dir, "FCS_matrices/Filtering_step1_singlets.rds")), singlets)


### filter live cells
singlets = readRDS(file = concat(c(in_dir, "FCS_matrices/Filtering_step1_singlets.rds")))

Get_live<-function(in_dir, mat, sample, gate_name){
  main = sample
  out_file_tmp = concat(c(in_dir,"FCS_matrices/tmp.txt"))
  cells1 = mat
  x =  "Live/dead"
  y = "SSC-A" #### change to CD19?
  cells1[which(cells1[,x]<=0),x] = 0.1
  plot(cells1[,x], cells1[,y], pch=20, col=cols[4],xlab="Live/dead", ylab="SSC-A", main=main,log="x",cex.axis = cex, cex.lab = cex, cex = 0.3, cex.main = cex, xlim = c(1, max(cells1[,x])))
  # coords <- locator() ### click points on graph then click escape 
  # write.table(cbind(coords$x, coords$y), file = out_file_tmp, append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = FALSE,col.names = TRUE, qmethod = c("escape", "double"),fileEncoding = "")
  poly1=polygon[which(polygon_header== gate_name),]
  ioall <- in.out(poly1,cells1[,c(x,y)])
  points(cells1[ioall,c(x,y)], pch=16, col="white") 
  points(cells1[ioall,c(x,y)], pch=16, col=cols[3], cex = 0.3) 
  polygon(poly1,lwd = 2) 
  perc_captured_1<-signif(length(which(ioall==TRUE))*100/length(ioall), digits=3)
  text(mean(poly1[,1])/1.5, mean(poly1[,2]), labels = concat(c(perc_captured_1 ," %\nLive cells")),offset =0.5, pos=1,cex = cex, font = 4)
  live_cells<-cells1[ioall,]
  return(live_cells)
}
list<-Get_polygons(in_dir)
polygon_header<-list[[1]]
polygon<-list[[2]]

fileout=concat(c(in_dir, "FCS_matrices/Filtering_step1_live.jpeg"))
w = 1800
jpeg(file=fileout, height=w*14, width=w*12 , res=600)
par(mfrow= c(17,12), mar = c(4, 4, 3, 0.1))

live = NULL

#### new gates 
new_gates1 =c(1,3,6)
new_gates2 = c(25:35)
new_gates3 = c(61:93)
new_gates4 = c(94:144)
new_gates5 = c(49:60)

new_gates = c(new_gates1, new_gates2, new_gates3)

for(ind in c(1:length(singlets))){
  sample = concat(c(ind,"  ", names(locations)[ind]))
  mat = singlets[[ind]]
  gate_name = "Live/dead"
  if(ind %in% new_gates1){gate_name = "Live/dead1"}
  if(ind %in% new_gates2){gate_name = "Live/dead2"}
  if(ind %in% new_gates3){gate_name = "Live/dead2"}
  if(ind %in% new_gates4){gate_name = "Live/dead3"}
  live_cells = Get_live(in_dir, mat, sample, gate_name)
  live = c(live, list(live_cells))
  print (concat(c(ind,"/", length(locations)," done")))
}
dev.off()
names(live) = names(locations)

saveRDS(file = concat(c(in_dir, "FCS_matrices/Filtering_step1_live.rds")), live)

### output to table for proliferation modelling in python
for(ind in c(1:length(live))){
  sample = names(locations)[ind]
  mat = live[[ind]]
  out_file = concat(c(in_dir,"FCS_matrices/MFI_matrix_",sample,".txt"))
  write.table(mat, file = out_file, append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = FALSE,col.names = TRUE, qmethod = c("escape", "double"),fileEncoding = "")
  print (concat(c(ind,"/", length(locations)," done")))
}

### get all other gates
live = readRDS(file = concat(c(in_dir, "FCS_matrices/Filtering_step1_live.rds")))

Get_gates<-function(in_dir,mat){
  thresholds <- as.matrix(read.csv(concat(c(in_dir, "MFI_thresholds.txt")), head=TRUE, sep="\t"))
  samples= thresholds[,"Day_sample"]
  colnames(thresholds) = gsub(".","/",colnames(thresholds), fixed = T)
  cols_use = colnames(mat)[which(colnames(mat) %in% colnames(thresholds)==T)]
  thresholds = apply(thresholds[,cols_use], 2, as.numeric) 
  rownames(thresholds) = samples
  colnames(thresholds) = cols_use
  return(thresholds)
}
thresholds =Get_gates(in_dir, live[[1]])

######
library(ggplot2)
library(MASS)

Get_gate1<-function(in_dir, mat, sample, thresholds, ind){
  main = concat(c(ind ," ",sample))
  out_file_tmp = concat(c(in_dir,"FCS_matrices/tmp.txt"))
  cells1 = mat
  cells1[which(cells1<=1)] = 1
  x =  "CD19"
  y = "CD14"
  
  smoothScatter(cells1[,x], (cells1[,y]), log="",colramp = colorRampPalette(c("blue", "green", "yellow", "red")), 
                xlab=x, ylab=concat(c("",y,"")), main=main)
  
  x_threshold = thresholds[sample, x]
  y_threshold = thresholds[sample, y]
  segments(0, y_threshold, 1e8, y_threshold,  col = "red", lwd = 2)
  segments( x_threshold, 1e8, x_threshold, 0, col = "red", lwd = 2)
  
  w2 = which(cells1[,y] <=thresholds[sample, y])
  mat1 = mat[w2,]
  # coords <- locator() ### click points on graph then click escape 
  # write.table(cbind(coords$x, coords$y), file = out_file_tmp, append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = FALSE,col.names = TRUE, qmethod = c("escape", "double"),fileEncoding = "")
  
  return(mat1)
}

names(live)[which(names(live) %in% rownames(thresholds)==F)]

w = which(names(live) %in% rownames(thresholds))
live = live[w]

fileout=concat(c(in_dir, "FCS_matrices/Filtering_step1_Gate1.jpeg"))
w = 1800
jpeg(file=fileout, height=w*12, width=w*12 , res=600)
par(mfrow= c(17,12), mar = c(4, 4, 3, 0.1))

CD19_pos_CD3_10_neg = NULL
for(ind in c(1:length(live))){
  sample = names(live)[ind]
  mat = live[[ind]]
  mat1 = Get_gate1(in_dir, mat, sample, thresholds, ind)
  CD19_pos_CD3_10_neg = c(CD19_pos_CD3_10_neg, list(mat1))
  print (concat(c(ind,"/", length(live)," done")))
}
names(CD19_pos_CD3_10_neg) = names(live)
dev.off()

saveRDS(file = concat(c(in_dir, "FCS_matrices/Filtering_step1_Gate1.rds")), CD19_pos_CD3_10_neg)

########
CD19_pos_CD3_10_neg = readRDS(file = concat(c(in_dir, "FCS_matrices/Filtering_step1_Gate1.rds")))


Get_gate2<-function(in_dir, mat, sample, thresholds, ind){
  main = concat(c(ind ," ",sample))
  out_file_tmp = concat(c(in_dir,"FCS_matrices/tmp.txt"))
  cells1 = mat
  cells1[which(cells1<=1)] = 1
  x =  "CD19"
  y = "CD3/CD10"
  smoothScatter(cells1[,x], cells1[,y], log="",colramp = colorRampPalette(c("blue", "green", "yellow", "red")), 
                xlab=x, ylab=y, main=main)
  
  x_threshold = thresholds[sample, x]
  y_threshold = thresholds[sample, y]
  segments(0, y_threshold, 1e8, y_threshold,  col = "red", lwd = 2)
  segments( x_threshold, 1e8, x_threshold, 0, col = "red", lwd = 2)
  
  w1 = which(cells1[,x] >=thresholds[sample, x])
  w2 = which(cells1[,y] <=thresholds[sample, y])
  mat1 = mat[intersect(w1,w2),]
  # coords <- locator() ### click points on graph then click escape 
  # write.table(cbind(coords$x, coords$y), file = out_file_tmp, append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = FALSE,col.names = TRUE, qmethod = c("escape", "double"),fileEncoding = "")
  
  return(mat1)
}


fileout=concat(c(in_dir, "FCS_matrices/Filtering_step1_Gate2.jpeg"))
w = 1800
jpeg(file=fileout, height=w*14, width=w*12 , res=600)
par(mfrow= c(17,12), mar = c(4, 4, 3, 0.1))

gate2_cells = NULL
for(ind in c(1:length(CD19_pos_CD3_10_neg))){
  sample = names(CD19_pos_CD3_10_neg)[ind]
  mat = CD19_pos_CD3_10_neg[[ind]]
  mat1 = Get_gate2(in_dir, mat, sample, thresholds, ind)
  gate2_cells = c(gate2_cells, list(mat1))
  print (concat(c(ind,"/", length(locations)," done")))
}
names(gate2_cells) = names(live)[c(1:length(gate2_cells))]
dev.off()

saveRDS(file = concat(c(in_dir, "FCS_matrices/Filtering_step1_Gate2.rds")), gate2_cells)

##########
gate2_cells = readRDS(file = concat(c(in_dir, "FCS_matrices/Filtering_step1_Gate2.rds")))

Get_gate3<-function(in_dir, mat, sample, thresholds, ind){
  main = concat(c(ind ," ",sample))
  out_file_tmp = concat(c(in_dir,"FCS_matrices/tmp.txt"))
  cells1 = mat
  cells1[which(cells1<=1)] = 1
  x =  "IgD"
  y = "CD27"
  smoothScatter(cells1[,x], (cells1[,y]), log="",colramp = colorRampPalette(c("blue", "green", "yellow", "red")), 
                xlab=x, ylab=y, main=main)
  
  x_threshold = thresholds[sample, x]
  y_threshold = thresholds[sample, y]
  segments(0, y_threshold, 1e8, y_threshold,  col = "red", lwd = 2)
  segments( x_threshold, 1e8, x_threshold, 0, col = "red", lwd = 2)
  
  w1 = which(cells1[,x] >=thresholds[sample, x])
  w2 = which(cells1[,y] <=thresholds[sample, y])
  mat1 = mat[intersect(w1,w2),]
  # coords <- locator() ### click points on graph then click escape 
  # write.table(cbind(coords$x, coords$y), file = out_file_tmp, append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = FALSE,col.names = TRUE, qmethod = c("escape", "double"),fileEncoding = "")
  
  return(mat1)
}

fileout=concat(c(in_dir, "FCS_matrices/Filtering_step1_Gate3.jpeg"))
w = 1800
jpeg(file=fileout, height=w*12, width=w*12 , res=600)
par(mfrow= c(17,12), mar = c(4, 4, 3, 0.1))

gate3_cells = NULL
for(ind in c(1:length(gate2_cells))){
  sample = names(gate2_cells)[ind]
  mat = gate2_cells[[ind]]
  mat = CD19_pos_CD3_10_neg[[ind]]
  mat1 = Get_gate3(in_dir, mat, sample, thresholds, ind)
  gate3_cells = c(gate3_cells, list(mat1))
  print (concat(c(ind,"/", length(locations)," done")))
}
names(gate3_cells) = names(live)
dev.off()



##########
gate2_cells = readRDS(file = concat(c(in_dir, "FCS_matrices/Filtering_step1_Gate2.rds")))

Get_gate4<-function(in_dir, mat, sample, thresholds, ind){
  main = concat(c(ind ," ",sample))
  out_file_tmp = concat(c(in_dir,"FCS_matrices/tmp.txt"))
  cells1 = mat
  cells1[which(cells1<=1)] = 1
  x =  "IgM"
  y = "CD38"
  smoothScatter(cells1[,x], cells1[,y], log="",colramp = colorRampPalette(c("blue", "green", "yellow", "red")), 
                xlab=x, ylab=y, main=main)
  
  x_threshold = thresholds[sample, x]
  y_threshold = thresholds[sample, y]
  segments(0, y_threshold, 1e8, y_threshold,  col = "red", lwd = 2)
  segments( x_threshold, 1e8, x_threshold, 0, col = "red", lwd = 2)
  
  w1 = which(cells1[,x] >=thresholds[sample, x])
  w2 = which(cells1[,y] <=thresholds[sample, y])
  mat1 = mat[intersect(w1,w2),]
  # coords <- locator() ### click points on graph then click escape 
  # write.table(cbind(coords$x, coords$y), file = out_file_tmp, append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = FALSE,col.names = TRUE, qmethod = c("escape", "double"),fileEncoding = "")
  
  return(mat1)
}

fileout=concat(c(in_dir, "FCS_matrices/Filtering_step1_Gate4.jpeg"))
w = 1800
jpeg(file=fileout, height=w*12, width=w*12 , res=600)
par(mfrow= c(17,12), mar = c(4, 4, 3, 0.1))

gate4_cells = NULL
for(ind in c(1:length(gate2_cells))){
  sample = names(gate2_cells)[ind]
  mat = gate2_cells[[ind]]
  mat = CD19_pos_CD3_10_neg[[ind]]
  mat1 = Get_gate4(in_dir, mat, sample, thresholds, ind)
  gate4_cells = c(gate4_cells, list(mat1))
  print (concat(c(ind,"/", length(locations)," done")))
}
names(gate4_cells) = names(live)
dev.off()

#########

Get_gate5<-function(in_dir, mat, sample, thresholds, ind){
  main = concat(c(ind ," ",sample))
  out_file_tmp = concat(c(in_dir,"FCS_matrices/tmp.txt"))
  cells1 = mat
  cells1[which(cells1<=1)] = 1
  x =  "CD21"
  y = "CD22"
  smoothScatter(cells1[,x], (cells1[,y]), log="",colramp = colorRampPalette(c("blue", "green", "yellow", "red")), 
                xlab=x, ylab=y, main=main)
  
  x_threshold = thresholds[sample, x]
  y_threshold = thresholds[sample, y]
  segments(0, y_threshold, 1e8, y_threshold,  col = "red", lwd = 2)
  segments( x_threshold, 1e8, x_threshold, 0, col = "red", lwd = 2)
  
  w1 = which(cells1[,x] >=thresholds[sample, x])
  w2 = which(cells1[,y] <=thresholds[sample, y])
  mat1 = mat[intersect(w1,w2),]
  # coords <- locator() ### click points on graph then click escape 
  # write.table(cbind(coords$x, coords$y), file = out_file_tmp, append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = FALSE,col.names = TRUE, qmethod = c("escape", "double"),fileEncoding = "")
  
  return(mat1)
}

fileout=concat(c(in_dir, "FCS_matrices/Filtering_step1_Gate5.jpeg"))
w = 1800
jpeg(file=fileout, height=w*12, width=w*12 , res=600)
par(mfrow= c(17,12), mar = c(4, 4, 3, 0.1))

gate5_cells = NULL
for(ind in c(1:length(gate2_cells))){
  sample = names(gate2_cells)[ind]
  mat = gate2_cells[[ind]]
  mat1 = Get_gate5(in_dir, mat, sample, thresholds, ind)
  gate5_cells = c(gate5_cells, list(mat1))
  print (concat(c(ind,"/", length(locations)," done")))
}
names(gate5_cells) = names(live)
dev.off()

#########

Get_gate6<-function(in_dir, mat, sample, thresholds, ind){
  main = concat(c(ind ," ",sample))
  out_file_tmp = concat(c(in_dir,"FCS_matrices/tmp.txt"))
  cells1 = mat
  cells1[which(cells1<=1)] = 1
  x =  "CD24"
  y = "CD27"
  smoothScatter(cells1[,x], (cells1[,y]), log="",colramp = colorRampPalette(c("blue", "green", "yellow", "red")), 
                xlab=x, ylab=y, main=main)
  
  x_threshold = thresholds[sample, x]
  y_threshold = thresholds[sample, y]
  segments(0, y_threshold, 1e8, y_threshold,  col = "red", lwd = 2)
  segments( x_threshold, 1e8, x_threshold, 0, col = "red", lwd = 2)
  
  w1 = which(cells1[,x] >=thresholds[sample, x])
  w2 = which(cells1[,y] <=thresholds[sample, y])
  mat1 = mat[intersect(w1,w2),]
  # coords <- locator() ### click points on graph then click escape 
  # write.table(cbind(coords$x, coords$y), file = out_file_tmp, append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = FALSE,col.names = TRUE, qmethod = c("escape", "double"),fileEncoding = "")
  
  return(mat1)
}

fileout=concat(c(in_dir, "FCS_matrices/Filtering_step1_Gate6.jpeg"))
w = 1800
jpeg(file=fileout, height=w*12, width=w*12 , res=600)
par(mfrow= c(17,12), mar = c(4, 4, 3, 0.1))

gate5_cells = NULL
for(ind in c(1:length(gate2_cells))){
  sample = names(gate2_cells)[ind]
  mat = gate2_cells[[ind]]
  mat1 = Get_gate6(in_dir, mat, sample, thresholds, ind)
  gate5_cells = c(gate5_cells, list(mat1))
  print (concat(c(ind,"/", length(locations)," done")))
}
names(gate5_cells) = names(live)
dev.off()

Get_gate7<-function(in_dir, mat, sample, thresholds, ind){
  main = concat(c(ind ," ",sample))
  out_file_tmp = concat(c(in_dir,"FCS_matrices/tmp.txt"))
  cells1 = mat
  cells1[which(cells1<=1)] = 1
  x =  "IgM"
  y = "CD38"
  smoothScatter(cells1[,x], (cells1[,y]), log="",colramp = colorRampPalette(c("blue", "green", "yellow", "red")), 
                xlab=x, ylab=y, main=main)
  
  x_threshold = thresholds[sample, x]
  y_threshold = thresholds[sample, y]
  segments(0, y_threshold, 1e8, y_threshold,  col = "red", lwd = 2)
  segments( x_threshold, 1e8, x_threshold, 0, col = "red", lwd = 2)
  
  w1 = which(cells1[,x] >=thresholds[sample, x])
  w2 = which(cells1[,y] <=thresholds[sample, y])
  mat1 = mat[intersect(w1,w2),]
  # coords <- locator() ### click points on graph then click escape 
  # write.table(cbind(coords$x, coords$y), file = out_file_tmp, append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = FALSE,col.names = TRUE, qmethod = c("escape", "double"),fileEncoding = "")
  
  return(mat1)
}

fileout=concat(c(in_dir, "FCS_matrices/Filtering_step1_Gate7.jpeg"))
w = 1800
jpeg(file=fileout, height=w*12, width=w*12 , res=600)
par(mfrow= c(17,12), mar = c(4, 4, 3, 0.1))

gate5_cells = NULL
for(ind in c(1:length(gate2_cells))){
  sample = names(gate2_cells)[ind]
  mat = gate2_cells[[ind]]
  mat1 = Get_gate7(in_dir, mat, sample, thresholds, ind)
  gate5_cells = c(gate5_cells, list(mat1))
  print (concat(c(ind,"/", length(locations)," done")))
}
names(gate5_cells) = names(live)
dev.off()

Get_gate8<-function(in_dir, mat, sample, thresholds, ind){
  main = concat(c(ind ," ",sample))
  out_file_tmp = concat(c(in_dir,"FCS_matrices/tmp.txt"))
  cells1 = mat
  cells1[which(cells1<=1)] = 1
  x =  "IgM"
  y = "IgD"
  smoothScatter(cells1[,x], (cells1[,y]), log="",colramp = colorRampPalette(c("blue", "green", "yellow", "red")), 
                xlab=x, ylab=y, main=main)
  
  x_threshold = thresholds[sample, x]
  y_threshold = thresholds[sample, y]
  segments(0, y_threshold, 1e8, y_threshold,  col = "red", lwd = 2)
  segments( x_threshold, 1e8, x_threshold, 0, col = "red", lwd = 2)
  
  w1 = which(cells1[,x] >=thresholds[sample, x])
  w2 = which(cells1[,y] <=thresholds[sample, y])
  mat1 = mat[intersect(w1,w2),]
  # coords <- locator() ### click points on graph then click escape 
  # write.table(cbind(coords$x, coords$y), file = out_file_tmp, append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = FALSE,col.names = TRUE, qmethod = c("escape", "double"),fileEncoding = "")
  
  return(mat1)
}

fileout=concat(c(in_dir, "FCS_matrices/Filtering_step1_Gate8.jpeg"))
w = 1800
jpeg(file=fileout, height=w*12, width=w*12 , res=600)
par(mfrow= c(17,12), mar = c(4, 4, 3, 0.1))

gate5_cells = NULL
for(ind in c(1:length(gate2_cells))){
  sample = names(gate2_cells)[ind]
  mat = gate2_cells[[ind]]
  mat1 = Get_gate8(in_dir, mat, sample, thresholds, ind)
  gate5_cells = c(gate5_cells, list(mat1))
  print (concat(c(ind,"/", length(locations)," done")))
}
names(gate5_cells) = names(live)
dev.off()
############### plot results

Combine_all_cell_frequencies_subsetted<-function(in_dir){
  live = readRDS(file = concat(c(in_dir, "FCS_matrices/Filtering_step1_live.rds")))
  singlets = readRDS(file = concat(c(in_dir, "FCS_matrices/Filtering_step1_singlets.rds")))
  
  Get_gates<-function(in_dir,mat){
    thresholds <- as.matrix(read.csv(concat(c(in_dir, "MFI_thresholds.txt")), head=TRUE, sep="\t"))
    samples= thresholds[,"Day_sample"]
    colnames(thresholds) = gsub(".","/",colnames(thresholds), fixed = T)
    cols_use = colnames(mat)[which(colnames(mat) %in% colnames(thresholds)==T)]
    thresholds = apply(thresholds[,cols_use], 2, as.numeric) 
    rownames(thresholds) = samples
    colnames(thresholds) = cols_use
    return(thresholds)
  }
  thresholds =Get_gates(in_dir, live[[1]])
  
  markers = setdiff(colnames(live[[1]]), c("FSC-A","FSC-H","FSC-W" ,"SSC-A","SSC-H","SSC-W" ,"Proliferation Dye" , "Live/dead" ,"Time" ))
  markers = c( "CD3/CD10", "CD14"   ,   "CD19"  ,  "IgD"  ,  "IgM"   ,   "CD27"   , "CD21"   , "CD38"    , "CD24"    ,   "CD22")
  
  all_populations = NULL
  for(i in c(1:length(markers))){
    if(i == 1){
      all_populations = c( concat(c(markers[i],"+")),  concat(c(markers[i],"-")))
    }else{
      new_pops = c(paste(all_populations, concat(c(markers[i],"+"))), paste(all_populations, concat(c(markers[i],"-"))))
      all_populations = new_pops
    }
    print(i)
    print(length(all_populations))
  }
  length(unique(all_populations))
  
  run = grep("Prolif_Check", names(live), invert = T)
  
  
  binary_matrix_list = NULL  
  percentage_live = NULL
  
  for(ind in run){
    sample = names(live)[ind]
    mat = live[[ind]]
    mat1 = mat[,markers]*0
    for(i in c(1:length(markers))){
      threshold = thresholds[sample, markers[i]]
      x = mat[,markers[i]]
      mat1[which(x>threshold),markers[i]] = 1
    }
    mat2 = mat1
    for(i in c(1:length(markers))){
      mat2[which(mat1[, markers[i]]>0), markers[i]] = concat(c(markers[i],"+"))
      mat2[which(mat1[, markers[i]]==0), markers[i]] = concat(c(markers[i],"-"))
    }
    binary_matrix_list = c(binary_matrix_list, list(mat2))
    print(ind) 
    live_perc = length(live[[sample]][,1])*100/length(singlets[[sample]][,1])
    percentage_live = c(percentage_live, live_perc)
    
  }
  names(binary_matrix_list) = names(live)[run]
  names(percentage_live)  = names(live)[run]
  
  level_1_names = NULL
  level_2a_names = NULL
  level_2b_names = NULL
  level_2c_names = NULL
  CD21 = NULL
  CD22 = NULL
  CD24 = NULL
  seq = ceiling(seq(from = 1, to = length(binary_matrix_list), length = 10))
  names(binary_matrix_list)[seq]
  for(ind in seq){
    print(ind)
    level_1_names = sort(unique(c(level_1_names, apply(unique(binary_matrix_list[[ind]][,c("CD3/CD10", "CD14"   ,   "CD19"  )]), 1, paste, collapse = " "))))
    level_2a_names = sort(unique(c(level_2a_names, apply(unique(binary_matrix_list[[ind]][,c("IgD"  ,  "CD27"  )]), 1, paste, collapse = " "))))
    level_2b_names = sort(unique(c(level_2b_names, apply(unique(binary_matrix_list[[ind]][,c("IgD"  ,  "CD27","IgM","CD38"  )]), 1, paste, collapse = " "))))
    level_2c_names = sort(unique(c(level_2c_names, apply(unique(binary_matrix_list[[ind]][,c("IgD"  ,  "CD27","IgM","CD38" ,"CD24" )]), 1, paste, collapse = " "))))
    
    CD21 = sort(unique(c(CD21, binary_matrix_list[[ind]][,"CD21"])))
    CD22 = sort(unique(c(CD22, binary_matrix_list[[ind]][,"CD22"])))
    CD24 = sort(unique(c(CD24, binary_matrix_list[[ind]][,"CD24"])))
    
  }
  
  level_2d_names = c(level_2b_names[grep("IgD- CD27+", level_2b_names, fixed = T)], level_2b_names[grep("IgD+ CD27+", level_2b_names, fixed = T)])
  level_2_names = level_2a_names
  level_3a_names = c(paste(level_2_names, "CD21+"), paste(level_2_names, "CD21-"))
  level_3b_names = c(paste(level_2_names, "CD22+"), paste(level_2_names, "CD22-"))

  level_1a_names = c(paste(level_2b_names, "CD21+"), paste(level_2b_names, "CD21-"))
  level_1b_names = c(paste(level_2b_names, "CD22+"), paste(level_2b_names, "CD22-"))
  
  p <- as.matrix(read.csv("~/Google_Drive/Projects/Fibromyalgia/Cell_culture/FMS/B_cell_phenotypes.txt", head=T, sep="\t"))
  cell_type =p[,1]
  l1 = apply(p[,c(2,3,4,5 )], 1, paste, collapse = " ")
  names(cell_type) = l1
  level_2b_names = sort(unique(cell_type[level_2b_names]))
  

  level_2d_names= sort(apply(expand.grid(level_2b_names, CD21), 1, paste, collapse = " : "))
  level_3a_names= sort(apply(expand.grid(level_2b_names, CD22), 1, paste, collapse = " : "))
  level_3b_names= sort(apply(expand.grid(level_2b_names, CD24), 1, paste, collapse = " : "))
  
  counts_level1 = matrix(data= 0, nrow = length(binary_matrix_list), ncol = length(level_1_names), dimnames = c(list(names(binary_matrix_list)), list(level_1_names)))
  counts_level1a = matrix(data= 0, nrow = length(binary_matrix_list), ncol = length(level_1a_names), dimnames = c(list(names(binary_matrix_list)), list(level_1a_names)))
  counts_level1b = matrix(data= 0, nrow = length(binary_matrix_list), ncol = length(level_1b_names), dimnames = c(list(names(binary_matrix_list)), list(level_1b_names)))
  counts_level2 = matrix(data= 0, nrow = length(binary_matrix_list), ncol = length(level_2_names), dimnames = c(list(names(binary_matrix_list)), list(level_2_names)))
  counts_level2b = matrix(data= 0, nrow = length(binary_matrix_list), ncol = length(level_2b_names), dimnames = c(list(names(binary_matrix_list)), list(level_2b_names)))
  counts_level2c = matrix(data= 0, nrow = length(binary_matrix_list), ncol = length(level_2c_names), dimnames = c(list(names(binary_matrix_list)), list(level_2c_names)))
  counts_level2d = matrix(data= 0, nrow = length(binary_matrix_list), ncol = length(level_2d_names), dimnames = c(list(names(binary_matrix_list)), list(level_2d_names)))
  counts_level3a = matrix(data= 0, nrow = length(binary_matrix_list), ncol = length(level_3a_names), dimnames = c(list(names(binary_matrix_list)), list(level_3a_names)))
  counts_level3b = matrix(data= 0, nrow = length(binary_matrix_list), ncol = length(level_3b_names), dimnames = c(list(names(binary_matrix_list)), list(level_3b_names)))
  lists_annotations = NULL
  
  for(ind in c(1:length(binary_matrix_list))){
    mat = binary_matrix_list[[ind]]
    l1 = apply(mat[,c("CD3/CD10", "CD14"   ,   "CD19"  )], 1, paste, collapse = " ")
    t = table(l1)
    counts_level1[names(binary_matrix_list)[ind], names(t)] = t
    
    w = which(l1 == "CD3/CD10- CD14- CD19+")
    l2 = apply(mat[w,c("IgD"  ,  "CD27")], 1, paste, collapse = " ")
    t = table(l2)
    t = t[which(names(t) %in% colnames(counts_level2))]
    counts_level2[names(binary_matrix_list)[ind], names(t)] = t
    
    l3b = apply(mat[w,c("IgD"  ,  "CD27","IgM","CD38"   )], 1, paste, collapse = " ")
    l3b = cell_type[l3b]
    l3a = apply(cbind(l3b, mat[w,c("CD22")]), 1, paste, collapse = " : ")
    t = table(l3a)
    t = t[which(names(t) %in% colnames(counts_level3a))]
    counts_level3a[names(binary_matrix_list)[ind], names(t)] = t
   
    l3b = apply(mat[w,c("IgD"  ,  "CD27","IgM","CD38"   )], 1, paste, collapse = " ")
    l3b = cell_type[l3b]
    l3a = apply(cbind(l3b, mat[w,c("CD24")]), 1, paste, collapse = " : ")
    t = table(l3a)
    t = t[which(names(t) %in% colnames(counts_level3b))]
    counts_level3b[names(binary_matrix_list)[ind], names(t)] = t
    
    l3b = apply(mat[w,c("IgD"  ,  "CD27","IgM","CD38"   )], 1, paste, collapse = " ")
    l3b = cell_type[l3b]
    t = table(l3b)
    t = t[which(names(t) %in% colnames(counts_level2b))]
    counts_level2b[names(binary_matrix_list)[ind], names(t)] = t
    lists_annotations = c(lists_annotations, list(l3b))
    
    l3b = apply(mat[w,c("IgD"  ,  "CD27", "IgM","CD38" ,"CD24"  )], 1, paste, collapse = " ")
    t = table(l3b)
    t = t[which(names(t) %in% colnames(counts_level2c))]
    counts_level2c[names(binary_matrix_list)[ind], names(t)] = t
    
    l3b = apply(mat[w,c("IgD"  ,  "CD27","IgM","CD38"   )], 1, paste, collapse = " ")
    l3b = cell_type[l3b]
    l3a = apply(cbind(l3b, mat[w,c("CD21")]), 1, paste, collapse = " : ")
    t = table(l3a)
    t = t[which(names(t) %in% colnames(counts_level2d))]
    counts_level2d[names(binary_matrix_list)[ind], names(t)] = t
    
    mat = binary_matrix_list[[ind]]
    l1 = apply(mat[w,c("IgD"  ,  "CD27", "IgM","CD38","CD21" )], 1, paste, collapse = " ")
    t = table(l1)
    counts_level1a[names(binary_matrix_list)[ind], names(t)] = t
    
    mat = binary_matrix_list[[ind]]
    l1 = apply(mat[w,c("IgD"  ,  "CD27", "IgM","CD38" ,"CD22" )], 1, paste, collapse = " ")
    t = table(l1)
    counts_level1b[names(binary_matrix_list)[ind], names(t)] = t
    
    print(ind)
  }
  
  names(lists_annotations) = names(binary_matrix_list)
  saveRDS(file = concat(c(in_dir, "FCS_matrices/Annotated_B_cells.rds")), lists_annotations)
  
  all_counts = c(list(counts_level1), list(counts_level2), list(counts_level3a), list(counts_level3b), list(counts_level1a), list(counts_level1b), list(counts_level2b), list(counts_level2c), list(counts_level2d))
  names(all_counts) = paste("Level", c("1","2", "3a", "3b","1a","1b","2b","2c","2d"))
  
  all_live = rowSums(counts_level1)
  all_counts_norm_within = all_counts
  all_counts_norm_live = all_counts
  for(i in c(1:length(all_counts))){
    mat = all_counts[[i]]
    mat1 = mat
    for(ind in c(1:length(mat[,1]))){
      mat[ind,] = mat[ind,]*100/sum(mat[ind,])
      mat1[ind,] = mat1[ind,]*100/all_live[ind]
    }
    all_counts_norm_within[[i]] = mat
    all_counts_norm_live[[i]] = mat1
    print(colSums(mat))
  }
  
  ##### get means per patient
  
  # Define a function to remove outliers using IQR
  remove_outliers <- function(x) {
    Q1 <- quantile(x, 0.35)
    Q3 <- quantile(x, 0.65)
    mean = median(x)
    IQR <- Q3 - Q1
    if(sd(x)!=0){
      x = x[which(x > (mean - 2 * IQR))]
      x = x[which(x < (mean + 2 * IQR))]
    }
    return(x)
  }
  
  samples = rownames(all_counts_norm_within[[1]])
  patient_sub = gsub("_A","", gsub("_B","", gsub("_C","", samples, fixed = T), fixed = T), fixed = T)
  patient_subs = sort(unique(patient_sub))
  
  per_patient_counts_norm_within = NULL
  per_patient_counts_norm_live = NULL
  per_patient_perc_live = NULL
  
  for(ind in c(1:length(all_counts_norm_within))){
    mat1 = all_counts_norm_within[[ind]]
    mat2 = all_counts_norm_live[[ind]]
    populations_use = colnames(mat1)
    all_counts_norm_within_per_patient = matrix(data= 0, nrow = length(patient_subs), ncol = length(populations_use), dimnames = c(list(patient_subs), list(populations_use)))
    all_counts_norm_live_per_patient = matrix(data= 0, nrow = length(patient_subs), ncol = length(populations_use), dimnames = c(list(patient_subs), list(populations_use)))
    for(p in c(1:length(patient_subs))){
      w = samples[which(patient_sub==patient_subs[p])]
      means = apply(mat1[w,], 2, function(x){mean(remove_outliers(x))})
      all_counts_norm_within_per_patient[patient_subs[p],]  = means
      means = apply(mat2[w,], 2, function(x){mean(remove_outliers(x))})
      all_counts_norm_live_per_patient[patient_subs[p],]  = means
      if(ind ==1){
        mean_perc_live = mean(remove_outliers(percentage_live[w]))
        per_patient_perc_live = c(per_patient_perc_live, mean_perc_live)
      }
    }
    per_patient_counts_norm_within = c(per_patient_counts_norm_within, list(all_counts_norm_within_per_patient))
    per_patient_counts_norm_live = c(per_patient_counts_norm_live, list(all_counts_norm_live_per_patient))
    print (ind)
  } 
  names(per_patient_counts_norm_within) = names(all_counts_norm_within)
  names(per_patient_counts_norm_live) = names(all_counts_norm_within)
  names(per_patient_perc_live) = patient_subs
  per_patient_counts_norm_live = c(per_patient_counts_norm_live, list(cbind(per_patient_perc_live)))
  names(per_patient_counts_norm_live)[length(names(per_patient_counts_norm_live))] = "% live of singlets"
  saveRDS(file = concat(c(in_dir, "FCS_matrices/Counts_per_patient_timepoint_within.rds")), per_patient_counts_norm_within)
  saveRDS(file = concat(c(in_dir, "FCS_matrices/Counts_per_patient_timepoint_live.rds")), per_patient_counts_norm_live)
}

Boxplot_custom<-function(groups, main, width_plot, colsx){
  factors = names(groups)
  max = max(c(unlist(groups), unlist(groups))*1.35)
  min = 0
  if(min(unlist(groups))<0){
    min = -1}
  b = (max-min)*0.034
  ylab = ""
  draw_signif_lines = TRUE
  y = max(c(unlist(groups), unlist(groups))*1)+b
  max_width = width_plot
  max_scale = min(c(max,100))
  range = max-min
  if(range>50){scale = c(-100:100)*20}
  if(range>200){scale = c(-100:100)*50}
  if(range<=50){scale = c(-100:100)*10}
  if(range<=30){scale = c(-100:100)*5}
  if(range <10){scale = c(-100:100)*2.5}
  if(range <5){scale = c(-100:100)*1}
  if(range <4){scale = c(-100:100)*0.5}
  if(range <1.5){scale = c(-100:1000)*0.2}
  if(range <0.5){scale = c(-100:100)*0.1}
  if(range <0.1){scale = c(-100:100)*0.01}
  if(range <0.01){scale = c(-100:100)*0.001}
  cex = 0.9
  Fun<-function(x){x}
  
  scale = scale[intersect(which(scale<= max_scale), which(scale>=min))]
  plot(c(0.5, max_width +0.5),c(min, max), pch=20, col="white",xlab="",ylab ="",cex=cex, cex.lab=cex+0.1,	cex.axis=cex,cex.main = cex, col.axis = "white",tck=0, mgp = c(2,0,0), main = main, axes = FALSE, ylim = c(min, max))
  mtext(side = 2, text = ylab, line = 2.8,cex= cex-0.1, las = 3, font = 1)
  mtext(side = 1, text = factors, line = 0.35,cex= cex-0.1,  at = c(1:length(factors)), las = 2, font = 1)
  segments(0.5,Fun(scale),length(groups)+0.5,Fun(scale),col = "grey",lwd = 1,lty = 3 )
  segments(0.5,Fun(scale),length(groups)+0.5,Fun(scale),col = "grey",lwd = 1,lty = 3 )
  mtext(side = 2, text = scale, line = 0.15,cex= cex-0.1,  at =Fun(scale), las = 2, font = 1)
  width = 0.38
  index = 1
  l = length(groups)
  l1 = length(groups[[1]])
  
  for(i in c(1:l)){
    points1=as.numeric(groups[[i]])
    box1<-c(as.numeric(quantile(points1))[3], as.numeric(quantile(points1, probs = c(0.1, 0.9))), as.numeric(quantile(points1))[c(2, 4)])	
    Draw_box_plot(box1,i,width,colsx[i],1, colsx1[i])
    points(rep(i, length(points1)),points1, pch =21, col=colsx[i],bg = colsx1[i], cex = 0.7)
  }
}

Draw_box_plot<-function(box,x,width,c,lwd,line_col){
  segments(x, box[2], x, box[3], col = line_col,lwd =lwd)
  segments(x-(width/2), box[2], x+(width/2), box[2], col = line_col,lwd =lwd)
  segments(x-(width/2), box[3], x+(width/2), box[3], col = line_col,lwd =lwd)
  rect(x-width, box[4], x+width, box[5], col = c,lwd =lwd, border = line_col)
  segments(x-width, box[1], x+width, box[1], col = line_col,lwd=2*lwd)}

Means_factor = function(factor, x){
  m = NULL
  for(i1 in c(1:length(levels(factor)))){
    x1 = x[which(factor==levels(factor)[i1])]
    x1 = x1[which(x1!=-1)]
    m = c(m, mean(x1))}
  return(m)}

Medians_factor = function(factor, x){
  m = NULL
  for(i1 in c(1:length(levels(factor)))){
    x1 = x[which(factor==levels(factor)[i1])]
    x1 = x1[which(x1!=-1)]
    m = c(m, median(x1))}
  return(m)}

Get_signficance_between_groups<-function(){
  
  ################################################################################################# clinical information
  file = "~/Google_Drive/Projects/Fibromyalgia/Data/BCR/Clinical_information.txt"
  p <- as.matrix(read.csv(file, head=TRUE, sep="\t"))
  
  sample_clin = p[,"SampleID"]
  age = as.numeric(p[,"Age"])
  gender = (p[,"Sex"])
  case_control = p[,"Disease"]
  ethnicity = p[,"Ethnicity"]
  Children.birthed = p[,"Children.birthed"]
  
  names(age) = sample_clin
  names(gender) = sample_clin
  names(case_control) = sample_clin
  names(ethnicity) = sample_clin
  names(Children.birthed) = sample_clin
  
  ##############
  per_patient_counts_norm_within =readRDS(file = concat(c(in_dir, "FCS_matrices/Counts_per_patient_timepoint_within.rds")))
  per_patient_counts_norm_live =readRDS(file = concat(c(in_dir, "FCS_matrices/Counts_per_patient_timepoint_live.rds")))
  
  per_patient_counts_norm_live = per_patient_counts_norm_live[which(names(per_patient_counts_norm_live)!="% live of singlets")]
  mat = per_patient_counts_norm_within[["Level 2b"]]
  mat = mat[,which(colnames(mat)!="Naive")]
  for(i in c(1:length(mat[,1]))){
    mat[i,] = mat[i,]*100/sum(mat[i,])
  }
  per_patient_counts_norm_within = c(per_patient_counts_norm_within, list(mat))
  names(per_patient_counts_norm_within)[length(per_patient_counts_norm_within)] = "Level 2b-non-naive"
  
  mat = per_patient_counts_norm_live[["Level 2b"]]
  mat = mat[,which(colnames(mat)!="Naive")]
  for(i in c(1:length(mat[,1]))){
    mat[i,] = mat[i,]*100/sum(mat[i,])
  }
  per_patient_counts_norm_live = c(per_patient_counts_norm_live, list(mat))
  names(per_patient_counts_norm_live)[length(per_patient_counts_norm_live)] = "Level 2b-non-naive"
  
  
  
  ### make list split by day and type
  list_proportions = NULL
  names_proportions = NULL
  for(ind in c(1:length(per_patient_counts_norm_within))){
    mat1 = per_patient_counts_norm_within[[ind]]
    mat2 = per_patient_counts_norm_live[[ind]]
    name = names(per_patient_counts_norm_within)[ind]
    
    #### for % within gate
    day0_mat1 = mat1[grep("Day_0", rownames(mat1)), ]
    day5_mat1 = mat1[grep("Day_0", rownames(mat1), invert = T), ]
    rownames(day0_mat1) = gsub("Day_0_", "", rownames(day0_mat1), fixed = T)
    rownames(day5_mat1) = gsub("Day_5_", "", rownames(day5_mat1), fixed = T)
    rownames(day5_mat1) = gsub("Day_6_", "", rownames(day5_mat1), fixed = T)
    
    w = which(rownames(day0_mat1) %in% sample_clin==F)
    rownames(day0_mat1)[w] = gsub("FMHP0","FMSH",rownames(day0_mat1)[w])
    w = which(rownames(day5_mat1) %in% sample_clin==F)
    rownames(day5_mat1)[w] = gsub("FMHP0","FMSH",rownames(day5_mat1)[w])
    
    day0_mat1 = day0_mat1[sort(rownames(day0_mat1)), ]
    day5_mat1 = day5_mat1[sort(rownames(day5_mat1)), ]
    rownames(day0_mat1)==rownames(day5_mat1)
    change_counts_for_all_populations_norm = (day5_mat1-day0_mat1)/(day5_mat1+day0_mat1)
    pseudocount <- 0.001
    change_counts_for_all_populations_norm1 = log10((day5_mat1+pseudocount)/(day0_mat1+pseudocount))
    
    
    list_proportions = c(list_proportions, list(day0_mat1), list(day5_mat1), list(change_counts_for_all_populations_norm), list(change_counts_for_all_populations_norm1))
    names_proportions = c(names_proportions, paste(concat(c(name," within group")), c("day 0", "day 5", "change", "FC")))
    
    #### for % of live
    day0_mat1 = mat2[grep("Day_0", rownames(mat2)), ]
    day5_mat1 = mat2[grep("Day_0", rownames(mat2), invert = T), ]
    rownames(day0_mat1) = gsub("Day_0_", "", rownames(day0_mat1), fixed = T)
    rownames(day5_mat1) = gsub("Day_5_", "", rownames(day5_mat1), fixed = T)
    rownames(day5_mat1) = gsub("Day_6_", "", rownames(day5_mat1), fixed = T)
    
    w = which(rownames(day0_mat1) %in% sample_clin==F)
    rownames(day0_mat1)[w] = gsub("FMHP0","FMSH",rownames(day0_mat1)[w])
    w = which(rownames(day5_mat1) %in% sample_clin==F)
    rownames(day5_mat1)[w] = gsub("FMHP0","FMSH",rownames(day5_mat1)[w])
    
    day0_mat1 = day0_mat1[sort(rownames(day0_mat1)), ]
    day5_mat1 = day5_mat1[sort(rownames(day5_mat1)), ]
    rownames(day0_mat1)==rownames(day5_mat1)
    change_counts_for_all_populations_norm = (day5_mat1-day0_mat1)/(day5_mat1+day0_mat1)
    pseudocount <- 0.001
    change_counts_for_all_populations_norm1 = log10((day5_mat1+pseudocount)/(day0_mat1+pseudocount))
    
    list_proportions = c(list_proportions, list(day0_mat1), list(day5_mat1), list(change_counts_for_all_populations_norm), list(change_counts_for_all_populations_norm1))
    names_proportions = c(names_proportions, paste(concat(c(name," of live")), c("day 0", "day 5", "change", "FC")))
    
  }
  names(list_proportions) = names_proportions
  
  per_patient_counts_norm_live =readRDS(file = concat(c(in_dir, "FCS_matrices/Counts_per_patient_timepoint_live.rds")))
  other = names(per_patient_counts_norm_live)[which(names(per_patient_counts_norm_live) %in% names(per_patient_counts_norm_within)==F)]
  mat1 = per_patient_counts_norm_live[[other]]
  
  day0_mat1 = cbind(mat1[grep("Day_0", rownames(mat1)), ])
  day5_mat1 = cbind(mat1[grep("Day_0", rownames(mat1), invert = T), ])
  rownames(day0_mat1) = gsub("Day_0_", "", rownames(day0_mat1), fixed = T)
  rownames(day5_mat1) = gsub("Day_5_", "", rownames(day5_mat1), fixed = T)
  rownames(day5_mat1) = gsub("Day_6_", "", rownames(day5_mat1), fixed = T)
  
  mat1 = cbind(day0_mat1, day5_mat1)
  mat2 = cbind((day5_mat1-day0_mat1)/(day5_mat1+day0_mat1))
  colnames(mat1) = c("Day 0","Day 5")
  colnames(mat2) = "proportional change"
  
  w = which(rownames(mat1) %in% sample_clin==F)
  rownames(mat1)[w] = gsub("FMHP0","FMSH",rownames(mat1)[w])
  w = which(rownames(mat2) %in% sample_clin==F)
  rownames(mat2)[w] = gsub("FMHP0","FMSH",rownames(mat2)[w])
  
  mat1 = mat1[sort(rownames(mat1)), ]

  list_proportions = c(list_proportions, list(mat1), list(mat2))
  names(list_proportions)[c(length(names(list_proportions))-1, length(names(list_proportions)))] = c(other, concat(c(other," change")))
  
  ### perc B cells
  mat1 = cbind(per_patient_counts_norm_within[["Level 1"]][,"CD3/CD10- CD14- CD19+"])
  
  day0_mat1 = cbind(mat1[grep("Day_0", rownames(mat1)), ])
  day5_mat1 = cbind(mat1[grep("Day_0", rownames(mat1), invert = T), ])
  rownames(day0_mat1) = gsub("Day_0_", "", rownames(day0_mat1), fixed = T)
  rownames(day5_mat1) = gsub("Day_5_", "", rownames(day5_mat1), fixed = T)
  rownames(day5_mat1) = gsub("Day_6_", "", rownames(day5_mat1), fixed = T)
  
  mat1 = cbind(day0_mat1, day5_mat1)
  mat2 = cbind((day5_mat1-day0_mat1)/(day5_mat1+day0_mat1))
  colnames(mat1) = c("Day 0","Day 5")
  colnames(mat2) = "proportional change"
  
  w = which(rownames(mat1) %in% sample_clin==F)
  rownames(mat1)[w] = gsub("FMHP0","FMSH",rownames(mat1)[w])
  w = which(rownames(mat2) %in% sample_clin==F)
  rownames(mat2)[w] = gsub("FMHP0","FMSH",rownames(mat2)[w])
  
  mat1 = mat1[sort(rownames(mat1)), ]
  mat2 = mat2[sort(rownames(mat2)), ]
  
  list_proportions = c(list_proportions, list(mat1), list(mat2))
  names(list_proportions)[c(length(names(list_proportions))-1, length(names(list_proportions)))] = c("% B cells of live singlets", concat(c(other," change")))
  saveRDS(file = concat(c(in_dir, "FCS_matrices/All_populations_boxplots_all.rds")), list_proportions)
  
  #####
  exp_groups = list_proportions
  
  ### 1. run tests between disease cohorts-max severity
  groups_ids = NULL
  Sources = sort(unique(case_control))
  for(i in c(1:length(Sources))){
    w = which(case_control==Sources[i])
    w = names(w)[which(names(w) %in% rownames(list_proportions[[1]]))]
    groups_ids = c(groups_ids, list((w)))
  }
  names(groups_ids) = Sources

  fileout1=concat(c(in_dir, "FCS_matrices/All_populations_boxplots_all.pdf"))
  w=4.5
  pdf(file=fileout1, height=w*1.4*6, width=w*2.5*6)
  par(mfrow= c(6,6), mar = c(35,4,4,0.5))
  summary_stats = NULL
  for(ind in c(1:length(exp_groups))){
    mat_stat = cbind(exp_groups[[ind]])
    use = which(apply(mat_stat, 2, function(x){length(unique(x))})>=10)
    mat_stat = cbind(mat_stat[,use])
    factor_sex = factor(gender[rownames(mat_stat)])
    factor_age = age[rownames(mat_stat)]
    factor_preg = as.numeric(Children.birthed[rownames(mat_stat)])
    factor_preg[which(factor_preg>0)] = 1
    factor = factor(case_control[rownames(mat_stat)])
    sources= levels(factor)
    
    #library(lmPerm)    ## non parametric manova
    #result <- lmp(mat_stat ~ factor +factor_age+factor_sex+factor_preg)
    #mancova_model <- manova(mat_stat1 ~ factor +factor_age+factor_sex+factor_preg)
    #summary_univariate <- summary.aov(mancova_model)
    if(dim(mat_stat)[2]>1){
      fit = manova(formula = mat_stat ~ factor +factor_age+factor_sex+factor_preg)
      p1 = summary.aov(fit)
      nam = gsub(" Response ","",names(p1))
      p_value = NULL
      means = NULL
      i1 = 0
      for(i in p1){
        i1 = i1+1
        p_value = c(p_value, i$'Pr(>F)'[1]) 
        if(length(mean)==0){means = Means_factor(factor, mat_stat[,i1])
        }else{means = rbind(means, Means_factor(factor, mat_stat[,i1]))}
      }
      p_value[which(is.na(p_value))] = 2
      names(p_value) = nam
      print(min(p_value))
      #print(length(which(p_value<0.05)))
      colnames(means) = paste("mean.group.", c(1:length(means[1,])))
      combined_p_value = cbind(p_value ,means)
      rownames(combined_p_value) = nam
      p.group = rep(names(exp_groups)[ind], length(nam))
      sublabel = gsub(concat(c(names(exp_groups)[ind],"..")), "", colnames(mat_stat), fixed = T)
      x = cbind(p.group, sublabel,combined_p_value)
      if(length(summary_stats)==0){summary_stats = x
      }else{summary_stats = rbind(summary_stats, x)}
    }else{
      fit = aov(formula = mat_stat ~ factor +factor_age+factor_sex+factor_preg)
      p1 = summary.aov(fit)
      p_value =  p1[[1]][1,4]
      print(min(p_value))
      
      combined_p_value = rbind(c(p_value, Means_factor(factor, mat_stat)))
      colnames(combined_p_value) = c("p-value", paste("mean.group.", c(1:length(levels(factor)))))
      rownames(combined_p_value) = "overall"
      sublabel = "overall"
      p.group = names(exp_groups)[ind]
      x = cbind(p.group, sublabel,combined_p_value)
      if(length(summary_stats)==0){summary_stats = x
      }else{summary_stats = rbind(summary_stats, x)}
    }

    
    mat = cbind(exp_groups[[ind]])
    p_value_all = rep(1, length(mat_stat[1,]))
    p_value_all[use] = p_value
    groups = NULL
    pvals =NULL
    include = NULL
    for(i in c(1:length(mat[1,]))){
      g1 = NULL
      lens= NULL
      for(p in c(1:length(groups_ids))){
        x = mat[groups_ids[[p]],i]
        x=x[which(is.na(x)==F)]
        g1 = c(g1, list( x))
        lens = c(lens, length(x))}
      names(g1) = names(groups_ids)
      if(min(lens)>=2){
        pvals = c(pvals, wilcox.test(g1[[1]], y=g1[[2]])$p.value)
      }else{pvals = c(pvals, 1)}
      groups = c(groups, list(g1))
      include = c(include, i)
      
    }
    names(groups)=colnames(mat)[include]
    names(pvals)=colnames(mat)[include]
    pvals = rep(1, length(colnames(exp_groups[[ind]])))
    names(pvals) = colnames(exp_groups[[ind]])
    pvals[names(p_value)] = p_value
    pvals[which(is.na(pvals))] = 1
    
    factors1 = names(groups_ids)
    factors = colnames(exp_groups[[ind]])
    if(is.null(factors)){factors = "overall"}
    main = names(exp_groups)[ind]
    max = max(c(unlist(groups), unlist(groups))*1.35)
    min = min(c(unlist(groups), unlist(groups))*1.15)
    min = 0
    if(min(unlist(groups))<0){
      min = -1
      max = 1.15}
    if(max(unlist(groups))>1.15){max = max(c(unlist(groups), unlist(groups))*1.35)}
    b = (max-min)*0.034
    ylab = ""
    draw_signif_lines = TRUE
    y = max(c(unlist(groups), unlist(groups))*1)+b
    max_width = 32
    max_scale = min(c(max,100))
    range = max-min
    if(range>50){scale = c(-100:100)*20}
    if(range<=50){scale = c(-100:100)*10}
    if(range<=30){scale = c(-100:100)*5}
    if(range <10){scale = c(-100:100)*2.5}
    if(range <5){scale = c(-100:100)*1}
    if(range <4){scale = c(-100:100)*0.5}
    if(range <1.5){scale = c(-100:1000)*0.2}
    if(range <0.5){scale = c(-100:100)*0.1}
    if(range <0.1){scale = c(-100:100)*0.01}
    if(range <0.01){scale = c(-100:100)*0.001}
    cex = 0.9
    Fun<-function(x){x}
    
    scale = scale[intersect(which(scale<= max_scale), which(scale>=min))]
    plot(c(1.5, max_width +0.5),c(min, max), pch=20, col="white",xlab="",ylab ="",cex=cex, cex.lab=cex+0.1,	cex.axis=cex,cex.main = cex, col.axis = "white",tck=0, mgp = c(2,0,0), main = main, axes = FALSE, ylim = c(min, max))
    mtext(side = 2, text = ylab, line = 2.8,cex= cex-0.1, las = 3, font = 1)
    mtext(side = 1, text = gsub("_","/",factors,fixed=T), line = 0.35,cex= cex-0.1,  at = c(1:length(factors)), las = 2, font = 1)
    segments(0.5,Fun(scale),length(groups)+0.5,Fun(scale),col = "grey",lwd = 1,lty = 3 )
    segments(0.5,Fun(scale),length(groups)+0.5,Fun(scale),col = "grey",lwd = 1,lty = 3 )
    mtext(side = 2, text = scale, line = 0.15,cex= cex-0.1,  at =Fun(scale), las = 2, font = 1)
    width = 0.17
    index = 1
    l = length(groups)
    l1 = length(groups[[1]])
    shift = c(1:l1)
    shift = (mean(shift)-shift)
    shift = shift*0.23/max(shift)
    
    library(RColorBrewer)
    cols1 =  add.alpha (brewer.pal(8, "Dark2")[c(3,6)], alpha = 0.95)
    cols =  add.alpha (cols1, alpha = 0.5)
    
    for(i in c(1:l)){
      for(i1 in c(1:l1)){
        points1=as.numeric(groups[[i]][[i1]])
        box1<-c(as.numeric(quantile(points1))[3], as.numeric(quantile(points1, probs = c(0.1, 0.9))), as.numeric(quantile(points1))[c(2, 4)])	
        Draw_box_plot(box1,i-shift[i1],width,cols[i1],1, cols1[i1])
        points(rep(i-shift[i1], length(points1)),points1, pch =21, col=cols[i1],bg = cols[i1], cex = 0.7)
      }}
    
    for(i in c(1:l)){	
      b = max*0.035
      signif_threshold = 0.05
      pval1= "NS"
      if(pvals[i]<signif_threshold){pval1 = "*"
        if(pvals[i]<signif_threshold/10){pval1 = "**"}
        y = max(unlist(groups[[i]]))
        y = y+3*b
        text(i, y+2*b, labels = pval1, cex = 1.3)
      }}
  }
  plot(c(0,1), c(0,1), pch = 21, col = "white", bg ="white", xlab = "", ylab = "", main = "", axes=F)
  
  legend("topleft", Sources, pch = 21,cex= 0.8, bty="n", pt.bg = cols, col = cols, pt.lwd = 2, text.font = 2)
  dev.off()

  out_file = concat(c(in_dir, "FCS_matrices/All_populations_boxplots_all.txt"))
  write.table(summary_stats, file = out_file, append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = TRUE,col.names = TRUE, qmethod = c("escape", "double"),fileEncoding = "")
  
}
Get_signficance_between_groups()

Combine_all_cell_frequencies_subsetted_proliferation<-function(in_dir){
  live = readRDS(file = concat(c(in_dir, "FCS_matrices/Filtering_step1_live.rds")))
  list_proliferation = readRDS(file= concat(c(in_dir, "FCS_matrices/Filtering_step1_proliferation_n_day5.rds")))
  
  Get_gates<-function(in_dir,mat){
    thresholds <- as.matrix(read.csv(concat(c(in_dir, "MFI_thresholds.txt")), head=TRUE, sep="\t"))
    samples= thresholds[,"Day_sample"]
    colnames(thresholds) = gsub(".","/",colnames(thresholds), fixed = T)
    cols_use = colnames(mat)[which(colnames(mat) %in% colnames(thresholds)==T)]
    thresholds = apply(thresholds[,cols_use], 2, as.numeric) 
    rownames(thresholds) = samples
    colnames(thresholds) = cols_use
    return(thresholds)
  }
  thresholds =Get_gates(in_dir, live[[1]])
  
  markers = setdiff(colnames(live[[1]]), c("FSC-A","FSC-H","FSC-W" ,"SSC-A","SSC-H","SSC-W" ,"Proliferation Dye" , "Live/dead" ,"Time" ))
  markers = c( "CD3/CD10", "CD14"   ,   "CD19"  ,  "IgD"  ,  "IgM"   ,   "CD27"   , "CD21"   , "CD38"    , "CD24"    ,   "CD22")

  p <- as.matrix(read.csv("~/Google_Drive/Projects/Fibromyalgia/Cell_culture/FMS/B_cell_phenotypes.txt", head=T, sep="\t"))
  cell_type =p[,1]
  l1 = apply(p[,c(2,3,4,5 )], 1, paste, collapse = " ")
  names(cell_type) = l1
  
  all_populations = NULL
  for(i in c(1:length(markers))){
    if(i == 1){
      all_populations = c( concat(c(markers[i],"+")),  concat(c(markers[i],"-")))
    }else{
      new_pops = c(paste(all_populations, concat(c(markers[i],"+"))), paste(all_populations, concat(c(markers[i],"-"))))
      all_populations = new_pops
    }
    print(i)
    print(length(all_populations))
  }
  length(unique(all_populations))
  
  run1 = grep("Prolif_Check", names(live), invert = T)
  run2 = grep("Day_0", names(live), invert = T)
  run = intersect(run1, run2)
  
  binary_matrix_list = NULL

  for(ind in run){
    sample = names(live)[ind]
    mat = live[[ind]]
    mat_sub = list_proliferation[[sample]]
    mat = mat[names(mat_sub),]
    mat1 = mat[,markers]*0
    for(i in c(1:length(markers))){
      threshold = thresholds[sample, markers[i]]
      x = mat[,markers[i]]
      mat1[which(x>threshold),markers[i]] = 1
    }
    mat2 = mat1
    for(i in c(1:length(markers))){
      mat2[which(mat1[, markers[i]]>0), markers[i]] = concat(c(markers[i],"+"))
      mat2[which(mat1[, markers[i]]==0), markers[i]] = concat(c(markers[i],"-"))
    }
    binary_matrix_list = c(binary_matrix_list, list(mat2))
    print(ind) 
  }
  names(binary_matrix_list) = names(live)[run]
  
  level_1_names = NULL
  level_2a_names = NULL
  level_2b_names = NULL
  level_2c_names = NULL
  seq = ceiling(seq(from = 1, to = length(binary_matrix_list), length = 20))
  names(binary_matrix_list)[seq]
  for(ind in seq){
    print(ind)
    level_1_names = sort(unique(c(level_1_names, apply(unique(binary_matrix_list[[ind]][,c("CD3/CD10", "CD14"   ,   "CD19"  )]), 1, paste, collapse = " "))))
    level_2a_names = sort(unique(c(level_2a_names, apply(unique(binary_matrix_list[[ind]][,c("IgD"  ,  "CD27"  )]), 1, paste, collapse = " "))))
    level_2b_names = sort(unique(c(level_2b_names, apply(unique(binary_matrix_list[[ind]][,c("IgD"  ,  "CD27","IgM","CD38"  )]), 1, paste, collapse = " "))))
    level_2c_names = sort(unique(c(level_2c_names, apply(unique(binary_matrix_list[[ind]][,c("CD21", "CD22", "CD24" )]), 1, paste, collapse = " "))))
  }
  level_2b_names = cell_type[level_2b_names]
  level_2c_names= apply(expand.grid(cell_type, level_2c_names), 1, paste, collapse = " : ")
    
  level_2d_names = c(level_2b_names[grep("IgD- CD27+", level_2b_names, fixed = T)], level_2b_names[grep("IgD+ CD27+", level_2b_names, fixed = T)])
  level_2_names = level_2a_names
  level_3a_names = c(paste(level_2_names, "CD21+"), paste(level_2_names, "CD21-"))
  level_3b_names = c(paste(level_2_names, "CD22+"), paste(level_2_names, "CD22-"))
  
  level_1a_names = c(paste(level_1_names, "CD21+"), paste(level_1_names, "CD21-"))
  level_1b_names = c(paste(level_1_names, "CD22+"), paste(level_1_names, "CD22-"))
  
  
  ### proliferation metrics
  DI_level1 = matrix(data= NA, nrow = length(binary_matrix_list), ncol = length(level_1_names), dimnames = c(list(names(binary_matrix_list)), list(level_1_names)))
  DI_level2 = matrix(data= NA, nrow = length(binary_matrix_list), ncol = length(level_2_names), dimnames = c(list(names(binary_matrix_list)), list(level_2_names)))
  DI_level2b = matrix(data= NA, nrow = length(binary_matrix_list), ncol = length(level_2b_names), dimnames = c(list(names(binary_matrix_list)), list(level_2b_names)))
  DI_level2c = matrix(data= NA, nrow = length(binary_matrix_list), ncol = length(level_2c_names), dimnames = c(list(names(binary_matrix_list)), list(level_2c_names)))
  DI_level2d = matrix(data= NA, nrow = length(binary_matrix_list), ncol = length(level_2d_names), dimnames = c(list(names(binary_matrix_list)), list(level_2d_names)))
  DI_level3a = matrix(data= NA, nrow = length(binary_matrix_list), ncol = length(level_3a_names), dimnames = c(list(names(binary_matrix_list)), list(level_3a_names)))
  DI_level3b = matrix(data= NA, nrow = length(binary_matrix_list), ncol = length(level_3b_names), dimnames = c(list(names(binary_matrix_list)), list(level_3b_names)))
  
  PI_level1 = matrix(data= NA, nrow = length(binary_matrix_list), ncol = length(level_1_names), dimnames = c(list(names(binary_matrix_list)), list(level_1_names)))
  PI_level2 = matrix(data= NA, nrow = length(binary_matrix_list), ncol = length(level_2_names), dimnames = c(list(names(binary_matrix_list)), list(level_2_names)))
  PI_level2b = matrix(data= NA, nrow = length(binary_matrix_list), ncol = length(level_2b_names), dimnames = c(list(names(binary_matrix_list)), list(level_2b_names)))
  PI_level2c = matrix(data= NA, nrow = length(binary_matrix_list), ncol = length(level_2c_names), dimnames = c(list(names(binary_matrix_list)), list(level_2c_names)))
  PI_level2d = matrix(data= NA, nrow = length(binary_matrix_list), ncol = length(level_2d_names), dimnames = c(list(names(binary_matrix_list)), list(level_2d_names)))
  PI_level3a = matrix(data= NA, nrow = length(binary_matrix_list), ncol = length(level_3a_names), dimnames = c(list(names(binary_matrix_list)), list(level_3a_names)))
  PI_level3b = matrix(data= NA, nrow = length(binary_matrix_list), ncol = length(level_3b_names), dimnames = c(list(names(binary_matrix_list)), list(level_3b_names)))
  
  perc_undiv_level1 = matrix(data= NA, nrow = length(binary_matrix_list), ncol = length(level_1_names), dimnames = c(list(names(binary_matrix_list)), list(level_1_names)))
  perc_undiv_level2 = matrix(data= NA, nrow = length(binary_matrix_list), ncol = length(level_2_names), dimnames = c(list(names(binary_matrix_list)), list(level_2_names)))
  perc_undiv_level2b = matrix(data= NA, nrow = length(binary_matrix_list), ncol = length(level_2b_names), dimnames = c(list(names(binary_matrix_list)), list(level_2b_names)))
  perc_undiv_level2c = matrix(data= NA, nrow = length(binary_matrix_list), ncol = length(level_2c_names), dimnames = c(list(names(binary_matrix_list)), list(level_2c_names)))
  perc_undiv_level2d = matrix(data= NA, nrow = length(binary_matrix_list), ncol = length(level_2d_names), dimnames = c(list(names(binary_matrix_list)), list(level_2d_names)))
  perc_undiv_level3a = matrix(data= NA, nrow = length(binary_matrix_list), ncol = length(level_3a_names), dimnames = c(list(names(binary_matrix_list)), list(level_3a_names)))
  perc_undiv_level3b = matrix(data= NA, nrow = length(binary_matrix_list), ncol = length(level_3b_names), dimnames = c(list(names(binary_matrix_list)), list(level_3b_names)))
  
  prop_divided_level1 = matrix(data= NA, nrow = length(binary_matrix_list), ncol = length(level_1_names), dimnames = c(list(names(binary_matrix_list)), list(level_1_names)))
  prop_divided_level2 = matrix(data= NA, nrow = length(binary_matrix_list), ncol = length(level_2_names), dimnames = c(list(names(binary_matrix_list)), list(level_2_names)))
  prop_divided_level2b = matrix(data= NA, nrow = length(binary_matrix_list), ncol = length(level_2b_names), dimnames = c(list(names(binary_matrix_list)), list(level_2b_names)))
  prop_divided_level2c = matrix(data= NA, nrow = length(binary_matrix_list), ncol = length(level_2c_names), dimnames = c(list(names(binary_matrix_list)), list(level_2c_names)))
  prop_divided_level2d = matrix(data= NA, nrow = length(binary_matrix_list), ncol = length(level_2d_names), dimnames = c(list(names(binary_matrix_list)), list(level_2d_names)))
  prop_divided_level3a = matrix(data= NA, nrow = length(binary_matrix_list), ncol = length(level_3a_names), dimnames = c(list(names(binary_matrix_list)), list(level_3a_names)))
  prop_divided_level3b = matrix(data= NA, nrow = length(binary_matrix_list), ncol = length(level_3b_names), dimnames = c(list(names(binary_matrix_list)), list(level_3b_names)))
  
  prop_undivided_level1 = matrix(data= NA, nrow = length(binary_matrix_list), ncol = length(level_1_names), dimnames = c(list(names(binary_matrix_list)), list(level_1_names)))
  prop_undivided_level2 = matrix(data= NA, nrow = length(binary_matrix_list), ncol = length(level_2_names), dimnames = c(list(names(binary_matrix_list)), list(level_2_names)))
  prop_undivided_level2b = matrix(data= NA, nrow = length(binary_matrix_list), ncol = length(level_2b_names), dimnames = c(list(names(binary_matrix_list)), list(level_2b_names)))
  prop_undivided_level2c = matrix(data= NA, nrow = length(binary_matrix_list), ncol = length(level_2c_names), dimnames = c(list(names(binary_matrix_list)), list(level_2c_names)))
  prop_undivided_level2d = matrix(data= NA, nrow = length(binary_matrix_list), ncol = length(level_2d_names), dimnames = c(list(names(binary_matrix_list)), list(level_2d_names)))
  prop_undivided_level3a = matrix(data= NA, nrow = length(binary_matrix_list), ncol = length(level_3a_names), dimnames = c(list(names(binary_matrix_list)), list(level_3a_names)))
  prop_undivided_level3b = matrix(data= NA, nrow = length(binary_matrix_list), ncol = length(level_3b_names), dimnames = c(list(names(binary_matrix_list)), list(level_3b_names)))
  
  for(ind in c(1:length(binary_matrix_list))){
    mat = binary_matrix_list[[ind]]
    proliferation = list_proliferation[[names(binary_matrix_list)[ind]]]
    divided  = which(proliferation!=0)
    undivided  = which(proliferation==0)
    
    l1 = apply(mat[,c("CD3/CD10", "CD14", "CD19"  )], 1, paste, collapse = " ")
    t = table(l1, proliferation)
    t1 = rbind(t[,which(as.numeric(colnames(t))!=0)])
    perc_undiv_level1[names(binary_matrix_list)[ind], rownames(t)] = t[,"0"]*100/rowSums(t)
    DI_level1[names(binary_matrix_list)[ind], rownames(t)] = apply(t, 1, function(x){mean(rep(as.numeric(colnames(t)), x))})
    PI_level1[names(binary_matrix_list)[ind], rownames(t)] = apply(t1, 1, function(x){mean(rep(as.numeric(colnames(t1)), x))})
    t = table(l1[divided])
    total1 = 1*100/table(l1)[names(t)]
    prop_divided_level1[names(binary_matrix_list)[ind], rownames(t)] = total1
    t = table(l1[undivided])
    total1 = 1*100/table(l1)[names(t)]
    prop_undivided_level1[names(binary_matrix_list)[ind], rownames(t)] = total1
    
    l2 = apply(mat[,c("IgD"  ,  "CD27")], 1, paste, collapse = " ")
    t = table(l2, proliferation)
    t = t[which(rowSums(t)>=10),]
    t1 = rbind(t[,which(as.numeric(colnames(t))!=0)])
    perc_undiv_level2[names(binary_matrix_list)[ind], rownames(t)] = t[,"0"]*100/rowSums(t)
    DI_level2[names(binary_matrix_list)[ind], rownames(t)] = apply(t, 1, function(x){mean(rep(as.numeric(colnames(t)), x))})
    PI_level2[names(binary_matrix_list)[ind], rownames(t)] = apply(t1, 1, function(x){mean(rep(as.numeric(colnames(t1)), x))})
    t = table(l2[divided])
    total1 = 1*100/table(l2)[names(t)]
    if(sum(t)>=10){prop_divided_level2[names(binary_matrix_list)[ind], rownames(t)] = total1}
    t = table(l2[undivided])
    total1 = 1*100/table(l2)[names(t)]
    if(sum(t)>=10){prop_undivided_level2[names(binary_matrix_list)[ind], rownames(t)] = total1}
    
    l3a = apply(mat[,c("IgD", "CD27", "CD21"  )], 1, paste, collapse = " ")
    t = table(l3a, proliferation)
    t = t[which(rowSums(t)>=10),]
    t1 = rbind(t[,which(as.numeric(colnames(t))!=0)])
    perc_undiv_level3a[names(binary_matrix_list)[ind], rownames(t)] = t[,"0"]*100/rowSums(t)
    DI_level3a[names(binary_matrix_list)[ind], rownames(t)] = apply(t, 1, function(x){mean(rep(as.numeric(colnames(t)), x))})
    PI_level3a[names(binary_matrix_list)[ind], rownames(t)] = apply(t1, 1, function(x){mean(rep(as.numeric(colnames(t1)), x))})
    t = table(l3a[divided])
    total1 = 1*100/table(l3a)[names(t)]
    if(sum(t)>=10){
      prop_divided_level3a[names(binary_matrix_list)[ind],] = 0
      prop_divided_level3a[names(binary_matrix_list)[ind], rownames(t)] = total1}
    t = table(l3a[undivided])
    total1 = 1*100/table(l3a)[names(t)]
    if(sum(t)>=10){
      prop_undivided_level3a[names(binary_matrix_list)[ind],] = 0
      prop_undivided_level3a[names(binary_matrix_list)[ind], rownames(t)] = total1}
    
    l3b = apply(mat[,c("IgD"  ,  "CD27", "CD22"  )], 1, paste, collapse = " ")
    t = table(l3b, proliferation)
    t = t[which(rowSums(t)>=10),]
    t1 = rbind(t[,which(as.numeric(colnames(t))!=0)])
    perc_undiv_level3b[names(binary_matrix_list)[ind], rownames(t)] = t[,"0"]*100/rowSums(t)
    DI_level3b[names(binary_matrix_list)[ind], rownames(t)] = apply(t, 1, function(x){mean(rep(as.numeric(colnames(t)), x))})
    PI_level3b[names(binary_matrix_list)[ind], rownames(t)] = apply(t1, 1, function(x){mean(rep(as.numeric(colnames(t1)), x))})
    t = table(l3b[divided])
    total1 = 1*100/table(l3b)[names(t)]
    if(sum(t)>=10){
      prop_divided_level3b[names(binary_matrix_list)[ind],] = 0
      prop_divided_level3b[names(binary_matrix_list)[ind], rownames(t)] = total1}
    t = table(l3b[undivided])
    total1 = 1*100/table(l3b)[names(t)]
    if(sum(t)>=10){
      prop_undivided_level3b[names(binary_matrix_list)[ind],] = 0
      prop_undivided_level3b[names(binary_matrix_list)[ind], rownames(t)] = total1}
    
    l3b = apply(mat[,c("IgD"  ,  "CD27","IgM","CD38"   )], 1, paste, collapse = " ")
    l3b = cell_type[l3b]
    t = table(l3b, proliferation)
    t = t[which(rowSums(t)>=10),]
    t1 = rbind(t[,which(as.numeric(colnames(t))!=0)])
    perc_undiv_level2b[names(binary_matrix_list)[ind], rownames(t)] = t[,"0"]*100/rowSums(t)
    DI_level2b[names(binary_matrix_list)[ind], rownames(t)] = apply(t, 1, function(x){mean(rep(as.numeric(colnames(t)), x))})
    PI_level2b[names(binary_matrix_list)[ind], rownames(t)] = apply(t1, 1, function(x){mean(rep(as.numeric(colnames(t1)), x))})
    t = table(l3b[divided])
    total1 = 1*100/table(l3b)[names(t)]
    if(sum(t)>=10){
      prop_divided_level2b[names(binary_matrix_list)[ind],] = 0
      prop_divided_level2b[names(binary_matrix_list)[ind], rownames(t)] = total1}
    t = table(l3b[undivided])
    total1 = 1*100/table(l3b)[names(t)]
    if(sum(t)>=10){
      prop_undivided_level2b[names(binary_matrix_list)[ind],] = 0
      prop_undivided_level2b[names(binary_matrix_list)[ind], rownames(t)] = total1}
    
    l3ba = apply(mat[,c("CD21", "CD22", "CD24"  )], 1, paste, collapse = " ")
    l3a = apply(cbind(l3b, l3ba), 1, paste, collapse = " : ")
    t = table(l3a, proliferation)
    t = t[which(rowSums(t)>=5),]
    t1 = rbind(t[,which(as.numeric(colnames(t))!=0)])
    perc_undiv_level2c[names(binary_matrix_list)[ind], rownames(t)] = t[,"0"]*100/rowSums(t)
    DI_level2c[names(binary_matrix_list)[ind], rownames(t)] = apply(t, 1, function(x){mean(rep(as.numeric(colnames(t)), x))})
    PI_level2c[names(binary_matrix_list)[ind], rownames(t)] = apply(t1, 1, function(x){mean(rep(as.numeric(colnames(t1)), x))})
    t = table(l3a[divided])
    total1 = 1*100/table(l3a)[names(t)]
    t = t[which(names(t) %in% colnames(prop_divided_level2c))]
    if(sum(t)>=10){
      prop_divided_level2c[names(binary_matrix_list)[ind],] = 0
      prop_divided_level2c[names(binary_matrix_list)[ind], rownames(t)] = total1}
    t = table(l3a[undivided])
    total1 = 1*100/table(l3a)[names(t)]
    if(sum(t)>=10){
      prop_undivided_level2c[names(binary_matrix_list)[ind],] = 0
      prop_undivided_level2c[names(binary_matrix_list)[ind], rownames(t)] = total1}
    
    print(ind)
  }
  #### proportion of divided and undivided
  
  all_perc_undiv = c(list(perc_undiv_level1), list(perc_undiv_level2), list(perc_undiv_level3a), list(perc_undiv_level3b), list(perc_undiv_level2b), list(perc_undiv_level2c))
  all_DI = c(list(DI_level1), list(DI_level2), list(DI_level3a), list(DI_level3b), list(DI_level2b), list(DI_level2c))
  all_PI = c(list(PI_level1), list(PI_level2), list(PI_level3a), list(PI_level3b), list(PI_level2b), list(PI_level2c))
  all_prop_undiv = c(list(prop_undivided_level1), list(prop_undivided_level2), list(prop_undivided_level3a), list(prop_undivided_level3b), list(prop_undivided_level2b), list(prop_undivided_level2c))  
  all_prop_div = c(list(prop_divided_level1), list(prop_divided_level2), list(prop_divided_level3a), list(prop_divided_level3b), list(prop_divided_level2b), list(prop_divided_level2c))  
  
  names(all_perc_undiv) = paste("Perc. undivided-Level", c("1","2", "3a", "3b","2b","2c"))
  names(all_DI) = paste("DI-Level", c("1","2", "3a", "3b","2b","2c"))
  names(all_PI) = paste("PI-Level", c("1","2", "3a", "3b","2b","2c"))
  names(all_prop_undiv) = paste("Prop. of undivided cells-Level", c("1","2", "3a", "3b","2b","2c"))
  names(all_prop_div) = paste("Prop. of divided cells-Level", c("1","2", "3a", "3b","2b","2c"))
  
  all_counts = c(all_DI, all_PI, all_perc_undiv,all_prop_undiv, all_prop_div)
  
  ##### get means per patient
  
  # Define a function to remove outliers using IQR
  remove_outliers <- function(x) {
    if(length(unique(x))>1){
      Q1 <- quantile(x, 0.35)
      Q3 <- quantile(x, 0.65)
      mean = median(x)
      IQR <- Q3 - Q1
      if(sd(x)!=0){
        x = x[which(x > (mean - 2 * IQR))]
        x = x[which(x < (mean + 2 * IQR))]
      }
      return(x)
    }else{return(mean(x))}
  }
  
  samples = rownames(all_counts[[1]])
  patient_sub = gsub("_A","", gsub("_B","", gsub("_C","", samples, fixed = T), fixed = T), fixed = T)
  patient_subs = sort(unique(patient_sub))
  
  per_patient_counts = NULL
  for(ind in c(1:length(all_counts))){
    mat1 = all_counts[[ind]]
    populations_use = colnames(mat1)
    all_counts_norm = matrix(data= 0, nrow = length(patient_subs), ncol = length(populations_use), dimnames = c(list(patient_subs), list(populations_use)))
    for(p in c(1:length(patient_subs))){
      w = samples[which(patient_sub==patient_subs[p])]
      means = apply(cbind(mat1[w,]), 2, function(x){mean(remove_outliers(x[which(is.na(x)==F)]))})
      all_counts_norm[patient_subs[p],]  = means
    }
    per_patient_counts = c(per_patient_counts, list(all_counts_norm))
    print (ind)
  } 
  names(per_patient_counts) = names(all_counts)
  saveRDS(file = concat(c(in_dir, "FCS_matrices/Proliferation_per_patient.rds")), per_patient_counts)
}
Combine_all_cell_frequencies_subsetted_proliferation(in_dir)

Get_signficance_between_groups_proliferation<-function(){
  
  ################################################################################################# clinical information
  file = "~/Google_Drive/Projects/Fibromyalgia/Data/BCR/Clinical_information.txt"
  p <- as.matrix(read.csv(file, head=TRUE, sep="\t"))
  
  sample_clin = p[,"SampleID"]
  age = as.numeric(p[,"Age"])
  gender = (p[,"Sex"])
  case_control = p[,"Disease"]
  ethnicity = p[,"Ethnicity"]
  Children.birthed = p[,"Children.birthed"]
  
  names(age) = sample_clin
  names(gender) = sample_clin
  names(case_control) = sample_clin
  names(ethnicity) = sample_clin
  names(Children.birthed) = sample_clin
  
  ##############
  per_patient_counts = readRDS(file = concat(c(in_dir, "FCS_matrices/Proliferation_per_patient.rds")))
  
  #####
  exp_groups = per_patient_counts
  ### split the Level 2c groups by cell type and by cD21CD22CD24 status
  split_terms = names(exp_groups)[grep("Level 2c", names(exp_groups))]
  exp_groups_names = names(exp_groups)
  for(ind in c(1:length(split_terms))){
    mat = exp_groups[[split_terms[ind]]]
    str = do.call(rbind, strsplit(colnames(mat)," : ", fixed = T))
    m1 = sort(unique(str[,1]))
    m2 = sort(unique(str[,2]))
    for(i in c(1:length(m1))){
      mat1 = mat[,which(str[,1]==m1[i])]
      exp_groups = c(exp_groups, list(mat1))
      exp_groups_names = c(exp_groups_names, concat(c(split_terms[ind]," : ", m1[i])))
      }
    for(i in c(1:length(m2))){
      mat1 = mat[,which(str[,2]==m2[i])]
      exp_groups = c(exp_groups, list(mat1))
      exp_groups_names = c(exp_groups_names, concat(c(split_terms[ind]," : ", m2[i])))}
  }
  
  names(exp_groups) = exp_groups_names
  
  a = exp_groups[["PI-Level 2b"]][,1]
  b = rowSums(exp_groups[["PI-Level 2b"]])
  
  plot(a,b)
  
  ### 1. run tests between disease cohorts-max severity
  names = rownames(exp_groups[[1]])
  names_pat = gsub("Day_5_","", gsub("Day_6_","", names, fixed = T), fixed = T)
  names_pat = gsub("FMHP0","FMSH", names_pat, fixed = T)
  names_pat[which(names_pat %in% names(case_control)==F)]
  names(names_pat) = names
  groups_ids = NULL
  Sources = sort(unique(case_control))
  for(i in c(1:length(Sources))){
    w = which(case_control==Sources[i])
    w = names(w)[which(names(w) %in% names_pat)]
    groups_ids = c(groups_ids, list((w)))
  }
  names(groups_ids) = Sources
  
  fileout1=concat(c(in_dir, "FCS_matrices/All_populations_boxplots_all_proliferation.pdf"))
  w=4.5
  pdf(file=fileout1, height=w*1.4*6, width=w*2.5*6)
  par(mfrow= c(6,6), mar = c(35,4,4,0.5))
  summary_stats = NULL
  for(ind in c(1:length(exp_groups))){
    mat_stat = cbind(exp_groups[[ind]])
    rownames(mat_stat) = gsub("Day_5_", "", gsub("Day_6_", "", rownames(mat_stat) , fixed = T), fixed = T)
    rownames(mat_stat) = gsub("FMHP0","FMSH", rownames(mat_stat), fixed = T)
    use = which(apply(mat_stat, 2, function(x){length(unique(x))})>=10)
    mat_stat = cbind(mat_stat[,use])
    use = which(apply(mat_stat, 2, function(x){length(which(is.na(x)==T))})<=5)
    mat_stat = cbind(mat_stat[,use])
    factor_sex = factor(gender[rownames(mat_stat)])
    factor_age = age[rownames(mat_stat)]
    factor_preg = as.numeric(Children.birthed[rownames(mat_stat)])
    factor_preg[which(factor_preg>0)] = 1
    factor = factor(case_control[rownames(mat_stat)])
    sources= levels(factor)
    
    #x = mat_stat[,"IgD- CD27+ IgM- CD38+"]
    #x[order(x)]
    #library(lmPerm)    ## non parametric manova
    #result <- lmp(mat_stat ~ factor +factor_age+factor_sex+factor_preg)
    if(dim(mat_stat)[2]>1){
      #fit = manova(formula = log10(mat_stat+1) ~ factor +factor_age+factor_sex+factor_preg)
      fit = manova(formula = log10(mat_stat+1)~ factor +factor_age+factor_sex)
      p1 = summary.aov(fit)
      nam = gsub(" Response ","",names(p1))
      p_value = NULL
      means = NULL
      i1 = 0
      for(i in p1){
        i1 = i1+1
        p_value = c(p_value, i$'Pr(>F)'[1]) 
        if(length(mean)==0){means = Means_factor(factor, mat_stat[,i1])
        }else{means = rbind(means, Means_factor(factor, mat_stat[,i1]))}
      }
      p_value[which(is.na(p_value))] = 2
      names(p_value) = nam
      print(min(p_value))
      colnames(means) = paste("mean.group.", c(1:length(means[1,])))
      combined_p_value = cbind(p_value ,means)
      rownames(combined_p_value) = nam
      p.group = rep(names(exp_groups)[ind], length(nam))
      sublabel = gsub(concat(c(names(exp_groups)[ind],"..")), "", colnames(mat_stat), fixed = T)
      x = cbind(p.group, sublabel,combined_p_value)
      if(length(summary_stats)==0){summary_stats = x
      }else{summary_stats = rbind(summary_stats, x)}
    }
    if(dim(mat_stat)[2]==1){
      fit = aov(formula = mat_stat ~ factor +factor_age+factor_sex)#+factor_preg)
      p1 = summary.aov(fit)
      p_value =  p1[[1]][1,4]
      print(min(p_value))
      
      combined_p_value = rbind(c(p_value, Means_factor(factor, mat_stat)))
      colnames(combined_p_value) = c("p-value", paste("mean.group.", c(1:length(levels(factor)))))
      rownames(combined_p_value) = "overall"
      sublabel = "overall"
      p.group = names(exp_groups)[ind]
      x = cbind(p.group, sublabel,combined_p_value)
      if(length(summary_stats)==0){summary_stats = x
      }else{summary_stats = rbind(summary_stats, x)}
    }
    if(dim(mat_stat)[2]>=1){
      
      
      #mat_stat = cbind(exp_groups[[ind]])
      #rownames(mat_stat) = gsub("Day_5_", "", gsub("Day_6_", "", rownames(mat_stat) , fixed = T), fixed = T)
      #rownames(mat_stat) = gsub("FMHP0","FMSH", rownames(mat_stat), fixed = T)
      mat = mat_stat
      p_value_all = rep(1, length(mat_stat[1,]))
      p_value_all[use] = p_value
      groups = NULL
      pvals =NULL
      all_data = NULL
      all_pats = NULL
      all_cell_types = NULL
      for(i in c(1:length(mat[1,]))){
        g1 = NULL
        for(p in c(1:length(groups_ids))){
          x = mat[groups_ids[[p]],i]
          x=x[which(x!=-1)]
          g1 = c(g1, list( x))
          all_data = c(all_data, x)
          all_pats = c(all_pats, names(x))
          all_cell_types = c(all_cell_types, rep(colnames(mat)[i], length(x)))
        }
        names(g1) = names(groups_ids)
        pvals = c(pvals, wilcox.test(g1[[1]], y=g1[[2]])$p.value)
        groups = c(groups, list(g1))
      }
      names(groups)=colnames(mat)
      names(pvals)=colnames(mat)
      
      factor_sex = factor(gender[all_pats])
      factor_age = age[all_pats]
      factor_preg = as.numeric(Children.birthed[all_pats])
      factor_preg[which(factor_preg>0)] = 1
      factor = factor(case_control[all_pats])
      sources= levels(factor)
      
      fit = aov(formula = all_data ~ factor+factor_sex+factor_age+factor_preg)
      pv = summary.lm(fit)$ coefficients["factorFMS",4]
      p_value_use = rep(1, length(groups))
      names(p_value_use) = names(groups)
      p_value_use[names(p_value)] = p_value

      factors1 = names(groups_ids)
      factors = names(groups)
      if(is.null(factors)){factors = "overall"}
      main = concat(c(names(exp_groups)[ind],"\nANOVA p-val:", signif(pv, digits = 3)))
      max = max(c(unlist(groups), unlist(groups))*1.35)
      min = min(c(unlist(groups), unlist(groups))*1.15)
      min = 0
      if(min(unlist(groups))<0){
        min = -1
        max = 1.15}
      b = (max-min)*0.034
      ylab = ""
      draw_signif_lines = TRUE
      y = max(c(unlist(groups), unlist(groups))*1)+b
      max_width = 32
      max_scale = min(c(max,100))
      range = max-min
      if(range>50){scale = c(-100:100)*20}
      if(range<=50){scale = c(-100:100)*10}
      if(range<=30){scale = c(-100:100)*5}
      if(range <10){scale = c(-100:100)*2.5}
      if(range <5){scale = c(-100:100)*1}
      if(range <4){scale = c(-100:100)*0.5}
      if(range <1.5){scale = c(-100:1000)*0.2}
      if(range <0.5){scale = c(-100:100)*0.1}
      if(range <0.1){scale = c(-100:100)*0.01}
      if(range <0.01){scale = c(-100:100)*0.001}
      cex = 0.9
      Fun<-function(x){x}
      
      scale = scale[intersect(which(scale<= max_scale), which(scale>=min))]
      plot(c(1.5, max_width +0.5),c(min, max), pch=20, col="white",xlab="",ylab ="",cex=cex, cex.lab=cex+0.1,	cex.axis=cex,cex.main = cex, col.axis = "white",tck=0, mgp = c(2,0,0), main = main, axes = FALSE, ylim = c(min, max))
      mtext(side = 2, text = ylab, line = 2.8,cex= cex-0.1, las = 3, font = 1)
      mtext(side = 1, text = gsub("_","/",factors,fixed=T), line = 0.35,cex= cex-0.1,  at = c(1:length(factors)), las = 2, font = 1)
      segments(0.5,Fun(scale),length(groups)+0.5,Fun(scale),col = "grey",lwd = 1,lty = 3 )
      segments(0.5,Fun(scale),length(groups)+0.5,Fun(scale),col = "grey",lwd = 1,lty = 3 )
      mtext(side = 2, text = scale, line = 0.15,cex= cex-0.1,  at =Fun(scale), las = 2, font = 1)
      width = 0.17
      index = 1
      l = length(groups)
      l1 = length(groups[[1]])
      shift = c(1:l1)
      shift = (mean(shift)-shift)
      shift = shift*0.23/max(shift)
      
      library(RColorBrewer)
      cols1 =  add.alpha (brewer.pal(8, "Dark2")[c(3,6)], alpha = 0.95)
      cols =  add.alpha (cols1, alpha = 0.5)
      
      for(i in c(1:l)){
        for(i1 in c(1:l1)){
          points1=as.numeric(groups[[i]][[i1]])
          box1<-c(as.numeric(quantile(points1))[3], as.numeric(quantile(points1, probs = c(0.1, 0.9))), as.numeric(quantile(points1))[c(2, 4)])	
          Draw_box_plot(box1,i-shift[i1],width,cols[i1],1, cols1[i1])
          points(rep(i-shift[i1], length(points1)),points1, pch =21, col=cols[i1],bg = cols[i1], cex = 0.7)
        }}
      
      for(i in c(1:l)){	
        b = max*0.035
        signif_threshold = 0.05
        pval1= "NS"
        if(p_value_use[i]<signif_threshold){
          pval1 = "*"
          if(p_value_use[i]<signif_threshold/10){pval1 = "**"}
          y = max(unlist(groups[[i]]))
          y = y+3*b
          text(i, y+2*b, labels = pval1, cex = 1.3)
        }
      }
    }}
  plot(c(0,1), c(0,1), pch = 21, col = "white", bg ="white", xlab = "", ylab = "", main = "", axes=F)
  
  legend("topleft", Sources, pch = 21,cex= 0.8, bty="n", pt.bg = cols, col = cols, pt.lwd = 2, text.font = 2)
  dev.off()
  
  out_file = concat(c(in_dir, "FCS_matrices/All_populations_boxplots_all_proliferation.txt"))
  write.table(summary_stats, file = out_file, append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = TRUE,col.names = TRUE, qmethod = c("escape", "double"),fileEncoding = "")
  
  ### heatmap of p-values
  
  
  Plot_heatmaps_of_differences<-function(p, main){
    group = p[,2]
    str = do.call(rbind, strsplit(group," : ", fixed = T))
    group = str[,1]
    class = str[,2]
    p_value = as.numeric(p[,3])
    mean1 = as.numeric(p[,4])
    mean2 = as.numeric(p[,5])
    diff = mean1-mean2
    mean = (mean1+mean2)/2
    
    classes = sort(unique(class))
    groups = sort(unique(group))
    
    m_DGE_genes_p_value = matrix(data = 2, nrow = length(groups),ncol = length(classes), dimnames = c(list(groups), list(classes)))
    m_DGE_genes_FC = matrix(data = 0, nrow = length(groups),ncol = length(classes), dimnames = c(list(groups), list(classes)))
    m_DGE_genes_mean = matrix(data = 0, nrow = length(groups),ncol = length(classes), dimnames = c(list(groups), list(classes)))
    for(g in c(1:length(groups))){
      w = which(group==groups[g])
      m_DGE_genes_p_value[groups[g], class[w]] = p_value[w]
      m_DGE_genes_FC[groups[g], class[w]] = diff[w]
      m_DGE_genes_mean[groups[g], class[w]] = mean[w]
    }
    
    matrix = m_DGE_genes_p_value
    matrix_mean = m_DGE_genes_mean
    matrix_FC = m_DGE_genes_FC
    
    m_col = matrix+2
    library(RColorBrewer)
    cols1 =  add.alpha (c(brewer.pal(8, "Dark2")[c(3,6)],"grey"), alpha = 0.95)
    cols =  add.alpha (cols1, alpha = 0.5)
    cols2 =  add.alpha (cols, alpha = 0.5)	
    
    ynames = rownames(matrix)
    xnames = colnames(matrix)
    ### plot heatmap
    x_range = c(1,length(xnames))
    y_range = c(1,length(ynames))
    cex = 0.9
    
    fileout1=concat(c(in_dir, "FCS_matrices/Heatmap", main,"_differences_proliferation.pdf"))
    w=6
    pdf(file=fileout1, height=w*0.8, width=w*1.85)
    y_range_max = 12
    x_start = 0.5
    par(mfrow= c(1,2), mar = c(5,12,8,5))
    
    shift = y_range_max - max(y_range)
    plot(c(0.5,max(x_range)), c(-1,(length(y_range))+0.5), pch=20, col="white",xlab="",ylab ="",cex=cex, cex.lab=cex+0.1, cex.axis=cex,cex.main = cex, col.axis = "white",tck=0, mgp = c(2,0,0), main = "", 
         ylim = c(0, max(y_range_max)-0),xlim = c(x_start,max(x_range)+5.5), axes = F)
    mtext(side = 3, text = gsub(".","-",xnames, fixed = T), line = 0.1,cex=0.7,  at =c(1:length(xnames)), las = 2, font = 1)
    label = gsub("_"," ",gsub("Percentage","%",gsub("mean","Mean",ynames, fixed = T), fixed = T), fixed = T)
    label = gsub("J gene freq by uniq VDJ ","",gsub(".","/",label, fixed = T), fixed = T)
    label = gsub("V gene freq by cluster ","",label, fixed = T)
    mtext(side = 2, text = label, line = 0.01,cex=0.7,  at =c(1:length(ynames))+shift, las = 1, font = 1)
    
    m_col = round(matrix, digits = 0)+1
    
    segments(0.5, c(1:length(ynames))+shift , max(x_range), c(1:length(ynames))+shift, col = "grey", lty = 3, lwd = 0.5)
    segments(c(1:length(xnames)),1+shift, c(1:length(xnames)) , max(y_range)+0.5+shift,  col = "grey", lty = 3, lwd = 0.5)
    size = matrix_mean^0.25
    size = 1.5*size/max(size)
    pches = c( 25,24)
    for(i in c(1:length(m_col[,1]))){
      for(j in c(1:length(m_col[1,]))){
        if(matrix[i,j]!=0){
          pch = 21
          col = cols1[3]
          if(matrix[i,j]<0.05){
            if(matrix_FC[i,j]>0){
              pch = pches[1]
              col = cols1[1]}
            if(matrix_FC[i,j]<0){
              pch = pches[2]
              col = cols1[2]}
          }
          points(j,i+shift, pch = pch, bg = col, col = col, cex = size[i,j])
        }}}
    
    plot(c(0,1), c(0,1), pch=20, col="white",xlab="",ylab ="",col.axis = "white",tck=0, mgp = c(2,0,0), main ='', axes = FALSE)
    legs = c("significantly lower in FMS", "significantly higher in FMS", "NS")
    legend("topleft", legs, pch = c((pches),21),cex= 0.8, bty="n", pt.bg = cols1[c(1:length(cols1))], col = NA, pt.cex = 1.2)
    
    dev.off()
    
  }
  p <- summary_stats
  w2 = grep("DI-Level 2c : ", p[,"p.group"])
  w1 = grep("CD21", p[,"p.group"], invert = T)
  w = intersect(w1,w2)
  main = "DI_Level_2c"
  p = p[w,]
  Plot_heatmaps_of_differences(p, main)
  
  p <- summary_stats
  w2 = grep("PI-Level 2c : ", p[,"p.group"])
  w1 = grep("CD21", p[,"p.group"], invert = T)
  w = intersect(w1,w2)
  main = "PI_Level_2c"
  p = p[w,]
  Plot_heatmaps_of_differences(p, main)
  
}
Get_signficance_between_groups_proliferation()

Combine_all_cell_frequencies_subsetted_perc_expressing<-function(in_dir){
  live = readRDS(file = concat(c(in_dir, "FCS_matrices/Filtering_step1_live.rds")))

  Get_gates<-function(in_dir,mat){
    thresholds <- as.matrix(read.csv(concat(c(in_dir, "MFI_thresholds.txt")), head=TRUE, sep="\t"))
    samples= thresholds[,"Day_sample"]
    colnames(thresholds) = gsub(".","/",colnames(thresholds), fixed = T)
    cols_use = colnames(mat)[which(colnames(mat) %in% colnames(thresholds)==T)]
    thresholds = apply(thresholds[,cols_use], 2, as.numeric) 
    rownames(thresholds) = samples
    colnames(thresholds) = cols_use
    return(thresholds)
  }
  thresholds =Get_gates(in_dir, live[[1]])
  
  markers = setdiff(colnames(live[[1]]), c("FSC-A","FSC-H","FSC-W" ,"SSC-A","SSC-H","SSC-W" ,"Proliferation Dye" , "Live/dead" ,"Time" ))
  markers = c( "CD3/CD10", "CD14"   ,   "CD19"  ,  "IgD"  ,  "IgM"   ,   "CD27"   , "CD21"   , "CD38"    , "CD24"    ,   "CD22")
  
  all_populations = NULL
  for(i in c(1:length(markers))){
    if(i == 1){
      all_populations = c( concat(c(markers[i],"+")),  concat(c(markers[i],"-")))
    }else{
      new_pops = c(paste(all_populations, concat(c(markers[i],"+"))), paste(all_populations, concat(c(markers[i],"-"))))
      all_populations = new_pops
    }
    print(i)
    print(length(all_populations))
  }
  length(unique(all_populations))
  
  run = grep("Prolif_Check", names(live), invert = T)
  
  binary_matrix_list = NULL  
  percentage_live = NULL
  
  for(ind in run){
    sample = names(live)[ind]
    mat = live[[ind]]
    mat1 = mat[,markers]*0
    for(i in c(1:length(markers))){
      threshold = thresholds[sample, markers[i]]
      x = mat[,markers[i]]
      mat1[which(x>threshold),markers[i]] = 1
    }
    mat2 = mat1
    for(i in c(1:length(markers))){
      mat2[which(mat1[, markers[i]]>0), markers[i]] = concat(c(markers[i],"+"))
      mat2[which(mat1[, markers[i]]==0), markers[i]] = concat(c(markers[i],"-"))
    }
    binary_matrix_list = c(binary_matrix_list, list(mat2))
    print(ind) 
  }
  names(binary_matrix_list) = names(live)[run]

  p <- as.matrix(read.csv("~/Google_Drive/Projects/Fibromyalgia/Cell_culture/FMS/B_cell_phenotypes.txt", head=T, sep="\t"))
  cell_type =p[,1]
  l1 = apply(p[,c(2,3,4,5 )], 1, paste, collapse = " ")
  names(cell_type) = l1
  
  level_1_names = NULL
  level_2_names = NULL
  level_3_names = NULL
  
  seq = ceiling(seq(from = 1, to = length(binary_matrix_list), length = 10))
  names(binary_matrix_list)[seq]
  for(ind in seq){
    print(ind)
    l3b = apply(binary_matrix_list[[ind]][,c("IgD"  ,  "CD27","IgM","CD38"   )], 1, paste, collapse = " ")
    l3b = cell_type[l3b]
    l3ba = apply(binary_matrix_list[[ind]][,c("CD21", "CD22", "CD24"  )], 1, paste, collapse = " ")
    combinations <- expand.grid(unique(l3b), unique(l3ba))
    l3a = apply(combinations, 1, paste, collapse = " : ")
    level_3_names = sort(unique(c(level_3_names, l3a)))
    print(length(level_3_names))
    
    l3ba = apply(binary_matrix_list[[ind]][,c("CD22", "CD24" )], 1, paste, collapse = " ")
    combinations <- expand.grid(unique(l3b), unique(l3ba))
    l3a = apply(combinations, 1, paste, collapse = " : ")
    level_1_names = sort(unique(c(level_1_names, l3a)))
    
    l3ba = apply(binary_matrix_list[[ind]][,c("CD21", "CD22" )], 1, paste, collapse = " ")
    combinations <- expand.grid(unique(l3b), unique(l3ba))
    l3a = apply(combinations, 1, paste, collapse = " : ")
    level_2_names = sort(unique(c(level_2_names, l3a)))
  }

  counts_level1 = matrix(data= NA, nrow = length(binary_matrix_list), ncol = length(level_1_names), dimnames = c(list(names(binary_matrix_list)), list(level_1_names)))
  counts_level2 = matrix(data= NA, nrow = length(binary_matrix_list), ncol = length(level_2_names), dimnames = c(list(names(binary_matrix_list)), list(level_2_names)))
  counts_level3 = matrix(data= NA, nrow = length(binary_matrix_list), ncol = length(level_3_names), dimnames = c(list(names(binary_matrix_list)), list(level_3_names)))

  for(ind in c(1:length(binary_matrix_list))){
    mat = binary_matrix_list[[ind]]
    
    l3b = apply(mat[,c("IgD"  ,  "CD27","IgM","CD38"   )], 1, paste, collapse = " ")
    l3b = cell_type[l3b]
    l3ba = apply(mat[,c("CD21", "CD22", "CD24"  )], 1, paste, collapse = " ")
    l3a = apply(cbind(l3b, l3ba), 1, paste, collapse = " : ")
    t = table(l3b, l3ba)
    for(i in c(1:length(t[,1]))){
      t[i,] = t[i,]*100/sum(t[i,])
      names = paste(rownames(t)[i], colnames(t), sep = " : ")
      counts_level3[names(binary_matrix_list)[ind], names] = t[i,]
    }
    print (ind)
    
    l3ba = apply(mat[,c( "CD22", "CD24"  )], 1, paste, collapse = " ")
    l3a = apply(cbind(l3b, l3ba), 1, paste, collapse = " : ")
    t = table(l3b, l3ba)
    for(i in c(1:length(t[,1]))){
      t[i,] = t[i,]*100/sum(t[i,])
      names = paste(rownames(t)[i], colnames(t), sep = " : ")
      counts_level1[names(binary_matrix_list)[ind], names] = t[i,]
    }
    
    l3ba = apply(mat[,c( "CD21", "CD22"  )], 1, paste, collapse = " ")
    l3a = apply(cbind(l3b, l3ba), 1, paste, collapse = " : ")
    t = table(l3b, l3ba)
    for(i in c(1:length(t[,1]))){
      t[i,] = t[i,]*100/sum(t[i,])
      names = paste(rownames(t)[i], colnames(t), sep = " : ")
      counts_level2[names(binary_matrix_list)[ind], names] = t[i,]
    }
    print (ind)
  }
  all_counts = c(list(counts_level1), list(counts_level2), list(counts_level3))
  names(all_counts) = c("level 1", "level 2", "level 3")
  ##### get means per patient
  
  # Define a function to remove outliers using IQR
  remove_outliers <- function(x) {
    x = x[which(is.na(x)==F)]
    if(length(unique(x))>1){
      Q1 <- quantile(x, 0.35)
      Q3 <- quantile(x, 0.65)
      mean = median(x)
      IQR <- Q3 - Q1
      if(sd(x)!=0){
        x = x[which(x > (mean - 2 * IQR))]
        x = x[which(x < (mean + 2 * IQR))]
      }
    }else{
      x = mean(x)
    }
    return(x)
  }
  
  samples = rownames(all_counts[[1]])
  patient_sub = gsub("_A","", gsub("_B","", gsub("_C","", samples, fixed = T), fixed = T), fixed = T)
  patient_subs = sort(unique(patient_sub))
  
  per_patient_counts_norm = NULL
  
  for(ind in c(1:length(all_counts))){
    mat1 = all_counts[[ind]]
    populations_use = colnames(mat1)
    all_counts_norm = matrix(data= 0, nrow = length(patient_subs), ncol = length(populations_use), dimnames = c(list(patient_subs), list(populations_use)))
    for(p in c(1:length(patient_subs))){
      w = samples[which(patient_sub==patient_subs[p])]
      means = apply(mat1[w,], 2, function(x){mean(remove_outliers(x))})
      all_counts_norm[patient_subs[p],]  = means
    }
    per_patient_counts_norm = c(per_patient_counts_norm, list(all_counts_norm))
    print (ind)
  } 
  names(per_patient_counts_norm) = names(all_counts)
  saveRDS(file = concat(c(in_dir, "FCS_matrices/Counts_per_patient_timepoint_perc_of_markers.rds")), per_patient_counts_norm)
}
Combine_all_cell_frequencies_subsetted_perc_expressing(in_dir)

Get_signficance_between_groups_perc_expressing<-function(){
  ################################################################################################# clinical information
  file = "~/Google_Drive/Projects/Fibromyalgia/Data/BCR/Clinical_information.txt"
  p <- as.matrix(read.csv(file, head=TRUE, sep="\t"))
  
  sample_clin = p[,"SampleID"]
  age = as.numeric(p[,"Age"])
  gender = (p[,"Sex"])
  case_control = p[,"Disease"]
  ethnicity = p[,"Ethnicity"]
  Children.birthed = p[,"Children.birthed"]
  
  names(age) = sample_clin
  names(gender) = sample_clin
  names(case_control) = sample_clin
  names(ethnicity) = sample_clin
  names(Children.birthed) = sample_clin
  
  ##############
  per_patient_counts = readRDS(file = concat(c(in_dir, "FCS_matrices/Counts_per_patient_timepoint_perc_of_markers.rds")))
  
  ### make list split by day and type
  list_proportions = NULL
  names_proportions = NULL
  for(ind in c(1:length(per_patient_counts))){
    mat1 = per_patient_counts[[ind]]
    name = names(per_patient_counts)[ind]
    
    #### for % within gate
    day0_mat1 = mat1[grep("Day_0", rownames(mat1)), ]
    day5_mat1 = mat1[grep("Day_0", rownames(mat1), invert = T), ]
    rownames(day0_mat1) = gsub("Day_0_", "", rownames(day0_mat1), fixed = T)
    rownames(day5_mat1) = gsub("Day_5_", "", rownames(day5_mat1), fixed = T)
    rownames(day5_mat1) = gsub("Day_6_", "", rownames(day5_mat1), fixed = T)
    
    w = which(rownames(day0_mat1) %in% sample_clin==F)
    rownames(day0_mat1)[w] = gsub("FMHP0","FMSH",rownames(day0_mat1)[w])
    w = which(rownames(day5_mat1) %in% sample_clin==F)
    rownames(day5_mat1)[w] = gsub("FMHP0","FMSH",rownames(day5_mat1)[w])
    
    day0_mat1 = day0_mat1[sort(rownames(day0_mat1)), ]
    day5_mat1 = day5_mat1[sort(rownames(day5_mat1)), ]
    rownames(day0_mat1)==rownames(day5_mat1)
    change_counts_for_all_populations_norm = (day5_mat1-day0_mat1)/(day5_mat1+day0_mat1)
    
    ## split by cell type: 
    colnames = colnames(day0_mat1)
    str = do.call(rbind, strsplit(colnames, " : ", fixed = T))
    unique_groups = sort(unique(str[,1]))
    for(i in c(1:length(unique_groups))){
      cols_use = colnames[which(str[,1] ==unique_groups[i])]
      list_proportions = c(list_proportions, list(day0_mat1[,cols_use]), list(day5_mat1[,cols_use]), list(change_counts_for_all_populations_norm[,cols_use]))
      names_proportions = c(names_proportions, paste(concat(c(name, unique_groups[i])), c("day 0", "day 5", "change")))
    }
  }
    
    
   
  names(list_proportions) = names_proportions
  
  saveRDS(file = concat(c(in_dir, "FCS_matrices/All_populations_boxplots_all_perc_expressing.rds")), list_proportions)
  
  #####
  exp_groups = list_proportions
  
  ### 1. run tests between disease cohorts-max severity
  names = rownames(exp_groups[[1]])
  names_pat = gsub("Day_5_","", gsub("Day_6_","", names, fixed = T), fixed = T)
  names_pat = gsub("FMHP0","FMSH", names_pat, fixed = T)
  names_pat[which(names_pat %in% names(case_control)==F)]
  names(names_pat) = names
  groups_ids = NULL
  Sources = sort(unique(case_control))
  for(i in c(1:length(Sources))){
    w = which(case_control==Sources[i])
    w = names(w)[which(names(w) %in% names_pat)]
    groups_ids = c(groups_ids, list((w)))
  }
  names(groups_ids) = Sources
  
  fileout1=concat(c(in_dir, "FCS_matrices/All_populations_boxplots_all_perc_expressing.pdf"))
  w=4.5
  pdf(file=fileout1, height=w*1.4*3, width=w*2.5*9)
  par(mfrow= c(3,9), mar = c(35,4,4,0.5))
  summary_stats = NULL
  for(ind in c(1:length(exp_groups))){
    mat_stat = cbind(exp_groups[[ind]])
    rownames(mat_stat) = gsub("Day_5_", "", gsub("Day_6_", "", rownames(mat_stat) , fixed = T), fixed = T)
    rownames(mat_stat) = gsub("FMHP0","FMSH", rownames(mat_stat), fixed = T)
    use = which(apply(mat_stat, 2, function(x){length(unique(x))})>=10)
    mat_stat = cbind(mat_stat[,use])
    factor_sex = factor(gender[rownames(mat_stat)])
    factor_age = age[rownames(mat_stat)]
    factor_preg = as.numeric(Children.birthed[rownames(mat_stat)])
    factor_preg[which(factor_preg>0)] = 1
    factor = factor(case_control[rownames(mat_stat)])
    sources= levels(factor)
    
    
    #library(lmPerm)    ## non parametric manova
    #result <- lmp(mat_stat ~ factor +factor_age+factor_sex+factor_preg)
    if(dim(mat_stat)[2]>1){
      fit = manova(formula = mat_stat ~ factor +factor_age+factor_sex+factor_preg)
      p1 = summary.aov(fit)
      nam = gsub(" Response ","",names(p1))
      p_value = NULL
      means = NULL
      i1 = 0
      for(i in p1){
        i1 = i1+1
        p_value = c(p_value, i$'Pr(>F)'[1]) 
        if(length(mean)==0){means = Means_factor(factor, mat_stat[,i1])
        }else{means = rbind(means, Means_factor(factor, mat_stat[,i1]))}
      }
      p_value[which(is.na(p_value))] = 2
      names(p_value) = nam
      print(min(p_value))
      colnames(means) = paste("mean.group.", c(1:length(means[1,])))
      combined_p_value = cbind(p_value ,means)
      rownames(combined_p_value) = nam
      p.group = rep(names(exp_groups)[ind], length(nam))
      sublabel = gsub(concat(c(names(exp_groups)[ind],"..")), "", colnames(mat_stat), fixed = T)
      x = cbind(p.group, sublabel,combined_p_value)
      if(length(summary_stats)==0){summary_stats = x
      }else{summary_stats = rbind(summary_stats, x)}
    
      }else{
      fit = aov(formula = mat_stat ~ factor +factor_age+factor_sex+factor_preg)
      p1 = summary.aov(fit)
      p_value =  p1[[1]][1,4]
      print(min(p_value))
      
      combined_p_value = rbind(c(p_value, Means_factor(factor, mat_stat)))
      colnames(combined_p_value) = c("p-value", paste("mean.group.", c(1:length(levels(factor)))))
      rownames(combined_p_value) = "overall"
      sublabel = "overall"
      p.group = names(exp_groups)[ind]
      x = cbind(p.group, sublabel,combined_p_value)
      if(length(summary_stats)==0){summary_stats = x
      }else{summary_stats = rbind(summary_stats, x)}
    }
    
    mat = mat_stat
    p_value_all = rep(1, length(mat_stat[1,]))
    p_value_all[use] = p_value
    groups = NULL
    pvals =NULL
    use = NULL
    all_data = NULL
    all_pats = NULL
    all_cell_types = NULL
    for(i in c(1:length(mat[1,]))){
      g1 = NULL
      for(p in c(1:length(groups_ids))){
        x = mat[groups_ids[[p]],i]
        x=x[which(x!=-1)]
        x=x[which(is.na(x)==F)]
        g1 = c(g1, list( x))
        all_data = c(all_data, x)
        all_pats = c(all_pats, names(x))
        all_cell_types = c(all_cell_types, rep(colnames(mat)[i], length(x)))
      }
      names(g1) = names(groups_ids)
      groups = c(groups, list(g1))
    }
    names(groups)=colnames(mat)

    factor_sex = factor(gender[all_pats])
    factor_age = age[all_pats]
    factor_preg = as.numeric(Children.birthed[all_pats])
    factor_preg[which(factor_preg>0)] = 1
    factor = factor(case_control[all_pats])
    sources= levels(factor)
    
    fit = aov(formula = all_data ~ factor+factor_sex+factor_age+factor_preg)
    pv = summary.lm(fit)$ coefficients["factorFMS",4]
    
    factors1 = names(groups_ids)
    factors = names(groups)
    if(is.null(factors)){factors = "overall"}
    main = concat(c(names(exp_groups)[ind],"\nANOVA p-val:", signif(pv, digits = 3)))
    max = max(c(unlist(groups), unlist(groups))*1.35)
    min = min(c(unlist(groups), unlist(groups))*1.15)
    min = 0
    if(min(unlist(groups))<0){
      min = -1
      max = 1.15}
    b = (max-min)*0.034
    ylab = ""
    draw_signif_lines = TRUE
    y = max(c(unlist(groups), unlist(groups))*1)+b
    max_width = 32
    max_scale = min(c(max,100))
    range = max-min
    if(range>50){scale = c(-100:100)*20}
    if(range<=50){scale = c(-100:100)*10}
    if(range<=30){scale = c(-100:100)*5}
    if(range <10){scale = c(-100:100)*2.5}
    if(range <5){scale = c(-100:100)*1}
    if(range <4){scale = c(-100:100)*0.5}
    if(range <1.5){scale = c(-100:1000)*0.2}
    if(range <0.5){scale = c(-100:100)*0.1}
    if(range <0.1){scale = c(-100:100)*0.01}
    if(range <0.01){scale = c(-100:100)*0.001}
    cex = 0.9
    Fun<-function(x){x}
    
    scale = scale[intersect(which(scale<= max_scale), which(scale>=min))]
    plot(c(1.5, max_width +0.5),c(min, max), pch=20, col="white",xlab="",ylab ="",cex=cex, cex.lab=cex+0.1,	cex.axis=cex,cex.main = cex, col.axis = "white",tck=0, mgp = c(2,0,0), main = main, axes = FALSE, ylim = c(min, max))
    mtext(side = 2, text = ylab, line = 2.8,cex= cex-0.1, las = 3, font = 1)
    mtext(side = 1, text = gsub("_","/",factors,fixed=T), line = 0.35,cex= cex-0.1,  at = c(1:length(factors)), las = 2, font = 1)
    segments(0.5,Fun(scale),length(groups)+0.5,Fun(scale),col = "grey",lwd = 1,lty = 3 )
    segments(0.5,Fun(scale),length(groups)+0.5,Fun(scale),col = "grey",lwd = 1,lty = 3 )
    mtext(side = 2, text = scale, line = 0.15,cex= cex-0.1,  at =Fun(scale), las = 2, font = 1)
    width = 0.17
    index = 1
    l = length(groups)
    l1 = length(groups[[1]])
    shift = c(1:l1)
    shift = (mean(shift)-shift)
    shift = shift*0.23/max(shift)
    
    library(RColorBrewer)
    cols1 =  add.alpha (brewer.pal(8, "Dark2")[c(3,6)], alpha = 0.95)
    cols =  add.alpha (cols1, alpha = 0.5)
    
    for(i in c(1:l)){
      for(i1 in c(1:l1)){
        points1=as.numeric(groups[[i]][[i1]])
        box1<-c(as.numeric(quantile(points1))[3], as.numeric(quantile(points1, probs = c(0.1, 0.9))), as.numeric(quantile(points1))[c(2, 4)])	
        Draw_box_plot(box1,i-shift[i1],width,cols[i1],1, cols1[i1])
        points(rep(i-shift[i1], length(points1)),points1, pch =21, col=cols[i1],bg = cols[i1], cex = 0.7)
      }}
    
    for(i in c(1:l)){	
      b = max*0.035
      signif_threshold = 0.05
      pval1= "NS"
      if(p_value_all[i]<signif_threshold){pval1 = "*"
      if(p_value_all[i]<signif_threshold/10){pval1 = "**"}
      y = max(unlist(groups[[i]]))
      y = y+3*b
      text(i, y+2*b, labels = pval1, cex = 1.3)
      }}
  }
  plot(c(0,1), c(0,1), pch = 21, col = "white", bg ="white", xlab = "", ylab = "", main = "", axes=F)
  
  legend("topleft", Sources, pch = 21,cex= 0.8, bty="n", pt.bg = cols1, col = NA, pt.lwd = 2, text.font = 2)
  dev.off()
  
  out_file = concat(c(in_dir, "FCS_matrices/All_populations_boxplots_all_perc_expressing.txt"))
  write.table(summary_stats, file = out_file, append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = TRUE,col.names = TRUE, qmethod = c("escape", "double"),fileEncoding = "")
  
  ### heatmap of p-values
  
  Plot_heatmaps_of_differences<-function(p, main){
    group = p[,2]
    str = do.call(rbind, strsplit(group," : ", fixed = T))
    group = str[,1]
    class = str[,2]
    p_value = as.numeric(p[,3])
    mean1 = as.numeric(p[,4])
    mean2 = as.numeric(p[,5])
    diff = mean1-mean2
    mean = (mean1+mean2)/2
    
    classes = sort(unique(class))
    groups = sort(unique(group))
    
    m_DGE_genes_p_value = matrix(data = 2, nrow = length(groups),ncol = length(classes), dimnames = c(list(groups), list(classes)))
    m_DGE_genes_FC = matrix(data = 0, nrow = length(groups),ncol = length(classes), dimnames = c(list(groups), list(classes)))
    m_DGE_genes_mean = matrix(data = 0, nrow = length(groups),ncol = length(classes), dimnames = c(list(groups), list(classes)))
    for(g in c(1:length(groups))){
      w = which(group==groups[g])
      m_DGE_genes_p_value[groups[g], class[w]] = p_value[w]
      m_DGE_genes_FC[groups[g], class[w]] = diff[w]
      m_DGE_genes_mean[groups[g], class[w]] = mean[w]
    }
    
    matrix = t(m_DGE_genes_p_value)
    matrix_mean = t(m_DGE_genes_mean)
    matrix_FC = t(m_DGE_genes_FC)
    
    m_col = matrix+2
    library(RColorBrewer)
    cols1 =  add.alpha (c(brewer.pal(8, "Dark2")[c(3,6)],"grey"), alpha = 0.95)
    cols =  add.alpha (cols1, alpha = 0.5)
    cols2 =  add.alpha (cols, alpha = 0.5)	
    
    ynames = rownames(matrix)
    xnames = colnames(matrix)
    ### plot heatmap
    x_range = c(1,length(xnames))
    y_range = c(1,length(ynames))
    cex = 0.9
    
    fileout1=concat(c(in_dir, "FCS_matrices/Heatmap", main,"_differences_expression.pdf"))
    w=6
    pdf(file=fileout1, height=w*0.7, width=w*1.95)
    y_range_max = 12
    x_start = 0.5
    par(mfrow= c(1,2), mar = c(3,12,10,5))
    
    shift = y_range_max - max(y_range)
    plot(c(0.5,max(x_range)), c(-1,(length(y_range))+0.5), pch=20, col="white",xlab="",ylab ="",cex=cex, cex.lab=cex+0.1, cex.axis=cex,cex.main = cex, col.axis = "white",tck=0, mgp = c(2,0,0), main = "", 
         ylim = c(0, max(y_range_max)+0.2),xlim = c(x_start,max(x_range)+5.5), axes = F)
    mtext(side = 3, text = gsub(".","-",xnames, fixed = T), line = 0.1,cex=0.7,  at =c(1:length(xnames)), las = 2, font = 1)
    label = gsub("_"," ",gsub("Percentage","%",gsub("mean","Mean",ynames, fixed = T), fixed = T), fixed = T)
    label = gsub("J gene freq by uniq VDJ ","",gsub(".","/",label, fixed = T), fixed = T)
    label = gsub("V gene freq by cluster ","",label, fixed = T)
    mtext(side = 2, text = label, line = 0.01,cex=0.7,  at =c(1:length(ynames))+shift, las = 1, font = 1)
    
    m_col = round(matrix, digits = 0)+1
    
    segments(0.5, c(1:length(ynames))+shift , max(x_range), c(1:length(ynames))+shift, col = "grey", lty = 3, lwd = 0.5)
    segments(c(1:length(xnames)),1+shift, c(1:length(xnames)) , max(y_range)+0.5+shift,  col = "grey", lty = 3, lwd = 0.5)
    size = matrix_mean^0.25
    size = 1.5*size/max(size)
    pches = c( 25,24)
    for(i in c(1:length(m_col[,1]))){
      for(j in c(1:length(m_col[1,]))){
        if(matrix[i,j]!=0){
          pch = 21
          col = cols1[3]
          if(matrix[i,j]<0.05){
            if(matrix_FC[i,j]>0){
              pch = pches[1]
              col = cols1[1]}
            if(matrix_FC[i,j]<0){
              pch = pches[2]
              col = cols1[2]}
          }
          points(j,i+shift, pch = pch, bg = col, col = col, cex = size[i,j])
        }}}
    
    plot(c(0,1), c(0,1), pch=20, col="white",xlab="",ylab ="",col.axis = "white",tck=0, mgp = c(2,0,0), main ='', axes = FALSE)
    legs = c("significantly lower in FMS", "significantly higher in FMS", "NS")
    legend("topleft", legs, pch = c((pches),21),cex= 0.8, bty="n", pt.bg = cols1[c(1:length(cols1))], col = NA, pt.cex = 1.2)
    
    dev.off()
    
  }
  
  levels = c("level 1", "level 2", "level 3")
  for(l in c(1:length(levels))){
    p= summary_stats[intersect(grep("day 0", summary_stats[,1]), grep(levels[l], summary_stats[,1])),]
    main = concat(c("level_",l,"_Day_0"))
    Plot_heatmaps_of_differences(p, main)
    
    p= summary_stats[intersect(grep("day 5", summary_stats[,1]), grep(levels[l], summary_stats[,1])),]
    main = concat(c("level_",l,"_Day_5"))
    min(as.numeric(p[,"p_value"]))
    Plot_heatmaps_of_differences(p, main)
  
    p= summary_stats[intersect(grep("change", summary_stats[,1]), grep(levels[l], summary_stats[,1])),]
    main = concat(c("level_",l,"_change"))
    min(as.numeric(p[,"p_value"]))
    head(sort(as.numeric(p[,"p_value"])))
    Plot_heatmaps_of_differences(p, main)
  }
}
Get_signficance_between_groups_perc_expressing()