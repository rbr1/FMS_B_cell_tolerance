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
  p=p[which(p[,"Experiment.type"] %in% c("FMS data", "Proliferation check")),]
  w = c(grep("Day_5",p[,"Day_sample"]), grep("Day_6",p[,"Day_sample"]), grep("Prolif_Check",p[,"Day_sample"]))
  p=p[w,]
  locations = gsub(".//",in_dir,p[,"Location"], fixed = T)
  samples= p[,"Day_sample"]
  names(locations) = samples
  patient = p[,"Patient"]
  return(cbind(locations, patient))
}
in_dir = "~/Google_Drive/Projects/Fibromyalgia/Cell_culture/FMS/"
file_table = concat(c(in_dir, "Files_information.txt"))

locations_patient = Get_files_data(in_dir, file_table)
locations = locations_patient[,1]
patient = locations_patient[,2]
patients = sort(unique(patient))

###############################################################

library(flowCore)

cex = 0.9
cols<-NULL
alpha<-0.5
cols<-c(cols, rgb(rbind(c(1,0,0)), alpha=alpha))
cols<-c(cols, rgb(rbind(c(0.4,0.9,0.4)), alpha=0.8))
cols<-c(cols, rgb(rbind(c(0,0.3,1)), alpha=alpha/2))
cols<-c(cols, rgb(rbind(c(0.5,0.5,0.5)), alpha=alpha))
cols<-c(cols, rgb(rbind(c(0,1,1)), alpha=alpha))

gate2_cells = readRDS(file = concat(c(in_dir, "FCS_matrices/Filtering_step1_Gate2.rds")))
live = readRDS(file = concat(c(in_dir, "FCS_matrices/Filtering_step1_live.rds")))


p <- as.matrix(read.csv("~/Google_Drive/Projects/Fibromyalgia/Cell_culture/FMS/Proliferation_thresholds.txt", head=TRUE, sep="\t"))
zero_prolif = as.numeric(p[,2])
separation = as.numeric(p[,3])
names(zero_prolif) = p[,1]
names(separation) = p[,1]

Get_proliferation<-function(list_exp_mat, main, midpoints_between_means1_sub){
  Get_hist_points<-function(hist){
    x1 = c(min(hist$mids), hist$mids, max(hist$mids))
    y1 = c(0,hist$density,0)
    return(cbind(x1,y1))
  }
  hists = NULL
  range = NULL
  rangey = NULL
  for(i in c(1:length(list_exp_mat))){
    mat = list_exp_mat[[i]]
    #mat[which(mat<0)] = 0
    x = log10(mat[,"Proliferation Dye"]+1)
    x = x[which(is.na(x)==F)]
    #hist(x, breaks = 100)
    h <- hist(x, breaks = 300, plot = FALSE)
    hists = c(hists, list(h))
    xy = Get_hist_points(h)
    spline_fit <- smooth.spline(xy[,1], xy[,2], spar = 0.25)
    spline_fit$ y
    range = range(c(range, x))
    rangey = range(c(rangey, spline_fit$ y))
  }
  
  col123 = add.alpha(rainbow(length(list_exp_mat)), alpha = 0.5)
  #par(mfrow= c(1,4), mar = c(4, 4, 3, 0.1))
  plot(Get_hist_points(hists[[1]]),type = "l",lwd = 2, col = "white", xlim = c(0,5.5), ylim = rangey, main = main)
  
  for(i in c(1:length(list_exp_mat))){
    xy = Get_hist_points(hists[[i]])
    spline_fit <- smooth.spline(xy[,1], xy[,2], spar = 0.25)
    lines(spline_fit, col = col123[i], lwd = 2)
    #points(,type = "l",lwd = 2, col = col123[i])
  }
  segments(midpoints_between_means1_sub, 0, midpoints_between_means1_sub, 1000, col = "grey", lwd = 3)
  b = (midpoints_between_means1_sub[1]-midpoints_between_means1_sub[2])/2
  for(i in c(1:length(midpoints_between_means1_sub))){
    text(midpoints_between_means1_sub[i]+b,y= max(rangey)*0.65, labels = i-1, col = "black", cex = 0.9, pos = 1)
  }
  
  # Add a legend
  leg = names(list_exp_mat)
  leg[grep("Prolif_Check",leg)] = "Prolif_Check"
  legend("topleft", legend = leg, fill = col123)
}

fileout=concat(c(in_dir, "FCS_matrices/Filtering_step1_proliferation_grouped.jpeg"))
w = 1800
jpeg(file=fileout, height=w*5, width=w*6 , res=600)
par(mfrow= c(5,6), mar = c(4, 4, 3, 0.1))


which(patients=="FMHP020")

overall = NULL
for(s in c(1:length(patients))){
  sample_group = locations[which(patient==patients[s])]
  sample_group = c(sample_group[grep("Day_5", sample_group)], sample_group[grep("Prolif_Check", sample_group)])
  sample_group2 = sample_group[which(names (sample_group) %in% names(gate2_cells)==F)]
  sample_group = sample_group[which(names (sample_group) %in% names(gate2_cells))]
  print(as.character(names(sample_group)))
  list_exp_mat = NULL
  x_all = NULL
  for(ind in c(1:length(sample_group))){
    sample = names(sample_group)[ind]
    list_exp_mat = c(list_exp_mat, list(gate2_cells[[sample]]))
    mat = gate2_cells[[sample]]
    x = log10(mat[,"Proliferation Dye"]+1)
    x = x[which(is.na(x)==F)]
    x_all = c(x_all, x)
  }
  if(length(sample_group2)!=0){
    for(ind in c(1:length(sample_group2))){
      sample = names(sample_group2)[ind]
      list_exp_mat = c(list_exp_mat, list(live[[sample]]))
      
      mat = live[[sample]]
      x = log10(mat[,"Proliferation Dye"]+1)
      x = x[which(is.na(x)==F)]
      x_all = c(x_all, x)
    }
    names(list_exp_mat) = c(names(sample_group), names(sample_group2))
  }else{names(list_exp_mat) = c(names(sample_group))}
  overall = c(overall, list(x))
  midpoints_between_means1_sub = seq(from = zero_prolif[patients[s]], to = zero_prolif[patients[s]]-(separation[patients[s]]*10), length = 10)
  main =  patients[s]
  Get_proliferation(list_exp_mat,main, midpoints_between_means1_sub)
}

dev.off()
