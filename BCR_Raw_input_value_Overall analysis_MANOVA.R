Boxplot_custom<-function(groups, main, width_plot, colsx){
  factors = names(groups)
  max = max(c(unlist(groups), unlist(groups))*1.35)
  min = 0
  if(min(unlist(groups))<0){
    min = min(c(unlist(groups), unlist(groups))*1.15)}
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
concat = function(v) {
  res = ""
  for (i in 1:length(v)){res = paste0(res,v[i])}
  res
}
add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha)) }

###########
out_dir = "~/Google_Drive/Projects/Fibromyalgia/Data/BCR/"
batch = "FMS1"
################################################################################################# clinical information
file = "~/Google_Drive/Projects/Fibromyalgia/Data/BCR/Clinical_information.txt"
p <- as.matrix(read.csv(file, head=TRUE, sep="\t"))
p=p[grep("Lupus", p[,"PMH"], invert = T),]
p=p[grep("carcinoma", p[,"PMH"], invert = T),]
sample_clin = p[,"SampleID"]
age = as.numeric(p[,"Age"])
gender = (p[,"Sex"])
case_control = p[,"Disease"]
ethnicity = p[,"Ethnicity"]
Children.birthed =as.numeric( p[,"Children.birthed"])

names(age) = sample_clin
names(gender) = sample_clin
names(case_control) = sample_clin
names(ethnicity) = sample_clin
names(Children.birthed) = sample_clin


#################################################################################################
file = "~/Google_Drive/Projects/Fibromyalgia/Data/BCR/All_raw_values_FMS1.txt"
p <- as.matrix(read.csv(file, head=TRUE, sep="\t"))
id = rownames(p)
intersect(id, sample_clin)
id[which(id %in% sample_clin==F)]

w = which(id %in% sample_clin==T)
p=p[w,]
headers = colnames(p)
headers = headers[grep("IGHEP",headers, invert = T)]
headers = headers[grep("IGHGP",headers, invert = T)]
headers = headers[grep("Reyni",headers, invert = T)]

p=p[,headers]

id = rownames(p)

Age = age[id]
Sex = gender[id]
Diagnosis = case_control[id]

################ get values
headers = colnames(p)

exp_group = strsplit(headers,"..",fixed = T)
exp_isotype = NULL
for(i in c(1:length(exp_group))){
  exp_isotype = c(exp_isotype, exp_group[[i]][2])
  exp_group[i] = exp_group[[i]][1]
}
exp_group = unlist(exp_group)
exp_group_full = cbind(headers,exp_isotype,exp_group)
exp_groups = unique(exp_group)
mat = p[, headers]
mat_numeric = matrix(data = -2, nrow = length(id), ncol = length(headers), dimnames = c(list(id),list(headers)))
for(i in c(1:length(headers))){
  mat_numeric[, headers[i]] = as.numeric(p[, headers[i]])
}

################################################
### 1. run tests between disease cohorts-max severity
groups_ids = NULL
Sources = sort(unique(Diagnosis))
for(i in c(1:length(Sources))){
  w = which(Diagnosis==Sources[i])
  groups_ids = c(groups_ids, list(id[w]))
}
names(groups_ids) =Sources

#################################################
summary_stats = NULL
prefix = "Disease"
library(RColorBrewer)
cols1 =  add.alpha (brewer.pal(8, "Dark2")[c(3,6)], alpha = 0.95)
cols =  add.alpha (cols1, alpha = 0.5)

fileout1=concat(c(out_dir, prefix,"_boxplots_all.pdf"))
w=4.5
pdf(file=fileout1, height=w*1.*3, width=w*1.4*4)
par(mfrow= c(4,4), mar = c(12,4,4,0.5))
all_BCR_values = c(list(exp_group_full))
for(ind in c(1:length(exp_groups))){
	w_analysis = exp_group_full[which(exp_group== exp_groups[ind]),]
	mat = mat_numeric[,w_analysis[,1]]
	if(length(grep("Gini_Index", exp_groups[ind]))==1){mat = mat[,grep("all", colnames(mat), invert = T)]}
	if(exp_groups[ind] %in% c("D5","D10","D50")){mat = mat[,grep("all", colnames(mat), invert = T)]}
	all_BCR_values = c(all_BCR_values, list(mat))
	if(length(mat[1,])<=25){
	  colus = colnames(mat)
	  n_per_group = matrix(data = 0, nrow = length(names(groups_ids)), ncol = length(colus), dimnames = c(list(names(groups_ids)),list(colus)))
	  var_per_group = matrix(data = 0, nrow = length(names(groups_ids)), ncol = length(colus), dimnames = c(list(names(groups_ids)),list(colus)))
	  for(i in c(1:length(Sources))){
	    for(j in c(1:length(colus))){
	      w = which(mat[groups_ids[[i]], colus[j]]!=-1)
	      n_per_group[i,j] = length(w)
	      var_per_group[i,j] = length(unique(mat[groups_ids[[i]], colus[j]]))
	    }}
	  w1 = which(apply(var_per_group, 2, min)>=3)
	  w2 = which(apply(n_per_group, 2, min)>=3)
	  run_headers = intersect(w1,w2)
	  mat = mat[, run_headers]
	  mat[which(mat==-1)] = NA
	  header_info = w_analysis[run_headers,]
	  
	  mat_stat = mat[unlist(groups_ids),]
	  factor_sex = factor(Sex[rownames(mat_stat)])
	  factor_age = Age[rownames(mat_stat)]
	  factor_preg = Children.birthed[rownames(mat_stat)]
	  factor = factor(Diagnosis)
	  
	  fit = manova(formula = mat_stat ~ factor +factor_age+factor_sex+factor_preg)
	  p1 = summary.aov(fit)
	  nam = gsub(" Response ","",names(p1))
	  p_value = NULL
	  means = NULL
	  i1 = 0
	  for(i in p1){
	    i1 = i1+1
	    p_value = c(p_value, i$'Pr(>F)'[1]) 
	    if(length(mean)==0){means = Means_factor(factor, mat[,i1])
	    }else{means = rbind(means, Means_factor(factor, mat[,i1]))}
	  }
	  p_value[which(is.na(p_value))] = 2
	  names(p_value) = nam
	  print(min(p_value))
	  
	  sublabel = gsub(concat(c(exp_groups[ind],"..")), "", colnames(mat_stat), fixed = T)
	  
	  colnames(means) = paste("mean.group.", c(1:length(means[1,])))
	  combined_p_value = cbind(p_value ,means)
	  rownames(combined_p_value) = nam
	  p.group = rep(exp_groups[ind], length(p_value))
	  x = cbind(p.group, sublabel,combined_p_value)
	  if(length(summary_stats)==0){summary_stats = x
	  }else{summary_stats = rbind(summary_stats, x)}
	  
	  mat = mat_stat
	  groups = NULL
	  pvals =NULL
	  for(i in c(1:length(mat[1,]))){
	    g1 = NULL
	    for(p in c(1:length(groups_ids))){
	      x = mat[groups_ids[[p]],i]
	      x=x[which(x!=-1)]
	      g1 = c(g1, list( x))}
	    names(g1) = names(groups_ids)
	    pvals = c(pvals, wilcox.test(g1[[1]], y=g1[[2]])$p.value)
	    groups = c(groups, list(g1))
	  }
	  names(groups)=colnames(mat)
	  names(pvals)=colnames(mat)
	  
	  factors1 = names(groups_ids)
	  factors = gsub("IGHD.M","IGHD/M", gsub(".IGH","-IGH", sublabel, fixed = T), fixed = T)
	  main = exp_groups[ind]
	  max = max(c(unlist(groups), unlist(groups))*1.35)
	  min = min(c(unlist(groups), unlist(groups))*1.15)
	  min = 0
	  b = (max-min)*0.034
	  ylab = ""
	  draw_signif_lines = TRUE
	  y = max(c(unlist(groups), unlist(groups))*1)+b
	  max_width = 25
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
	  plot(c(0.5, max_width +0.5),c(min, max), pch=20, col="white",xlab="",ylab ="",cex=cex, cex.lab=cex+0.1,	cex.axis=cex,cex.main = cex, col.axis = "white",tck=0, mgp = c(2,0,0), main = main, axes = FALSE, ylim = c(min, max))
	  mtext(side = 2, text = ylab, line = 2.8,cex= cex-0.1, las = 3, font = 1)
	  mtext(side = 1, text = gsub("_","/",factors,fixed=T), line = 0.35,cex= cex-0.1,  at = c(1:length(factors)), las = 2, font = 1)
	  segments(0.5,Fun(scale),length(groups)+0.5,Fun(scale),col = "grey",lwd = 1,lty = 3 )
	  segments(0.5,Fun(scale),length(groups)+0.5,Fun(scale),col = "grey",lwd = 1,lty = 3 )
	  mtext(side = 2, text = scale, line = 0.15,cex= cex-0.1,  at =Fun(scale), las = 2, font = 1)
	  width = 0.18
	  index = 1
	  l = length(groups)
	  l1 = length(groups[[1]])
	  shift = c(1:l1)
	  shift = (mean(shift)-shift)
	  shift = shift*0.25/max(shift)
	  
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
	    if(p_value[i]<signif_threshold){pval1 = "*"
	    y = max(unlist(groups[[i]]))
	    y = y+3*b
	    text(i, y+2*b, labels = pval1, cex = 1.3)
	    }}
	}
}
plot(c(0,1), c(0,1), pch = 21, col = "white", bg ="white", xlab = "", ylab = "", main = "", axes=F)

legend("topleft", Sources, pch = 21,cex= 0.8, bty="n", pt.bg = cols, col = cols, pt.lwd = 2, text.font = 2)

dev.off()
names(all_BCR_values) = c( "all", exp_groups)

saveRDS(file = concat(c(out_dir, prefix,"_boxplots_all.rds")), all_BCR_values)



fileout1=concat(c(out_dir, prefix,"_boxplots_all_large.pdf"))
w=4.5
pdf(file=fileout1, height=w*1.*3, width=w*4*1)
par(mfrow= c(4,1), mar = c(12,4,4,0.5))

for(ind in c(1:length(exp_groups))){
  w_analysis = exp_group_full[which(exp_group== exp_groups[ind]),]
  mat = mat_numeric[,w_analysis[,1]]
  if(length(mat[1,])>20){
    mat[which(mat==-1)] = NA
    a = apply(mat, 2, function(x){length(which(is.na(x)))})
    mat = mat[,which(a<10)]
    if(length(unique(sort(mat)))>5){
      w_analysis = w_analysis[which(a<10),]
      colus = colnames(mat)
      n_per_group = matrix(data = 0, nrow = length(names(groups_ids)), ncol = length(colus), dimnames = c(list(names(groups_ids)),list(colus)))
      var_per_group = matrix(data = 0, nrow = length(names(groups_ids)), ncol = length(colus), dimnames = c(list(names(groups_ids)),list(colus)))
      for(i in c(1:length(Sources))){
        for(j in c(1:length(colus))){
          w = which(mat[groups_ids[[i]], colus[j]]!=-1)
          n_per_group[i,j] = length(w)
          var_per_group[i,j] = length(unique(mat[groups_ids[[i]], colus[j]]))
        }}
      w1 = which(apply(var_per_group, 2, min)>=3)
      w2 = which(apply(n_per_group, 2, min)>=3)
      run_headers = intersect(w1,w2)
      if(length(run_headers)>5){
        mat = mat[, run_headers]
      
        header_info = w_analysis[run_headers,]
        
        mat_stat = mat[unlist(groups_ids),]
        factor_sex = factor(Sex[rownames(mat_stat)])
        factor_age = Age[rownames(mat_stat)]
        factor_preg = Children.birthed[rownames(mat_stat)]
        factor = factor(Diagnosis)
        
        fit = manova(formula = mat_stat ~ factor +factor_age+factor_sex+factor_preg)
        p1 = summary.aov(fit)
        nam = gsub(" Response ","",names(p1))
        p_value = NULL
        means = NULL
        i1 = 0
        for(i in p1){
          i1 = i1+1
          p_value = c(p_value, i$'Pr(>F)'[1]) 
          if(length(mean)==0){means = Means_factor(factor, mat[,i1])
          }else{means = rbind(means, Means_factor(factor, mat[,i1]))}
        }
        p_value[which(is.na(p_value))] = 2
        names(p_value) = nam
        print(min(p_value))
        
        sublabel = gsub(concat(c(exp_groups[ind],"..")), "", colnames(mat_stat), fixed = T)
        
        colnames(means) = paste("mean.group.", c(1:length(means[1,])))
        combined_p_value = cbind(p_value ,means)
        rownames(combined_p_value) = nam
        p.group = rep(exp_groups[ind], length(p_value))
        x = cbind(p.group, sublabel,combined_p_value)
        if(length(summary_stats)==0){summary_stats = x
        }else{summary_stats = rbind(summary_stats, x)}
        
        
        
        mat = mat_stat
        groups = NULL
        for(i in c(1:length(mat[1,]))){
          g1 = NULL
          for(p in c(1:length(groups_ids))){
            x = mat[groups_ids[[p]],i]
            x=x[which(x!=-1)]
            g1 = c(g1, list( x))}
          names(g1) = names(groups_ids)
          groups = c(groups, list(g1))
        }
        names(groups)=colnames(mat)
        
        factors1 = names(groups_ids)
        factors = sublabel
        main = exp_groups[ind]
        max = max(c(unlist(groups), unlist(groups))*1.35)
        min = min(c(unlist(groups), unlist(groups))*1.15)
        min = 0
        b = (max-min)*0.034
        ylab = ""
        draw_signif_lines = TRUE
        y = max(c(unlist(groups), unlist(groups))*1)+b
        max_width = 80
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
        plot(c(0.5, max_width +0.5),c(min, max), pch=20, col="white",xlab="",ylab ="",cex=cex, cex.lab=cex+0.1,	cex.axis=cex,cex.main = cex, col.axis = "white",tck=0, mgp = c(2,0,0), main = main, axes = FALSE, ylim = c(min, max))
        mtext(side = 2, text = ylab, line = 2.8,cex= cex-0.1, las = 3, font = 1)
        mtext(side = 1, text = gsub("/Non Debris/Singlets/"," ",factors,fixed=T), line = 0.35,cex= cex-0.1,  at = c(1:length(factors)), las = 2, font = 1)
        segments(0.5,Fun(scale),length(groups)+0.5,Fun(scale),col = "grey",lwd = 1,lty = 3 )
        segments(0.5,Fun(scale),length(groups)+0.5,Fun(scale),col = "grey",lwd = 1,lty = 3 )
        mtext(side = 2, text = scale, line = 0.15,cex= cex-0.1,  at =Fun(scale), las = 2, font = 1)
        width = 0.18
        index = 1
        l = length(groups)
        l1 = length(groups[[1]])
        shift = c(1:l1)
        shift = (mean(shift)-shift)
        shift = shift*0.25/max(shift)
        
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
          if(p_value[i]<signif_threshold){pval1 = "*"
          y = max(unlist(groups[[i]]))
          y = y+3*b
          text(i, y+2*b, labels = pval1, cex = 1.3)
          }}
      }}
}}
plot(c(0,1), c(0,1), pch = 21, col = "white", bg ="white", xlab = "", ylab = "", main = "", axes=F)

legend("topleft", Sources, pch = 21,cex= 0.8, bty="n", pt.bg = cols, col = cols, pt.lwd = 2, text.font = 2)

dev.off()

out_file_table = concat(c(out_dir, prefix,"_statistical_summary_age_sex_adjusted.txt"))
write.table(summary_stats, file = out_file_table, append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = FALSE,col.names = TRUE, qmethod = c("escape", "double"),fileEncoding = "")


#### 
prefix = "Disease"
file = concat(c(out_dir, prefix,"_statistical_summary_age_sex_adjusted.txt"))
p <- as.matrix(read.csv(file, head=TRUE, sep="\t"))
p=p[which(as.numeric(p[,3])<0.05),]
group = p[,1]
class = p[,2]
p_value = as.numeric(p[,3])
mean1 = as.numeric(p[,4])
mean2 = as.numeric(p[,5])
diff = mean1-mean2
mean = (mean1+mean2)/2
exp_groups_sub = sort(unique(group))
classes = sort(unique(class))

fileout1=concat(c(out_dir, prefix,"_boxplots_all_large_signif_only.pdf"))
w=4.5
pdf(file=fileout1, height=w*1.1*3, width=w*4*1)
par(mfrow= c(4,1), mar = c(12,4,4,0.5))

for(ind in c(1:length(exp_groups_sub))){
  w_analysis = exp_group_full[which(exp_group== exp_groups_sub[ind]),]
  w_analysis = w_analysis[which(w_analysis[,2] %in% classes),]
  mat = mat_numeric[,w_analysis[,1]]
  if(length(mat[1,])>20){
    mat[which(mat==-1)] = NA
    a = apply(mat, 2, function(x){length(which(is.na(x)))})
    mat = mat[,which(a<10)]
    if(length(unique(sort(mat)))>5){
      w_analysis = w_analysis[which(a<10),]
      colus = colnames(mat)
      n_per_group = matrix(data = 0, nrow = length(names(groups_ids)), ncol = length(colus), dimnames = c(list(names(groups_ids)),list(colus)))
      var_per_group = matrix(data = 0, nrow = length(names(groups_ids)), ncol = length(colus), dimnames = c(list(names(groups_ids)),list(colus)))
      for(i in c(1:length(Sources))){
        for(j in c(1:length(colus))){
          w = which(mat[groups_ids[[i]], colus[j]]!=-1)
          n_per_group[i,j] = length(w)
          var_per_group[i,j] = length(unique(mat[groups_ids[[i]], colus[j]]))
        }}
      w1 = which(apply(var_per_group, 2, min)>=3)
      w2 = which(apply(n_per_group, 2, min)>=3)
      run_headers = intersect(w1,w2)
      if(length(run_headers)>5){
        mat = mat[, run_headers]
        
        header_info = w_analysis[run_headers,]
        
        mat_stat = mat[unlist(groups_ids),]
        factor_sex = factor(Sex[rownames(mat_stat)])
        factor_age = Age[rownames(mat_stat)]
        factor = factor(Diagnosis)
        
        fit = manova(formula = mat_stat ~ factor +factor_age+factor_sex)
        p1 = summary.aov(fit)
        nam = gsub(" Response ","",names(p1))
        p_value = NULL
        means = NULL
        i1 = 0
        for(i in p1){
          i1 = i1+1
          p_value = c(p_value, i$'Pr(>F)'[1]) 
          if(length(mean)==0){means = Means_factor(factor, mat[,i1])
          }else{means = rbind(means, Means_factor(factor, mat[,i1]))}
        }
        p_value[which(is.na(p_value))] = 2
        names(p_value) = nam
        print(min(p_value))
        
        sublabel = gsub(concat(c(exp_groups_sub[ind],"..")), "", colnames(mat_stat), fixed = T)
        
        colnames(means) = paste("mean.group.", c(1:length(means[1,])))
        combined_p_value = cbind(p_value ,means)
        rownames(combined_p_value) = nam
        p.group = rep(exp_groups_sub[ind], length(p_value))
        x = cbind(p.group, sublabel,combined_p_value)
        if(length(summary_stats)==0){summary_stats = x
        }else{summary_stats = rbind(summary_stats, x)}
        
        mat = mat_stat
        groups = NULL
        for(i in c(1:length(mat[1,]))){
          g1 = NULL
          for(p in c(1:length(groups_ids))){
            x = mat[groups_ids[[p]],i]
            x=x[which(x!=-1)]
            g1 = c(g1, list( x))}
          names(g1) = names(groups_ids)
          groups = c(groups, list(g1))
        }
        names(groups)=colnames(mat)
        
        factors1 = names(groups_ids)
        factors = sublabel
        main = exp_groups_sub[ind]
        max = max(c(unlist(groups), unlist(groups))*1.35)
        min = min(c(unlist(groups), unlist(groups))*1.15)
        min = 0
        b = (max-min)*0.034
        ylab = ""
        draw_signif_lines = TRUE
        y = max(c(unlist(groups), unlist(groups))*1)+b
        max_width = 80
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
        plot(c(0.5, max_width +0.5),c(min, max), pch=20, col="white",xlab="",ylab ="",cex=cex, cex.lab=cex+0.1,	cex.axis=cex,cex.main = cex, col.axis = "white",tck=0, mgp = c(2,0,0), main = main, axes = FALSE, ylim = c(min, max))
        mtext(side = 2, text = ylab, line = 2.8,cex= cex-0.1, las = 3, font = 1)
        mtext(side = 1, text = gsub("/Non Debris/Singlets/"," ",factors,fixed=T), line = 0.35,cex= cex-0.1,  at = c(1:length(factors)), las = 2, font = 1)
        segments(0.5,Fun(scale),length(groups)+0.5,Fun(scale),col = "grey",lwd = 1,lty = 3 )
        segments(0.5,Fun(scale),length(groups)+0.5,Fun(scale),col = "grey",lwd = 1,lty = 3 )
        mtext(side = 2, text = scale, line = 0.15,cex= cex-0.1,  at =Fun(scale), las = 2, font = 1)
        width = 0.18
        index = 1
        l = length(groups)
        l1 = length(groups[[1]])
        shift = c(1:l1)
        shift = (mean(shift)-shift)
        shift = shift*0.25/max(shift)
        
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
          if(p_value[i]<signif_threshold){pval1 = "*"
          y = max(unlist(groups[[i]]))
          y = y+3*b
          text(i, y+2*b, labels = pval1, cex = 1.3)
          }}
      }}
  }}
plot(c(0,1), c(0,1), pch = 21, col = "white", bg ="white", xlab = "", ylab = "", main = "", axes=F)

legend("topleft", Sources, pch = 21,cex= 0.8, bty="n", pt.bg = cols, col = cols, pt.lwd = 2, text.font = 2)

dev.off()

################################################ heatmap of changes
Plot_heatmaps_of_differences<-function(out_dir){
  prefix = "Disease"
  file = concat(c(out_dir, prefix,"_statistical_summary_age_sex_adjusted.txt"))
  p <- as.matrix(read.csv(file, head=TRUE, sep="\t"))
  group = p[,1]
  class = p[,2]
  p_value = as.numeric(p[,3])
  mean1 = as.numeric(p[,4])
  mean2 = as.numeric(p[,5])
  diff = mean1-mean2
  mean = (mean1+mean2)/2
  
  classes = sort(unique(class))
  groups = sort(unique(group))
  
  classes_uses = c(list(c("IGHA1" ,"IGHA2" ,"IGHD","IGHE" ,"IGHG1"  ,"IGHG2" ,"IGHG3","IGHG4","IGHM")), 
                   list(c("IGHJ1","IGHJ2","IGHJ3" ,"IGHJ4","IGHJ5","IGHJ6" )), 
                   list(classes[grep("IGHV", classes)]))
  names(classes_uses) = c("isotypes_raw","IGHJ", "IGHV")
  
  classes_uses[[3]] = classes_uses[[3]] [grep("IGHJ",  classes_uses[[3]] , invert = T)]
  
  for(c in c(1:length(classes_uses))){
    g1 = NULL
    for(g in c(1:length(groups))){
      w = which(group==groups[g])
      cl = class[w]
      if(length(which(classes_uses[[c]] %in% cl))>= min(c(length(classes_uses[[c]]), 15))){
        print (groups[g])
        g1 = c(g1, groups[g])
      }
      if(c==1){if(groups[g] =="Percentage_unique_BCRs_in_switched"){g1 = c(g1, groups[g])}}
    }
    
    m_DGE_genes_p_value = matrix(data = 2, nrow = length(g1),ncol = length((classes_uses[[c]])), dimnames = c(list(g1), list((classes_uses[[c]]))))
    m_DGE_genes_FC = matrix(data = 0, nrow = length(g1),ncol = length((classes_uses[[c]])), dimnames = c(list(g1), list((classes_uses[[c]]))))
    m_DGE_genes_mean = matrix(data = 0, nrow = length(g1),ncol = length((classes_uses[[c]])), dimnames = c(list(g1), list((classes_uses[[c]]))))
    for(g in c(1:length(g1))){
      w = intersect(which(group==g1[g]), which(class %in% classes_uses[[c]]))
      cl = class[w]
      m_DGE_genes_p_value[g1[g], cl] = p_value[w]
      m_DGE_genes_FC[g1[g], cl] = diff[w]
      m_DGE_genes_mean[g1[g], cl] = mean[w]*100/sum(mean[w]) ## normalised
    }
    order = c(1:length(rownames(m_DGE_genes_p_value)))
    if(c==1){
      order = rev(rownames(m_DGE_genes_p_value)[c(grep("unique",rownames(m_DGE_genes_p_value)),
                                                  grep("Gini",rownames(m_DGE_genes_p_value)), grep("D5",rownames(m_DGE_genes_p_value)), grep("D10",rownames(m_DGE_genes_p_value)), 
                                                  grep("mutations",rownames(m_DGE_genes_p_value)),grep("CDR3",rownames(m_DGE_genes_p_value)))])
      order[c(8,9)] = order[c(9,8)]
    }
    if(c==2){
      order = rev(rownames(m_DGE_genes_p_value)[c(grep("unmutated",rownames(m_DGE_genes_p_value)),
                                                  grep("_mutated",rownames(m_DGE_genes_p_value)),
                                                  grep("switched",rownames(m_DGE_genes_p_value)))])
    }
    if(c==3){
      order = rownames(m_DGE_genes_p_value)[c(grep("Mean_cluster_size_",rownames(m_DGE_genes_p_value), invert = T))]
      order1 = c(order[c(grep("Class",order, invert = F))], order[c(grep("_mutated",order, invert = F))], order[c(grep("_unmutated",order, invert = F))])
      order1 = order1[c(grep("_singleton",order1, invert = T))]
      order2 = setdiff(order, order1)
      order2 = order2[c(grep("expand",order2, invert = T))]
      order2 = order2[c(grep("_singleton",order2, invert = T))]
      order = rev(c(order1, order2))
      order = order[grep("by_cluster", order)]
    }
    
    matrix = m_DGE_genes_p_value[order, ]
    matrix_mean = m_DGE_genes_mean[order, ]
    matrix_FC = m_DGE_genes_FC[order, ]
    
    main = names(classes_uses)[c]
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
    
    fileout1=concat(c(out_dir,"Heatmap", main,"_differences_", batch,".pdf"))
    if(c %in% c(1,2)){
      w=6
      pdf(file=fileout1, height=w*0.8, width=w*1.85)
      y_range_max = 12
      x_start = 0.5
    }
    if(c %in% c(3)){
      w=9
      pdf(file=fileout1, height=w*0.8, width=w*2.5)
      y_range_max = 34
      x_start = 2.5
    }
    par(mfrow= c(1,2), mar = c(5,12,6,5))
    
    shift = y_range_max - max(y_range)
    plot(c(0.5,max(x_range)), c(-1,(length(y_range))+0.5), pch=20, col="white",xlab="",ylab ="",cex=cex, cex.lab=cex+0.1, cex.axis=cex,cex.main = cex, col.axis = "white",tck=0, mgp = c(2,0,0), main = "", 
         ylim = c(0, max(y_range_max)-0),xlim = c(x_start,max(x_range)+5.5), axes = F)
    mtext(side = 3, text = gsub(".","-",xnames, fixed = T), line = 0.01,cex=0.7,  at =c(1:length(xnames)), las = 2, font = 1)
    label = gsub("_"," ",gsub("Percentage","%",gsub("mean","Mean",ynames, fixed = T), fixed = T), fixed = T)
    label = gsub("J gene freq by uniq VDJ ","",gsub(".","/",label, fixed = T), fixed = T)
    label = gsub("V gene freq by cluster ","",label, fixed = T)
    label = gsub(" unique BCRs per isotype group"," isotype usage (overall)",label, fixed = T)
    label = gsub(" unique BCRs in switched"," isotype usage (switched B cells)",label, fixed = T)
    label = gsub("Cluster Gini Index","Clonal Diversification Index",label, fixed = T)
    label = gsub("Vertex Gini Index","Clonal Expansion Index",label, fixed = T)
    mtext(side = 2, text = label, line = 0.2,cex=0.7,  at =c(1:length(ynames))+shift, las = 1, font = 1)
    
    m_col = round(matrix, digits = 0)+1
    
    segments(0.5, c(1:length(ynames))+shift , max(x_range), c(1:length(ynames))+shift, col = "grey", lty = 3, lwd = 0.5)
    segments(c(1:length(xnames)),1+shift, c(1:length(xnames)) , max(y_range)+0.5+shift,  col = "grey", lty = 3, lwd = 0.5)
    size = matrix_mean^0.25
    size = 1.8*size/max(size)
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

}
Plot_heatmaps_of_differences(out_dir)

################################################ LDA to separate out patients from controls
Plot_overall_projection_on_significant_only<-function(out_dir){
  prefix = "Disease"
  file = concat(c(out_dir, prefix,"_statistical_summary_age_sex_adjusted.txt"))
  p <- as.matrix(read.csv(file, head=TRUE, sep="\t"))
  p=p[which(as.numeric(p[,3])<0.05), ]
  group = p[,1]
  class = p[,2]
  
  groups = sort(unique(group))
  names = apply(cbind(group, class), 1, paste, collapse = "..")
  which(names %in% colnames(mat_numeric))
  #names = exp_group_full[which(exp_group %in% groups),1]
  
  large_matrix = mat_numeric[,names]
  a = apply(large_matrix, 2, function(x){length(which(is.na(x)))})
  large_matrix = large_matrix[,which(a==0)]
  a = apply(large_matrix, 2, function(x){length(unique(x))})
  large_matrix = large_matrix[,which(a>=10)]
  a = apply(large_matrix, 2, max)
  large_matrix = large_matrix[,which(a>min(a))]
  
  large_matrix = large_matrix[,grep("V_gene_freq", colnames(large_matrix), invert = T)]
  large_matrix = large_matrix[,grep("Mean_cluster_size", colnames(large_matrix), invert = T)]
  large_matrix = large_matrix[,grep("VJ_gene", colnames(large_matrix), invert = T)]
  keep = colnames(large_matrix)[grep("J_gene_freq_by_uniq_VDJ_IGHD.IGHM", colnames(large_matrix))]
  remove = setdiff(colnames(large_matrix)[grep("J_gene_freq_by_uniq_VDJ", colnames(large_matrix))], keep)
  large_matrix = large_matrix[,which(colnames(large_matrix) %in% remove ==F)]
  large_matrix = large_matrix[,grep("V_gene_replacement_frequency", colnames(large_matrix), invert = T)]

  
  factor_sex = factor(Sex[rownames(large_matrix)])
  factor_age = Age[rownames(large_matrix)]
  factor = factor(Diagnosis)
  
  library(umap)
  library(labdsv)
  x <- pca(large_matrix,dim=10)
  pca_info = x
  dimred = x$scores
  loading = x$loading
  
  loading_rank = loading
  for(i in c(1:length(loading_rank[1,]))){
    loading_rank[, i] = rank(abs(loading[,i]))}
  
  overall_rank = apply(loading_rank, 1, min)
  out = cbind(loading, loading_rank, overall_rank)
  
  out_file_table = concat(c(out_dir, "PCA_loadings_",batch,".txt"))
  write.table(out, file = out_file_table, append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = TRUE,col.names = TRUE, qmethod = c("escape", "double"),fileEncoding = "")
  
  par(mfrow= c(1,2), mar = c(4,4,4,3))
  plot(dimred, pch = 21, bg = factor)
  
  pvals = NULL
  for(i in c(1:length(dimred[1,]))){
    groups = c(list(dimred[which(factor ==levels(factor)[1]),i]), list(dimred[which(factor ==levels(factor)[2]),i]))
    names(groups) = levels(factor)
    pval = wilcox.test(groups[[1]], groups[[2]])$p.value
    pvals = c(pvals, pval)}
  
  w_plot = which(pvals<0.03)
  w_plot1 = which(pvals<0.05)
  
  library(RColorBrewer)
  cols1 =  add.alpha (c(brewer.pal(8, "Dark2")[c(3,6)],"grey"), alpha = 0.95)
  cols =  add.alpha (cols1, alpha = 0.5)
  cols2 =  add.alpha (cols, alpha = 0.5)	
  
  ### plot PCA values between groups
  groups = NULL
  pvals =NULL
  for(i in c(1:length(dimred[1,]))){
    g1 = c(list(dimred[which(factor ==levels(factor)[1]),i]), list(dimred[which(factor ==levels(factor)[2]),i]))
    names(g1) = levels(factor)
    pvals = c(pvals, wilcox.test(g1[[1]], y=g1[[2]])$p.value)
    groups = c(groups, list(g1))
  }
  names(groups)=colnames(dimred)
  names(pvals)=colnames(dimred)
  
  sort(pvals)
  
  fileout1=concat(c(out_dir,"PCA_differences2_", batch,".pdf"))
  w=4
  pdf(file=fileout1, height=w*1, width=w*2.7)
  par(mfrow= c(1,2), mar = c(5,5,5,5))
  
  factors1 = names(groups[1])
  factors = names(groups)
  main = "PCAs"
  max = max(c(unlist(groups), unlist(groups))*1.35)
  min = min(c(unlist(groups), unlist(groups))*1.15)
  b = (max-min)*0.034
  ylab = ""
  draw_signif_lines = TRUE
  y = max(c(unlist(groups), unlist(groups))*1)+b
  max_width = length(groups)
  max_scale = min(c(max,100))
  range = max-min
  if(range>50){scale = c(-100:100)*20}
  if(range>100){scale = c(-100:100)*50}
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
  mtext(side = 1, text = gsub("_","/",factors,fixed=T), line = 0.35,cex= cex-0.1,  at = c(1:length(factors)), las = 2, font = 1)
  segments(0.5,Fun(scale),length(groups)+0.5,Fun(scale),col = "grey",lwd = 1,lty = 3 )
  segments(0.5,Fun(scale),length(groups)+0.5,Fun(scale),col = "grey",lwd = 1,lty = 3 )
  mtext(side = 2, text = scale, line = 0.15,cex= cex-0.1,  at =Fun(scale), las = 2, font = 1)
  width = 0.18
  index = 1
  l = length(groups)
  l1 = length(groups[[1]])
  shift = c(1:l1)
  shift = (mean(shift)-shift)
  shift = shift*0.25/max(shift)
  
  for(i in c(1:l)){
    for(i1 in c(1:l1)){
      points1=as.numeric(groups[[i]][[i1]])
      box1<-c(as.numeric(quantile(points1))[3], as.numeric(quantile(points1, probs = c(0.1, 0.9))), as.numeric(quantile(points1))[c(2, 4)])	
      Draw_box_plot(box1,i-shift[i1],width,cols[i1],1, cols1[i1])
      points(rep(i-shift[i1], length(points1)),points1, pch =21, col=cols[i1],bg = cols[i1], cex = 0.7)
    }}
  p_value = pvals
  for(i in c(1:l)){	
    b = max*0.035
    signif_threshold = 0.05
    pval1= "NS"
    if(p_value[i]<signif_threshold){pval1 = "*"
    y = max(unlist(groups[[i]]))
    y = y+3*b
    text(i, y+2*b, labels = pval1, cex = 1.3)
    }}
  
  plot(c(0,1), c(0,1), pch = 21, col = "white", bg ="white", xlab = "", ylab = "", main = "", axes=F)
  legend("topleft", as.character(Sources), pch = 21,cex= 0.8, bty="n", pt.bg = cols, col = cols, pt.lwd = 2, text.font = 2)
  
  dev.off()
  
  
  
  fileout1=concat(c(out_dir,"PCA_differences_", batch,".pdf"))
  w=6
  pdf(file=fileout1, height=w*1, width=w*1)
  par(mfrow= c(2,2), mar = c(5,5,5,5))
  
  w_plot = c(1,2)
  i = w_plot[1]
  groups = c(list(dimred[which(factor ==levels(factor)[1]),i]), list(dimred[which(factor ==levels(factor)[2]),i]))
  names(groups) = levels(factor)
  pval = wilcox.test(groups[[1]], groups[[2]])$p.value
  main = concat(c("PCA",i,"\np-value:",signif(pval, digits = 3)))
  colsx1 = cols1
  colsx =  cols
  width_plot = 5
  
  Boxplot_custom(groups, main, width_plot, colsx)
  
  
  Get_confidence_interval_plot<-function(x,y){
    data <- as.data.frame(x = x, y = y)
    p <- ggplot(data, aes(x, y)) +
      geom_point(color = "blue", alpha = 0.5) +
      stat_ellipse(level = 0.85, color = "red", size = 1.2) +  # 95% CI
      theme_minimal()
    # Extract the ellipse coordinates
    ellipse_data <- ggplot_build(p)$data[[2]]  # For 95% CI
    contour_data = ellipse_data
    return(contour_data)
  }  
  dimred1 = x$scores[,c(w_plot[1], w_plot[2])]
  plot(dimred1, pch = 21,  col = cols1[factor], bg = cols1[factor])
  
  for(i in c(1:length(levels(factor)))){
    w = which(factor == levels(factor)[i])
    x = dimred1[w,1]
    y = dimred1[w,2]
    contour_data = Get_confidence_interval_plot(x,y)
    polygon(contour_data$x, contour_data$y, col = add.alpha(cols1[i], 0.25), border = cols1[i], lwd = 2)
  }
  
  for(i in c(1:length(levels(factor)))){
    w = which(factor == levels(factor)[i])
    x = dimred1[w,1]
    y = dimred1[w,2]
    dataEllipse(x, y, levels=c(0.45), add = TRUE, col = cols1[i], center.cex = 0.0001, fill =T, fill.alpha =0.1 , plot.points = FALSE)
  }
  
  i = w_plot1[2]
  groups = c(list(dimred[which(factor ==levels(factor)[1]),i]), list(dimred[which(factor ==levels(factor)[2]),i]))
  names(groups) = levels(factor)
  pval = wilcox.test(groups[[1]], groups[[2]])$p.value
  main = concat(c("PCA",i,"\np-value:",signif(pval, digits = 3)))
  colsx1 = cols1
  colsx =  cols
  width_plot = 5
  
  Boxplot_custom(groups, main, width_plot, colsx)
  
  plot(summary(pca_info)[2,], xlab = "PCA dimension", ylab = "Proportion of \nvariance explained", pch = 21, col = add.alpha("red", 0.5),
       bg = add.alpha("red", 0.5), xlim = c(0,length(summary(pca_info)[2,])))
  points(summary(pca_info)[2,],  type= "l", lwd = 2, col = add.alpha("red", 0.5))
  
  dev.off()
  
  ### plot contribution of features to PCA
  library(pheatmap)
  t = rev(sort(abs(apply(loading, 1, max))))#[40]
  
  w1 = which(apply(loading, 1, max)>=t)
  plot(sort(abs(apply(loading[w1,], 1, max))))
  plot(sort(abs(apply(loading, 1, max))))
  
  fileout1=concat(c(out_dir,"PCA_contributions_", batch,".pdf"))
  w=5
  pdf(file=fileout1, height=w*0.8, width=w*1.6)
  par(mfrow= c(1,1), mar = c(5,5,5,5))
  pheatmap(loading[w1,], 
           cluster_rows = TRUE, 
           cluster_cols = TRUE, 
           width = 10,                 # Adjust the width of the heatmap
           height = 10,    
           main = "Contribution of Features to PCA Dimensions",
           display_numbers = TRUE)
  dev.off()
  
  ### linear classification
  Linear_classification_CV<-function(){
    data <- as.data.frame(large_matrix)
    data = as.data.frame(dimred[,c(1:2)])
    data <- scale(data)
    data <- as.data.frame(as.matrix(data[,-length(data[1,])]))
    data$Class <- as.factor(factor)
    
    n <- nrow(data)
    predicted <- rep(0, n)
    actual <- data$Class
    
    # LOOCV loop
    library(randomForest)
    for (i in 1:n) {
      train_data <- data[-i, ]
      test_data <- data[i, ]
      model <- glm(Class ~ ., data = train_data, family = binomial)
      prob <- predict(model, newdata = test_data, type = "response")
      predicted[i] <- ifelse(prob > 0.5, 1, 0)
      #model <- randomForest(Class ~ ., data = train_data, ntree = 500)
      #prob <- predict(model, newdata = test_data, type = "response")
      #predicted[i] <- as.character(prob)
    }
    
    # Convert to factors for confusion matrix
    predicted <- factor(levels(actual)[factor(predicted, levels = c(0, 1))])
    #predicted <- factor(predicted)
    # Calculate confusion matrix
    library(caret)
    cm <- confusionMatrix(predicted, actual)
    
    # Extract performance metrics
    accuracy <- cm$overall["Accuracy"]
    sensitivity <- cm$byClass["Sensitivity"]
    specificity <- cm$byClass["Specificity"]
    
    # Display results
    cat("Accuracy:", accuracy, "\n")
    cat("Sensitivity:", sensitivity, "\n")
    cat("Specificity:", specificity, "\n")
    
  }
  
  #Supervised Partial Least Squares Discriminant Analysis (PLS-DA)
  PLS_DA<-function(){
    library(caret)      # For model training and evaluation
    library(pls)         # PLS regression
    X <- large_matrix
    Y <- as.factor(factor)
    
    # PLS-DA model
    plsda_model <- plsda(X, Y, ncomp = 2)  # ncomp = number of components
    
    library(ggplot2)    # For visualization
    library(dplyr)      # For data manipulation
    
    # Extract PLS scores
    scores <- as.data.frame(as.matrix(plsda_model$scores)[,c(1:2)])
    colnames(scores) <- c("PLS1", "PLS2")  # Rename for clarity
    scores$Class <- Y  # Add class labels
    
    ggplot(scores, aes(x = PLS1, y = PLS2, color = Class)) +
      geom_point(size = 4, alpha = 0.8) +
      labs(title = "PLS-DA Score Plot", x = "PLS Component 1", y = "PLS Component 2") +
      theme_minimal() +
      scale_color_manual(values = c("orange", "purple")) +
      theme(legend.position = "right",
            plot.title = element_text(hjust = 0.5, face = "bold", size = 16))
    
    pred <- predict(plsda_model, X)
    
    # Confusion matrix
    conf_matrix <- confusionMatrix(as.factor(pred), Y)
    print(conf_matrix)
    
    # Extract accuracy, sensitivity, and specificity
    accuracy <- conf_matrix$overall["Accuracy"]
    sensitivity <- conf_matrix$byClass["Sensitivity"]
    specificity <- conf_matrix$byClass["Specificity"]
    
    cat("Accuracy:", accuracy, "\n")
    cat("Sensitivity:", sensitivity, "\n")
    cat("Specificity:", specificity, "\n")
    
    }
  
  data <- as.data.frame(large_matrix)
  data$Class <- as.factor(factor)
  # Fit logistic regression
  model <- glm(Class ~ ., data = data, family = binomial)

  # Model summary
  summary(model)
  
  # Make predictions
  predictions <- predict(model, data, type = "response")
  # Convert to binary class labels
  class_pred <- ifelse(predictions > 0.5, 1, 0)
  
  # Confusion matrix
  table(Predicted = class_pred, Actual = data$Class)
  
  # Make predictions
  
  # Plot predicted probabilities
  ggplot(data, aes(x = prob, fill = class)) +
    geom_histogram(bins = 10, alpha = 0.7, position = "identity") +
    #scale_fill_manual(values = c("orange", "purple")) +
    labs(title = "Predicted Probabilities by Class",
         x = "Predicted Probability",
         y = "Count") +
    theme_minimal()
  
  data$prob <- predict(model, type = "response")
  
  lda_df <- data.frame(LD1 = data$prob, Species = factor)
  
  fileout1=concat(c(out_dir,"GLM_separation_", batch,".pdf"))
  w=2.5
  pdf(file=fileout1, height=w*0.7, width=w*1.2)
  par(mfrow= c(1,1), mar = c(5,5,5,5))
  # Plot using ggplot2
  ggplot(lda_df, aes(x = LD1, fill = factor, color = factor)) +
    geom_density(binwidth = 0.3,alpha = 0.5) +
    scale_fill_manual(values = c("FMS" = "orange", "Control" = "purple")) +
    scale_color_manual(values = c("FMS" = "orange", "Control" = "purple")) +
    labs(title = "Density Plot of LDA Dimension 1",
         x = "LDA Dimension 1",
         y = "Density") +
    theme_minimal()
  dev.off()
  
  # confusion matrix
  fileout1=concat(c(out_dir,"GLM_confusion_matrix_", batch,".pdf"))
  w=1.4
  pdf(file=fileout1, height=w*1.1, width=w*1.7)
  par(mfrow= c(1,1), mar = c(5,5,5,5))
  class_pred[which(class_pred==1)] = "FMS"
  class_pred[which(class_pred==0)] = "Control"
  cm_df <- as.data.frame(  table(Predicted = class_pred, Actual = data$Class))
  ggplot(data = cm_df, aes(x = Predicted, y = Actual, fill = Freq)) +
    geom_tile() +
    geom_text(aes(label = Freq), color = "black", size = 5) +
    scale_fill_gradient(low = "white", high = "orange") +
    labs(title = "Confusion Matrix",
         x = "Predicted Label",
         y = "True Label") +
    theme_minimal()
  dev.off()
  
  # Order by the absolute value of LD1
  ## get variable importance
  z_values <- abs(summary(model)$coefficients[, "z value"])
  importance <- sort(z_values, decreasing = TRUE)
  lda_df <- as.data.frame(importance[which(names(importance)!="(Intercept)")])
  colnames(lda_df) = "importance"
  lda_df$feature = rownames(lda_df)
  lda_df$feature = gsub("J_gene_freq_by_uniq_VDJ_", "% J gene ", lda_df$feature, fixed = T)
  lda_df$feature = gsub("Percentage_unique_BCRs_in_switched", "% isotype (switched B cells)", lda_df$feature, fixed = T)
  lda_df$feature = gsub("Percentage_unique_BCRs_per_isotype_group", "% isotype (overall)", lda_df$feature, fixed = T)
  lda_df$feature = gsub("Mean_mutations_per_BCR", "SHM", lda_df$feature, fixed = T)
  lda_df$feature = gsub("Percentage_unmutated", "% unmutated BCRs", lda_df$feature, fixed = T)
  lda_df$feature = gsub("..", " ", lda_df$feature, fixed = T)
  lda_df$feature = gsub("_", " ", lda_df$feature, fixed = T)
  lda_df$feature = gsub("IGHD.IGHM", "IGHD/M", lda_df$feature, fixed = T)
  lda_df = lda_df[order(lda_df$importance),]
  fileout1=concat(c(out_dir,"GLM_contributions_", batch,".pdf"))
  w=4
  pdf(file=fileout1, height=w*1.2, width=w*1.2)
  par(mfrow= c(1,1), mar = c(5,5,5,5))
  ggplot(lda_df, aes(x = reorder(feature, abs(importance)), y = importance)) +
    geom_bar(stat = "identity", fill = "skyblue") +
    coord_flip() +  # Flip coordinates to make it horizontal
    labs(title = "Contribution to GLM",
         x = "Feature",
         y = "GLM feature importance") +
    theme(plot.margin = margin(2,2,2,2, "cm"))+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  dev.off()
  
  
  # Generate ROC curve
  roc_curve <- roc(data$Class, data$prob)
  
  # Plot ROC curve
  plot(roc_curve, col = "blue", lwd = 2, main = "ROC Curve")
  abline(a = 0, b = 1, lty = 2, col = "red")
  auc(roc_curve)
 
  
  
  
  
  
  library(MASS)
  lda_model <- lda(factor ~ large_matrix)
  print(lda_model)
  plot(rev(sort(abs(lda_model$ scaling)))[1:100])
  
  filtered_large_matrix = large_matrix[,which(abs(lda_model$ scaling)> 0.75)]
  dim(filtered_large_matrix)
  colnames(filtered_large_matrix)
  lda_model <- lda(factor ~ filtered_large_matrix)
  
  
  predictions <- predict(lda_model, factor)
  lda_values <- predictions$x
  # Confusion matrix
  confusion_matrix = table(Predicted = predictions$class, Actual = factor)
  confusion_matrix
  # Calculate the accuracy
  accuracy <- mean(predictions$class == factor)
  print(paste("Accuracy:", accuracy))
  
  colors <- c("FMS" = "red", "Control" = "darkgreen")
  
  # Plot the LDA values
  plot(lda_values[,1], lda_values[,1], col=colors[factor],
       pch=19, xlab="Linear Discriminant 1", ylab="Linear Discriminant 2",
       main="LDA of FMS Dataset")
  
  # Add a legend
  legend("bottomright", legend=levels(factor), col=colors, pch=19)
  
  ##### plot histograms
  library(RColorBrewer)
  cols1 =  add.alpha (brewer.pal(8, "Dark2")[c(3,6)], alpha = 0.95)
  cols =  add.alpha (cols1, alpha = 0.5)
  
  library(ggplot2)
 
  
 
  lda_coefficients = lda_model$ scaling
  rownames(lda_coefficients) = gsub("filtered_large_matrix","",gsub("..",": ",gsub("Percentage","%",rownames(lda_coefficients), fixed = T), fixed = T), fixed = T)
  rownames(lda_coefficients) = gsub("_"," ",rownames(lda_coefficients), fixed = T)
  rownames(lda_coefficients) = gsub("IGHD.M","IGHD/M",rownames(lda_coefficients), fixed = T)
  rownames(lda_coefficients) = gsub(".","-",rownames(lda_coefficients), fixed = T)
  rownames(lda_coefficients) = gsub("IGHD IGHM","IGHD/M",rownames(lda_coefficients), fixed = T)
  # Convert to data frame
  lda_df <- as.data.frame(lda_coefficients)
  
  # Add feature names
  lda_df$Feature <- rownames(lda_df)
  
  # Select only the first LDA dimension
  lda_df <- lda_df[, c("LD1", "Feature")]
  
  
  
}

Plot_overall_projection<-function(out_dir){
  large_matrix = NULL
  for(ind in c(1:length(groups))){
    w_analysis = exp_group_full[which(exp_group== groups[ind]),]
    mat = mat_numeric[,w_analysis[,1]]
    mat[which(mat==-1)] = NA
    a = apply(mat, 2, function(x){length(which(is.na(x)))})
    mat = mat[,which(a<10)]
    if(length(unique(sort(mat)))>5){
      mat = mat[unlist(groups_ids),]
      if(length(large_matrix)==0){large_matrix = mat
      }else{
        large_matrix = cbind(large_matrix, mat)
      }
    }
  }
  
  
  a = apply(large_matrix, 2, function(x){length(which(is.na(x)))})
  large_matrix = large_matrix[,which(a==0)]
  a = apply(large_matrix, 2, function(x){length(unique(x))})
  large_matrix = large_matrix[,which(a>=10)]
  a = apply(large_matrix, 2, max)
  large_matrix = large_matrix[,which(a>min(a))]
  
  large_matrix = large_matrix[,grep("V_gene_freq", colnames(large_matrix), invert = T)]
  large_matrix = large_matrix[,grep("Mean_cluster_size", colnames(large_matrix), invert = T)]
  large_matrix = large_matrix[,grep("VJ_gene", colnames(large_matrix), invert = T)]
  #large_matrix = large_matrix[,grep("J_gene_freq_by_uniq_VDJ", colnames(large_matrix), invert = T)]
  large_matrix = large_matrix[,grep("J_gene_freq_by_uniq_VDJ_IGH", colnames(large_matrix), invert = T)]
  large_matrix = large_matrix[,grep("mean_CDR3_charge", colnames(large_matrix), invert = T)]
  large_matrix = large_matrix[,grep("V4_34_AVY_unmut", colnames(large_matrix), invert = T)]
  large_matrix = large_matrix[,grep("V4_34_NHS_unmut", colnames(large_matrix), invert = T)]
  
    
  factor_sex = factor(Sex[rownames(large_matrix)])
  factor_age = Age[rownames(large_matrix)]
  factor = factor(Diagnosis)
  
  library(umap)
  library(labdsv)
  x <- pca(large_matrix,dim=10)
  pca_info = x
  dimred = x$scores
  loading = x$loading
  
  loading_rank = loading
  for(i in c(1:length(loading_rank[1,]))){
    loading_rank[, i] = rank(abs(loading[,i]))}
  
  overall_rank = apply(loading_rank, 1, min)
  out = cbind(loading, loading_rank, overall_rank)
  
  out_file_table = concat(c(out_dir, "PCA_loadings_",batch,".txt"))
  write.table(out, file = out_file_table, append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = TRUE,col.names = TRUE, qmethod = c("escape", "double"),fileEncoding = "")
  
  
  
  par(mfrow= c(1,2), mar = c(4,4,4,3))
  plot(dimred, pch = 21, bg = factor)
  
  pvals = NULL
  for(i in c(1:length(dimred[1,]))){
    groups = c(list(dimred[which(factor ==levels(factor)[1]),i]), list(dimred[which(factor ==levels(factor)[2]),i]))
    names(groups) = levels(factor)
    pval = wilcox.test(groups[[1]], groups[[2]])$p.value
    pvals = c(pvals, pval)}
  
  w_plot = which(pvals<0.03)
  w_plot1 = which(pvals<0.05)
  
  library(RColorBrewer)
  cols1 =  add.alpha (c(brewer.pal(8, "Dark2")[c(3,6)],"grey"), alpha = 0.95)
  cols =  add.alpha (cols1, alpha = 0.5)
  cols2 =  add.alpha (cols, alpha = 0.5)	
  
  ### plot PCA values between groups
  groups = NULL
  pvals =NULL
  for(i in c(1:length(dimred[1,]))){
    g1 = c(list(dimred[which(factor ==levels(factor)[1]),i]), list(dimred[which(factor ==levels(factor)[2]),i]))
    names(g1) = levels(factor)
    pvals = c(pvals, wilcox.test(g1[[1]], y=g1[[2]])$p.value)
    groups = c(groups, list(g1))
  }
  names(groups)=colnames(dimred)
  names(pvals)=colnames(dimred)
  
  sort(pvals)
  
  fileout1=concat(c(out_dir,"PCA_differences2_", batch,".pdf"))
  w=4
  pdf(file=fileout1, height=w*1, width=w*2.7)
  par(mfrow= c(1,2), mar = c(5,5,5,5))
  
  factors1 = names(groups[1])
  factors = names(groups)
  main = "PCAs"
  max = max(c(unlist(groups), unlist(groups))*1.35)
  min = min(c(unlist(groups), unlist(groups))*1.15)
  b = (max-min)*0.034
  ylab = ""
  draw_signif_lines = TRUE
  y = max(c(unlist(groups), unlist(groups))*1)+b
  max_width = length(groups)
  max_scale = min(c(max,100))
  range = max-min
  if(range>50){scale = c(-100:100)*20}
  if(range>100){scale = c(-100:100)*50}
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
  mtext(side = 1, text = gsub("_","/",factors,fixed=T), line = 0.35,cex= cex-0.1,  at = c(1:length(factors)), las = 2, font = 1)
  segments(0.5,Fun(scale),length(groups)+0.5,Fun(scale),col = "grey",lwd = 1,lty = 3 )
  segments(0.5,Fun(scale),length(groups)+0.5,Fun(scale),col = "grey",lwd = 1,lty = 3 )
  mtext(side = 2, text = scale, line = 0.15,cex= cex-0.1,  at =Fun(scale), las = 2, font = 1)
  width = 0.18
  index = 1
  l = length(groups)
  l1 = length(groups[[1]])
  shift = c(1:l1)
  shift = (mean(shift)-shift)
  shift = shift*0.25/max(shift)
  
  for(i in c(1:l)){
    for(i1 in c(1:l1)){
      points1=as.numeric(groups[[i]][[i1]])
      box1<-c(as.numeric(quantile(points1))[3], as.numeric(quantile(points1, probs = c(0.1, 0.9))), as.numeric(quantile(points1))[c(2, 4)])	
      Draw_box_plot(box1,i-shift[i1],width,cols[i1],1, cols1[i1])
      points(rep(i-shift[i1], length(points1)),points1, pch =21, col=cols[i1],bg = cols[i1], cex = 0.7)
    }}
  p_value = pvals
  for(i in c(1:l)){	
    b = max*0.035
    signif_threshold = 0.05
    pval1= "NS"
    if(p_value[i]<signif_threshold){pval1 = "*"
    y = max(unlist(groups[[i]]))
    y = y+3*b
    text(i, y+2*b, labels = pval1, cex = 1.3)
    }}
  
  plot(c(0,1), c(0,1), pch = 21, col = "white", bg ="white", xlab = "", ylab = "", main = "", axes=F)
  legend("topleft", as.character(Sources), pch = 21,cex= 0.8, bty="n", pt.bg = cols, col = cols, pt.lwd = 2, text.font = 2)
  
  dev.off()
  
  
  
  fileout1=concat(c(out_dir,"PCA_differences_", batch,".pdf"))
  w=6
  pdf(file=fileout1, height=w*1, width=w*1)
  par(mfrow= c(2,2), mar = c(5,5,5,5))
  
  i = w_plot[1]
  groups = c(list(dimred[which(factor ==levels(factor)[1]),i]), list(dimred[which(factor ==levels(factor)[2]),i]))
  names(groups) = levels(factor)
  pval = wilcox.test(groups[[1]], groups[[2]])$p.value
  main = concat(c("PCA",i,"\np-value:",signif(pval, digits = 3)))
  colsx1 = cols1
  colsx =  cols
  width_plot = 5
  
  Boxplot_custom(groups, main, width_plot, colsx)
  
  dimred1 = x$scores[,c(w_plot[1], w_plot1[2])]
  plot(dimred1, pch = 21,  col = cols1[factor], bg = cols1[factor])
  
  library(ellipse)
  for(i in c(1:length(levels(factor)))){
    w = which(factor == levels(factor)[i])
    x = dimred1[w,1]
    y = dimred1[w,2]
    dataEllipse(x, y, levels=c(0.45), add = TRUE, col = cols1[i], center.cex = 0.0001, fill =T, fill.alpha =0.1 , plot.points = FALSE)
  }
  
  library(car)
  for(i in c(1:length(levels(factor)))){
    w = which(factor == levels(factor)[i])
    x = dimred1[w,1]
    y = dimred1[w,2]
    dataEllipse(x, y, levels=c(0.45), add = TRUE, col = cols1[i], center.cex = 0.0001, fill =T, fill.alpha =0.1 , plot.points = FALSE)
  }
  
  i = w_plot1[2]
  groups = c(list(dimred[which(factor ==levels(factor)[1]),i]), list(dimred[which(factor ==levels(factor)[2]),i]))
  names(groups) = levels(factor)
  pval = wilcox.test(groups[[1]], groups[[2]])$p.value
  main = concat(c("PCA",i,"\np-value:",signif(pval, digits = 3)))
  colsx1 = cols1
  colsx =  cols
  width_plot = 5
  
  Boxplot_custom(groups, main, width_plot, colsx)
  
  plot(summary(pca_info)[2,], xlab = "PCA dimension", ylab = "Proportion of \nvariance explained", pch = 21, col = add.alpha("red", 0.5),
       bg = add.alpha("red", 0.5), xlim = c(0,length(summary(pca_info)[2,])))
  points(summary(pca_info)[2,],  type= "l", lwd = 2, col = add.alpha("red", 0.5))
  
  dev.off()
  
  ### plot contribution of features to PCA
  library(pheatmap)
  t = rev(sort(abs(apply(loading, 1, max))))[40]
  
  w1 = which(apply(loading, 1, max)>=t)
  plot(sort(abs(apply(loading[w1,], 1, max))))
  plot(sort(abs(apply(loading, 1, max))))
  
  fileout1=concat(c(out_dir,"PCA_contributions_", batch,".pdf"))
  w=6
  pdf(file=fileout1, height=w*1.1, width=w*1.6)
  par(mfrow= c(1,1), mar = c(5,5,5,5))
  pheatmap(loading[w1,], 
           cluster_rows = TRUE, 
           cluster_cols = TRUE, 
           width = 10,                 # Adjust the width of the heatmap
           height = 10,    
           main = "Contribution of Features to PCA Dimensions",
           display_numbers = TRUE)
  dev.off()
  
  ## run linear discriminant analysis
  
  library(MASS)
  lda_model <- lda(factor ~ large_matrix)
  print(lda_model)
  plot(rev(sort(abs(lda_model$ scaling)))[1:100])
  
  filtered_large_matrix = large_matrix[,which(abs(lda_model$ scaling)> 0.75)]
  dim(filtered_large_matrix)
  colnames(filtered_large_matrix)
  lda_model <- lda(factor ~ filtered_large_matrix)
  
  
  predictions <- predict(lda_model, factor)
  lda_values <- predictions$x
  # Confusion matrix
  confusion_matrix = table(Predicted = predictions$class, Actual = factor)
  confusion_matrix
  # Calculate the accuracy
  accuracy <- mean(predictions$class == factor)
  print(paste("Accuracy:", accuracy))
  
  colors <- c("FMS" = "red", "Control" = "darkgreen")
  
  # Plot the LDA values
  plot(lda_values[,1], lda_values[,1], col=colors[factor],
       pch=19, xlab="Linear Discriminant 1", ylab="Linear Discriminant 2",
       main="LDA of FMS Dataset")
  
  # Add a legend
  legend("bottomright", legend=levels(factor), col=colors, pch=19)
  
  ##### plot histograms
  library(RColorBrewer)
  cols1 =  add.alpha (brewer.pal(8, "Dark2")[c(3,6)], alpha = 0.95)
  cols =  add.alpha (cols1, alpha = 0.5)
  
  library(ggplot2)
  lda_df <- data.frame(LD1 = lda_values[, 1], Species = factor)
  
  fileout1=concat(c(out_dir,"LDA_separation_", batch,".pdf"))
  w=2.5
  pdf(file=fileout1, height=w*0.7, width=w*1.2)
  par(mfrow= c(1,1), mar = c(5,5,5,5))
  # Plot using ggplot2
  ggplot(lda_df, aes(x = LD1, fill = factor, color = factor)) +
    geom_density(binwidth = 0.3,alpha = 0.5) +
    scale_fill_manual(values = c("FMS" = "orange", "Control" = "purple")) +
    scale_color_manual(values = c("FMS" = "orange", "Control" = "purple")) +
    labs(title = "Density Plot of LDA Dimension 1",
         x = "LDA Dimension 1",
         y = "Density") +
    theme_minimal()
  dev.off()
  
  # confusion matrix
  fileout1=concat(c(out_dir,"LDA_confusion_matrix_", batch,".pdf"))
  w=1.7
  pdf(file=fileout1, height=w*1.1, width=w*1.6)
  par(mfrow= c(1,1), mar = c(5,5,5,5))
  cm_df <- as.data.frame(confusion_matrix)
  ggplot(data = cm_df, aes(x = Predicted, y = Actual, fill = Freq)) +
    geom_tile() +
    geom_text(aes(label = Freq), color = "black", size = 5) +
    scale_fill_gradient(low = "white", high = "orange") +
    labs(title = "Confusion Matrix",
         x = "Predicted Label",
         y = "True Label") +
    theme_minimal()
  dev.off()
  
  lda_coefficients = lda_model$ scaling
  rownames(lda_coefficients) = gsub("filtered_large_matrix","",gsub("..",": ",gsub("Percentage","%",rownames(lda_coefficients), fixed = T), fixed = T), fixed = T)
  rownames(lda_coefficients) = gsub("_"," ",rownames(lda_coefficients), fixed = T)
  rownames(lda_coefficients) = gsub("IGHD.M","IGHD/M",rownames(lda_coefficients), fixed = T)
  rownames(lda_coefficients) = gsub(".","-",rownames(lda_coefficients), fixed = T)
  rownames(lda_coefficients) = gsub("IGHD IGHM","IGHD/M",rownames(lda_coefficients), fixed = T)
  # Convert to data frame
  lda_df <- as.data.frame(lda_coefficients)
  
  # Add feature names
  lda_df$Feature <- rownames(lda_df)
  
  # Select only the first LDA dimension
  lda_df <- lda_df[, c("LD1", "Feature")]
  
  # Order by the absolute value of LD1
  lda_df <- lda_df[order(abs(lda_df$LD1), decreasing = TRUE), ]
  
  fileout1=concat(c(out_dir,"LDA_contributions_", batch,".pdf"))
  w=3.5
  pdf(file=fileout1, height=w*0.9, width=w*1.6)
  par(mfrow= c(1,1), mar = c(5,5,5,5))
  ggplot(lda_df, aes(x = reorder(Feature, abs(LD1)), y = LD1)) +
    geom_bar(stat = "identity", fill = "skyblue") +
    coord_flip() +  # Flip coordinates to make it horizontal
    labs(title = "Contribution to LDA 1",
         x = "Feature",
         y = "LDA 1 Coefficient") +
    theme_minimal()
  dev.off()
  
}

SVM<-function(){
  # Run an SVM on this using LOOCV
  
  library(gmodels)
  
  library(e1071)
  model <-svm(dimred1,factor, type = "C-classification",kernel = "linear")
  pred <- predict(model, dimred1)
  
  err_metric=function(CM){
    TN =CM[1,1]
    TP =CM[2,2]
    FP =CM[1,2]
    FN =CM[2,1]
    precision =(TP)/(TP+FP)
    recall_score =(FP)/(FP+TN)
    f1_score=2*((precision*recall_score)/(precision+recall_score))
    accuracy_model  =(TP+TN)/(TP+TN+FP+FN)
    False_positive_rate =(FP)/(FP+TN)
    False_negative_rate =(FN)/(FN+TP)
    sensitivity = (TP)/(TP+FN)
    print(paste("Precision value of the model: ",round(precision,2)))
    print(paste("Accuracy of the model: ",round(accuracy_model,2)))
    print(paste("Recall value of the model: ",round(recall_score,2)))
    print(paste("False Positive rate of the model: ",round(False_positive_rate,2)))
    print(paste("False Negative rate of the model: ",round(False_negative_rate,2)))
    print(paste("f1 score of the model: ",round(f1_score,2)))
    a = c(precision, accuracy_model, sensitivity, recall_score, False_positive_rate, False_negative_rate,f1_score )
    names(a) = c("precision", "accuracy_model", "sensitivity","recall_score", "False_positive_rate", "False_negative_rate","f1_score")
    return(a)
  }
  
  y.svm <- rep(NA, nrow(dimred1))
  model_fitness = NULL
  for (i in 1:nrow(dimred1)*2) {
    # 1. Logistic regression
    prop = 0.5
    sample_test = sample.int(length(dimred1[,1]), length(dimred1[,1])*prop)
    sample_train = setdiff(c(1:length(dimred1[,1])), sample_test)
    testset <- dimred1[sample_test,]
    trainset <- dimred1[sample_train,]
    
    logit_m <-svm(trainset,factor[sample_train], type = "C-classification",kernel = "linear")
    logit_P = predict(logit_m , newdata = testset ,type = 'response' )
    CM= table(factor[sample_test] , logit_P)
    mf = err_metric(CM)
    if(length(model_fitness)==0){model_fitness = mf
    }else{model_fitness= rbind(model_fitness, mf)}
  }
  #### compare to shuffled
  model_fitness1 = NULL
  for (i in 1:nrow(dimred1)*2) {
    # 1. Logistic regression
    prop = 0.5
    sample_test = sample.int(length(dimred1[,1]), length(dimred1[,1])*prop)
    sample_train = setdiff(c(1:length(dimred1[,1])), sample_test)
    testset <- dimred1[sample_test,]
    trainset <- dimred1[sample_train,]
    
    logit_m <-svm(trainset,sample(factor[sample_train], length(sample_train), replace = F), type = "C-classification",kernel = "linear")
    logit_P = predict(logit_m , newdata = testset ,type = 'response' )
    CM= table(factor[sample_test] , logit_P)
    mf = err_metric(CM)
    if(length(model_fitness1)==0){model_fitness1 = mf
    }else{model_fitness1= rbind(model_fitness1, mf)}
  }

  w = c(1:3)
  model_fitness2 = cbind(model_fitness[,w], model_fitness1[,w])
  groups = NULL
  for(i in c(1:length(model_fitness2[1,]))){
    x = model_fitness2[,i]
    x = x[which(is.na(x)==F)]
    groups = c(groups, list(x))
  }
  names(groups) = colnames(model_fitness2)
  boxplot(groups)
  
    
  for (i in 1:nrow(dimred1)) {
    testset <- dimred1[i,]
    trainset <- dimred1[-i,]
    model <-svm(trainset,factor[-i], type = "C-classification",kernel = "linear")
    y.svm[i] <- as.character(predict(model, rbind(testset)))
  }
  
  #basic R solution
  table(y.svm,factor)
  
  #Output similar to what users of SPSS or SAS expects
  CrossTable(y.svm,factor)
  
  ## draw ROC curve
  testset <- dimred1[i,]
  trainset <- dimred1[-i,]
  
  library(e1071)
  library(readxl)
  library(caret)
  
  class1.svm.model <- svm(Class ~ ., data = class1.trainset,cost=20, cross=10,type="C-classification",kernel="radial",na.action=na.omit)
  class1.svm.pred <- predict(class1.svm.model, class1.testset)
  finalmatrix<-data.matrix(class1.svm.pred, rownames.force = F)
  test<-table(pred = class1.svm.pred, true = class1.testset[,c(15768)])
  
  confusionMatrix(test)
  
  prop = 0.25
  sample_test = sample.int(length(dimred1[,1]), length(dimred1[,1])*prop)
  sample_train = setdiff(c(1:length(dimred1[,1])), sample_test)
  testset <- dimred1[sample_test,]
  trainset <- dimred1[sample_train,]
  
  library(pROC)
  
  roc_svm_test <- roc(response = factor[sample_test], predictor =as.numeric(class1.svm.pred))
  plot(roc_svm_test, add = TRUE,col = "red", print.auc=TRUE, print.auc.x = 0.5, print.auc.y = 0.3)
  legend(0.3, 0.2, legend = c("test-svm"), lty = c(1), col = c("blue"))
  
  
  
  
  ####
 
  # 1. Logistic regression
  prop = 0.5
  sample_test = sample.int(length(dimred1[,1]), length(dimred1[,1])*prop)
  sample_train = setdiff(c(1:length(dimred1[,1])), sample_test)
  testset <- dimred1[sample_test,]
  trainset <- dimred1[sample_train,]
  
  logit_m <-svm(trainset,factor[sample_train], type = "C-classification",kernel = "linear")
  #logit_m =glm(formula = default~. ,data =train_data ,family='binomial')
  summary(logit_m)
  logit_P = predict(logit_m , newdata = testset ,type = 'response' )
  CM= table(factor[sample_test] , logit_P)
  print(CM)
  err_metric(CM)
  
  

  
}
  
  ############
  
  
prop_test = 0.3
test_n = sample.int(length(factor), ceiling(length(factor)*prop_test))
train_n = setdiff(c(1:length(factor)), test_n)
table(factor[test_n])
table(factor[train_n])

test = large_matrix[test_n,]
train= large_matrix[train_n,]

library(MASS)
z <- lda(train, factor[train_n])
pred = as.character(predict(z, train)$class)
table(pred, factor[train_n])

pred = as.character(predict(z, test)$class)
table(pred, factor[test_n])

### plot LDA
z <- lda(large_matrix, factor, CV = T)
pred = as.character(predict(z, large_matrix)$class)
table(pred, factor)

z <- lda(large_matrix, factor)
x = z$ scaling

plot(sort(x))
n = 20
lower = colnames(large_matrix)[which(x<=sort(x)[n])]
upper = colnames(large_matrix)[which(x>=rev(sort(x))[n])]


large_matrix_small = large_matrix[,c(lower, upper)]
z <- lda(large_matrix_small, factor)
pred = as.character(predict(z, large_matrix_small)$class)
table(pred, factor)


z <- lda(large_matrix_small, factor, CV = T)
table(z$class, factor)
concat(c(length(which(z$class==factor))*100/length(factor)," %"))



library(umap)
library(labdsv)
x <- pca(large_matrix,dim=10)
dimred = x$scores
par(mfrow= c(1,2), mar = c(4,4,4,3))
plot(dimred, pch = 21, bg = factor)

groups = c(list(dimred[which(factor ==levels(factor)[1]),1]), list(dimred[which(factor ==levels(factor)[2]),1]))
names(groups) = levels(factor)
pval = wilcox.test(groups[[1]], groups[[2]])$p.value
boxplot(groups, main = signif(pval, digits = 3))

pvals = NULL
for(i in c(1:length(dimred[1,]))){
  groups = c(list(dimred[which(factor ==levels(factor)[1]),i]), list(dimred[which(factor ==levels(factor)[2]),i]))
  names(groups) = levels(factor)
  pval = wilcox.test(groups[[1]], groups[[2]])$p.value
  pvals = c(pvals, pval)}
  
i = 6
groups = c(list(dimred[which(factor ==levels(factor)[1]),i]), list(dimred[which(factor ==levels(factor)[2]),i]))
names(groups) = levels(factor)
pval = wilcox.test(groups[[1]], groups[[2]])$p.value
boxplot(groups, main = signif(pval, digits = 3))

loadings= x$ loadings

uppers = NULL
lowers = NULL
for(i in c(1:length(dimred[1,]))){
  load = loadings[,i]
  plot(sort(load))
  n = 10
  lower = names(load)[which(load<=sort(load)[n])]
  upper = names(load)[which(load>=rev(sort(load))[n])]
  uppers = c(uppers,list(upper))
  lowers = c(lowers,list(lower))
}
  
dimred = dimred[,which(pvals<=sort(pvals)[2])]

x <- pca(large_matrix,dim=10)
dimred = x$scores[,c(1,6)]
plot(dimred, pch = 21, bg = factor)

library(car)
for(i in c(1:length(levels(factor)))){
  w = which(factor == levels(factor)[i])
  x = dimred[w,1]
  y = dimred[w,2]
  dataEllipse(x, y, levels=c(0.55), add = TRUE, col = i, center.cex = 0.0001, fill =T, fill.alpha =0.1 , plot.points = FALSE)
}

###### with LDA enriched features


library(umap)
library(labdsv)
large_matrix_small = large_matrix[,sort(unique(c(uppers[[1]], uppers[[6]], lowers[[1]], lowers[[6]])))]
x <- pca(large_matrix_small,dim=10)
dimred = x$scores
par(mfrow= c(1,2), mar = c(4,4,4,3))
plot(dimred, pch = 21, bg = factor)

groups = c(list(dimred[which(factor ==levels(factor)[1]),1]), list(dimred[which(factor ==levels(factor)[2]),1]))
names(groups) = levels(factor)
pval = wilcox.test(groups[[1]], groups[[2]])$p.value
boxplot(groups, main = signif(pval, digits = 3))

pvals = NULL
for(i in c(1:length(dimred[1,]))){
  groups = c(list(dimred[which(factor ==levels(factor)[1]),i]), list(dimred[which(factor ==levels(factor)[2]),i]))
  names(groups) = levels(factor)
  pval = wilcox.test(groups[[1]], groups[[2]])$p.value
  pvals = c(pvals, pval)}

i = 5
groups = c(list(dimred[which(factor ==levels(factor)[1]),i]), list(dimred[which(factor ==levels(factor)[2]),i]))
names(groups) = levels(factor)
pval = wilcox.test(groups[[1]], groups[[2]])$p.value
boxplot(groups, main = signif(pval, digits = 3))


plot(dimred[,c(1,5)], pch = 21, bg = factor)

library(car)
for(i in c(1:length(levels(factor)))){
  w = which(factor == levels(factor)[i])
  x = dimred[w,1]
  y = dimred[w,5]
  dataEllipse(x, y, levels=c(0.45), add = TRUE, col = i, center.cex = 0.0001, fill =T, fill.alpha =0.1 , plot.points = FALSE)
}






