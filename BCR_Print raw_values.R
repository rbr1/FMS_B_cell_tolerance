# #### Run on Oxford compute set up
# # srun -p short --cpus-per-task 1 --pty bash

##### run 
module purge
module load Anaconda3/2024.02-1
conda init
source ~/.bashrc
conda --version
conda activate seurat_env
R

### run jaccard analysis
#Rscript AnalysisJaccard.R -o /well/immune-rep/shared/MISEQ/FMS/ -r FMS1 -g BCR -b /well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE/Samples_FMS1_post.txt -t 1


concat = function(v) {
	res = ""
	for (i in 1:length(v)){res = paste(res,v[i],sep="")}
	res
}
###########################################
in_dir = "/well/immune-rep/shared/MISEQ/FMS/ORIENTATED_SEQUENCES/"
out_dir = "/well/immune-rep/shared/MISEQ/FMS/ORIENTATED_SEQUENCES/COMBINED_OUTPUTS/"
batch = "FMS1"
###########################################
files = c(
  concat(c(in_dir, "Filtering_report*")), 
  concat(c(in_dir, "ISOTYPER/Cluster_per_cluster_network_parameters_SUMMARY_SUBSAMPLED*")),
  concat(c(in_dir, "ISOTYPER/Cluster_per_sequence_network_parameters*")), 
  concat(c(in_dir, "ISOTYPER/CDR3_lengths_overall*")), 
  concat(c(in_dir, "ISOTYPER/SHM_Unmutated_sequences*")), 
  concat(c(in_dir, "ISOTYPER/Isotye_normalised_overlap_frequencies_uniq*")), 
  concat(c(in_dir, "ISOTYPER/SHM_Mutation_summmary_selection*")), 
  concat(c(in_dir, "ISOTYPER/Secondary_rearrangements_*")), 
  concat(c(in_dir, "ISOTYPER/Cluster_expansion_isotype*")), 
  concat(c(in_dir, "ISOTYPER/V_gene_IGHV4_34_quantification*")), 
  concat(c(in_dir, "ISOTYPER/Secondary_rearrangements_clone_sizes*")),
  concat(c(in_dir, "ISOTYPER/Isotype_overlapping_frequencies_*")),
  concat(c(in_dir, "ISOTYPER/CDR3_charge*")),
  concat(c(in_dir, "ISOTYPER/J_gene_grouped_isotype_frequency*")),
  concat(c(in_dir, "ISOTYPER/V_gene_summary_cluster_file*")),
  concat(c(in_dir, "ISOTYPER/V_gene_per_cluster_VJ_gene_usage_by_cluster_classification*"))
)

# check files
for(i in c(1:length(files))){
  command = concat(c("ls -l ", files[i]))
  a =system(command, intern = TRUE)
  print (files[i])
  print (length(a))
}

# collapse into single files
outfiles = gsub(concat(c(in_dir, "ISOTYPER/")), concat(c(out_dir, "All_")), files)
outfiles = gsub(concat(c(in_dir, "")), concat(c(out_dir, "All_")), outfiles)
outfiles = gsub("/All_COMBINED_OUTPUTS", "", outfiles, fixed = T)
outfiles = gsub("*", ".txt", outfiles, fixed = T)
outfiles = gsub("_F.txt", ".txt", outfiles, fixed = T)

for(i in c(1:length(files))){
  command = concat(c("cat ", files[i]," > ", outfiles[i]))
  a =system(command, intern = TRUE)
  print (a)
}


#######################################################################
file = concat(c(out_dir, "All_Filtering_report.txt"))
p <- as.matrix(read.csv(file, head=TRUE, sep="\t"))
p=p[which(p[,1]!="Directory"),]
reads = as.numeric(p[,"N.BCR.filtered..post.ORF.filtering."]) 
w = which(reads>1000)
id_use = as.character(p[w,"Sample"])
reads_total = as.numeric(p[w,"N.BCR.filtered..post.ORF.filtering."])
names(reads_total) = id_use
ids_all = id_use
#############################################
Function1<-function(out_dir){
  file = concat(c(out_dir, "All_Cluster_per_cluster_network_parameters_SUMMARY_SUBSAMPLED.txt"))
  p <- as.matrix(read.csv(file, head=TRUE, sep="\t"))
  p=p[which(as.character(p[,"X.Id"]) %in% ids_all),]
  p=p[setdiff(c(1:length(p[,1])), grep("P", as.character(p[,"Isotype"]))),]
  
  id = as.character(p[,"X.Id"])
  ids = sort(unique(id))
  class = as.character(p[,"Isotype"])
  classes = sort(unique(class))
  reads = as.numeric(p[,"N.reads"])
  vgini = as.numeric(p[,"Vertex.Gini.Index"])
  N.reads = as.numeric(p[,"N.reads"])
  vertices = as.numeric(p[,"N.vertices"])
  N_clusters = as.numeric(p[,"N_clusters"])
  cgini =  as.numeric(p[,"Cluster.Gini.Index"])
  # cgini_reads =  as.numeric(p[,"cgini_clone_sampling"])
  mean_vertex_size =  as.numeric(p[,"mean_vertex_size"])
  max_clust_size =  as.numeric(p[,"max_clust_pop"])
  max_vertex_size =  as.numeric(p[,"max_vertex_pop"])
  
  # vrenyi = 1-(as.numeric(p[,"Vertex.Reyni"])/log(as.numeric(p[,"N.vertices"])))
  # crenyi =  1-(as.numeric(p[,"Cluster_Renyi"])/log(as.numeric(p[,"N_clusters"])))
  # vrenyi[which(vrenyi<0)] = -1
  # crenyi[which(crenyi<0)] = -1
  vrenyi = as.numeric(p[,"Vertex.Reyni"])
  crenyi = as.numeric(p[,"Cluster_Renyi"])
  
  class_reads = matrix(data = 0, nrow = length(ids_all),ncol = length(classes), dimnames=c(list(ids_all), list(classes)))
  v_gini = matrix(data = -1, nrow = length(ids_all),ncol = length(classes), dimnames=c(list(ids_all), list(classes)))
  c_gini = matrix(data = -1, nrow = length(ids_all),ncol = length(classes), dimnames=c(list(ids_all), list(classes)))
  c_renyi = matrix(data = -1, nrow = length(ids_all),ncol = length(classes), dimnames=c(list(ids_all), list(classes)))
  v_renyi = matrix(data = -1, nrow = length(ids_all),ncol = length(classes), dimnames=c(list(ids_all), list(classes)))
  
  m_mean_vertex_size = matrix(data = -1, nrow = length(ids_all),ncol = length(classes), dimnames=c(list(ids_all), list(classes)))
  m_max_clust_size = matrix(data = 0, nrow = length(ids_all),ncol = length(classes), dimnames=c(list(ids_all), list(classes)))
  m_max_vertex_size = matrix(data = 0, nrow = length(ids_all),ncol = length(classes), dimnames=c(list(ids_all), list(classes)))
  
  for(i in c(1:length(ids))){
    for (c in c(1:length(classes))){
      w = intersect(which(id==ids[i]), which(class==classes[c]))
      if(length(w)>=1){
        class_reads[i,c] = mean(reads[w])
        v_gini[i,c] = mean(vgini[w])
        c_gini[i,c] = mean(cgini[w])
        m_mean_vertex_size[i,c] = mean(mean_vertex_size[w])
        m_max_clust_size[i,c] = mean(max_clust_size[w])
        m_max_vertex_size[i,c] = mean(max_vertex_size[w])
        c_renyi[i,c] = mean(crenyi[w])
        v_renyi[i,c] = mean(vrenyi[w])
      }
    }
  }
  
  analysis_names = c("Vertex Gini Index","Cluster Gini Index","mean vertex size","Percentage max cluster size","Percentage max vertex size","Vertex Reyni", "Cluster Reyni")
  analysis_matrices = c(list(v_gini), list(c_gini), list(m_mean_vertex_size), list(m_max_clust_size), list(m_max_vertex_size), list(c_renyi),list(v_renyi))
  
  for(ind in c(1:length(analysis_names))){
    m = analysis_matrices[[ind]]
    cn = colnames(m)
    cn1= apply(cbind(analysis_names[ind], cn),1,paste, collapse = "..")
    colnames(m) = cn1
    analysis_matrices[[ind]]= m
  }
  analysis_matrices1 = analysis_matrices
  return(analysis_matrices1)
}
analysis_matrices1<-Function1(out_dir)

#############################################
Function2<-function(out_dir){
  file = concat(c(out_dir, "All_Cluster_per_sequence_network_parameters.txt"))
  p <- as.matrix(read.csv(file, head=TRUE, sep="\t"))
  p=p[which(as.character(p[,"X.Id"]) %in% ids_all),]
  p=p[setdiff(c(1:length(p[,1])), grep("P", as.character(p[,"Isotype"]))),]
  id = as.character(p[,"X.Id"])
  ids = sort(unique(id))
  class = as.character(p[,"Isotype"])
  classes = sort(unique(class))
  reads_per_isotype = as.numeric(p[,"N.reads"])
  unique_reads_per_isotype = as.numeric(p[,"N.vertices"])
  
  m_reads_per_isotype = matrix(data = 0, nrow = length(ids_all),ncol = length(classes), dimnames=c(list(ids_all), list(classes)))
  m_unique_reads_per_isotype = matrix(data = 0, nrow = length(ids_all),ncol = length(classes), dimnames=c(list(ids_all), list(classes)))
  
  for(i in c(1:length(ids))){
    for (c in c(1:length(classes))){
      w = intersect(which(id==ids[i]), which(class==classes[c]))
      if(length(w)>=1){
        m_reads_per_isotype[i,c] = mean(reads_per_isotype[w])
        m_unique_reads_per_isotype[i,c] = mean(unique_reads_per_isotype[w])
      }
    }
  }
  
  c1 = c("Class_switched","IGHD,IGHM_mutated","IGHD,IGHM_unmutated")
  c2 = c( "IGHA1","IGHA2","IGHD","IGHE","IGHG1","IGHG2","IGHG3","IGHG4","IGHM"  )
  c2 = c2[which(c2 %in% classes)]
  c3 = c( "IGHA1","IGHA2","IGHE","IGHG1","IGHG2","IGHG3","IGHG4" )
  
  m_reads_per_isotype_group = m_reads_per_isotype[,c1]
  m_unique_reads_per_isotype_group = m_unique_reads_per_isotype[,c1]
  m_reads_per_isotype_single = m_reads_per_isotype[,c2]
  m_unique_reads_per_isotype_single = m_unique_reads_per_isotype[,c2]
  m_unique_reads_per_isotype_switched = m_unique_reads_per_isotype[,c3]
  
  for(i in c(1:length(ids))){
    m_reads_per_isotype_group[i,] = m_reads_per_isotype_group[i,]*100/sum(m_reads_per_isotype_group[i,])
    m_unique_reads_per_isotype_group[i,] = m_unique_reads_per_isotype_group[i,]*100/sum(m_unique_reads_per_isotype_group[i,])
    m_reads_per_isotype_single[i,] = m_reads_per_isotype_single[i,]*100/sum(m_reads_per_isotype_single[i,])
    m_unique_reads_per_isotype_single[i,] = m_unique_reads_per_isotype_single[i,]*100/sum(m_unique_reads_per_isotype_single[i,])
    m_unique_reads_per_isotype_switched[i,] = m_unique_reads_per_isotype_switched[i,]*100/sum(m_unique_reads_per_isotype_switched[i,])
  }
  
  analysis_names = c( "Percentage unique_BCRs_per_isotype", "Percentage unique_BCRs_per_isotype_group", "Percentage unique_BCRs_in_switched")
  analysis_matrices = c(list(m_unique_reads_per_isotype_group), list(m_unique_reads_per_isotype_single), list(m_unique_reads_per_isotype_switched))
  
  for(ind in c(1:length(analysis_names))){
    m = analysis_matrices[[ind]]
    cn = colnames(m)
    cn1= apply(cbind(analysis_names[ind], cn),1,paste, collapse = "..")
    colnames(m) = cn1
    analysis_matrices[[ind]]= m
  }
  analysis_matrices2 = analysis_matrices
  return(analysis_matrices2)
}
analysis_matrices2<-Function2(out_dir)

#############################################
Function3<-function(out_dir){
  file = concat(c(out_dir, "All_CDR3_lengths_overall.txt"))
  p <- as.matrix(read.csv(file, head=TRUE, sep="\t"))
  p=p[which(as.character(p[,"X.sample"]) %in% ids_all),]
  p=p[setdiff(c(1:length(p[,1])), grep("P", as.character(p[,"isotype"]))),]
  p=p[which(as.numeric(p[,"Number_of_BCRs"])>10),]
  id = as.character(p[,"X.sample"])
  ids = sort(unique(id))
  class = as.character(p[,"isotype"])
  chains = sort(unique(class))
  value = as.numeric(p[,"mean_CDR3_length"])
  reads = as.numeric(p[,"Number_of_BCRs"])
  
  values = matrix(data = -1, nrow = length(ids_all),ncol = length(chains), dimnames=c(list(ids_all), list(chains)))
  print_info = list()
  for(i in c(1:length(id))){
    values[id[i], class[i]] = value[i]
  }
  
  analysis_matrices = list(values)
  analysis_names = c("Mean CDR3 lengths")
  for(ind in c(1:length(analysis_names))){
    m = analysis_matrices[[ind]]
    cn = colnames(m)
    cn1= apply(cbind(analysis_names[ind], cn),1,paste, collapse = "..")
    colnames(m) = cn1
    analysis_matrices[[ind]]= m
  }
  analysis_matrices3 = analysis_matrices
  return(analysis_matrices3)
}
analysis_matrices3<-Function3(out_dir)

#############################################
Function4<-function(out_dir){
  file = concat(c(out_dir, "All_SHM_Unmutated_sequences.txt"))
  p <- as.matrix(read.csv(file, head=TRUE, sep="\t"))
  p=p[which(as.character(p[,"X.Sample"]) %in% ids_all),]
  p=p[setdiff(c(1:length(p[,1])), grep("P", as.character(p[,"isotype"]))),]
  p=p[which(is.na(as.numeric(p[,"mean.mutations"]))==F),]
  p=p[which(as.numeric(p[,"number.of.unique.sequences"])>10),]
  id = as.character(p[,"X.Sample"])
  ids = sort(unique(id))
  class = as.character(p[,"isotype"])
  chains = sort(unique(class))
  value = as.numeric(p[,"mean.mutations"])
  reads = as.numeric(p[,"number.of.unique.sequences"])
  perc_unmutated = as.numeric(p[,"perc_unumtated"])
  
  values = matrix(data = -1, nrow = length(ids_all),ncol = length(chains), dimnames=c(list(ids_all), list(chains)))
  m_perc_unmutated = matrix(data = -1, nrow = length(ids_all),ncol = length(chains), dimnames=c(list(ids_all), list(chains)))
  print_info = list()
  for(i in c(1:length(id))){
    values[id[i], class[i]] = value[i]
    m_perc_unmutated[id[i], class[i]] = perc_unmutated[i]
  }
  
  analysis_matrices = c(list(values),list(m_perc_unmutated))
  analysis_names = c("Mean SHM per BCR","Percentage unmutated")
  for(ind in c(1:length(analysis_names))){
    m = analysis_matrices[[ind]]
    cn = colnames(m)
    cn1= apply(cbind(analysis_names[ind], cn),1,paste, collapse = "..")
    colnames(m) = cn1
    analysis_matrices[[ind]]= m
  }
  analysis_matrices4 = analysis_matrices
  return(analysis_matrices4)
}
analysis_matrices4<-Function4(out_dir)

#############################################
Function5<-function(out_dir){
  file = concat(c(out_dir, "All_Isotye_normalised_overlap_frequencies_uniq.txt"))
  p <- as.matrix(read.csv(file, head=TRUE, sep="\t"))
  p=p[which(as.character(p[,"X.sample"]) %in% ids_all),]
  # p=p[setdiff(c(1:length(p[,1])), grep("IGHG4",as.character(p[,"iso1"]))),]
  # p=p[setdiff(c(1:length(p[,1])), grep("IGHG4",as.character(p[,"iso2"]))),]
  id = as.character(p[,"X.sample"])
  ids = sort(unique(id))
  class1 = as.character(p[,"iso1"])
  class2 = as.character(p[,"iso2"])
  classes = sort(unique(class1))
  classes = classes[grep("P", classes, invert = T)]
  classes = classes[grep("IGH", classes, invert = F)]
  w = intersect(which(class1 %in% classes), which(class2 %in% classes))
  p=p[w,]
  id = as.character(p[,"X.sample"])
  ids = sort(unique(id))
  class1 = as.character(p[,"iso1"])
  class2 = as.character(p[,"iso2"])
  classes = sort(unique(class1))
  classes = classes[grep("P", classes, invert = T)]
  classes = classes[grep("IGH", classes, invert = F)]
  mean_overlap_proportion1 = as.numeric(p[,"mean_overlap"])
  class12 = apply(cbind(class1, class2),1,paste,collapse = "-")
  class12s = sort(unique(class12))
  
  overlap = matrix(data = -1, nrow = length(ids_all),ncol = length(class12s), dimnames=c(list(ids_all), list(class12s)))
  for(i in c(1:length(ids))){
    w = which(id==ids[i])
    if(length(w)>0){
      overlap[ids[i],]=0
      overlap[ids[i], class12[w]]= mean_overlap_proportion1[w]
    }
  }
  
  analysis_names = "Relative class switching absolute"
  analysis_matrices = list(overlap)
  for(ind in c(1:length(analysis_names))){
    m = analysis_matrices[[ind]]
    cn = colnames(m)
    cn1= apply(cbind(analysis_names[ind], cn),1,paste, collapse = "..")
    colnames(m) = cn1
    analysis_matrices[[ind]]= m
  }
  analysis_matrices5 = analysis_matrices
  return(analysis_matrices5)
}
analysis_matrices5<-Function5(out_dir)

#############################################
Function6<-function(out_dir){
  file = concat(c(out_dir, "All_SHM_Mutation_summmary_selection.txt"))
  p <- as.matrix(read.csv(file, head=TRUE, sep="\t"))
  p=p[which(as.character(p[,"X.sample"]) %in% ids_all),]
  p=p[setdiff(c(1:length(p[,1])), grep("P", as.character(p[,"chain"]))),]
  p=p[which(as.numeric(p[,"total_BCRs_counted"])>10),]
  id = as.character(p[,"X.sample"])
  ids = sort(unique(id))
  class = as.character(p[,"chain"])
  type = as.character(p[,"count_type"])
  types = sort(unique(type))
  chains = sort(unique(class))
  value = as.numeric(p[,"mean.value"])
  
  value_list = list()
  for(t in c(1:length(types))){
    values = matrix(data = -1, nrow = length(ids_all),ncol = length(chains), dimnames=c(list(ids_all), list(chains)))
    w = which(type== types[t])
    for(i in c(1:length(w))){
      values[id[w[i]], class[w[i]]] = value[w[i]]}
    value_list = c(value_list, list(values))
  }
  a = value_list [[which(types=="mean CDR_mm per BCR")]]+ value_list [[which(types=="mean FWR_mm per BCR")]]
  b = value_list [[which(types=="mean nonsilent per BCR" )]]+ value_list [[which(types=="mean silent per BCR")]]
  
  types = c(types, "Mean mutations per BCR")
  value_list = c(value_list, list(a))
  
  w = which(types %in% c("mean CDR_FWR_ratio", "Mean mutations per BCR" ,"FWR3_mm", "mean CDR_mm per BCR","mean FWR_mm per BCR"  ))
  analysis_matrices = value_list[w]
  analysis_names = types[w]
  for(ind in c(1:length(analysis_matrices))){
    m = analysis_matrices[[ind]]
    cn = colnames(m)
    cn1= apply(cbind(analysis_names[ind], cn),1,paste, collapse = "..")
    colnames(m) = cn1
    analysis_matrices[[ind]]= m
  }
  analysis_matrices6 = analysis_matrices
  return(analysis_matrices6)
}
analysis_matrices6<-Function6(out_dir)

#############################################
Function7<-function(out_dir){
  file = concat(c(out_dir, "All_Secondary_rearrangements.txt"))
  p <- as.matrix(read.csv(file, head=TRUE, sep="\t"))
  p=p[which(as.character(p[,"X.sample"]) %in% ids_all),]
  p=p[grep("IGH",as.character(p[,"chain"])),]
  w = setdiff(p[,"chain"], p[,"chain"][grep("IGHV",as.character(p[,"chain"]))])
  p=p[which(as.character(p[,"chain"]) %in% w),]
  w = setdiff(p[,"chain"], p[,"chain"][grep("mut",as.character(p[,"chain"]))])
  p=p[which(as.character(p[,"chain"]) %in% w),]
  
  id = as.character(p[,"X.sample"])
  ids = sort(unique(id))
  class = as.character(p[,"chain"])
  id_replacement_freq =as.numeric(p[,"percentage"])
  total_isotype = as.numeric(p[,"total"])
  w1 = which(total_isotype>50)
  
  classes = sort(unique(class))
  types = c(list(id_replacement_freq))
  
  analysis_name = c("V gene replacement frequency")
  means = matrix(data = -1, nrow = length(ids_all),ncol = length(classes), dimnames=c(list(ids_all), list(classes)))
  for(i in c(1:length(ids))){
    w = intersect(which(id==ids[i]), w1)
    means[ids[i], class[w]] = id_replacement_freq[w]
  }
  
  analysis_matrices = list(means)
  analysis_names = c("V gene replacement frequency")
  for(ind in c(1:length(analysis_matrices))){
    m = analysis_matrices[[ind]]
    cn = colnames(m)
    cn1= apply(cbind(analysis_names[ind], cn),1,paste, collapse = "..")
    colnames(m) = cn1
    analysis_matrices[[ind]]= m
  }
  analysis_matrices7 = analysis_matrices
  return(analysis_matrices7)
}
analysis_matrices7<-Function7(out_dir)

#############################################
Function8<-function(out_dir){
  file = concat(c(out_dir, "All_Cluster_expansion_isotype.txt"))
  p <- as.matrix(read.csv(file, head=TRUE, sep="\t"))
  p=p[which(as.character(p[,"X.sample"]) %in% ids_all),]
  p=p[which(as.numeric(p[,"d5"]) !=-1),]
  id = as.character(p[,"X.sample"])
  ids = sort(unique(id))
  isotype = as.character(p[,"isotype"])
  d5 = as.numeric(p[,"d5"])
  d10 = as.numeric(p[,"d10"])
  d50 = as.numeric(p[,"d50"])
  classes = sort(unique(isotype))
  
  m_d5 = matrix(data = -1, nrow = length(ids_all),ncol = length(classes), dimnames=c(list(ids_all), list(classes)))
  m_d10 = matrix(data = -1, nrow = length(ids_all),ncol = length(classes), dimnames=c(list(ids_all), list(classes)))
  m_d50 = matrix(data = -1, nrow = length(ids_all),ncol = length(classes), dimnames=c(list(ids_all), list(classes)))
  
  for(i in c(1:length(ids))){
    w = which(id==ids[i])
    m_d5[ids[i],isotype[w]]= d5[w]
    m_d10[ids[i],isotype[w]]= d10[w]
    m_d50[ids[i],isotype[w]]= d50[w]
  }
  analysis_matrices = c(list(m_d5),list(m_d10),list(m_d50))
  analysis_names = c("D5","D10","D50")
  for(ind in c(1:length(analysis_matrices))){
    m = analysis_matrices[[ind]]
    cn = colnames(m)
    cn1= apply(cbind(analysis_names[ind], cn),1,paste, collapse = "..")
    colnames(m) = cn1
    analysis_matrices[[ind]]= m
  }
  
  analysis_matrices8 = analysis_matrices
  return(analysis_matrices8)
}
analysis_matrices8<-Function8(out_dir)

#############################################
Function9<-function(out_dir){
  file = concat(c(out_dir, "All_V_gene_IGHV4_34_quantification.txt"))
  p <- as.matrix(read.csv(file, head=TRUE, sep="\t"))
  p=p[which(as.character(p[,"X.sample"]) %in% ids_all),]
  p=p[which(as.numeric(p[,"total_all"]) >10),]
  id = as.character(p[,"X.sample"])
  ids = sort(unique(id))
  isotype = as.character(p[,"isotype"])
  V4_34_AVY_total_unmut = as.numeric(p[,"V4_34_AVY_total_unmut"])*100/as.numeric(p[,"total_all"])
  V4_34_NHS_total_unmut = as.numeric(p[,"V4_34_NHS_total_unmut"])*100/as.numeric(p[,"total_all"])
  V4_34_AVY_NHS_total_unmut = as.numeric(p[,"V4_34_AVY_NHS_total_unmut"])*100/as.numeric(p[,"total_all"])
  classes = sort(unique(isotype))
  V4_34_AVY_NHS_total_unmut_count = as.numeric(p[,"V4_34_AVY_NHS_total_unmut"])
  all_counts = as.numeric(p[,"total_all"])
  
  m_V4_34_AVY_total_unmut = matrix(data = -1, nrow = length(ids_all),ncol = length(classes), dimnames=c(list(ids_all), list(classes)))
  m_V4_34_NHS_total_unmut = matrix(data = -1, nrow = length(ids_all),ncol = length(classes), dimnames=c(list(ids_all), list(classes)))
  m_V4_34_AVY_NHS_total_unmut = matrix(data = -1, nrow = length(ids_all),ncol = length(classes), dimnames=c(list(ids_all), list(classes)))
  m_V4_34_AVY_NHS_total_unmut_count = matrix(data = -1, nrow = length(ids_all),ncol = length(classes), dimnames=c(list(ids_all), list(classes)))
  m_all_count = matrix(data = -1, nrow = length(ids_all),ncol = length(classes), dimnames=c(list(ids_all), list(classes)))
  
  for(i in c(1:length(ids))){
    w = which(id==ids[i])
    m_V4_34_AVY_total_unmut[ids[i],isotype[w]]= V4_34_AVY_total_unmut[w]
    m_V4_34_NHS_total_unmut[ids[i],isotype[w]]= V4_34_NHS_total_unmut[w]
    m_V4_34_AVY_NHS_total_unmut[ids[i],isotype[w]]= V4_34_AVY_NHS_total_unmut[w]
    m_V4_34_AVY_NHS_total_unmut_count[ids[i],isotype[w]]= V4_34_AVY_NHS_total_unmut_count[w]
    m_all_count[ids[i],isotype[w]]= all_counts[w]
  }
  IGHDM = rep(-1, length(ids_all))
  names(IGHDM)= ids_all
  class_switched = IGHDM
  for(i in c(1:length(ids))){
    w1 = which(classes %in% c("IGHA1", "IGHA2" , "IGHE",  "IGHG1" ,"IGHG2" ,"IGHG3" ,"IGHG4"))
    w2 = which(classes %in% c("IGHD", "IGHM"))
    w = which(m_V4_34_AVY_NHS_total_unmut[ids[i],]!=-1)
    IGHDM[ids[i]] = sum(m_V4_34_AVY_NHS_total_unmut_count[ids[i],intersect(w,w2)])*100/sum(m_all_count[ids[i],intersect(w,w2)])
    class_switched[ids[i]] = sum(m_V4_34_AVY_NHS_total_unmut_count[ids[i],intersect(w,w1)])*100/sum(m_all_count[ids[i],intersect(w,w1)])
    
  }
  m_V4_34_AVY_NHS_total_unmut = cbind(m_V4_34_AVY_NHS_total_unmut, IGHDM, class_switched)
  
  analysis_matrices = c(list(m_V4_34_AVY_total_unmut),list(m_V4_34_NHS_total_unmut),list(m_V4_34_AVY_NHS_total_unmut))
  analysis_names = c("V4_34_AVY_unmut","V4_34_NHS_unmut","V4_34_AVY_NHS_unmut")
  for(ind in c(1:length(analysis_matrices))){
    m = analysis_matrices[[ind]]
    cn = colnames(m)
    cn1= apply(cbind(analysis_names[ind], cn),1,paste, collapse = "..")
    colnames(m) = cn1
    analysis_matrices[[ind]]= m
  }
  
  analysis_matrices9 = analysis_matrices
  return(analysis_matrices9)
}
analysis_matrices9<-Function9(out_dir)

#############################################
Function10<-function(out_dir){
  file = concat(c(out_dir, "All_Secondary_rearrangements_clone_sizes.txt"))
  p <- as.matrix(read.csv(file, head=TRUE, sep="\t"))
  p=p[which(as.character(p[,"X.sample"]) %in% ids_all),]
  
  id = as.character(p[,"X.sample"])
  ids = sort(unique(id))
  d5_norm =as.numeric(p[,"d5_norm"])
  d5_secondary =as.numeric(p[,"d5_secondary"])
  mean_clone_size_norm =as.numeric(p[,"mean_clone_size_norm"])
  mean_clone_size_secondary =as.numeric(p[,"X_mean_clone_size_secondary"])
  w1 = which(as.numeric(p[,"n_secondary"])>5)
  groups = c(list(mean_clone_size_norm[w1]), list(mean_clone_size_secondary[w1] ))
  boxplot(groups)
  
  groups = c(list(d5_norm[w1]), list(d5_secondary[w1] ))
  boxplot(groups)
  
  headers = c("d5_norm","d5_secondary","mean_clone_size_norm", "mean_clone_size_secondary")
  m_all_count = matrix(data = -1, nrow = length(ids_all),ncol = length(headers), dimnames=c(list(ids_all), list(headers)))
  x=  cbind(d5_norm,d5_secondary,mean_clone_size_norm,mean_clone_size_secondary)
  m_all_count[id,] = x
  
  
  analysis_matrices = list(m_all_count)
  analysis_names = c("V gene replacement clonal expansion")
  for(ind in c(1:length(analysis_matrices))){
    m = analysis_matrices[[ind]]
    cn = colnames(m)
    cn1= apply(cbind(analysis_names[ind], cn),1,paste, collapse = "..")
    colnames(m) = cn1
    analysis_matrices[[ind]]= m
  }
  analysis_matrices10 = analysis_matrices
  return(analysis_matrices10)
}
analysis_matrices10<-Function10(out_dir)

#############################################
Function11<-function(out_dir){
  file = concat(c(out_dir, "All_Isotype_overlapping_frequencies_.txt"))
  p <- as.matrix(read.csv(file, head=TRUE, sep="\t"))
  p=p[which(as.character(p[,"X.ID"]) %in% ids_all),]
  table(p[,"sample_depth"])
  p=p[which(as.numeric(p[,"sample_depth"])==50),] ## use the 50 sampling depth
  id = as.character(p[,"X.ID"])
  ids = sort(unique(id))
  class1 = as.character(p[,"isotype1"])
  class2 = as.character(p[,"isotype2"])
  classes = sort(unique(class1))
  classes = classes[grep("P", classes, invert = T)]
  classes = classes[grep("IGH", classes, invert = F)]
  w = intersect(which(class1 %in% classes), which(class2 %in% classes))
  p=p[w,]
  id = as.character(p[,"X.ID"])
  ids = sort(unique(id))
  class1 = as.character(p[,"isotype1"])
  class2 = as.character(p[,"isotype2"])
  classes = sort(unique(class1))
  classes = classes[grep("P", classes, invert = T)]
  classes = classes[grep("IGH", classes, invert = F)]
  mean_overlap_proportion1 = as.numeric(p[,"mean_overlap"])
  class12 = apply(cbind(class1, class2),1,paste,collapse = "-")
  class12s = sort(unique(class12))
  
  overlap = matrix(data = -1, nrow = length(ids_all),ncol = length(class12s), dimnames=c(list(ids_all), list(class12s)))
  for(i in c(1:length(ids))){
    w = which(id==ids[i])
    if(length(w)>0){
      overlap[ids[i],]=0
      overlap[ids[i], class12[w]]= mean_overlap_proportion1[w]
    }
  }
  
  analysis_names = "Relative class switching normalised"
  analysis_matrices = list(overlap)
  for(ind in c(1:length(analysis_names))){
    m = analysis_matrices[[ind]]
    cn = colnames(m)
    cn1= apply(cbind(analysis_names[ind], cn),1,paste, collapse = "..")
    colnames(m) = cn1
    analysis_matrices[[ind]]= m
  }
  analysis_matrices11 = analysis_matrices
  return(analysis_matrices11)
}
analysis_matrices11<-Function11(out_dir)

#############################################
Function12<-function(out_dir){
  file = concat(c(out_dir, "All_CDR3_charge.txt"))
  p <- as.matrix(read.csv(file, head=TRUE, sep="\t"))
  p=p[which(as.character(p[,"X.sample"]) %in% ids_all),]
  p=p[which(as.numeric(p[,"Number_of_BCRs"]) >=5),]
  id = as.character(p[,"X.sample"])
  ids = sort(unique(id))
  isotype = as.character(p[,"isotype"])
  mean_CDR3_R_K_residues = as.numeric(p[,"mean_CDR3_R_K_residues"])
  classes = sort(unique(isotype))
  
  m_mean_CDR3_R_K_residues = matrix(data = -1, nrow = length(ids_all),ncol = length(classes), dimnames=c(list(ids_all), list(classes)))
  
  for(i in c(1:length(ids))){
    w = which(id==ids[i])
    m_mean_CDR3_R_K_residues[ids[i],isotype[w]]= mean_CDR3_R_K_residues[w]
    }
  analysis_matrices = c(list(m_mean_CDR3_R_K_residues))
  analysis_names = c("mean_CDR3_charge")
  for(ind in c(1:length(analysis_matrices))){
    m = analysis_matrices[[ind]]
    cn = colnames(m)
    cn1= apply(cbind(analysis_names[ind], cn),1,paste, collapse = "..")
    colnames(m) = cn1
    analysis_matrices[[ind]]= m
  }
  
  analysis_matrices12 = analysis_matrices
  return(analysis_matrices12)
}
analysis_matrices12<-Function12(out_dir)

#############################################
Function13<-function(out_dir){
  file = concat(c(out_dir, "All_J_gene_grouped_isotype_frequency.txt"))
  p <- as.matrix(read.csv(file, head=TRUE, sep="\t"))
  p=p[which(as.character(p[,"X.sample"]) %in% ids_all),]
  id = as.character(p[,"X.sample"])
  ids = sort(unique(id))
  isotype = as.character(p[,"class"])
  J.gene = as.character(p[,"J.gene"])
  count = as.numeric(p[,"uniq_read_freq"])
  classes = sort(unique(isotype))
  J.genes = sort(unique(J.gene))
  
  all_classes = NULL
  for(c in c(1:length(classes))){
    m_J_gene = matrix(data = 0, nrow = length(ids_all),ncol = length(J.genes), dimnames=c(list(ids_all), list(J.genes)))
    for(j in c(1:length(J.genes))){
      w = intersect(which(J.gene==J.genes[j]), which(isotype==classes[c]))
      for(i in c(1:length(ids))){
        w1 = intersect(w, which(id==ids[i]))
        m_J_gene[id[w1], J.gene[w1]]= count[w1]
      }
    }
    for(i in c(1:length(ids))){
      if(sum(m_J_gene[ids[i],])>20){
        m_J_gene[ids[i],] = m_J_gene[ids[i],]*100/sum(m_J_gene[ids[i],])
      }}
    all_classes = c(all_classes, list(m_J_gene))
  }
  names(all_classes) = classes
  analysis_matrices = all_classes
  analysis_names = paste("J_gene_freq_by_uniq_VDJ", names(all_classes))
  for(ind in c(1:length(analysis_matrices))){
    m = analysis_matrices[[ind]]
    cn = colnames(m)
    cn1= apply(cbind(analysis_names[ind], cn),1,paste, collapse = "..")
    colnames(m) = cn1
    analysis_matrices[[ind]]= m
  }
  
  analysis_matrices13 = analysis_matrices
  return(analysis_matrices13)
}
analysis_matrices13<-Function13(out_dir)

#############################################
Function14<-function(out_dir){
  file = concat(c(out_dir, "All_V_gene_summary_cluster_file.txt"))
  p <- as.matrix(read.csv(file, head=TRUE, sep="\t"))
  p=p[which(as.character(p[,"X.sample"]) %in% ids_all),]
  id = as.character(p[,"X.sample"])
  ids = sort(unique(id))
  isotype = as.character(p[,"class"])
  gene = as.character(p[,"V.gene"])
  count = as.numeric(p[,"n_seq"])
  count_per_cluster = as.numeric(p[,"n_clusters"])
  mean_cluster_size =count/count_per_cluster
  classes = sort(unique(isotype))
  genes = sort(unique(gene))
  
  all_classes_per_VDJ = NULL
  for(c in c(1:length(classes))){
    m_J_gene = matrix(data = 0, nrow = length(ids_all),ncol = length(genes), dimnames=c(list(ids_all), list(genes)))
    for(j in c(1:length(genes))){
      w = intersect(which(gene==genes[j]), which(isotype==classes[c]))
      for(i in c(1:length(ids))){
        w1 = intersect(w, which(id==ids[i]))
        m_J_gene[id[w1], gene[w1]]= count[w1]
      }
    }
    for(i in c(1:length(ids))){
      if(sum(m_J_gene[ids[i],])>20){
        m_J_gene[ids[i],] = m_J_gene[ids[i],]*100/sum(m_J_gene[ids[i],])
      }else{m_J_gene[ids[i],] = -1}
      }
    all_classes_per_VDJ = c(all_classes_per_VDJ, list(m_J_gene))
  }
  names(all_classes_per_VDJ) = classes
  
  all_classes_per_cluster = NULL
  for(c in c(1:length(classes))){
    m_J_gene = matrix(data = 0, nrow = length(ids_all),ncol = length(genes), dimnames=c(list(ids_all), list(genes)))
    for(j in c(1:length(genes))){
      w = intersect(which(gene==genes[j]), which(isotype==classes[c]))
      for(i in c(1:length(ids))){
        w1 = intersect(w, which(id==ids[i]))
        m_J_gene[id[w1], gene[w1]]= count_per_cluster[w1]
      }
    }
    for(i in c(1:length(ids))){
      if(sum(m_J_gene[ids[i],])>20){
        m_J_gene[ids[i],] = m_J_gene[ids[i],]*100/sum(m_J_gene[ids[i],])
      }else{m_J_gene[ids[i],] = -1}}
    all_classes_per_cluster = c(all_classes_per_cluster, list(m_J_gene))
  }
  names(all_classes_per_cluster) = classes
  
  
  all_cluster_size= NULL
  for(c in c(1:length(classes))){
    m_J_gene = matrix(data = -1, nrow = length(ids_all),ncol = length(genes), dimnames=c(list(ids_all), list(genes)))
    for(j in c(1:length(genes))){
      w = intersect(which(gene==genes[j]), which(isotype==classes[c]))
      for(i in c(1:length(ids))){
        w1 = intersect(w, which(id==ids[i]))
        w1 = intersect(which(count>=10), w1)
        m_J_gene[id[w1], gene[w1]]= mean_cluster_size[w1]
      }
    }
    all_cluster_size = c(all_cluster_size, list(m_J_gene))
  }
  names(all_cluster_size) = classes
  
  #names(all_classes_per_VDJ) = paste("V_gene_freq_by_uniq_VDJ",  names(all_classes_per_VDJ))
  #names(all_classes_per_cluster) = paste("V_gene_freq_by_cluster",  names(all_classes_per_cluster))
  #names(all_cluster_size) = paste("Mean_cluster_size",  names(all_cluster_size))
  
  analysis_names = c(paste("V_gene_freq_by_uniq_VDJ",  names(all_classes_per_VDJ)), 
                     paste("V_gene_freq_by_cluster",  names(all_classes_per_cluster)), 
                     paste("Mean_cluster_size",  names(all_cluster_size)))
  analysis_matrices = c(all_classes_per_VDJ, all_classes_per_cluster, all_cluster_size)
  for(ind in c(1:length(analysis_matrices))){
    m = analysis_matrices[[ind]]
    cn = colnames(m)
    cn1= apply(cbind(analysis_names[ind], cn),1,paste, collapse = "..")
    colnames(m) = cn1
    analysis_matrices[[ind]]= m
  }
  
  analysis_matrices14 = analysis_matrices
  return(analysis_matrices14)
}
analysis_matrices14<-Function14(out_dir)

#############################################
Function15<-function(out_dir){
  file = concat(c(out_dir, "All_V_gene_per_cluster_VJ_gene_usage_by_cluster_classification.txt"))
  p <- as.matrix(read.csv(file, head=TRUE, sep="\t"))
  p=p[which(as.character(p[,"X.ID"]) %in% ids_all),]
  id = as.character(p[,"X.ID"])
  ids = sort(unique(id))
  isotype = as.character(p[,"classification"])
  gene = as.character(p[,"VJ"])
  count = as.numeric(p[,"number.of.sequences"])
  t = table(isotype)
  classes = names(t)[which(t>=100)]
  genes = sort(unique(gene))
  
  all_classes_per_VDJ = NULL
  for(c in c(1:length(classes))){
    m_J_gene = matrix(data = 0, nrow = length(ids_all),ncol = length(genes), dimnames=c(list(ids_all), list(genes)))
    for(j in c(1:length(genes))){
      w = intersect(which(gene==genes[j]), which(isotype==classes[c]))
      for(i in c(1:length(ids))){
        w1 = intersect(w, which(id==ids[i]))
        m_J_gene[id[w1], gene[w1]]= count[w1]
      }
    }
    for(i in c(1:length(ids))){
      if(sum(m_J_gene[ids[i],])>20){
        m_J_gene[ids[i],] = m_J_gene[ids[i],]*100/sum(m_J_gene[ids[i],])
      }}
    all_classes_per_VDJ = c(all_classes_per_VDJ, list(m_J_gene))
    print(c)
  }
  names(all_classes_per_VDJ) = classes
  
  analysis_names = c(paste("VJ_gene_freq_by_uniq_VDJ",  names(all_classes_per_VDJ)))
  analysis_matrices = c(all_classes_per_VDJ)
  for(ind in c(1:length(analysis_matrices))){
    m = analysis_matrices[[ind]]
    cn = colnames(m)
    cn1= apply(cbind(analysis_names[ind], cn),1,paste, collapse = "..")
    colnames(m) = cn1
    analysis_matrices[[ind]]= m
  }
  
  analysis_matrices15 = analysis_matrices
  return(analysis_matrices15)
}
analysis_matrices15<-Function15(out_dir)

###################################################################################################
print_info = c(analysis_matrices1, analysis_matrices2,  analysis_matrices3, analysis_matrices4, analysis_matrices5, 
               analysis_matrices6, analysis_matrices7, analysis_matrices8, analysis_matrices9, analysis_matrices10,
               analysis_matrices11, analysis_matrices12, analysis_matrices13, analysis_matrices14, analysis_matrices15)

overall_matrix = NULL
for(i in c(1:length(print_info))){
	if(length(overall_matrix)==0){
		overall_matrix = print_info[[i]]
	}else{overall_matrix = cbind(overall_matrix ,print_info[[i]][ids_all,])}
}
colnames(overall_matrix) = gsub(" ","_", colnames(overall_matrix), fixed= T)

out_file_table = concat(c(out_dir, "All_raw_values_",batch,".txt"))
write.table(overall_matrix, file = out_file_table, append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = T,col.names = TRUE, qmethod = c("escape", "double"),fileEncoding = "")


