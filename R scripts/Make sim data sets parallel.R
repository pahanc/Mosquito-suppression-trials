
#Define a function to make a simulated trial data set given the number of clusters per arm,
#suppression effect (G) and simulation number (runno)
make_sim_counts<-function(cluster_per_group, G, runno){
  
  num_months=12 #Number of months sampled per year
  num_group=2 #Number of arms
  num_period=3 #Number of trial periods, including a baseline year
  
  #Period
  #Make a vector of the period for each cluster of length cluster_per_group-by-num_period
  period<-rep(1,cluster_per_group)
  for (i in 2:num_period){
    period<-c(period,rep(i, cluster_per_group))
  }
  #Replicate period vector for trial arm
  period<-rep(period,num_group)

  #Cluster
  #Make a vector of cluster indexes of length cluster_per_group-by-num_period*num_group
  cluster<-1:cluster_per_group
  cluster2<-cluster
  for (i in 1:(num_period*num_group-1)){
    cluster2<-c(cluster2, cluster)
  }
  #Add an offset to place clusters into groups (trial arms)
  offset<-matrix(rep(0:(num_group-1), num_period*cluster_per_group)*cluster_per_group,
               byrow=T, nc=num_group)
  offset<-as.vector(offset)
  cluster2<-cluster2+offset

  #Intervention effects
  #Make a vector indicating whether each cluster receives the intervention
  #according to a parallel design with a baseline year
  treatment<-c(rep('control',cluster_per_group), rep('intervention', (num_period-1)*cluster_per_group), 
             rep('control', num_period*cluster_per_group))

  specs<-data.frame(period,cluster2,treatment)
  #order the data by cluster index (then by period)
  specs<-specs[order(specs$cluster2),]


  # Make the simulated data sets
  cluster_names<-c("Bana village","Bana market", "Pala", "Souroukoudingan")
  year_names<-c("2012","2013", "2014", "2017", "2018", "2019")
  for (rep in runno){
    counts_df<-data.frame(grp=NA, cluster=NA, period=NA, mth=NA,treatment=NA,count=NA,sel_cluster=NA,sel_year=NA,cluster_rf=NA,samp_no=NA, samp_index=NA, count_orig=NA )[numeric(0), ]
    index=1
    index2=1
    for (grp in 1:num_group){
      for (cluster in 1:cluster_per_group){
        cluster_rf<-rnorm(n=1,mean=0,sd=1/sqrt(obs.noise.village[1]))#Draw the cluster-level random effect
        year_ind_vec<-NULL
        for (period1 in 1:num_period){
          if (period1 ==1){
            sel_cluster_ind<-sample(1:4,1)#select a cluster at random to represent each cluster
            sel_cluster<-cluster_names[sel_cluster_ind]
            year_ind<-sample(x=seq(1,6),size=1)#select a year at random
            year_ind_vec<-c(year_ind_vec,year_ind)#store the selected years used to make each 3-year sequence
            sel_year<-year_names[year_ind]
          }
        if (period1>1){#Use a different year to generate cluster data for each period
          year_ind<-sample(x=seq(1,6)[-year_ind_vec],size=1)
          sel_year<-year_names[year_ind]
        }
        for (mth in 1:12){
          #Get the sampling locations for the selected cluster, year and month
          samp_inds<-which(metadata$samp_village==sel_cluster & metadata$samp_years==sel_year & metadata$samp_months==mth)
          #Get the simulated mosquito counts for the selected indices and
          #apply the suppression effect to the treatment clusters
          counts<-lp2_i[samp_inds,(grp-1)*cluster_per_group+cluster] + cluster_rf +log(G)*(specs$treatment[index]=='intervention')
          #Store the mosquito counts without applying the suppression effect
          counts_orig<-lp2_i[samp_inds,(grp-1)*cluster_per_group+cluster] + cluster_rf
          for (samp in 1:length(samp_inds)){
            if (length(counts)>0){
              counts_df[index2,"grp"]=grp
              counts_df[index2,"cluster"]=specs$cluster2[index]
              counts_df[index2,"period"]=specs$period[index]
              counts_df[index2,"mth"]=mth
              counts_df[index2,"treatment"]=specs$treatment[index]
              counts_df[index2,"count"]<-rpois(1,exp(counts[samp]))
              counts_df[index2,"sel_cluster"]<-sel_cluster
              counts_df[index2,"sel_year"]<-sel_year
              counts_df[index2,"cluster_rf"]<-cluster_rf
              counts_df[index2, "samp_no"]<-(grp-1)*cluster_per_group+cluster
              counts_df[index2, "samp_index"]<-samp_inds[samp]
              counts_df[index2,"count_orig"]<-rpois(1,exp(counts_orig[samp]))
              index2=index2+1
            }
          }
        }
        index=index+1
      }
    }
  }
  for (i in 1:nrow(counts_df)){
    if (counts_df$treatment[i]=='intervention') counts_df$treatment[i]=1
    if (counts_df$treatment[i]=='control') counts_df$treatment[i]=0
  }
  return(counts_df)
}

}
