
reproduction <- function(Diploid_ind, Tetraploid_ind, h, L, Big_U, selfing_rate, N){
  
  
  if(nrow(Diploid_ind) < (2*N)){
    vector_proba_sampling = c(rep(1, times = nrow(Diploid_ind)/2), rep(2, times = nrow(Tetraploid_ind)/4))}
  else{vector_proba_sampling = c(rep(1, times = nrow(Diploid_ind)/2))}
  vector_proba_unred_gam = proba_unreduced_gametes(h, L,N,Diploid_ind, Tetraploid_ind)
  
  ## Reproduction
  
  Nb_offsprings = 0
  
  diplo_ind_temp = c(NULL)
  tetra_ind_temp = c(NULL)
  
  repeat{
    
    Parent_1 = sample(1:N, 1, replace = F, prob = vector_proba_sampling)
    if(Parent_1 <= (nrow(Diploid_ind)/2)){Ploidy_parent1 = "diploid"}
    else{Ploidy_parent1 = "tetraploid"}
    
    if(Ploidy_parent1 == "diploid" && vector_proba_unred_gam[Parent_1] <= runif(1, 0, 1)){Gamete_parent1 = "haplo"}
    else if(Ploidy_parent1 == "diploid" && vector_proba_unred_gam[Parent_1] >= runif(1, 0, 1)){Gamete_parent1 = "diplo"}
    else if(Ploidy_parent1 == "tetraploid" && vector_proba_unred_gam[Parent_1] <= runif(1, 0, 1)){Gamete_parent1 = "diplo"}
    else{Gamete_parent1 = "tetra"}
    
    ### Selfing or outcrossing
    
    if(selfing_rate <= runif(1, 0, 1)){
      
      repeat{Parent_2 = sample(1:N, 1, replace = F, prob = vector_proba_sampling)
      if(Parent_2 != Parent_1){break}
      }
      
      if(Parent_2 <= (nrow(Diploid_ind)/2)){Ploidy_parent2 = "diploid"}
      else{Ploidy_parent2 = "tetraploid"}
      
      if(Ploidy_parent2 == "diploid" && vector_proba_unred_gam[Parent_2] <= runif(1, 0, 1)){Gamete_parent2 = "haplo"}
      else if(Ploidy_parent2 == "diploid" && vector_proba_unred_gam[Parent_2] >= runif(1, 0, 1)){Gamete_parent2 = "diplo"}
      else if(Ploidy_parent2 == "tetraploid" && vector_proba_unred_gam[Parent_2] <= runif(1, 0, 1)){Gamete_parent2 = "diplo"}
      else{Gamete_parent2 = "tetra"}
    }
    else{
      Parent_2 = Parent_1
      Ploidy_parent2 = Ploidy_parent1
      
      if(Ploidy_parent1 == "diploid" && vector_proba_unred_gam[Parent_2] <= runif(1, 0, 1)){Gamete_parent2 = "haplo"}
      else if(Ploidy_parent1 == "diploid" && vector_proba_unred_gam[Parent_2] >= runif(1, 0, 1)){Gamete_parent2 = "diplo"}
      else if(Ploidy_parent1 == "tetraploid" && vector_proba_unred_gam[Parent_2] <= runif(1, 0, 1)){Gamete_parent2 = "diplo"}
      else{Gamete_parent2 = "tetra"}
    }
    
    if(Ploidy_parent1 == "diploid" && Ploidy_parent2 == "diploid" && Gamete_parent1 == "haplo" && Gamete_parent2 == "haplo"){
      Nb_offsprings = Nb_offsprings + 1
      
      # Gametes 1 & 2
      
      Gamete_1_temp = c(NULL)
      Gamete_2_temp = c(NULL)
      
      for(l in 1:L){
        
        temp_ind_1 = c(Diploid_ind[2*Parent_1,l],Diploid_ind[(2*Parent_1)-1,l])
        temp_ind_2 = c(Diploid_ind[2*Parent_2,l],Diploid_ind[(2*Parent_2)-1,l])
        
        Gamete_1_temp[l] = sample(temp_ind_1, 1, replace = F)
        Gamete_2_temp[l] = sample(temp_ind_2, 1, replace = F)
        
      }
      
      Gamete_1_mut = mutation(Big_U, L, Gamete_1_temp)
      Gamete_2_mut = mutation(Big_U, L, Gamete_2_temp)
      
      ploidy = c(rep("Diploid", times = 2))
      
      Offspirng_temp_fin = rbind(Gamete_1_mut,Gamete_2_mut)
      Offspirng_fin = cbind(Offspirng_temp_fin, ploidy)
      
      if(length(diplo_ind_temp) == 0){diplo_ind_temp = Offspirng_fin}
      else{diplo_ind_temp = rbind(diplo_ind_temp,Offspirng_fin)}
    }
    else if(Ploidy_parent1 == "diploid" && Ploidy_parent2 == "diploid" && Gamete_parent1 == "diplo" && Gamete_parent2 == "diplo"){
      Nb_offsprings = Nb_offsprings + 1
      
      # Gametes 1 & 2
      
      Gamete_1.1_temp = Diploid_ind[((2*Parent_1)-1),1:L]
      Gamete_1.2_temp = Diploid_ind[(2*Parent_1),1:L]
      Gamete_2.1_temp = Diploid_ind[((2*Parent_2)-1),1:L]
      Gamete_2.2_temp = Diploid_ind[(2*Parent_2),1:L]
      
      Gamete_1.1_mut = mutation(Big_U, L, Gamete_1.1_temp)
      Gamete_1.2_mut = mutation(Big_U, L, Gamete_1.2_temp)
      Gamete_2.1_mut = mutation(Big_U, L, Gamete_2.1_temp)
      Gamete_2.2_mut = mutation(Big_U, L, Gamete_2.2_temp)
      
      ploidy = c(rep("Tetraploid", times = 4))
      
      Offspirng_temp_fin = rbind(Gamete_1.1_mut,Gamete_1.2_mut,Gamete_2.1_mut,Gamete_2.2_mut)
      Offspirng_fin = cbind(Offspirng_temp_fin, ploidy)
      
      if(length(tetra_ind_temp) == 0){tetra_ind_temp = Offspirng_fin}
      else{tetra_ind_temp = rbind(tetra_ind_temp,Offspirng_fin)}
    }
    else if(Ploidy_parent1 != Ploidy_parent2 && Gamete_parent1 == Gamete_parent2){
      Nb_offsprings = Nb_offsprings + 1
      
      if(Ploidy_parent1 == "diploid"){
        
        Gamete_1.1_temp = Diploid_ind[((2*Parent_1)-1),1:L]
        Gamete_1.2_temp = Diploid_ind[(2*Parent_1),1:L]
        
        Gamete_2.1_temp = c(NULL)
        Gamete_2.2_temp = c(NULL)
        
        for(l in 1:L){
          
          temp_ind_2 = c(Tetraploid_ind[4*(Parent_2-(nrow(Diploid_ind)/2)),l],Tetraploid_ind[(4*(Parent_2-(nrow(Diploid_ind)/2)))-3,l],Tetraploid_ind[(4*(Parent_2-(nrow(Diploid_ind)/2)))-2,l],Tetraploid_ind[(4*(Parent_2-(nrow(Diploid_ind)/2)))-1,l])
          
          Gamete_2_temp = c(sample(temp_ind_2, 2, replace = F))
          
          Gamete_2.1_temp[l] = Gamete_2_temp[1]
          Gamete_2.2_temp[l] = Gamete_2_temp[2]
        }
        
        Gamete_1.1_mut = mutation(Big_U, L, Gamete_1.1_temp)
        Gamete_1.2_mut = mutation(Big_U, L, Gamete_1.2_temp)
        Gamete_2.1_mut = mutation(Big_U, L, Gamete_2.1_temp)
        Gamete_2.2_mut = mutation(Big_U, L, Gamete_2.2_temp)
        
        ploidy = c(rep("Tetraploid", times = 4))
        
        Offspirng_temp_fin = rbind(Gamete_1.1_mut,Gamete_1.2_mut,Gamete_2.1_mut,Gamete_2.2_mut)
        Offspirng_fin = cbind(Offspirng_temp_fin, ploidy)
        
        if(length(tetra_ind_temp) == 0){tetra_ind_temp = Offspirng_fin}
        else{tetra_ind_temp = rbind(tetra_ind_temp,Offspirng_fin)}
      }
      else{
        Gamete_2.1_temp = Diploid_ind[((2*Parent_2)-1),1:L]
        Gamete_2.2_temp = Diploid_ind[(2*Parent_2),1:L]
        
        Gamete_1.1_temp = c(NULL)
        Gamete_1.2_temp = c(NULL)
        
        for(l in 1:L){
          
          temp_ind_1 = c(Tetraploid_ind[4*(Parent_1-(nrow(Diploid_ind)/2)),l],Tetraploid_ind[(4*(Parent_1-(nrow(Diploid_ind)/2)))-3,l],Tetraploid_ind[(4*(Parent_1-(nrow(Diploid_ind)/2)))-2,l],Tetraploid_ind[(4*(Parent_1-(nrow(Diploid_ind)/2)))-1,l])
          
          Gamete_1_temp = c(sample(temp_ind_1, 2, replace = F))
          
          Gamete_1.1_temp[l] = Gamete_1_temp[1]
          Gamete_1.2_temp[l] = Gamete_1_temp[2]
        }
        
        Gamete_1.1_mut = mutation(Big_U, L, Gamete_1.1_temp)
        Gamete_1.2_mut = mutation(Big_U, L, Gamete_1.2_temp)
        Gamete_2.1_mut = mutation(Big_U, L, Gamete_2.1_temp)
        Gamete_2.2_mut = mutation(Big_U, L, Gamete_2.2_temp)
        
        ploidy = c(rep("Tetraploid", times = 4))
        
        Offspirng_temp_fin = rbind(Gamete_1.1_mut,Gamete_1.2_mut,Gamete_2.1_mut,Gamete_2.2_mut)
        Offspirng_fin = cbind(Offspirng_temp_fin, ploidy)
        
        if(length(tetra_ind_temp) == 0){tetra_ind_temp = Offspirng_fin}
        else{tetra_ind_temp = rbind(tetra_ind_temp,Offspirng_fin)}}
    }
    else if(Ploidy_parent1 == "tetraploid" && Ploidy_parent2 == "tetraploid" && Gamete_parent1 == "diplo" && Gamete_parent2 == "diplo"){
      Nb_offsprings = Nb_offsprings + 1
      
      Gamete_1.1_temp = c(NULL)
      Gamete_1.2_temp = c(NULL)
      
      Gamete_2.1_temp = c(NULL)
      Gamete_2.2_temp = c(NULL)
      
      for(l in 1:L){
        
        temp_ind_1 = c(Tetraploid_ind[4*(Parent_1-(nrow(Diploid_ind)/2)),l],Tetraploid_ind[(4*(Parent_1-(nrow(Diploid_ind)/2)))-3,l],Tetraploid_ind[(4*(Parent_1-(nrow(Diploid_ind)/2)))-2,l],Tetraploid_ind[(4*(Parent_1-(nrow(Diploid_ind)/2)))-1,l])
        temp_ind_2 = c(Tetraploid_ind[4*(Parent_2-(nrow(Diploid_ind)/2)),l],Tetraploid_ind[(4*(Parent_2-(nrow(Diploid_ind)/2)))-3,l],Tetraploid_ind[(4*(Parent_2-(nrow(Diploid_ind)/2)))-2,l],Tetraploid_ind[(4*(Parent_2-(nrow(Diploid_ind)/2)))-1,l])
        
        Gamete_1_temp = c(sample(temp_ind_1, 2, replace = F))
        Gamete_2_temp = c(sample(temp_ind_2, 2, replace = F))
        
        Gamete_1.1_temp[l] = Gamete_1_temp[1]
        Gamete_1.2_temp[l] = Gamete_1_temp[2]
        
        Gamete_2.1_temp[l] = Gamete_2_temp[1]
        Gamete_2.2_temp[l] = Gamete_2_temp[2]
      }
      
      Gamete_1.1_mut = mutation(Big_U, L, Gamete_1.1_temp)
      Gamete_1.2_mut = mutation(Big_U, L, Gamete_1.2_temp)
      Gamete_2.1_mut = mutation(Big_U, L, Gamete_2.1_temp)
      Gamete_2.2_mut = mutation(Big_U, L, Gamete_2.2_temp)
      
      ploidy = c(rep("Tetraploid", times = 4))
      
      Offspirng_temp_fin = rbind(Gamete_1.1_mut,Gamete_1.2_mut,Gamete_2.1_mut,Gamete_2.2_mut)
      Offspirng_fin = cbind(Offspirng_temp_fin, ploidy)
      
      if(length(tetra_ind_temp) == 0){tetra_ind_temp = Offspirng_fin}
      else{tetra_ind_temp = rbind(tetra_ind_temp,Offspirng_fin)}
    }
    
    if(Nb_offsprings == N){break}   
  }
  
  All_ind_temp = rbind(diplo_ind_temp, tetra_ind_temp)  
  
  return(All_ind_temp)
  
}

proba_unreduced_gametes <-function(h, L, N,Diploid_ind, Tetraploid_ind){
  
  list_proba_unreduced = c(NULL)
  
  k = 0
  if(nrow(Diploid_ind) !=0){
    for(i in 1:(nrow(Diploid_ind)/2)){
      k = k + 1
      
      val_temp = 0
      
      for(j in 1:L){
        if(Diploid_ind[2*i, j] == Diploid_ind [(2*i)-1, j]){val_temp = val_temp + (as.numeric(Diploid_ind[2*i, j])*(1/L))}
        else{val_temp = val_temp + h*(1/L)}
      }
      
      list_proba_unreduced[k] = val_temp
    }
  }
  
  if(k < N && is.null(nrow(Tetraploid_ind)) == F){
    for(i in 1:(nrow(Tetraploid_ind)/4)){
      k = k + 1
      
      val_temp = 0
      
      for(j in 1:L){
        if(Tetraploid_ind[4*i, j] == Tetraploid_ind[(4*i)-1, j] && Tetraploid_ind[4*i, j] == Tetraploid_ind[(4*i)-2, j] && Tetraploid_ind[4*i, j] == Tetraploid_ind[(4*i)-3, j]){val_temp = val_temp + (as.numeric(Tetraploid_ind[4*i, j])*(1/L))}
        else{val_temp = val_temp + h*(1/L)}
      }
      
      list_proba_unreduced[k] = val_temp
    }
  }
  
  return(list_proba_unreduced)
  
}

mutation <- function(Big_U, L, haplotype){
  Nb_mut = rpois(1, Big_U)
  
  if(Nb_mut > L){Nb_mut = L}
  
  position_mut = c(sample(1:L, Nb_mut, replace = F))
  if(Nb_mut != 0){
    for(i in 1:length(position_mut)){
      if(haplotype[position_mut[i]] != "0"){haplotype[position_mut[i]] = "0"}
      else{haplotype[position_mut[i]] = "1"}
    }
  }
  return(haplotype)
}

count_allele_pop_1 <- function(All_ind, L){
  
  freq_per_loc = c(NULL)
  for(i in 1:L){
    freq_per_loc[i] = sum(as.numeric(All_ind[,i]))/nrow(All_ind)
  }
  return(mean(freq_per_loc))
}

count_allele_diplo_1 <- function(Diploid_ind, L){
  
  freq_per_loc = c(NULL)
  for(i in 1:L){
    freq_per_loc[i] = sum(as.numeric(Diploid_ind[,i]))/nrow(Diploid_ind)
  }
  return(mean(freq_per_loc))
}

count_tetraploids <- function(Tetraploid_ind, L,N){
  
  nb_tetra = (nrow(Tetraploid_ind)/4)/N
  
  return(nb_tetra)
}

#L=5
#mut_loc=0.05
#N=10
#selfing_rate=0
#h=0.5
#Nb_gen = 100
#Nb_rep = 15

simulation_model <-function(L, mut_loc, N, Nb_gen, selfing_rate, h, Nb_rep){
  
  rep = 0
  
  freq_all_1 = c(NULL)
  freq_diplo_1 = c(NULL)
  freq_tetra = c(NULL)
  
  for (r in 1:Nb_rep){
    
    rep = rep + 1
    
    # Initialization 
    
    ## genomic mutation rate
    
    Big_U = mut_loc
    
    ## At t = 0, the population is diploid, and only produce haploid gametes
    
    haplotype_t0 = rep(0, times = L)
    
    for(m in 1:(2*N)){
      if(m == 1){All_ind = haplotype_t0}
      else{All_ind = rbind(All_ind,haplotype_t0)}
    }
    
    Ploidy = c(rep("Diploid", times = 2*N))
    All_ind = cbind(All_ind, Ploidy)
    
    Diploid_ind = subset(All_ind, All_ind[,L+1]=="Diploid")
    Tetraploid_ind = c(NULL)
    
    # Vector of mean frequency of allele 1, and number of tetraploids
    
    freq_all_pop_0 = c(NULL)
    freq_all_diplo_0 = c(NULL)
    freq_tetraploids = c(NULL)
    
    # Functions
    
    for(t in 1:Nb_gen){
      
      # Reproduction
      
      All_ind = reproduction(Diploid_ind, Tetraploid_ind, h, L, Big_U, selfing_rate, N)
      
      # Seperation of diploids and tetraploids individuals
      
      Diploid_ind = subset(All_ind, All_ind[,L+1]=="Diploid")
      Tetraploid_ind = subset(All_ind, All_ind[,L+1]=="Tetraploid")
      nrow(Diploid_ind)/2
      nrow(Tetraploid_ind)/4
      
      
      
      # Count of alleles 1 and tetraploids
      
      freq_all_pop_0[t] = count_allele_pop_1(All_ind, L)
      freq_all_diplo_0[t] = count_allele_diplo_1(Diploid_ind, L)
      freq_tetraploids[t] = count_tetraploids(Tetraploid_ind, L, N)
      
    }
    
    if(rep == 1){freq_all_1 = freq_all_pop_0
    freq_diplo_1 = freq_all_diplo_0
    freq_tetra = freq_tetraploids}
    else{freq_all_1 = rbind(freq_all_1,freq_all_pop_0)
    freq_diplo_1 = rbind(freq_diplo_1,freq_all_diplo_0)
    freq_tetra = rbind(freq_tetra,freq_tetraploids)}
    
  }
  
  write.table(freq_all_1, file=paste("freq_all1_s",selfing_rate,"_h",h,"_N",N,"_u",mut_loc,"_L",L,".txt",sep=""),row.names = F, dec = ".")
  write.table(freq_diplo_1, file=paste("freq_diplo1_s",selfing_rate,"_h",h,"_N",N,"_u",mut_loc,"_L",L,".txt",sep=""),row.names = F, dec = ".")
  write.table(freq_tetra, file=paste("freq_tetra_s",selfing_rate,"_h",h,"_N",N,"_u",mut_loc,"_L",L,".txt",sep=""),row.names = F, dec = ".")
}

simulation_model(2, 0.1, 100, 20000, 0, 0.5, 10)