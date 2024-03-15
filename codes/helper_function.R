# fig setetings -----------------------------------------------------------
fig_setting <- fig_theme <- theme_classic() +
  theme(legend.position = 'top')

select <- dplyr::select

result_processor <- function(results,
                             baseline_year=2001,
                             start_year = 1,
                             end_year = 18,
                             start_age = 3,
                             end_age = 18,
                             max_age = 111,
                             chosen_province='BC',
                             growth_type = "M3",
                             pop_years = c(2001:2018),
                             plot_dims = c(3,4),
                             save_fig=F,
                             scenario=NULL){
  start_age <- start_age + 1
  end_age <- end_age + 1
  save_plot <- function(gg_obj,nam,fig.width=10,fig.height=7,fig.dpi=600,save=save_fig){
    if(save){
      ggsave(paste0(fig_dir,'/',nam),gg_obj,width = fig.width,height=fig.height,dpi = fig.dpi)
    } else{
      "plot not saved"
    }
  }
  
  longer_df <- function(type){
    tmp_df <- results$outcome_matrix[type][[1]]
    tmp_df_male <- tmp_df[,,2] %>% 
      as.data.frame() %>% 
      mutate(year = row_number()+baseline_year-1) %>% 
      pivot_longer(-year,names_to="age",values_to="n") %>% 
      mutate(age = as.numeric(str_remove(age,"V"))-1) %>% 
      mutate(sex="male")
    tmp_df_female <- tmp_df[,,1] %>% 
      as.data.frame() %>% 
      mutate(year = row_number()+baseline_year-1) %>% 
      pivot_longer(-year,names_to="age",values_to="n") %>% 
      mutate(age = as.numeric(str_remove(age,"V"))-1) %>% 
      mutate(sex="female")  
    
    return(rbind(tmp_df_male,tmp_df_female))
  }
  
  extract_df <- function(type,sex='both',year='all',age='all'){
    if(!is.numeric(sex)){
      sex = ifelse(sex=='male',2,
                   ifelse(sex=='female',1,
                          sex))
    }
    if(year[1]=='all'){
      year <- 1:max_year
    }
    if(age[1]=='all'){
      age <- 1:max_age
    }
    
    if(type %in% c("exacerbation_severity","exacerbation_by_severity")){
      tmp_result <- c()
      for (j in 1:4){
        if (sex=="both"){
          tmp_result[[j]] <- results[[2]][type][[1]][year,age,1,j] + results[[2]][type][[1]][year,age,2,j]
        } else{
          tmp_result[[j]] <- results[[2]][type][[1]][year,age,sex,j]
        }
      }
      
      tmp_result <- tmp_result %>% lapply(.,function(x){
        dim_x <- dim(x)
        matrix(sapply(x,identity),nrow =dim_x[1],dim_x[2])
      })
      
      return(tmp_result)
    } else if(type %in% c("control")){
      tmp_result <- c()
      for (j in 1:3){
        if (sex=="both"){
          tmp_result[[j]] <- results[[2]][type][[1]][year,age,1,j] + results[[2]][type][[1]][year,age,2,j]
        } else{
          tmp_result[[j]] <- results[[2]][type][[1]][year,age,sex,j]
        }
      }
      
      tmp_result <- tmp_result %>% lapply(.,function(x){
        dim_x <- dim(x)
        matrix(sapply(x,identity),nrow =dim_x[1],dim_x[2])
      })
      
      return(tmp_result)
    }
    else{
      if(sex=='both'){
        return(results[[2]][type][[1]][,,1][year,age]+results[[2]][type][[1]][,,2][year,age]) 
      } else{
        return(results[[2]][type][[1]][,,sex][year,age])
      }
    }
    
    
  }
  
  alive <- sum(results$outcome_matrix$alive[start_year:end_year,start_age:end_age,1:2])
  death <- sum(results$outcome_matrix$death[start_year:end_year,start_age:end_age,1:2])
  emigration <- sum(results$outcome_matrix$emigration[start_year:end_year,start_age:end_age,1:2])
  total_years <- alive + death + emigration
  
  asthma_inc <- sum(results$outcome_matrix$asthma_incidence[start_year:end_year,start_age:end_age,1:2])
  asthma_prev <- sum(results$outcome_matrix$asthma_prevalence[start_year:end_year,start_age:end_age,1:2])
  asthma_prev_2018 <- sum(results$outcome_matrix$asthma_prevalence[end_year,start_age:end_age,1:2])
  mild <- sum(results$outcome_matrix$exacerbation_by_severity[start_year:end_year,start_age:end_age,1:2,1])
  moderate <- sum(results$outcome_matrix$exacerbation_by_severity[start_year:end_year,start_age:end_age,1:2,2])
  severe <- sum(results$outcome_matrix$exacerbation_by_severity[start_year:end_year,start_age:end_age,1:2,3])
  very_severe <- sum(results$outcome_matrix$exacerbation_by_severity[start_year:end_year,start_age:end_age,1:2,4])

  max_year <- nrow(results[[1]])
  
  df_ABE <- rbind(data.frame(Year=1:max_year,
                             n=extract_df(type = 'antibiotic_exposure',sex='male',age=1),
                             N = extract_df(type = 'alive',sex='male',age=1)+
                               extract_df(type = 'death',sex='male',age=1)+
                               extract_df(type = 'emigration',sex='male',age=1)) %>% 
                    mutate(sex=1),
                  data.frame(Year=1:max_year,
                             n=extract_df(type = 'antibiotic_exposure',sex='female',age=1),
                             N = extract_df(type = 'alive',sex='female',age=1)+
                               extract_df(type = 'death',sex='female',age=1)+
                               extract_df(type = 'emigration',sex='female',age=1)) %>% 
                    mutate(sex=0))%>% 
    mutate(rate=n/N*1000) %>% 
    mutate(sex=ifelse(sex==0,"Female","Male")) %>% 
    mutate(Year=Year+baseline_year-1)
  
  df_inc <- data.frame(Year=1:max_year,
                       as.data.frame(results$outcome_matrix$asthma_incidence[start_year:end_year,start_age:end_age,1] +
                                       results$outcome_matrix$asthma_incidence[start_year:end_year,start_age:end_age,2])) %>% 
    pivot_longer(-Year,names_to="age",values_to="value") %>% 
    mutate(age = str_remove(age,"V"),
           age = as.numeric(age) + start_age - 2,
           Year = Year+baseline_year-1)
  
  df_prev <- data.frame(Year=1:max_year,
                       as.data.frame(results$outcome_matrix$asthma_prevalence[start_year:end_year,start_age:end_age,1] +
                                       results$outcome_matrix$asthma_prevalence[start_year:end_year,start_age:end_age,2])) %>% 
    pivot_longer(-Year,names_to="age",values_to="value") %>% 
    mutate(age = str_remove(age,"V"),
           age = as.numeric(age) + start_age - 2,
           Year = Year+baseline_year-1)
  
  df_exac <- data.frame(Year=1:max_year,
                        as.data.frame(results$outcome_matrix$exacerbation[start_year:end_year,start_age:end_age,1] +
                                        results$outcome_matrix$exacerbation[start_year:end_year,start_age:end_age,2])) %>% 
    pivot_longer(-Year,names_to="age",values_to="value") %>% 
    mutate(age = str_remove(age,"V"),
           age = as.numeric(age) + start_age - 2,
           Year = Year+baseline_year-1)
  
  return(list(outcomes = 
                c(num_year = total_years,
              num_asthma_inc = asthma_inc,
              num_asthma_total_years= asthma_prev,
              num_prev_2018 = asthma_prev_2018,
              num_mild_exac = mild,
              num_moderate_exac = moderate,
              num_severe_exac = severe,
              num_very_severe_exac = very_severe,
              scenario= scenario),
              df_inc = df_inc %>% 
                mutate(scenario= scenario),
              df_prev =df_prev %>% 
                mutate(scenario= scenario),
              df_abx = df_ABE %>% 
                mutate(scenario = scenario),
              df_exac = df_exac %>% 
                mutate(scenario = scenario)))
}

compute_var <- function(p,sample_size){
  p*(1-p)/sample_size
}



