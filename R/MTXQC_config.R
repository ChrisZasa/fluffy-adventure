### Import of config_files

#conversion list    
con_se = read.csv("config_mtx/conversion_metabolite.csv", TRUE)

if (ncol(con_se) == 1) {
  con_se = read.csv("config_mtx/conversion_metabolite.csv", T, sep = ";")
  message("Colon separator detected in : conversion_metabolite.csv")
}

file_spec = read.csv("config_files/conv_filenames.csv", TRUE)
file_spec_manVal = read.csv("config_files/conv_filenames_manVal.csv", TRUE)

#Metabolic profile: relation of intermediates and pathways
pathway_profile  =  read.csv('config_mtx/letter_pathway_complete.csv', T)
pathway_profile = pathway_profile[, 1:7]


#reference file summarising mass of interest
mass_li = read.csv('config_mtx/incorp_calc_masses.csv',  T)
backup_mids = read.csv('config_mtx/mid_backups.csv', T)


