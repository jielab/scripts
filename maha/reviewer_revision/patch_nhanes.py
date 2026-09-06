from pathlib import Path
p=Path(__file__).resolve().parents[1]/'f/MAHA_nhanes.R'
s=p.read_text(encoding='utf-8')
def sub(a,b):
    global s
    assert a in s,a[:90]
    s=s.replace(a,b)
sub('source(file.path(.this_dir,"comm.f.R"))','source(file.path(.this_dir,"comm.f.R"))\nsource(file.path(.this_dir,"nhanes_meat_sensitivity.R"))')
sub('red_processed_meat = str_detect(fd, "beef|pork|lamb|bacon|sausage|ham|hot dog|frankfurter|pepperoni|salami"),', 'red_processed_meat = str_detect(fd, "beef|pork|lamb|bacon|sausage|ham|hot dog|frankfurter|pepperoni|salami"),\n\t\t\tprocessed_meat_proxy = ifelse(is.na(food_desc), NA, maha_processed_meat_flag(food_desc)),')
sub('"poultry","red_processed_meat","sweets_pastries"','"poultry","red_processed_meat","processed_meat_proxy","sweets_pastries"')
sub('upf_proxy = refined_grain + sweets_pastries + ssb + fried_fast + red_processed_meat,', 'upf_proxy = refined_grain + sweets_pastries + ssb + fried_fast + red_processed_meat,\n\t\t\tupf_processed_proxy = refined_grain + sweets_pastries + ssb + fried_fast + processed_meat_proxy,\n\t\t\tupf_nomeat_proxy = refined_grain + sweets_pastries + ssb + fried_fast,')
sub('maha_c_upf = qscore(upf_proxy, reverse = TRUE),','maha_c_upf = qscore(upf_proxy, reverse = TRUE),\n\t\t\tmaha_c_upf_processed = qscore(upf_processed_proxy, reverse = TRUE),\n\t\t\tmaha_c_upf_nomeat = qscore(upf_nomeat_proxy, reverse = TRUE),')
original='diet.maha.sum = rowmean_min(., c("maha_c_protein","maha_c_dairy","maha_c_veg","maha_c_fruit","maha_c_wholegrain","maha_c_fat","maha_c_upf","maha_c_alcohol","maha_c_sodium"), 0.60),'
sub(original, original+'\n\t\t\t'+original.replace('diet.maha.sum','diet.maha_processed.sum').replace('"maha_c_upf"','"maha_c_upf_processed"')+'\n\t\t\t'+original.replace('diet.maha.sum','diet.maha_nomeat.sum').replace('"maha_c_upf"','"maha_c_upf_nomeat"'))
sub('nhanes_construct <- nhanes_construct_profile(dat1_obs)', 'nhanes_construct <- nhanes_construct_profile(dat1_obs)\nwrite_xlsx(nhanes_meat_sensitivity(dat_mort_obs, covs.mort), "Reviewer.meat_proxy_sensitivity.xlsx")')
sub('spearman_matrix = as.data.frame(cor_mat, check.names = FALSE)', '''spearman_matrix = as.data.frame(cor_mat, check.names = FALSE),
    pairwise_N = as.data.frame(crossprod(!is.na(score_mat))),
    metadata = data.frame(N = nrow(score_mat), sample = "dat1_obs: observed repeated-recall mean, age >=20, eligible dietary weights; descriptive main cycles",
      cycles = paste(main_cycles, collapse = "; "), method = "Unweighted Spearman, pairwise complete observations",
      shared_figure = "Shared FigS1 and NHANES construct figure read this exact same correlation matrix")''')
p.write_text(s,encoding='utf-8')
