source("code/plot_dellog2fc_transcript_protein.R")


# averaging rows with the same gene symbol
avg_del_log2_fc <- protein_vs_transcript_coord %>% 
  group_by(Symbol) %>% 
  summarize_all(funs(mean), na.rm=T)

# function for creating matrix of upregulation (1)/downregulation (-1)/not detected (0)
upreg_or_downreg <- function(del_log2_fc_col) {
  ifelse(del_log2_fc_col < 0, -1, 1)
}
# checking to make sure the function works (it does)
upreg_or_downreg(avg_del_log2_fc$del_log2_fc_wcl2)

# applying the function across all of the cols of del_log2_fc
upreg_or_downreg_matrix <- as.tibble(apply(avg_del_log2_fc[2:length(colnames(avg_del_log2_fc))], 2, FUN = upreg_or_downreg))
# averaging the cols together to get a single metric
cumulative_upreg_or_downreg <- as.tibble(apply(upreg_or_downreg_matrix, 1, FUN = mean, na.rm = T))

# making a dataframe of changes in flux
cumulative_exp_change <- avg_del_log2_fc %>%
  select(Symbol) %>% 
  mutate(exp_change = cumulative_upreg_or_downreg$value) # need to include colname even though there's only one

#View(cumulative_exp_change)
#write_delim(cumulative_exp_change, "data/cumulative_exp_change.txt", delim = "\t")



# selecting only genes that had a cumulative change in expression != 0
#cumulative_exp_change_nonzero <- cumulative_exp_change %>% 
  filter(exp_change != 0)

#View(cumulative_exp_change_nonzero)

# writing to a file to be used for manual data annotation (pos/neg regulator, which pathways, references)
#write_delim(cumulative_exp_change_nonzero, "data/cumulative_exp_change_nonzero.txt", delim = "\t")
