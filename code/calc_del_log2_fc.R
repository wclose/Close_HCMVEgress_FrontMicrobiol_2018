# import variables from sam probes and weekes proteins
source("code/import_pellett_sam_probe.R")
source("code/import_weekes_protein.R")

#####

# calculating dellog2fc for transcripts

# calculate averages from triplacte samples
sam_probe_data$mock_avg <- rowMeans(sam_probe_data[c("A1.AVG_Signal", "B1.AVG_Signal", "C1.AVG_Signal")], na.rm = T)
sam_probe_data$towne_12_avg <- rowMeans(sam_probe_data[c("A2.AVG_Signal", "B2.AVG_Signal", "C2.AVG_Signal")], na.rm = T)
sam_probe_data$towne_96_avg <- rowMeans(sam_probe_data[c("A8.AVG_Signal", "B8.AVG_Signal", "C8.AVG_Signal")], na.rm = T)
sam_probe_data$ad169_12_avg <- rowMeans(sam_probe_data[c("A6.AVG_Signal", "B6.AVG_Signal", "C6.AVG_Signal")], na.rm = T)
sam_probe_data$ad169_96_avg <- rowMeans(sam_probe_data[c("A11.AVG_Signal", "B11.AVG_Signal", "C11.AVG_Signal")], na.rm = T)

# normalize data relative to mock and log2 transform
sam_probe_data$log2_towne_12_mock <- log2(sam_probe_data$towne_12_avg / sam_probe_data$mock_avg)
sam_probe_data$log2_towne_96_mock <- log2(sam_probe_data$towne_96_avg / sam_probe_data$mock_avg)
sam_probe_data$log2_ad169_12_mock <- log2(sam_probe_data$ad169_12_avg / sam_probe_data$mock_avg)
sam_probe_data$log2_ad169_96_mock <- log2(sam_probe_data$ad169_96_avg / sam_probe_data$mock_avg)

# subtract 12 hpi from 96 hpi to get dellog2fc for each strain
sam_probe_data$del_log2_fc_towne <- sam_probe_data$log2_towne_96_mock - sam_probe_data$log2_towne_12_mock
sam_probe_data$del_log2_fc_ad169 <- sam_probe_data$log2_ad169_96_mock - sam_probe_data$log2_ad169_12_mock

# check to make sure cols/rows look correct
#head(sam_probe_data)
#tail(sam_probe_data)

#####

# calculating dellog2fc for proteins

# average replicates (only mock has more than one sample)
weekes_data_host$pm2_mock_avg <- rowMeans(weekes_data_host[c("pm2_Mock_1", "pm2_Mock_2")], na.rm = T)
weekes_data_host$wcl2_mock_avg <- rowMeans(weekes_data_host[c("wcl2_Mock_1", "wcl2_Mock_2")], na.rm = T)

# normalize data relative to mock and log 2 transform
# only need to do for 12 and 96 hpi samples
weekes_data_host$log2_pm2_12_mock <- log2(weekes_data_host$pm2_12h / weekes_data_host$pm2_mock_avg)
weekes_data_host$log2_pm2_96_mock <- log2(weekes_data_host$pm2_96h / weekes_data_host$pm2_mock_avg)
weekes_data_host$log2_wcl2_12_mock <- log2(weekes_data_host$wcl2_12h / weekes_data_host$wcl2_mock_avg)
weekes_data_host$log2_wcl2_96_mock <- log2(weekes_data_host$wcl2_96h / weekes_data_host$wcl2_mock_avg)

# subtract 12 hpi from 96 hpi to get dellog2fc for each strain
weekes_data_host$del_log2_fc_pm2 <- weekes_data_host$log2_pm2_96_mock - weekes_data_host$log2_pm2_12_mock
weekes_data_host$del_log2_fc_wcl2 <- weekes_data_host$log2_wcl2_96_mock - weekes_data_host$log2_wcl2_12_mock

# check to make sure cols/rows look correct
#head(weekes_data_host)
#tail(weekes_data_host)
