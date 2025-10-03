library(tidyverse)
library(scales)
library(RMariaDB)
library(DBI)

# Connect
mysql_password = '***'
con <- dbConnect(
  drv = RMariaDB::MariaDB(),
  dbname = 'chembl_36',
  username = 'root',
  password = mysql_password,
  host = NULL, 
  port = 3306
)
#Check the connection
dbListTables(con)

#Get the data
ic50__query <- dbSendQuery(con, "SELECT a.activity_id, md.chembl_id, a.molregno, td.chembl_id, td.pref_name, assa.tid, a.standard_value, a.standard_relation, a.standard_units, a.activity_comment, a.data_validity_comment
												FROM activities a JOIN assays assa JOIN target_dictionary td JOIN molecule_dictionary md
												WHERE a.molregno = md.molregno AND
												a.assay_id = assa.assay_id AND
												assa.tid = td.tid AND
												a.standard_type = 'IC50' AND
												a.standard_relation = '='")
ic50_rslt_raw <- dbFetch(ic50__query)
#Prepare names and delete the most suspicious records
ic50_rslt <- ic50_rslt_raw |> repair_names() |> rename(target_ID = chembl_id1, compound_ID = chembl_id, compound_molregno = molregno, target_name = pref_name) |>
															mutate(general_data_validity_comment = if_else(is.na(data_validity_comment) | data_validity_comment == "Manually validated", "OK", "notOK")) |>
															filter(standard_value > 0 & !is.na(standard_value) & !is.nan(standard_value))
dbClearResult(ic50__query)
#Assess numbers of measurements for target-compound pairs
pairs <- ic50_rslt |> select(target_ID, compound_ID) |>
						group_by(target_ID, compound_ID) |>
						summarise(n_of_measurements=n(), .groups = 'drop') |>
						group_by(n_of_measurements) |>
						summarise(n_of_n = sum(n_of_measurements), .groups = 'drop') |>
						mutate(category = if_else(n_of_measurements > 10, '>10', as.factor(n_of_measurements))) |>
						select(category, n_of_n) |>
						group_by(category) |>
						summarise(n = sum(n_of_n), .groups = 'drop')
sum_measurements <- pairs |> pull(n) |> sum()
#Convert to pcts and arrange records for plotting using factor relveling
pairs_pct <- pairs |> mutate(pct = (n/sum_measurements)*100) |>
						select(category, pct) |>
						mutate(category = fct_relevel(category, '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '>10'))
#Plot
basic_summary_plot <- ggplot(pairs_pct, aes(x=category, y=pct)) + 
  								geom_bar(stat = "identity", fill = "#515095FF") +
  								scale_y_continuous(name="Процент записей", breaks=c(0,10,20,30,40,50,60,70,80,90,100), limits=c(0,100)) +
  								scale_x_discrete(name="Количество экспериментов") +
  								labs(title = "Данные ChEMBL v36", subtitle = str_glue("{format(sum_measurements, big.mark = '\u2009')} записи IC\u2085\u2080")) +
  								theme_classic() +
  								theme(plot.title = element_text(hjust = .95), plot.subtitle = element_text(hjust = .95))
basic_summary_plot
sqrt_summary_plot <- ggplot(pairs_pct, aes(x=category, y=pct)) + 
  								geom_bar(stat = "identity", fill = "#515095FF") +
  								scale_y_sqrt(name="Процент записей", breaks=c(0,1,5,10,20,30,40,50,60,70,80,90,100), limits=c(0,100)) +
  								scale_x_discrete(name="Количество экспериментов") +
  								labs(title = "Данные ChEMBL v36", subtitle = str_glue("{format(sum_measurements, big.mark = '\u2009')} записи IC\u2085\u2080")) +
  								theme_classic() +
  								theme(plot.title = element_text(hjust = .95), plot.subtitle = element_text(hjust = .95))
sqrt_summary_plot

#Stat on units
units_stat <- ic50_rslt |> select(standard_units) |> group_by(standard_units) |>
															summarize(n = n()) |>
															arrange(desc(n))

#IC50 total distribution, the usage of histogram is problematic, since the reasonable binsize does not allow to cover the whole range due to the technical limitations: The number of histogram bins must be less than 1,000,000.
ic50_vals <- ic50_rslt |> select(standard_value, standard_units) |> filter(!is.na(standard_value) & (standard_units=="nM")) |>
														mutate(pIC50 = log10( 1 / (standard_value*10^-9) ) )
ic50_vals_median <- ic50_vals |> pull(pIC50) |> median() |> round(digits = 1)
ic50_vals_mean <- ic50_vals |> pull(pIC50) |> mean() |> round(digits = 1)
#Prepare pIC50 intervals
ic50_intervals <- ic50_vals |> pull(pIC50) |> cut_width(width = 1, center = 0)
interval_stat <- ic50_intervals |> as_tibble() |> group_by(value) |> summarise(n = n()) |> rename(interval = value)
#Plot IC50
base_ic50_nm_plot <- ggplot(interval_stat, aes(x=interval, y=n)) +
										geom_bar(stat = "identity", fill = "#515095FF") +
										scale_y_continuous(name="Количество измерений", breaks=c(2000, 100000, 200000, 400000, 600000), labels = label_comma(big.mark = '\u2009')) +
										scale_x_discrete(name="Интервал pIC\u2085\u2080") +
										labs(title = "Данные ChEMBL v36", subtitle = str_glue("{format(interval_stat |> pull(n) |> sum(), big.mark = '\u2009')} записей pIC\u2085\u2080\u000A медиана = {ic50_vals_median}, среднее арифметическое = {ic50_vals_mean}")) +
										coord_flip() +
										theme_classic() +
										theme(plot.title = element_text(hjust = .95), plot.subtitle = element_text(hjust = .95))
base_ic50_nm_plot
sqrt_ic50_nm_plot <- ggplot(interval_stat, aes(x=interval, y=n)) +
										geom_bar(stat = "identity", fill = "#515095FF") +
										scale_y_sqrt(name="Количество измерений", breaks=c(1, 2000, 100000, 200000, 400000, 600000), labels = label_comma(big.mark = '\u2009')) +
										scale_x_discrete(name="Интервал pIC\u2085\u2080") +
										labs(title = "Данные ChEMBL v36", subtitle = str_glue("{format(interval_stat |> pull(n) |> sum(), big.mark = '\u2009')} записей pIC\u2085\u2080\u000A медиана = {ic50_vals_median}, среднее арифметическое = {ic50_vals_mean}")) +
										coord_flip() +
										theme_classic() +
										theme(plot.title = element_text(hjust = .95), plot.subtitle = element_text(hjust = .95))
sqrt_ic50_nm_plot

#Calculate mean and median values
ic50_avg <- ic50_rslt |> group_by(compound_ID, target_ID, standard_units, general_data_validity_comment) |>
												 summarize(n = n(),
												 						median_ic50 = median(standard_value |> na.omit()),
												 						mean_ic50 = mean(standard_value |> na.omit()), 
												 						target_name = unique(target_name),
												 						activity_comments = str_c( activity_comment  |> na.omit() |> unique(), collapse = ", ")) |>
												 ungroup() |>
												 mutate(median_to_mean = median_ic50 / mean_ic50)

#Visualize median to mean relation
ic50_medmean <- ic50_avg |> mutate(n = if_else(n > 10, as.factor(">10"), as.factor(n))) |>
														mutate(n = fct_relevel(n, '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '>10')) |>
														select(n, median_to_mean) |>
														filter(!is.na(median_to_mean)) |>
														filter(median_to_mean != 'NaN')
x <- ic50_medmean |> select(n) |> distinct()
y <- ic50_medmean |> select(median_to_mean) |> distinct()
medmean_plot <- ggplot(ic50_medmean, aes(x=n, y=median_to_mean, fill = n)) +
    						scale_x_discrete(name="Количество измерений IC\u2085\u2080", breaks=c('1','2','3','4','5','6','7','8','9','10','>10')) +
    						scale_fill_manual(values = c("#8281AFFF", "#FFEC9DFF", "#FAC881FF", "#F4A464FF", "#E87444FF", "#D9402AFF", "#BF2729FF", "#912534FF", "#64243EFF", "#3D1B28FF", "#161212FF")) +
    						scale_y_continuous(name="Отношение медианы к арифметическому среднему") +
    						labs(title = "Данные ChEMBL v36") +
    						geom_boxplot() +
    						theme_classic() +
    						theme(plot.title = element_text(hjust = .95)) +
    						theme(legend.position = "none")
medmean_plot

#Export the results, tables
write_tsv(ic50_rslt, ".../output_2025/IC50_values.tsv")
write_tsv(pairs_pct |> arrange(category) |> rename(`Number of measurements` = category, `%` = pct), ".../output_2025/Number_of_measurements.tsv")
write_tsv(interval_stat |> rename(`Interval of pIC50` = interval, `Number of measurements` = n), ".../output_2025/IC50_nM_intervals.tsv")
write_tsv(ic50_avg, ".../output_2025/IC50_averaged.tsv")
#Export the results, plots
ggsave(".../output_2025/1_Number_of_measurements__ChEMBL_v36.png", plot = basic_summary_plot, width = 7, height = 5, units = "in", dpi = 300)
ggsave(".../output_2025/1_sqrt_Number_of_measurements__ChEMBL_v36.png", plot = sqrt_summary_plot, width = 7, height = 5, units = "in", dpi = 300)
ggsave(".../output_2025/2_Number_of_measurements_in_intervals__ChEMBL_v36.png", plot = base_ic50_nm_plot, width = 7, height = 5, units = "in", dpi = 300)
ggsave(".../output_2025/2_sqrt_Number_of_measurements_in_intervals__ChEMBL_v36.png", plot = sqrt_ic50_nm_plot, width = 7, height = 5, units = "in", dpi = 300)
ggsave(".../output_2025/3_Median_to_Mean__ChEMBL_v36.png", plot = medmean_plot, width = 7, height = 5, units = "in", dpi = 300)