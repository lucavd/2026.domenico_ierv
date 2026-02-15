##############################################################################
# SCRIPT COMPLETO PER SALVARE E INVIARE IMMAGINI E .DOCX SU TELEGRAM
##############################################################################

library(tidyverse)
library(gtsummary)
library(flextable)
library(telegram.bot)

# Imposta il bot e la chat_id
readRenviron(".env")
bot <- telegram.bot::Bot(token = Sys.getenv("TELEGRAM_BOT_TOKEN"))
chat_id <- Sys.getenv("TELEGRAM_CHAT_ID")

# Percorsi per salvare i file
img_path  <- "~/Projects/MMM/simulazioni_tesi/ARTICOLO/img"
docx_path <- "~/Projects/MMM/simulazioni_tesi/ARTICOLO"

##############################################################################
# 1. Creazione e manipolazione dei dataset come da tuo script
##############################################################################

db <- bind_rows(final_res_1000, .id = 'id') |>
  rownames_to_column() |>
  mutate(mape = 100*abs((est - beta)/est)) |>
  mutate(beta = ifelse(beta < 0.2, "Beta = 0.18", "Beta = 0.40")) |>
  mutate(across(c(ntime, id, beta), as.factor)) |>
  mutate(method = case_match(method,
                             'glmer' ~ 'GLMM',
                             'gee'   ~ 'GEE',
                             'mmm'   ~ 'MMM'),
         method = factor(method, levels = c('GLMM', 'GEE', 'MMM')))

db_time <- db |>
  filter(grepl("^time", rowname)) |>
  group_by(method, ntime, N, beta) |>
  summarise(
    mapemean   = mean(mape, na.rm = TRUE),
    mapesd     = sd(mape, na.rm = TRUE),
    mapeiqr    = IQR(mape, na.rm = TRUE),
    sigtrue    = 100 - (sum(sig, na.rm = TRUE)/length(sig))*100,
    mapemedian = median(mape, na.rm = TRUE)
  ) |>
  as.data.frame()

##############################################################################
# 2. Esempio di tabella con gtsummary --> .docx --> invio Telegram
##############################################################################

# tbl_time_filtered <- db_time |>
#   filter(beta == 'Beta = 0.18' & ntime == 20) |>
#   tbl_strata(
#     strata = N,
#     .tbl_fun = ~ .x |>
#       tbl_summary(
#         by       = method,
#         missing  = "no",
#         statistic = all_continuous() ~ '{mean}',
#         type = c(mapemean, mapesd, mapemedian, mapeiqr, sigtrue) ~ 'continuous'
#       ),
#     .header = "N = **{strata}**"
#   ) |>
#   modify_header(
#     label ~ "**Method**",
#     all_stat_cols() ~ "**{level}**"
#   ) |>
#   modify_caption("**MAPE Beta = 0.18, Clusters = 20**")
#
# # Salvo la tabella in formato docx
# docx_file <- file.path(docx_path, "mape_beta018_nclust20.docx")
# flextable::save_as_docx(tbl_time_filtered, path = docx_file)
#
# # Invio il file docx su Telegram
# bot$sendDocument(chat_id = chat_id, document = docx_file)

##############################################################################
# 3. Creazione e salvataggio dei vari grafici, invio su Telegram
##############################################################################

# Plot 1: Mean MAPE
p1 <- ggplot(db_time, aes(x = N, y = mapemean)) +
  geom_point(aes(colour = method)) +
  geom_line(aes(colour = method)) +
  facet_grid(beta ~ ntime) +
  labs(
    y = 'MAPE',
    colour = "Method",
    x = 'Samples number (ID)',
    subtitle = 'Number of clusters'
  ) +
  theme(plot.subtitle = element_text(hjust = 0.5))

img_file1 <- file.path(img_path, "1_mean_mape.png")
ggsave(filename = img_file1, plot = p1, device = "png", dpi = 600)
bot$sendPhoto(chat_id = chat_id, photo = img_file1)

# Plot 2: Standard Deviation of MAPE
p2 <- ggplot(db_time, aes(x = N, y = mapesd)) +
  geom_point(aes(colour = method)) +
  geom_line(aes(colour = method)) +
  facet_grid(beta ~ ntime) +
  labs(
    y = 'APE sd',
    colour = "Method",
    x = 'Samples number (ID)',
    subtitle = 'Number of clusters'
  ) +
  theme(plot.subtitle = element_text(hjust = 0.5))

img_file2 <- file.path(img_path, "2_sd_mape.png")
ggsave(filename = img_file2, plot = p2, device = "png", dpi = 600)
bot$sendPhoto(chat_id = chat_id, photo = img_file2)

# Plot 3: Median MAPE
p3 <- ggplot(db_time, aes(x = N, y = mapemedian)) +
  geom_point(aes(colour = method)) +
  geom_line(aes(colour = method)) +
  facet_grid(beta ~ ntime) +
  labs(
    y = 'Median APE',
    colour = "Method",
    x = 'Samples number (ID)',
    subtitle = 'Number of clusters'
  ) +
  theme(plot.subtitle = element_text(hjust = 0.5))

img_file3 <- file.path(img_path, "3_median_mape.png")
ggsave(filename = img_file3, plot = p3, device = "png", dpi = 600)
bot$sendPhoto(chat_id = chat_id, photo = img_file3)

# Plot 4: IQR of MAPE
p4 <- ggplot(db_time, aes(x = N, y = mapeiqr)) +
  geom_point(aes(colour = method)) +
  geom_line(aes(colour = method)) +
  facet_grid(beta ~ ntime) +
  labs(
    y = 'APE IQR',
    colour = "Method",
    x = 'Samples number (ID)',
    subtitle = 'Number of clusters'
  ) +
  theme(plot.subtitle = element_text(hjust = 0.5))

img_file4 <- file.path(img_path, "4_iqr_mape.png")
ggsave(filename = img_file4, plot = p4, device = "png", dpi = 600)
bot$sendPhoto(chat_id = chat_id, photo = img_file4)

# Plot 5: Type II (false negative) error
p5 <- ggplot(db_time, aes(x = N, y = sigtrue)) +
  geom_point(aes(colour = method)) +
  geom_line(aes(colour = method)) +
  facet_grid(beta ~ ntime) +
  labs(
    title = 'Type II (false negative) error rate',
    y = '% of p-values > 0.05',
    colour = "Method",
    x = 'Samples number (ID)',
    subtitle = 'Number of clusters'
  ) +
  theme(plot.subtitle = element_text(hjust = 0.5))

img_file5 <- file.path(img_path, "5_type2err.png")
ggsave(filename = img_file5, plot = p5, device = "png", dpi = 600)
bot$sendPhoto(chat_id = chat_id, photo = img_file5)

##############################################################################
# 4. Grafici per la variabile "trt" (db_xe)
##############################################################################

db_xe <- db |>
  filter(grepl("trt", rowname)) |>
  group_by(method, ntime, N, beta) |>
  summarise(
    mapemean = mean(mape, na.rm = TRUE),
    mapesd   = sd(mape, na.rm = TRUE)
  ) |>
  mutate(mapemean = ifelse(mapemean > 10, 10, mapemean))

# Plot 6: MAPE per "trt"
p6 <- ggplot(db_xe, aes(x = N, y = mapemean)) +
  geom_point(aes(colour = method)) +
  geom_line(aes(colour = method)) +
  facet_grid(beta ~ ntime)

img_file6 <- file.path(img_path, "6_mean_mape_trt.png")
ggsave(filename = img_file6, plot = p6, device = "png", dpi = 600)
bot$sendPhoto(chat_id = chat_id, photo = img_file6)

##############################################################################
# 5. Grafici "SHRINK"
##############################################################################

db_shrink_1 <- bind_rows(final_res_1000, .id = 'id') |>
  rownames_to_column() |>
  mutate(mape = 100*abs((est - beta)/est)) |>
  mutate(method = factor(method, levels = c('glmer', 'gee', 'mmm'))) |>
  mutate(mape = case_when(
    method == 'glmer' ~ 100*abs((est*inv_shrink - beta)/(est*inv_shrink)),
    method == 'gee'   ~ 100*abs((est - beta)/est),
    method == 'mmm'   ~ 100*abs((est - beta)/est)
  )) |>
  mutate(beta = ifelse(beta < 0.2, "Beta = 0.18", "Beta = 0.40")) |>
  mutate(method = case_when(
    method == 'glmer' ~ 'GLMM shrink',
    method == 'gee'   ~ 'GEE',
    method == 'mmm'   ~ 'MMM'
  )) |>
  filter(method == 'GLMM shrink')

db_shrink <- rbind(db, db_shrink_1)

db_time_shrink <- db_shrink |>
  filter(grepl("^time", rowname)) |>
  group_by(method, ntime, N, beta) |>
  summarise(
    mapemean   = mean(mape, na.rm = TRUE),
    mapesd     = sd(mape, na.rm = TRUE),
    mapeiqr    = IQR(mape, na.rm = TRUE),
    sigtrue    = 100 - (sum(sig, na.rm = TRUE)/length(sig))*100,
    mapemedian = median(mape, na.rm = TRUE)
  ) |>
  as.data.frame()

# Plot 7: Mean MAPE con "SHRINK"
p7 <- ggplot(db_time_shrink, aes(x = N, y = mapemean)) +
  geom_point(aes(colour = method)) +
  geom_line(aes(colour = method)) +
  facet_grid(beta ~ ntime) +
  labs(
    y = 'MAPE',
    colour = "Method",
    x = 'Samples number (ID)',
    subtitle = 'Number of clusters'
  ) +
  theme(plot.subtitle = element_text(hjust = 0.5))

img_file7 <- file.path(img_path, "7_mean_mape_shrink.png")
ggsave(filename = img_file7, plot = p7, device = "png", dpi = 600)
bot$sendPhoto(chat_id = chat_id, photo = img_file7)

# Plot 8: Median MAPE con "SHRINK"
p8 <- ggplot(db_time_shrink, aes(x = N, y = mapemedian)) +
  geom_point(aes(colour = method)) +
  geom_line(aes(colour = method)) +
  facet_grid(beta ~ ntime) +
  labs(
    y = 'Median APE',
    colour = "Method",
    x = 'Samples number (ID)',
    subtitle = 'Number of clusters'
  ) +
  theme(plot.subtitle = element_text(hjust = 0.5))

img_file8 <- file.path(img_path, "8_median_mape_shrink.png")
ggsave(filename = img_file8, plot = p8, device = "png", dpi = 600)
bot$sendPhoto(chat_id = chat_id, photo = img_file8)

##############################################################################
# 6. Grafici Type I error
##############################################################################

db_type1 <- bind_rows(final_res_1000, .id = 'id') |>
  rownames_to_column() |>
  mutate(mape = 100*abs((est - beta)/est)) |>
  mutate(beta = ifelse(beta < 0.2, "0.18", "0.40")) |>
  mutate(across(c(ntime, id, beta), as.factor)) |>
  mutate(method = case_match(method,
                             'glmer' ~ 'GLMM',
                             'gee'   ~ 'GEE',
                             'mmm'   ~ 'MMM'
  ),
  method = factor(method, levels = c('GLMM', 'GEE', 'MMM')))

db_time_1 <- db_type1 |>
  filter(grepl("^time", rowname)) |>
  group_by(method, ntime, N, beta) |>
  summarise(
    mapemean = mean(mape, na.rm = TRUE),
    mapesd   = sd(mape, na.rm = TRUE),
    mapeiqr  = IQR(mape, na.rm = TRUE),
    sigtrue  = (sum(sig, na.rm = TRUE)/length(sig))*100
  ) |>
  as.data.frame()

# Plot 9: Type I (false positive) error rate
p9 <- ggplot(db_time_1, aes(x = N, y = sigtrue)) +
  geom_point(aes(colour = method)) +
  geom_line(aes(colour = method)) +
  facet_grid(~ ntime)  +
  labs(
    title = 'Type I (false positive) error rate',
    y = '% of p-values < 0.05',
    colour = "Method",
    x = 'Samples number (ID)',
    subtitle = 'Number of clusters'
  ) +
  theme(plot.subtitle = element_text(hjust = 0.5))

img_file9 <- file.path(img_path, "6_type1_mape.png")
ggsave(filename = img_file9, plot = p9, device = "png", dpi = 600)
bot$sendPhoto(chat_id = chat_id, photo = img_file9)

# Grafico Type I con 'glmer_shrink'
db_shrink_2 <- bind_rows(final_res_1000_type1, .id = 'id') |>
  rownames_to_column() |>
  mutate(mape = abs((est - beta)/est)) |>
  mutate(method = factor(method, levels = c('glmer', 'gee', 'mmm'))) |>
  mutate(mape = case_when(
    method == 'glmer' ~ abs((est*inv_shrink - beta)/(est*inv_shrink)),
    method == 'gee'   ~ abs((est - beta)/est),
    method == 'mmm'   ~ abs((est - beta)/est)
  )) |>
  mutate(beta = ifelse(beta < 0.2, "0.18", "0.40")) |>
  mutate(method = case_when(
    method == 'glmer' ~ 'glmer_shrink',
    method == 'gee'   ~ 'gee',
    method == 'mmm'   ~ 'mmm'
  )) |>
  filter(method == 'glmer_shrink')

db_shrink_type1 <- rbind(db_type1, db_shrink_2) |>
  filter(grepl("^time", rowname)) |>
  group_by(method, ntime, N, beta) |>
  summarise(
    mapemean = mean(mape, na.rm = TRUE),
    mapesd   = sd(mape, na.rm = TRUE),
    mapeiqr  = IQR(mape, na.rm = TRUE),
    sigtrue  = (sum(sig, na.rm = TRUE)/length(sig))*100
  ) |>
  as.data.frame()

# Plot 10: Type I (false positive) error rate con shrink
p10 <- ggplot(db_shrink_type1, aes(x = N, y = sigtrue)) +
  geom_point(aes(colour = method), position = position_dodge(width = 0.5)) +
  geom_line(aes(colour = method), position = position_dodge(width = 0.5)) +
  facet_grid(~ ntime) +
  labs(
    title = 'Type I (false positive) error rate',
    y = '% of p-values < 0.05',
    caption = 'glmer e glmer_shrink hanno gli stessi valori',
    colour = "Methods",
    x = 'Samples number (ID)',
    subtitle = 'Number of clusters'
  ) +
  theme(plot.subtitle = element_text(hjust = 0.5))

img_file10 <- file.path(img_path, "type1err.png")
ggsave(filename = img_file10, plot = p10, device = "png", dpi = 600)
bot$sendPhoto(chat_id = chat_id, photo = img_file10)

##############################################################################
# FINE
##############################################################################

