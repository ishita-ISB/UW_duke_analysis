library(tidyverse)
library(ComplexHeatmap)
library(jsonlite)
library(colorRamp2)
library(magick)

# Load model files
overExpressedMembersMatrix <- read_csv("/users/imukherj/omics4tb2/ishitaM/GBM_R01_stateMap/data/run2/miner_network_results_1/overExpressedMembers.csv") |>
  column_to_rownames("...1")

underExpressedMembersMatrix <- read_csv("/users/imukherj/omics4tb2/ishitaM/GBM_R01_stateMap/data/run2/miner_network_results_1/underExpressedMembers.csv") |>
  column_to_rownames("...1")
# get the matrix file
dfr = (overExpressedMembersMatrix-underExpressedMembersMatrix)


## Process Programs for annotating heatmap rows
## We only show disease relevant programs in the heatmap
###disease_programs <- read_csv("../GBM-Model-052022/CoxProportionalHazardsPrograms.csv") |>
  ###filter(`p-value` <= 0.05) |>
  ###rename("Program_ID" = "...1")

# Load programs ans states from model
programs_json = read_json(simplifyVector = T,"/users/imukherj/omics4tb2/ishitaM/GBM_R01_stateMap/data/run2/miner_network_results_1/transcriptional_programs.json")
states_json = read_json(simplifyVector = T, "/users/imukherj/omics4tb2/ishitaM/GBM_R01_stateMap/data/run2/miner_network_results_1/transcriptional_states.json")

# get ordering of the programs
df2 <- arrange(dfr, factor(rownames(dfr), levels=as.vector(unlist(programs_json))))

# get ordering of states
df3 <- df2 |>
  select(as.vector(unlist(states_json)))

# unlist prgrams dictionary
prog_df <- data.frame()
for(i in 1:length(programs_json)){
  tmp1 <- data.frame(programs_json[i]) |>
    rename("regulon" = 1) |>
    mutate(program = names(programs_json[i]))
  prog_df <- bind_rows(prog_df, tmp1)
}


# unlist states dictionary
state_df <- data.frame()
for(i in 1:length(states_json)){
  tmp1 <- data.frame(states_json[i]) |>
    rename("regulon" = 1) |>
    mutate(state = names(states_json[i]))
  state_df <- bind_rows(state_df, tmp1)
}



prog_df <- prog_df |>
  column_to_rownames("regulon")

prog_df_upT <- prog_df %>% mutate(idT = as.numeric(prog_df$program) + 1) %>% select(idT)
colnames(prog_df_upT) <- c("program")
prog_df <- prog_df_upT
# filter for disease specific prgrams
prog_df_disease <- prog_df 

# get row names
prog_df_rows <- prog_df |>
  rownames_to_column("regulon")

# get rownames for disease programs
prog_df_disease_rows <- prog_df |>
  rownames_to_column("regulon")



# reformat for annotations
prog_df_disease_final <- left_join(prog_df_rows,prog_df_disease_rows, by="regulon") |>
  select(regulon, program.y) |>
  rename("program" = "program.y") |>
  column_to_rownames("regulon") |>
  mutate(program = as.character(program)) |>
  replace_na(list(program = "")) |>
  mutate(color = if_else(program != "", "#dd99ff", "#ffffff"))

# Final formatting fr annotations
prog_df_annot <- prog_df_disease_final |>
  rownames_to_column("regulon") |>
  rownames_to_column("row_id") |>
  group_by(program) |>
  mutate(median_col = round(median(as.numeric(row_id)),digits = 0)) |>
  mutate(median_col = if_else(program == "", "", paste0(median_col))) |>
  filter(program != "") |>
  select(program,median_col) |>
  unique()

### State annotations
state_df <- state_df |>
  column_to_rownames("regulon")

state_df_upT <- state_df %>% mutate(idT = as.numeric(state_df$state) + 1) %>% select(idT)
colnames(state_df_upT) <- c("state")
state_df$state <- as.character(state_df_upT$state)
#state_df <- state_df_upT

# Final formatting fr annotations
state_df_annot <- state_df |>
  rownames_to_column("regulon") |>
  rownames_to_column("row_id") |>
  select(state)

state_list = list()
for(i in unique(state_df_annot$state)){
  state_list[i] = list(which(state_df_annot$state == i))
}

state_starts = vector()
for(i in unique(state_df_annot$state)){
  state_start = which(state_df_annot$state == i)[1]
  state_starts <- append(state_starts, state_start)
}



# color functions for bwr
col_fun = colorRamp2(c(-1, 0, 1), c("#2166ac", "#f7f7f7", "#b2182b"))

program_annotations <- 
  anno_mark(
    link_gp = gpar(fontsize = 2),
    extend = F,
    labels_gp = gpar(fontsize = 15),
    padding = 0.1,
    side = "left",
    which = "row",
    at = as.numeric(unlist(prog_df_annot$median_col)),
    labels = as.character(prog_df_annot$program))
  
# funtion for drawing state annotations
panel_fun = function(index, names) {
  grid.rect(gp = gpar(fill = "orange", col = "white"))
  grid.text(paste0(names), 0.5, 0.5, gp=gpar(fontsize = 14))
}

state_annotations <- HeatmapAnnotation(name = "States",  
  State = anno_block(
    align_to = state_list,
    panel_fun = panel_fun,
    gp = gpar(fill = "orange"),
    which = "column",
   ),
  gap = unit(6, "mm"),
  which="column")


# create patient count anotations
# change colun names for patients
colnames(df3) <- seq(1,  length(colnames(df3)))
ticks <- seq(0, length(colnames(df3)), 350)
ticks_df <- as.data.frame(ticks)# data_frame(patient_no = colnames(df3)) |>
  #mutate(label = if_else(patient_no %in% ticks, patient_no, ""))


patient_annotations <- HeatmapAnnotation(
  Patients = anno_mark(
    at = as.numeric(ticks_df$ticks),
    labels = as.character(ticks_df$ticks),
    side = "bottom",
    
  #  gp = gpar(fill = "#333333"),
    which = "column",
  ),
  which="column")



ht1 <- Heatmap(as.matrix(df3),#left_annotation=program_annotations,
        use_raster = T,
        name="Activity",
       #column_labels = ticks_df$label,
        cluster_row_slices = F,cluster_column_slices = F,
        row_title = "Programs",
        column_title = "Cells",
        row_title_side = "left",
        column_title_side = "bottom",
        row_dend_reorder = F,
        show_row_dend = F,
        raster_resize_mat = T,
        raster_magick_filter = F,
        show_column_dend = F,
        show_column_names = F,
        show_row_names = F,
        cluster_columns = F,
        cluster_rows = F,
        column_split = state_df_annot,
        column_gap = unit(1, "mm"),
       column_names_gp = grid::gpar(fontsize = 18),
       row_names_gp = grid::gpar(fontsize = 18),
       column_title_gp = gpar(fontsize = 20, fontface = "bold"),
       row_title_gp = gpar(fontsize = 20, fontface = "bold"),
        #split = prog_df_annot,
        #split=prog_df_disease_final,cluster_row_slices = F,
        col=colorRampPalette(c("blue", "white", "red"))(100),
        left_annotation = rowAnnotation(mark=program_annotations),
        top_annotation = state_annotations, #columnAnnotation(df=state_df_annot,  gp = gpar(color = "yellow"))
        bottom_annotation = patient_annotations
        )

pdf(file = "GBM_Model_Heatmapv2.pdf", width = 14, height = 10,)
ht1
dev.off()





