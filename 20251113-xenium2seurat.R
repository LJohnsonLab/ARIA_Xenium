library(Seurat)
library(tidyverse)
library(duckplyr)
library(qs2)


# Load Xenium data
# Raw data path in a T9 Samsung external drive
xenium_dir <- "/Volumes/T9/Xenium_ARIA/"

##################################
## Slide1_0021991 ###############
################################
xenium_21991 <- LoadXenium(paste0(xenium_dir,"Slide1_output-XETG00118__0021991__Region_1__20250606__202953"))

# Read the transcripts parquet file, filter out UNASSIGNED cell IDs,
# and calculate the number of transcripts and median QV per cell.
transcripts <- read_parquet_duckdb("/Volumes/T9/Xenium_ARIA/Slide1_output-XETG00118__0021991__Region_1__20250606__202953/transcripts.parquet") |> 
  # filter(cell_id!="UNASSIGNED") |> 
  summarize(nCount_raw = n(),
            nFeature_raw = n_distinct(feature_name),
            median_qv = median(qv, na.rm=TRUE), 
            .by = cell_id) |> 
  collect() # Finally, collect the results into a data frame.



# Read sample ROIs from CSV
# The CSV file should contain the polygon coordinates for each ROI
sample1 <- list(
  "KK4_465" = read_csv("spatial_data/KK4_465_cells_stats.csv",comment = "#" ),
  "KK4_492" = read_csv("spatial_data/KK4_492_cells_stats.csv",comment = "#" ),
  "KK4_504" = read_csv("spatial_data/KK4_504_cells_stats.csv",comment = "#" )) |> 
  bind_rows(.id = "sample_id") |> 
  select(sample_id, cell_id = "Cell ID", cell_area = "Area (µm^2)" )

# include sample_id and median QV in the seourat object metadata
meta_data1 <- xenium_21991@meta.data |> 
  rownames_to_column("cell_id") |>
  left_join(sample1, by = "cell_id") |>
  left_join(transcripts, by = "cell_id") |>
  column_to_rownames("cell_id")|> 
    mutate(slide = "Slide1_0021991")

# Update the Seurat object metadata
xenium_21991@meta.data <- meta_data1

# Filter out cells that:
# 1) don't fall within any ROI
# 2)have zero counts
# 3) median QV value < 20
xenium_21991_filtered <- subset(xenium_21991, 
  subset = !is.na(sample_id) & nCount_Xenium > 0 & median_qv >= 20)


  ##################################
  ## Slide2_0021998 ###############
  ################################
  xenium_0021998 <- LoadXenium(paste0(xenium_dir,"Slide2_output-XETG00118__0021998__Region_1__20250606__202953"))
  
  # Read the transcripts parquet file, filter out UNASSIGNED cell IDs,
  # and calculate the number of transcripts and median QV per cell.
  transcripts2 <- read_parquet_duckdb("/Volumes/T9/Xenium_ARIA/Slide2_output-XETG00118__0021998__Region_1__20250606__202953/transcripts.parquet") |> 
    # filter(cell_id!="UNASSIGNED") |> 
    summarize(nCount_raw = n(),
              nFeature_raw = n_distinct(feature_name),
              median_qv = median(qv, na.rm=TRUE), 
              .by = cell_id) |> 
    collect() # Finally, collect the results into a data frame.
  
  
  
  # Read sample ROIs from CSV
  # The CSV file should contain the polygon coordinates for each ROI
  sample2 <- list(
    "KK4_496" = read_csv("spatial_data/KK4_496_cells_stats.csv",comment = "#" ),
    "KK4_502" = read_csv("spatial_data/KK4_502_cells_stats.csv",comment = "#" ),
    "KK4_464" = read_csv("spatial_data/KK4_464_cells_stats.csv",comment = "#" )) |> 
    bind_rows(.id = "sample_id") |> 
    select(sample_id, cell_id = "Cell ID", cell_area = "Area (µm^2)" )
  
  # include sample_id and median QV in the seourat object metadata
  meta_data2 <- xenium_0021998@meta.data |> 
    rownames_to_column("cell_id") |>
    left_join(sample2, by = "cell_id") |>
    left_join(transcripts2, by = "cell_id") |>
    column_to_rownames("cell_id") |> 
    mutate(slide = "Slide2_0021998")
  
  # Update the Seurat object metadata
  xenium_0021998@meta.data <- meta_data2
  
  # Filter out cells that:
  # 1) don't fall within any ROI
  # 2)have zero counts
  # 3) median QV value < 20
  xenium_0021998_filtered <- subset(xenium_0021998, 
    subset = !is.na(sample_id) & nCount_Xenium > 0 & median_qv >= 20)
  
###############################################
# Merge the two filtered Seurat objects #######
###############################################
aria <- merge(xenium_21991_filtered, xenium_0021998_filtered,
               add.cell.ids = c("slide1", "slide2"),
              project ="ARIA")

aria$treatment <- ifelse(aria$sample_id %in% c("KK4_465","KK4_504","KK4_496"), "IgG", "Adu")

aria <- JoinLayers(aria, assay = "Xenium")

# Save the filtered Seurat object
## Use the qs2 package for faster saving/loading. It is much faster than RDS.

# qs_save(aria, file = "data/20251113_aria.qs2")
