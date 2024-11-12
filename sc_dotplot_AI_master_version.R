input_gene <- c("CXCR4", "LYN", "C3", "P2RX4", "NCF4", "FMOD", "DENND3")


# 在R中设置环境变量
Sys.setenv(R_LIBS_USER = "~/R/x86_64-pc-linux-gnu-library/4.4")
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
# 建议设置为总核心数的75%左右，保留一些核心给系统和其他进程
Sys.setenv(OPENBLAS_NUM_THREADS = "60")

# 加载parallel包
library(parallel)
# 设置全局默认核心数（留一个核心给系统）
options(mc.cores = detectCores() - 4)



library(data.table)
library(stringr)
library(dplyr)
library(ggplot2)
library(viridis)
library(tidyverse)
library(stringr)

input_location <- file.path("data", "Ref_Gene_location.txt")

input_data <- file.path("data", "merged_ref_data.csv")

proper_cell_names <- c(
  "B cell", "capillary endothelial cell", "conventional dendritic cell", "endothelial cell of artery", "endothelial cell of lymphatic vessel",
  "endothelial cell of venule", "fibroblast", "gamma-delta T cell", "helper T cell", "keratinocyte", "Langerhans cell", "macrophage", "mast cell",
  "mature NK T cell", "melanocyte", "mucosal invariant T cell", "natural killer cell", "neutrophil", "pericyte", "plasma cell",
  "plasmacytoid dendritic cell", "regulatory T cell", "Schwann cell", "T-helper 17 cell", "T cell", "vascular associated smooth muscle cell"
)

# proper_cell_names主要是应位有些名字比如CD4+会报错，其二是这个变量决定了最后绘制图的细胞类型
# 定义sc_dotplot函数
sc_dotplot <- function(Group, genelist) {
  Gene_location <- read.csv(input_location,
    header = FALSE
  )

  templist <- c("cell_type", "disease", genelist)
  templist <- tolower(templist)
  first_row <- tolower(Gene_location[1, ])
  match_index <- match(templist, first_row)

  # 使用fread函数来加载csv文件中的特定列
  data <- fread(input_data, select = match_index)
  genelist <- colnames(data)[c(-1, -2)]
  data <- as.data.frame(data)
  data <- data %>% filter(disease == Group)
  data_logical <- data
  data_logical[, genelist] <- lapply(data[, genelist], function(x) {
    x != 0
  })
  data_filtered <- data[rowSums(data[, genelist] != 0) > 0, ]

  # 初始化图
  plot_data <- data.frame()

  for (gene in genelist) {
    percent_exp <- data_logical %>%
      group_by(cell_type) %>%
      summarise(pct = mean(!!sym(gene)))

    avg_exp <- data_filtered %>%
      group_by(cell_type) %>%
      summarise(avg = mean(!!sym(gene)))

    plot_gene <- left_join(percent_exp, avg_exp)
    plot_gene$gene <- gene

    plot_data <- rbind(plot_data, plot_gene)
  }

  # 修改pct列值
  plot_data$pct <- plot_data$pct * 100

  # 保证基因名顺序与输入一致
  plot_data$gene <- factor(plot_data$gene, levels = input_gene)

  return(plot_data)
}

# 调用sc_dotplot函数生成p1和p2
p1 <- sc_dotplot("normal", input_gene)
p2 <- sc_dotplot("periodontitis", input_gene)

# 定义pvalue_counts函数
pvalue_counts <- function(genelist) {
  Gene_location <- read.csv(input_location,
    header = FALSE
  )

  genelist <- c("cell_type", "disease", genelist)
  genelist <- tolower(genelist)
  first_row <- tolower(Gene_location[1, ])
  match_index <- match(genelist, first_row)

  # 使用fread函数来加载csv文件中的特定列
  raw <- fread(input_data, select = match_index)

  # 选择分析的细胞类型
  SR <- raw$cell_type
  table(raw$cell_type)

  # 细胞名称
  cells <- proper_cell_names

  # 循环每个细胞类型
  purrr::map(cells, function(x) {
    # 筛选该细胞类型
    SR <- ifelse(str_detect(raw$cell_type, x), TRUE, FALSE)
    selcted_cell <- raw[SR, ]

    # 设置参考水平
    group_list <- ifelse(str_detect(selcted_cell$disease, "periodontitis"), "periodontitis", "normal")
    group_list <- factor(group_list, levels = c("normal", "periodontitis"))

    exp <- as.data.frame(t(selcted_cell))
    colnames(exp) <- exp[2, ]
    exp <- exp[c(-1, -2), ]

    for (i in 1:ncol(exp)) {
      if (is.character(exp[, i])) {
        exp[, i] <- as.numeric(exp[, i])
      }
    }

    genes_of_interest <- rownames(exp)
    group_info <- group_list
    genes <- genes_of_interest

    # 创建绘图数据框
    plot_data <- data.frame(
      Gene = rep(genes_of_interest, each = length(group_info)),
      Group = factor(rep(group_info, length(genes_of_interest))),
      Expression = as.vector(t(exp[genes_of_interest, ]))
    )

    results <- data.frame(
      Gene = character(),
      pvalue = numeric(),
      stars = character(),
      stringsAsFactors = FALSE
    )

    for (gene in genes) {
      gene_data <- subset(plot_data, Gene == gene)
      wilcox_result <- wilcox.test(Expression ~ Group, data = gene_data)
      pvalue <- wilcox_result$p.value

      stars <- ifelse(pvalue < 0.001, "3", ifelse(pvalue < 0.01, "2", ifelse(pvalue < 0.05, "1", "0")))

      newrow <- data.frame(Gene = gene, pvalue = pvalue, stars = stars)
      results <- rbind(results, newrow)
    }

    unique(gene_data$Group)
    return(results)
  })
}

res <- pvalue_counts(input_gene)

names(res) <- proper_cell_names

# 将list中的data.frame合并为一个大的data.frame
res <- do.call(rbind, Map(cbind, dataframe = names(res), res))

rownames(res) <- NULL
colnames(res)[1] <- "cell_type"
colnames(res)[2] <- "gene"

# 合并两个表格
p1 <- left_join(res, p1, by = c("cell_type", "gene"))
p2 <- left_join(res, p2, by = c("cell_type", "gene"))






generate_plots <- function(p1, p2, threshold, title, output_file) {
  # 将NA值设为0，但不显示它们
  p1_filtered <- p1[!(p1$pct == 0 | is.na(p1$pct)), ]
  p2_filtered <- p2[!(p2$pct == 0 | is.na(p2$pct)), ]

  # 根据阈值筛选数据
  if (threshold == 0) {
    p1_filtered <- p1_filtered[p1_filtered$stars != threshold, ]
    p2_filtered <- p2_filtered[p2_filtered$stars != threshold, ]
  } else {
    p1_filtered <- p1_filtered[p1_filtered$stars == threshold, ]
    p2_filtered <- p2_filtered[p2_filtered$stars == threshold, ]
  }

  p <- ggplot() +
    geom_point(
      data = p1_filtered,
      aes(
        x = factor(gene, levels = input_gene),
        y = cell_type,
        size = pct,
        color = avg
      ),
      position = position_nudge(x = -0.15)
    ) +
    geom_point(
      data = p2_filtered,
      aes(
        x = factor(gene, levels = input_gene),
        y = cell_type,
        size = pct,
        color = avg
      ),
      position = position_nudge(x = 0.15)
    ) +
    scale_color_viridis_c(
      option = "plasma",
      limits = c(0, 3.5),
      alpha = 0.6
    ) +
    labs(title = title) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      plot.title = element_text(hjust = 0.5),
      panel.grid.major = element_line(color = "grey90", linetype = "dotted"),
      panel.grid.minor = element_line(color = "grey90", linetype = "dotted")
    ) +
    guides(
      size = guide_legend(
        order = 2,
        title = "Percent Expressed",
        title.theme = element_text(size = 10)
      ),
      color = guide_colourbar(
        order = 1,
        title = "Average Expression",
        title.theme = element_text(size = 10)
      )
    ) +
    scale_size(
      range = c(0, 12),
      breaks = seq(10, 80, by = 10),
      limits = c(0, 100)
    )

  ggsave(output_file, plot = p, width = 12, height = 8)
}




# 生成不同阈值的图形
generate_plots(p1, p2, threshold = 0, title = "normal vs periodontitis(P<0.05)", output_file = "2024_10_29_all.PDf")
generate_plots(p1, p2, threshold = 3, title = "normal vs periodontitis(P<0.001)", output_file = "2024_10_29_3star.PDf")
generate_plots(p1, p2, threshold = 2, title = "normal vs periodontitis(0.001<P<0.01)", output_file = "2024_10_29_2star.PDf")
generate_plots(p1, p2, threshold = 1, title = "normal vs periodontitis(0.01<P<0.05)", output_file = "2024_10_29_1star.PDf")







generate_plots_all <- function(p1, p2, title, output_file) {
  # 只过滤掉NA值和0值，不进行p值筛选
  p1_filtered <- p1[!(p1$pct == 0 | is.na(p1$pct)), ]
  p2_filtered <- p2[!(p2$pct == 0 | is.na(p2$pct)), ]

  p <- ggplot() +
    geom_point(
      data = p1_filtered,
      aes(
        x = factor(gene, levels = input_gene),
        y = cell_type,
        size = pct,
        color = avg
      ),
      position = position_nudge(x = -0.15)
    ) +
    geom_point(
      data = p2_filtered,
      aes(
        x = factor(gene, levels = input_gene),
        y = cell_type,
        size = pct,
        color = avg
      ),
      position = position_nudge(x = 0.15)
    ) +
    scale_color_viridis_c(
      option = "plasma",
      limits = c(0, 3.5),
      alpha = 0.6
    ) +
    labs(title = title) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      plot.title = element_text(hjust = 0.5),
      panel.grid.major = element_line(color = "grey90", linetype = "dotted"),
      panel.grid.minor = element_line(color = "grey90", linetype = "dotted")
    ) +
    guides(
      size = guide_legend(
        order = 2,
        title = "Percent Expressed",
        title.theme = element_text(size = 10)
      ),
      color = guide_colourbar(
        order = 1,
        title = "Average Expression",
        title.theme = element_text(size = 10)
      )
    ) +
    scale_size(
      range = c(0, 12),
      breaks = seq(10, 80, by = 10),
      limits = c(0, 100)
    )

  ggsave(output_file, plot = p, width = 12, height = 8)
}


generate_plots_all(p1, p2,
  title = "normal vs periodontitis",
  output_file = "all-nofilter.PDF"
)








generate_heatmap <- function(p1, p2, title, output_file) {
  # 创建差异数据框
  diff_data <- p1 %>%
    select(cell_type, gene, avg, stars) %>%
    rename(avg_normal = avg) %>%
    left_join(
      p2 %>% select(cell_type, gene, avg) %>% rename(avg_perio = avg),
      by = c("cell_type", "gene")
    ) %>%
    mutate(
      # 计算log2 fold change
      # 添加小数以避免除0错误
      log2fc = log2((avg_perio + 0.01) / (avg_normal + 0.01)),
      # 创建显著性标记
      sig_stars = case_when(
        stars == "3" ~ "***",
        stars == "2" ~ "**",
        stars == "1" ~ "*",
        TRUE ~ ""
      )
    )

  # 计算log2fc的范围，用于对称的色标范围
  max_abs_fc <- max(abs(diff_data$log2fc), na.rm = TRUE)

  # 创建热图
  p <- ggplot(diff_data, aes(
    x = factor(gene, levels = input_gene),
    y = cell_type,
    fill = log2fc
  )) +
    geom_tile() +
    geom_text(aes(label = sig_stars),
      color = "black",
      size = 3
    ) +
    scale_fill_gradient2(
      low = "blue",
      mid = "white",
      high = "red",
      midpoint = 0,
      limits = c(-max_abs_fc, max_abs_fc),
      name = "Log2 Fold Change\n(Perio/Normal)"
    ) +
    labs(
      title = title,
      x = "Gene",
      y = "Cell Type"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      plot.title = element_text(hjust = 0.5),
      panel.grid = element_blank()
    )

  # 保存图片
  ggsave(output_file, plot = p, width = 12, height = 8)

  # 返回差异数据用于检查
  return(diff_data)
}



# 生成热图
result_data <- generate_heatmap(p1, p2,
  title = "Gene Expression Changes in Periodontitis vs Normal",
  output_file = "log2fc_heatmap.PDF"
)










generate_heatmap_filtered <- function(p1, p2, threshold, title, output_file) {
  # 创建差异数据框
  diff_data <- p1 %>%
    select(cell_type, gene, avg, stars) %>%
    rename(avg_normal = avg) %>%
    left_join(
      p2 %>% select(cell_type, gene, avg) %>% rename(avg_perio = avg),
      by = c("cell_type", "gene")
    ) %>%
    mutate(
      # 计算log2 fold change
      log2fc = log2((avg_perio + 0.01) / (avg_normal + 0.01)),
      # 创建显著性标记
      sig_stars = case_when(
        stars == "3" ~ "***",
        stars == "2" ~ "**",
        stars == "1" ~ "*",
        TRUE ~ ""
      )
    )

  # 根据阈值筛选数据
  if (threshold == 0) {
    # 显示所有显著的结果 (p < 0.05)
    diff_data <- diff_data %>% filter(stars != "0")
  } else {
    # 显示特定显著性水平的结果
    diff_data <- diff_data %>% filter(stars == as.character(threshold))
  }

  # 如果没有数据符合条件，返回NULL
  if (nrow(diff_data) == 0) {
    warning("No data meets the significance threshold")
    return(NULL)
  }

  # 计算log2fc的范围，用于对称的色标范围
  max_abs_fc <- max(abs(diff_data$log2fc), na.rm = TRUE)

  # 创建热图
  p <- ggplot(diff_data, aes(
    x = factor(gene, levels = input_gene),
    y = cell_type,
    fill = log2fc
  )) +
    geom_tile() +
    geom_text(aes(label = sig_stars),
      color = "black",
      size = 3
    ) +
    scale_fill_gradient2(
      low = "blue",
      mid = "white",
      high = "red",
      midpoint = 0,
      limits = c(-max_abs_fc, max_abs_fc),
      name = "Log2 Fold Change\n(Perio/Normal)"
    ) +
    labs(
      title = title,
      x = "Gene",
      y = "Cell Type"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      plot.title = element_text(hjust = 0.5),
      panel.grid = element_blank()
    )

  # 保存图片
  ggsave(output_file, plot = p, width = 12, height = 8)

  return(diff_data)
}

# 生成四张不同显著性水平的热图
result_all <- generate_heatmap_filtered(
  p1, p2,
  threshold = 0,
  title = "Expression Changes (P<0.05)",
  output_file = "heatmap_all_sig.PDF"
)

result_3star <- generate_heatmap_filtered(
  p1, p2,
  threshold = 3,
  title = "Expression Changes (P<0.001)",
  output_file = "heatmap_3star.PDF"
)

result_2star <- generate_heatmap_filtered(
  p1, p2,
  threshold = 2,
  title = "Expression Changes (0.001<P<0.01)",
  output_file = "heatmap_2star.PDF"
)

result_1star <- generate_heatmap_filtered(
  p1, p2,
  threshold = 1,
  title = "Expression Changes (0.01<P<0.05)",
  output_file = "heatmap_1star.PDF"
)

