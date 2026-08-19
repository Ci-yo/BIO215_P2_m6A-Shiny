# ---- app.R (CLEAN, final) ----
options(shiny.usecairo = TRUE)  # 让文字更清晰（Windows 上尤其有效）
library(shiny)
library(shinydashboard)
library(DT)
library(ggplot2)
library(dplyr)
suppressWarnings(library(randomForest))
library(fresh)

brand  <- "#0EA5E7"   # 主色（teal）
accent <- "#2563EB"   # 次强调（蓝）

theme <- create_theme(
  adminlte_sidebar(
    width          = "240px",
    light_bg       = "#FFFFFF",
    light_hover_bg = "#F1F5F9",
    light_color    = "#0F172A"
  ),
  adminlte_global(
    content_bg  = "transparent",  # ← 这里改成透明
    box_bg      = "#FFFFFF",
    info_box_bg = "#FFFFFF"
  )
)

# 只加载现成文件（不要训练！）
m6a_data     <- read.csv("m6A_imbalanced.csv", stringsAsFactors = FALSE)
rf_fit       <- readRDS("rf_fit.rds")
example_data <- read.csv("m6A_input_example.csv", stringsAsFactors = FALSE)

# ---- Helpers ----
`%||%` <- function(a, b) if (is.null(a)) b else a
is_valid_5mer <- function(x){
  x <- toupper(trimws(x %||% ""))
  nchar(x) == 5 & !grepl("[^ATCG]", x)
}

dna_encoding <- function(dna_strings){
  nn <- nchar(dna_strings[1])
  seq_m <- matrix(unlist(strsplit(dna_strings, "")), ncol = nn, byrow = TRUE)
  colnames(seq_m) <- paste0("nt_pos", 1:nn)
  seq_df <- as.data.frame(seq_m)
  seq_df[] <- lapply(seq_df, factor, levels = c("A","T","C","G"))
  seq_df
}

prediction_multiple <- function(ml_fit, feature_df, positive_threshold = 0.5){
  stopifnot(all(c("gc_content","RNA_type","RNA_region","exon_length",
                  "distance_to_junction","evolutionary_conservation","DNA_5mer")
                %in% colnames(feature_df)))
  
  # 因子与校验
  allowed_type   <- c("mRNA","lincRNA","lncRNA","pseudogene")
  allowed_region <- c("CDS","intron","3'UTR","5'UTR")
  bad_type   <- setdiff(unique(as.character(feature_df$RNA_type)),   allowed_type)
  bad_region <- setdiff(unique(as.character(feature_df$RNA_region)), allowed_region)
  if (length(bad_type))
    stop("Unknown RNA_type level(s): ", paste(bad_type, collapse = ", "))
  if (length(bad_region))
    stop("Unknown RNA_region level(s): ", paste(bad_region, collapse = ", "))
  
  feature_df$RNA_type   <- factor(feature_df$RNA_type,   levels = allowed_type)
  feature_df$RNA_region <- factor(feature_df$RNA_region, levels = allowed_region)
  
  # 5-mer 严格校验
  ok_5mer <- is_valid_5mer(feature_df$DNA_5mer)
  if (any(!ok_5mer)) {
    idx <- which(!ok_5mer)[1]
    stop("Invalid DNA_5mer at row ", idx,
         ": '", feature_df$DNA_5mer[idx], "'. Must be exactly 5 chars from A/T/C/G.")
  }
  
  enc <- dna_encoding(toupper(feature_df$DNA_5mer))
  pred_input <- cbind(
    feature_df[, c("gc_content","RNA_type","RNA_region","exon_length",
                   "distance_to_junction","evolutionary_conservation")],
    enc
  )
  
  p <- predict(ml_fit, newdata = pred_input, type = "prob")[, "Positive"]
  feature_df$predicted_m6A_prob   <- round(p, 3)
  feature_df$predicted_m6A_status <- ifelse(p > positive_threshold, "Positive", "Negative")
  feature_df
}

prediction_single <- function(ml_fit, gc_content, RNA_type, RNA_region, exon_length,
                              distance_to_junction, evolutionary_conservation, DNA_5mer,
                              positive_threshold = 0.5){
  single_df <- data.frame(
    gc_content = as.numeric(gc_content),
    RNA_type = as.character(RNA_type),
    RNA_region = as.character(RNA_region),
    exon_length = as.numeric(exon_length),
    distance_to_junction = as.numeric(distance_to_junction),
    evolutionary_conservation = as.numeric(evolutionary_conservation),
    DNA_5mer = toupper(as.character(DNA_5mer)),
    stringsAsFactors = FALSE
  )
  res <- prediction_multiple(ml_fit, single_df, positive_threshold)
  c(predicted_m6A_prob   = as.character(res$predicted_m6A_prob[1]),
    predicted_m6A_status = as.character(res$predicted_m6A_status[1]))
}

# ---- UI ----
ui <- dashboardPage(
  skin = "blue",                  # 也可 "black","green","purple","red","yellow"
  dashboardHeader(title = "Mini m6A Database & Predictor"),
  dashboardSidebar(
    collapsed = FALSE,            # TRUE = 默认收起
    sidebarMenu(
      id = "tabs",
      # 分组标题（见第3节更完整示例）
      tags$li(class = "header", "DATA & VISUALIZATION"),
      menuItem("Home", tabName="home", icon=icon("home")),
      menuItem("Data", tabName="data", icon=icon("table"),
               badgeLabel = nrow(m6a_data), badgeColor = "green"),
      menuItem("Visualize", tabName="viz", icon=icon("chart-pie")),
      menuItem("Download", tabName="down", icon=icon("download")),
      tags$li(class = "header", "PREDICTION"),
      menuItem("m6A Prediction Multiple", tabName="pred_multi", icon=icon("list")),
      menuItem("m6A Prediction Single",   tabName="pred_single", icon=icon("vial")),
      tags$li(class = "header", "INFO"),
      menuItem("Help", tabName = "about", icon = icon("circle-info"))
      
    )
  ),
  dashboardBody(
    # 轻量 CSS
    tags$head(tags$style(HTML("

@media (max-width: 992px){
  #workflow_img img{ max-height:300px; }  /* 小屏稍微收一点 */
}
/* Workflow：放大“背景板”，并等高 */
#wf_card .box-body{ padding:10px 12px; }
#wf_side .box-body{ min-height:520px; }     /* 右侧卡片等高 */

#wf_crop{
  overflow:hidden;                           /* 裁边以减少四周留白 */
  border-radius:12px;
  border:1px solid #E5E7EB;
  background:#ffffff;
  box-shadow:0 1px 2px rgba(0,0,0,.04);
  padding:6px;
}
#wf_crop img{
  width:100% !important;
  height:auto !important;
  max-height:520px;                          /* 上限（可调 500–560） */
  object-fit:contain;
  transform:scale(1.08);                     /* 轻微放大 8%，让“板”更大 */
  transform-origin:center top;
  display:block;
  margin:0 auto;
}

@media (max-width: 992px){
  #wf_side .box-body{ min-height:unset; }
  #wf_crop img{ max-height:420px; transform:scale(1.04); }
}


/* --- Single prediction top summary --- */
.pred-badges{display:flex;gap:10px;flex-wrap:wrap;align-items:center;}
.prog { width:100%; height:14px; background:#E5E7EB; border-radius:9999px; overflow:hidden; margin-top:8px; }
.prog > div { height:100%; }

/* 缩小右侧 Prediction 卡片正文的默认行高，让信息更紧凑 */
#single-details .box-body h5 { margin: 6px 0; }

/* 输入区标签更醒目 */
#single-inputs .control-label { font-weight: 600; }

/* 让顶部摘要卡片的徽章字号稍大一点 */
#single-summary .pill { font-size: 13px; padding: 4px 12px; }

/* —— Home: 流程图容器与图片 —— */
.home-figure{ display:flex; justify-content:center; }
#home_img{ width:100%; }
#home_img img{
  display:block; margin:0 auto;
  width:100%; height:auto;         /* 列多宽就多宽 */
  max-width:none;                   /* 取消上限 */
  background:#ffffff;
  border:1px solid #E5E7EB; border-radius:12px; padding:10px;
  box-shadow:0 1px 2px rgba(0,0,0,.04);
}

    /* —— 侧栏：浅色主题 + 高亮色 —— */
.skin-blue .main-sidebar, .skin-blue .left-side { background:#ffffff; }
.skin-blue .sidebar a { color:#334155; }                      /* 普通文字色 */
.skin-blue .sidebar-menu>li>a { color:#334155; }
.skin-blue .sidebar-menu>li.active>a,
.skin-blue .sidebar-menu>li:hover>a {
  background:#F1F5F9; color:#0f172a; border-left-color:#14B8A6;
}
.skin-blue .sidebar-menu>li.active>a { font-weight:600; }
.skin-blue .sidebar-menu>li.header {
  color:#64748B; background:#ffffff; font-weight:700; letter-spacing:.04em;
}

/* —— 宽度 240px（默认 230px），内容区同步左边距 —— */
.main-sidebar { width:240px !important; }
.content-wrapper, .right-side, .main-footer { margin-left:240px !important; }
@media (max-width:767px){
  .content-wrapper, .right-side, .main-footer { margin-left:0 !important; }
}

/* —— 行距 / 字号更紧凑 —— */
.sidebar-menu>li>a { padding:10px 14px; font-size:14px; }
.sidebar-menu>li { margin:2px 0; }

      .stepper {display:flex;gap:12px;align-items:center;margin:6px 0 16px;}
      .step {display:flex;align-items:center;gap:8px;}
      .step .dot {width:26px;height:26px;border-radius:50%;background:#3B82F6;color:#fff;
                  display:flex;align-items:center;justify-content:center;font-weight:700;}
      .step .txt {font-weight:600;}
      .sep {flex:1;height:2px;background:#E5E7EB;}
      .pill{display:inline-block;padding:3px 10px;border-radius:9999px;font-size:12px;font-weight:600;}
      .pill-green{background:#DCFCE7;color:#166534;}
      .pill-red{background:#FEE2E2;color:#991B1B;}
      .pill-blue{background:#DBEAFE;color:#1E40AF;}
      .pill-gray{background:#E5E7EB;color:#111827;}
      .btn-row{display:flex;gap:10px;flex-wrap:wrap;}
      .tight{margin-top:6px;margin-bottom:6px;}
      /* —— Home 视觉统一 —— */
.home-figure{ display:flex; justify-content:center; }
.home-figure img{
  display:block; margin:0 auto;
  height: 300px;        /* ← 想更长就调大：680/720/800 px 都行 */
  width: auto;          /* 保持等比；不会横向被拉伸 */
  max-width: 100%;      /* 宽度不超过列宽 */
  object-fit: contain;  /* 等比缩放，必要时左右留白 */
  background:#ffffff;
  border:1px solid #E5E7EB; border-radius:12px; padding:10px;
  box-shadow:0 1px 2px rgba(0,0,0,.04);
}

/* --- Multiple: 空态占位 --- */
.empty-state{
  display:flex; flex-direction:column; align-items:center; justify-content:center;
  gap:10px; padding:24px; height:360px;
  border:2px dashed #E5E7EB; border-radius:12px;
  background:#F8FAFC; color:#64748B; text-align:center;
}
.empty-state i{ font-size:38px; opacity:.7 }
/* 上传行：Browse + 两个按钮 同行；Browse 框固定较短 */
#upload-row{
  display:flex; align-items:center; gap:10px; flex-wrap:wrap;   /* 宽度不足时自动换行 */
}
#upload-row .upload-file{ flex:0 0 380px; max-width:380px; min-width:320px; }
#upload-row .upload-file .shiny-input-container{ margin:0; }   /* 去掉 fileInput 底部空隙 */
#upload-row .btn, #upload-row .btn-file{
  height:34px; padding:6px 12px; line-height:1.428; white-space:nowrap;
}
/* 小屏时：让 Browse 独占一行、按钮掉到下一行 */
@media (max-width:1200px){
  #upload-row .upload-file{ flex:1 1 100%; max-width:100%; }
}



/* --- Results：按钮与表格的间距 --- */
.results-actions{ margin-bottom:12px; }           /* Predict/Download 下方留白 */
#tbl_pred_wrapper{ margin-top:12px; }             /* DataTable 顶部整体下沉 */
#tbl_pred_wrapper .dataTables_length,
#tbl_pred_wrapper .dataTables_filter{ margin-top:6px; }


#multi-panel{ margin-top:8px; }


/* 小屏幕再允许换行，避免挤压 */
@media (max-width: 1024px){
  #upload-row{ flex-wrap:wrap; }
}

/* Results 顶部按钮行与内容的间距（强制生效） */
.box .results-actions{ margin-bottom:12px !important; }
#multi-panel{ margin-top:8px !important; }

/* Results 顶部按钮行与下面的内容拉开间距 */
.results-actions{ margin-bottom:12px; }

/* 统一让 Results 面板里，按钮行下方的任何内容都往下错开一点 */
#multi-panel{ margin-top:8px; }

/* 当已经预测并显示 DataTable 时：把“Show … entries / Search”再下沉一点 */
#tbl_pred_wrapper{ margin-top:10px; }
#tbl_pred_wrapper .dataTables_length,
#tbl_pred_wrapper .dataTables_filter{ margin-top:6px; }

/* ===== Home: Welcome(Hero) ===== */
.hero {
  padding: 18px 22px;
  background: #F8FAFC;
  border: 1px solid #E5E7EB;
  border-radius: 14px;
  box-shadow: 0 1px 2px rgba(0,0,0,.04);
}
.hero h1 {
  margin: 0 0 6px;
  font-weight: 800;
  letter-spacing: .2px;
}
.hero p.sub {
  color: #475569;
  margin: 0 0 12px;
  font-size: 16px;
}
.hero .btn-row a.btn {
  margin-right: 8px;
}

/* Home：右侧“使用步骤”列表更紧凑 */
.use-steps li { margin-bottom: 6px; }


/* box 头部与说明文字更干净 */
.box.solidHeader > .box-header { background:#F8FAFC; border-bottom:1px solid #E5E7EB; }
.box .help-block{ text-align:center; color:#6B7280; margin-top:8px; }

/* Required columns: wrap neatly without overlap */
#reqcols{
  display:flex;
  flex-wrap:wrap;
  gap: 8px 12px;          /* 行间/列间间距 */
  margin: 6px 0;          /* 与上下内容拉开 */
}
#reqcols .pill{           /* 只作用于这一块的 pill，避免影响其他地方 */
  line-height: 1.25;
  margin: 0;              /* 用 gap 控制间距 */
}

/* Single: center the percentage inside the progress bar */
#single-summary .progress{
  height: 18px;                 /* 比 14px 更不拥挤，可按需改 */
  margin-top: 6px;
}
#single-summary .progress .progress-bar{
  display: flex;                 /* 让内容用 flex 居中 */
  align-items: center;           /* 垂直居中 */
  justify-content: center;       /* 水平居中 */
  padding: 0;                    /* 去掉内边距，避免文字贴底 */
  line-height: 1;                /* 避免基线偏移 */
  font-weight: 700;
  font-size: 12px;               /* 可调：11–13 */
}
/* App title (top bar) — make it bolder */
.main-header .logo,
.main-header .logo .logo-lg{
  font-weight: 800 !important;   /* 700→900 更粗 */
}

/* Sidebar section headers (e.g., DATA & VISUALIZATION / PREDICTION / INFO) */
.sidebar-menu > li.header{
  font-weight: 700 !important;   /* 700→900 */
  letter-spacing: .06em;         /* 可选：略微加字距，更稳重 */
  color: #0f172a;                /* 可选：更深的文字色 */
}


/* 标题留点呼吸感 */
h2{ margin-top:6px; margin-bottom:14px; }

/* Make section titles bolder */
.content h3 { font-weight: 800; }                 /* e.g., 'Workflow' */
.box .box-header .box-title { font-weight: 800; } /* e.g., 'm6A background (brief)' */
.content h2 { font-weight: 800; }    

/* === Global theme tokens === */
:root{
  --brand: #0EA5B7;         /* 主品牌色（teal）*/
  --brand-dk:#0D9488;       /* 深一档 */
  --accent:#2563EB;         /* 次强调（蓝）*/
  --ink:#0F172A;            /* 文字主色 */
  --muted:#64748B;          /* 次要文字 */
  --bg1:#F6FAFD;            /* 背景渐变 1 */
  --bg2:#EEF7FB;            /* 背景渐变 2 */
  --panel:#FFFFFF;          /* 卡片底色 */
  --panel-border:#E5E7EB;   /* 卡片描边 */
}

/* === 整站背景（内容区） === */
.content-wrapper, .right-side, .main-footer{
  background: linear-gradient(180deg, var(--bg1) 0%, var(--bg2) 100%) !important;
}

/* 卡片（box）更现代：圆角 + 细描边 + 软阴影 */
.box{
  background: var(--panel) !important;
  border: 1px solid var(--panel-border) !important;
  border-radius: 14px !important;
  box-shadow: 0 1px 2px rgba(0,0,0,.04) !important;
}
.box.solidHeader > .box-header{
  background: rgba(255,255,255,.65) !important;
  border-bottom: 1px solid var(--panel-border) !important;
}


/* === 侧栏整体（白底 + 轻描边 + 毛玻璃感，可选） === */
.main-sidebar{
  background: rgba(255,255,255,.96) !important;
  border-right: 1px solid #E5E7EB !important;
  backdrop-filter: blur(3px);
}

/* 侧栏分组标题 */
.sidebar-menu > li.header{
  margin: 10px 12px 6px;
  color: var(--muted) !important;
  font-weight: 700 !important;
  letter-spacing: .04em;
}

/* === 侧栏菜单“按键”外观（软按钮） === */
.sidebar-menu > li > a{
  margin: 4px 10px;                 /* 与左侧边框拉开距离 */
  padding: 10px 12px;
  border-radius: 10px;              /* 圆角形成“按钮感” */
  color: var(--ink) !important;
  font-weight: 600;
}
.skin-blue .sidebar-menu > li > a i{
  color: #475569 !important;        /* 图标默认色 */
}

/* hover / active 高亮：淡色条+加深图标色 */
.skin-blue .sidebar-menu > li.active > a,
.skin-blue .sidebar-menu > li:hover > a{
  background: linear-gradient(90deg, rgba(14,165,183,.12), rgba(14,165,183,.18)) !important;
  color: var(--ink) !important;
  border-left: 3px solid var(--brand) !important;   /* 左侧高亮条 */
}
.skin-blue .sidebar-menu > li.active > a i,
.skin-blue .sidebar-menu > li:hover  > a i{
  color: var(--brand-dk) !important;
}

/* 侧栏右侧徽章（Data 数字）做成小药丸 */
.sidebar-menu .badge{
  background: #10B981 !important;
  color: #fff !important;
  border-radius: 9999px !important;
  padding: 4px 8px !important;
  font-weight: 700 !important;
}

/* 小屏时侧栏铺白即可 */
@media (max-width: 767px){
  .main-sidebar{ background:#fff !important; }
}


    "))),
    tabItems(
      tabItem(
        "home",
        
        # ——— Welcome / Hero ———
        box(width = 12, title = NULL, solidHeader = FALSE,
            div(class = "hero",
                h1("Welcome — Mini m6A Database & Prediction Server"),
                p(class = "sub",
                  HTML("An academic web resource for <b>curated m<sup>6</sup>A data exploration</b> and <b>site-level prediction</b>.")
                )
                # (per your request: no buttons here)
            )
        ),
        
        # ——— Workflow ———
        h3("Workflow"),
        fluidRow(
          # Left: figure only (no inner title)
          column(
            width = 8,
            box(width = 12, title = NULL, solidHeader = TRUE,
                div(class = "home-figure",
                    imageOutput("home_img", height = "auto")
                )
            )
          ),
          
          # Right: Background (top) + How to use (bottom)
          column(
            width = 4,
            
            # Background (brief) — moved here
            box(
              width = 12, title = "Background", solidHeader = TRUE,
              HTML("
          <ul>
  <li><b>Writers.</b> The METTL3–METTL14 catalytic core, with WTAP, VIRMA/KIAA1429, RBM15/15B and ZC3H13, deposits m<sup>6</sup>A co-transcriptionally on mRNA.</li>
  <li><b>Erasers.</b> FTO and ALKBH5 remove the mark via Fe<sup>2+</sup>/α-ketoglutarate–dependent demethylation.</li>
  <li><b>Readers.</b> YTH (YTHDF1–3; YTHDC1–2), IGF2BPs and hnRNPs interpret m<sup>6</sup>A to control splicing, export, stability and translation.</li>
  <li><b>Genomic context.</b> m<sup>6</sup>A is enriched around stop codons, in long internal exons and 3′UTRs, and depleted in many introns.</li>
  <li><b>Sequence context.</b> Sites often match an RRACH/DRACH-like motif; our <code>DNA_5mer</code> feature is a 5-nt proxy for this signal.</li>
  <li><b>This dataset.</b> Each row represents a candidate site with a label (<code>m6A_status</code>) and a quantitative signal (<code>m6A_level</code>). Predictive covariates capture sequence (<code>DNA_5mer</code>) and transcript context (<code>RNA_region</code>, <code>RNA_type</code>, <code>exon_length</code>, <code>distance_to_junction</code>, <code>gc_content</code>, <code>evolutionary_conservation</code>).</li>
</ul>
        ")
            ),
            
            # How to use
            box(
              width = 12, title = "How to use", solidHeader = TRUE,
              HTML("
          <ol class='use-steps'>
            <li><b>Data.</b> Browse the full table with per-column filters and horizontal scrolling; export filtered subsets.</li>
            <li><b>Visualize.</b> Inspect categorical proportions and histograms of <code>m6A_level</code>.</li>
            <li><b>Prediction → Multiple.</b> Download the template, prepare a CSV, set the decision threshold, then run batch prediction.</li>
            <li><b>Prediction → Single.</b> Enter one record (including <code>DNA_5mer</code>) to obtain probability and status in real time.</li>
          </ol>
          <p class='tight'><span class='pill pill-blue'>Note</span> The threshold governs <i>Positive/Negative</i> calls (≈0.50 balanced; ≥0.70 precision-oriented; ≤0.30 recall-oriented).</p>
        ")
            )
          )
        )
      )
      ,
      
      
      tabItem("data",
              h2("Data Table (DT)"),
              fluidRow(
                valueBoxOutput("dt_rows", width = 4),
                valueBoxOutput("dt_posrate", width = 4),
                valueBoxOutput("dt_cols", width = 4)
              ),
              fluidRow(
                box(width = 12, title = tagList(icon("table"), "Interactive table"),
                    solidHeader = TRUE,
                    div(
                      style = "display:flex;align-items:center;gap:10px;flex-wrap:wrap; margin-bottom:6px;",
                      tags$span(icon("circle-info"), " Tip: use per-column filters (row under headers). "),
                      tags$span("Use ", tags$code("Columns"), " to show/hide; drag headers to reorder; header is fixed while scrolling.")
                    ),
                    div(
                      style = "display:flex;gap:8px;align-items:center;flex-wrap:wrap;margin-bottom:8px;",
                      radioButtons("dt_dl_scope", "Download scope", inline = TRUE,
                                   choices = c("Filtered rows", "Selected rows"), selected = "Filtered rows"),
                      downloadButton("dt_dl_csv", "Download CSV"),
                      actionButton("dt_clear_sel", "Clear selection", icon = icon("xmark"))
                    ),
                    DTOutput("tbl")
                )
              )
      ),
      
      tabItem("viz",
              fluidRow(
                box(width = 12, status = "primary", solidHeader = TRUE,
                    title = tagList(icon("chart-pie"), "Pie chart"),
                    div(
                      style = "display:flex; align-items:center; gap:12px; margin-bottom:8px;",
                      tags$label("Categorical column:", `for` = "pie_col",
                                 style = "margin:0; font-weight:600;"),
                      selectInput("pie_col", NULL,
                                  choices = c("m6A_status","RNA_region","RNA_type"),
                                  selected = "m6A_status", width = "260px")
                    ),
                    plotOutput("pie_plot", height = 400)
                )
              ),
              fluidRow(
                box(width = 12, status = "info", solidHeader = TRUE,
                    title = tagList(icon("chart-bar"), "Histogram of m6A_level"),
                    fluidRow(
                      column(3, selectInput("f_status","m6A_status",
                                            choices = c("No filter", sort(unique(m6a_data$m6A_status))),
                                            selected = "No filter")),
                      column(3, selectInput("f_region","RNA_region",
                                            choices = c("No filter", sort(unique(m6a_data$RNA_region))),
                                            selected = "No filter")),
                      column(3, selectInput("f_type","RNA_type",
                                            choices = c("No filter", sort(unique(m6a_data$RNA_type))),
                                            selected = "No filter")),
                      column(3, sliderInput("bins","Bins:", min = 10, max = 100, value = 30, step = 1))
                    ),
                    plotOutput("hist_plot", height = 400)
                )
              )
      ),
      
      tabItem("down",
              h2("Download data"),
              fluidRow(
                valueBoxOutput("vb_nrows", width = 4),
                valueBoxOutput("vb_ncols", width = 4),
                valueBoxOutput("vb_types", width = 4)
              ),
              fluidRow(
                box(width = 7, title = "Build & download a custom subset", solidHeader = TRUE,
                    p("Choose filters and columns. The preview updates in real time; the download uses the same subset."),
                    fluidRow(
                      column(6, selectizeInput("down_type", "RNA_type (multi)", multiple = TRUE,
                                               choices  = sort(unique(m6a_data$RNA_type)),
                                               selected = sort(unique(m6a_data$RNA_type)))),
                      column(6, selectizeInput("down_region", "RNA_region (multi)", multiple = TRUE,
                                               choices  = sort(unique(m6a_data$RNA_region)),
                                               selected = sort(unique(m6a_data$RNA_region))))
                    ),
                    fluidRow(
                      column(6, selectizeInput("down_status", "m6A_status (multi)", multiple = TRUE,
                                               choices  = sort(unique(m6a_data$m6A_status)),
                                               selected = sort(unique(m6a_data$m6A_status)))),
                      column(6, sliderInput("down_level", "m6A_level range",
                                            min = min(m6a_data$m6A_level, na.rm = TRUE),
                                            max = max(m6a_data$m6A_level, na.rm = TRUE),
                                            value = range(m6a_data$m6A_level, na.rm = TRUE), step = 0.01))
                    ),
                    selectizeInput("down_cols", "Columns to include",
                                   choices  = names(m6a_data),
                                   selected = names(m6a_data),
                                   multiple = TRUE,
                                   options  = list(plugins = list("remove_button"))),
                    fluidRow(
                      column(4, radioButtons("down_fmt", "Format", inline = TRUE,
                                             choices = c("CSV", "TSV", "RDS"), selected = "CSV")),
                      column(4, radioButtons("down_comp", "Compression", inline = TRUE,
                                             choices = c("none", "gzip"), selected = "none")),
                      column(4, textInput("down_name", "File name (no extension)", value = "m6A_subset"))
                    ),
                    downloadButton("dl_subset", "Download subset"),
                    hr(),
                    h4("Preview (first 20 rows)"),
                    DTOutput("tbl_down_preview")
                ),
                box(width = 5, title = "Quick downloads", solidHeader = TRUE,
                    p("One-click exports without changing any filters (CSV only)."),
                    tags$div(style = "display:grid; grid-template-columns:1fr; gap:10px;",
                             downloadButton("dl_all", "All rows (CSV)", style = "width:100%;")
                    ),
                    tags$div(style = "height:8px;"),
                    h4("By RNA_type"),
                    tags$div(style = "display:grid; grid-template-columns:1fr 1fr; gap:8px;",
                             downloadButton("dl_q_mrna",  "mRNA (CSV)",      style = "width:100%;"),
                             downloadButton("dl_q_linc",  "lincRNA (CSV)",   style = "width:100%;"),
                             downloadButton("dl_q_lnc",   "lncRNA (CSV)",    style = "width:100%;"),
                             downloadButton("dl_q_pseudo","pseudogene (CSV)",style = "width:100%;")
                    ),
                    tags$hr(),
                    helpText(icon("circle-info"), "Tip: need TSV/RDS or gzip? Use the custom subset on the left.")
                )
              )
      ),
      
      tabItem("pred_multi",
              h3(tagList(icon("flask"), " m6A Prediction — Multiple (CSV)")),
              div(class = "stepper",
                  div(class="step", div(class="dot","1"), div(class="txt","Upload CSV")),
                  div(class="sep"),
                  div(class="step", div(class="dot","2"), div(class="txt","Set threshold")),
                  div(class="sep"),
                  div(class="step", div(class="dot","3"), div(class="txt","Predict & Download"))
              ),
              fluidRow(
                box(width = 5, title = tagList(icon("file-upload"), "Upload & Settings"), solidHeader = TRUE,
                    
                    ## —— 保留“Choose CSV”这一行 ——（用 label 绑定到 fileInput）
                    ## —— 保留“Choose CSV”这一行 ——（用 label 绑定到 fileInput）
                    tags$label("Choose CSV", `for` = "file", class = "upload-title"),
                    
                    ## —— 同行：Browse + Template + Example ——  ⬅️ 关键：外层用 id="upload-row"
                    ## —— 同行：Browse + Template + Example —— 关键是外层 id="upload-row"
                    div(
                      id = "upload-row",
                      div(
                        class = "upload-file",
                        style = "flex:0 0 300px; max-width:360px; min-width:300px;",  # 想更短改这里
                        fileInput("file", label = NULL, accept = ".csv",
                                  buttonLabel = "Browse…", width = "100%")
                      ),
                      downloadButton("btn_template","Template CSV", class = "btn btn-default"),
                      actionButton("btn_example","Example data",    class = "btn btn-default")
                    ),
                    tags$hr(),
                    h4("Positive threshold"),
                    div(class="btn-row",
                        sliderInput("thresh","", min=0, max=1, value=0.5, step=0.01, width = "100%"),
                        uiOutput("thr_badge")
                    ),
                    tags$hr(),
                    h4("Required columns"),
                    uiOutput("col_check")
                ),
                box(width = 7, title = tagList(icon("table"), "Results"), solidHeader = TRUE,
                    div(class="btn-row results-actions", style = "margin-bottom:12px;",
                        actionButton("btn_predict","Predict", icon = icon("play")),
                        downloadButton("btn_download","Download predicted.csv")
                    ),
                    div(id = "multi-panel", uiOutput("multi_panel"))
                )
                
              )
      ),
      
      tabItem("pred_single",
              h3(tagList(icon("vial"), " m6A Prediction — Single")),
              
              # 顶部：预测摘要（Status / Prob / Thr + 进度条）
              box(width = 12, id = "single-summary", title = "Prediction Summary", solidHeader = TRUE,
                  div(class = "pred-badges",
                      uiOutput("single_status_badge"),
                      uiOutput("single_prob_badge"),
                      uiOutput("single_thr_badge")
                  ),
                  uiOutput("single_prob_bar")   # 横向概率进度条
              ),
              
              # 左大右小：Inputs(7) | Prediction(5)
              fluidRow(
                column(7,
                       box(width=12, id = "single-inputs", title = "Inputs", solidHeader = TRUE,
                           selectInput("s_RNA_type","RNA_type", choices=c("mRNA","lincRNA","lncRNA","pseudogene"), selected="mRNA"),
                           selectInput("s_RNA_region","RNA_region", choices=c("CDS","intron","3'UTR","5'UTR"), selected="CDS"),
                           sliderInput("s_gc","gc_content", 0,1,0.5, step=0.01),
                           sliderInput("s_exon","exon_length", 0,15,10, step=0.1),
                           sliderInput("s_dist","distance_to_junction", 0,15,8, step=0.1),
                           sliderInput("s_cons","evolutionary_conservation", 0,1,0.5, step=0.01),
                           textInput("s_5mer","DNA_5mer (A/T/C/G, length 5)","GGACA"),
                           sliderInput("s_thr","Positive threshold", 0,1,0.5, step=0.01)
                       )
                ),
                column(5,
                       # 上：预测详情
                       box(width=12, id = "single-details", title = "Prediction Details", solidHeader = TRUE,
                           uiOutput("single_out"),
                           tags$hr(),
                           div(style="display:flex; gap:8px; flex-wrap:wrap;",
                               actionButton("single_add","Add to compare", icon=icon("plus")),
                               actionButton("single_clear","Clear list", icon=icon("trash")),
                               downloadButton("single_dl","Download compare.csv")
                           )
                       ),
                       # 下：对比表（与上面对齐）
                       box(width = 12, id = "single-compare", title = "Compare multiple tries", solidHeader = TRUE,
                           DTOutput("single_compare_tbl", height = "260px")
                       )
                )
              )     
      ),
      tabItem("about",
              h2("Help"),
              fluidRow(
                box(width = 6, title = "Quick Start", solidHeader = TRUE,
                    HTML("
        <ol>
          <li><b>Data</b>: Browse the full table; use per-column filters and horizontal scroll.</li>
          <li><b>Visualize</b>: Inspect categorical proportions and histograms of <code>m6A_level</code>.</li>
          <li><b>Prediction → Multiple</b>: Download the <i>Template CSV</i>, fill required columns, upload, set the threshold, then click <i>Predict</i>.</li>
          <li><b>Prediction → Single</b>: Enter one record in the form to obtain probability and status in real time.</li>
          <li><b>Download</b>: Export all rows or filtered subsets.</li>
        </ol>
        ")
                ),
                box(width = 6, title = "Model & Intended Use", solidHeader = TRUE,
                    HTML("
        <p>This app uses a <b>Random Forest</b> classifier to estimate the probability that a site is m<sup>6</sup>A-positive. 
        Features: <code>gc_content</code>, <code>RNA_type</code>, <code>RNA_region</code>, <code>exon_length</code>, 
        <code>distance_to_junction</code>, <code>evolutionary_conservation</code>, and <code>DNA_5mer</code> (one-hot encoded to <code>nt_pos1–5</code>).</p>
        <p><b>Intended use</b>: instructional and exploratory analyses. <b>Limitations</b>: not experimentally validated in all contexts; 
        performance can vary across datasets; results should not be used for clinical decision-making.</p>
        ")
                )
              ),
              fluidRow(
                box(width = 12, title = "Data Dictionary", solidHeader = TRUE,
                    HTML("
        <div style='overflow:auto'>
        <table style='border-collapse:collapse;width:100%'>
          <thead>
            <tr>
              <th style='border:1px solid #e5e7eb;padding:6px'>Field</th>
              <th style='border:1px solid #e5e7eb;padding:6px'>Type</th>
              <th style='border:1px solid #e5e7eb;padding:6px'>Description</th>
            </tr>
          </thead>
          <tbody>
            <tr><td style='border:1px solid #e5e7eb;padding:6px'>m6A_level</td><td style='border:1px solid #e5e7eb;padding:6px'>numeric</td><td style='border:1px solid #e5e7eb;padding:6px'>Quantitative signal used in visualization.</td></tr>
            <tr><td style='border:1px solid #e5e7eb;padding:6px'>m6A_status</td><td style='border:1px solid #e5e7eb;padding:6px'>factor</td><td style='border:1px solid #e5e7eb;padding:6px'>Label (Positive/Negative) in the dataset.</td></tr>
            <tr><td style='border:1px solid #e5e7eb;padding:6px'>RNA_type</td><td style='border:1px solid #e5e7eb;padding:6px'>factor</td><td style='border:1px solid #e5e7eb;padding:6px'>mRNA / lincRNA / lncRNA / pseudogene.</td></tr>
            <tr><td style='border:1px solid #e5e7eb;padding:6px'>RNA_region</td><td style='border:1px solid #e5e7eb;padding:6px'>factor</td><td style='border:1px solid #e5e7eb;padding:6px'>CDS / intron / 3'UTR / 5'UTR.</td></tr>
            <tr><td style='border:1px solid #e5e7eb;padding:6px'>gc_content</td><td style='border:1px solid #e5e7eb;padding:6px'>numeric (0–1)</td><td style='border:1px solid #e5e7eb;padding:6px'>GC fraction.</td></tr>
            <tr><td style='border:1px solid #e5e7eb;padding:6px'>exon_length</td><td style='border:1px solid #e5e7eb;padding:6px'>numeric</td><td style='border:1px solid #e5e7eb;padding:6px'>Exon length (scaled/log as in training).</td></tr>
            <tr><td style='border:1px solid #e5e7eb;padding:6px'>distance_to_junction</td><td style='border:1px solid #e5e7eb;padding:6px'>numeric</td><td style='border:1px solid #e5e7eb;padding:6px'>Distance to splice junction (scaled).</td></tr>
            <tr><td style='border:1px solid #e5e7eb;padding:6px'>evolutionary_conservation</td><td style='border:1px solid #e5e7eb;padding:6px'>numeric (0–1)</td><td style='border:1px solid #e5e7eb;padding:6px'>Conservation score.</td></tr>
            <tr><td style='border:1px solid #e5e7eb;padding:6px'>DNA_5mer</td><td style='border:1px solid #e5e7eb;padding:6px'>string (len=5)</td><td style='border:1px solid #e5e7eb;padding:6px'>Window of five bases (A/T/C/G only, uppercase).</td></tr>
            <tr><td style='border:1px solid #e5e7eb;padding:6px'>nt_pos1–5</td><td style='border:1px solid #e5e7eb;padding:6px'>factor (A/T/C/G)</td><td style='border:1px solid #e5e7eb;padding:6px'>One-hot expansion of DNA_5mer (not required in uploads).</td></tr>
            <tr><td style='border:1px solid #e5e7eb;padding:6px'>predicted_m6A_prob</td><td style='border:1px solid #e5e7eb;padding:6px'>numeric (0–1)</td><td style='border:1px solid #e5e7eb;padding:6px'>Positive-class probability from the model.</td></tr>
            <tr><td style='border:1px solid #e5e7eb;padding:6px'>predicted_m6A_status</td><td style='border:1px solid #e5e7eb;padding:6px'>factor</td><td style='border:1px solid #e5e7eb;padding:6px'>Thresholded call (Positive/Negative).</td></tr>
          </tbody>
        </table></div>
        ")
                )
              ),
              fluidRow(
                box(width = 7, title = "Input File Format (CSV)", solidHeader = TRUE,
                    HTML("
        <p><b>Required columns</b>: <code>gc_content, RNA_type, RNA_region, exon_length, distance_to_junction, evolutionary_conservation, DNA_5mer</code></p>
        <ul>
          <li><b>RNA_type</b> ∈ {mRNA, lincRNA, lncRNA, pseudogene}</li>
          <li><b>RNA_region</b> ∈ {CDS, intron, 3'UTR, 5'UTR}</li>
          <li><b>DNA_5mer</b>: exactly 5 uppercase bases from A/T/C/G.</li>
          <li>Numeric fields must use the dot decimal separator; UTF-8 CSV; comma delimiter.</li>
        </ul>
        <p>Use the <i>Template CSV</i> to avoid column-name or order problems.</p>
        ")
                ),
                box(width = 5, title = "Prediction Outputs & Threshold", solidHeader = TRUE,
                    HTML("
        <ul>
          <li><b>Probability</b> (<code>predicted_m6A_prob</code>): 0–1.</li>
          <li><b>Status</b> (<code>predicted_m6A_status</code>): <i>Positive</i> if probability &gt; threshold; else <i>Negative</i>.</li>
          <li><b>Threshold tips</b>: 0.50 (balanced), ≥0.70 (precision↑), ≤0.30 (recall↑).</li>
        </ul>
        <p>The downloaded <code>predicted.csv</code> appends both fields to the original columns.</p>
        ")
                )
              ),
              fluidRow(
                box(width = 12, title = "FAQ / Troubleshooting", solidHeader = TRUE,
                    HTML("
        <ul>
          <li><b>Missing required columns</b>: Re-create the file using the Template CSV.</li>
          <li><b>DNA_5mer validation failed</b>: Not exactly five A/T/C/G bases.</li>
          <li><b>Unknown RNA_type / RNA_region</b>: Values not in allowed sets.</li>
          <li><b>Encoding / delimiter issues</b>: Use UTF-8, comma-separated.</li>
        </ul>
        ")
                )
              ),
              fluidRow(
                box(width = 6, title = "Citation & Credits", solidHeader = TRUE,
                    HTML("<p>Please cite as: “Mini m6A Database & Predictor (Shiny App), author: Xiyao Chen, 2025.”</p>")
                ),
                box(width = 6, title = "Version", solidHeader = TRUE,
                    HTML("<ul style='margin-top:0'><li>App version: v1.0</li><li>Last updated: 2025-11-02</li></ul>")
                )
              )
      )
      
    )
  )
)

# ---- Server ----
server <- function(input, output, session){
  
  # Home illustration 渲染（把图片放在 app 根目录或 www/ 下）
  # Home illustration
  #output$home_img <- renderImage({
  # paths <- c("www/Mech_figure2.jpg", "Mech_figure2.jpg", "www/m6A_mechanism.png")
  #  p <- paths[file.exists(paths)][1]
  # if (is.na(p)) return(list(src = NULL))
  # ctype <- if (grepl("\\.png$", p, ignore.case = TRUE)) "image/png" else "image/jpeg"
  # list(src = p, contentType = ctype)  # 尺寸交给 CSS 控制
  #  }, deleteFile = FALSE)
  output$home_img <- renderImage({
    # 优先找 PNG，其次 JPG；可放在 app 根目录或 www/
    paths <- c("www/Mech_figure2.png", "Mech_figure2.png")
    p <- paths[file.exists(paths)][1]
    validate(need(length(p) && !is.na(p), "Please put Mech_figure2.png (or .jpg) into app/ or www/"))
    list(
      src = p,
      contentType = if (grepl("\\.png$", p, ignore.case = TRUE)) "image/png" else "image/jpeg",
      alt = "Workflow of the Mini m6A App"
      # 不设 width/height，让 CSS 控制，避免被固定高度裁剪
    )
  }, deleteFile = FALSE)
  
  
  
  
  
  # ---------- Data tab ----------
  num_cols <- intersect(names(m6a_data),
                        c("m6A_level","gc_content","evolutionary_conservation",
                          "m6A_coverage","total_coverage","distance_to_junction","exon_length"))
  
  output$tbl <- DT::renderDT({
    DT::datatable(
      m6a_data,
      filter    = "top",
      rownames  = FALSE,
      selection = "multiple",
      class     = "display nowrap compact",
      extensions = c("Buttons","ColReorder","FixedHeader","Scroller"),
      options = list(
        dom = "Bfrtip",
        buttons = list(
          list(extend = "colvis", text = "Columns"),
          list(extend = "copy",  text = "Copy"),
          list(extend = "csv",   text = "CSV"),
          list(extend = "excel", text = "Excel"),
          list(extend = "print", text = "Print")
        ),
        colReorder   = TRUE,
        fixedHeader  = TRUE,
        deferRender  = TRUE,
        scroller     = TRUE,
        scrollX      = TRUE,
        scrollY      = 420,
        pageLength   = 15,
        lengthMenu   = c(15,30,50,100),
        autoWidth    = TRUE
      )
    ) |>
      DT::formatRound(num_cols, digits = 3) |>
      DT::formatStyle(
        "m6A_status",
        color = DT::styleEqual(c("Positive","Negative"), c("#0D9488","#DC2626")),
        fontWeight = "bold"
      ) |>
      DT::formatStyle(
        "RNA_region",
        backgroundColor = DT::styleEqual(
          c("CDS","intron","3'UTR","5'UTR"),
          c("#DBEAFE","#DCFCE7","#FEF3C7","#FEE2E2")
        )
      )
  }, server = TRUE)
  
  data_filtered_now <- reactive({
    idx <- input$tbl_rows_all
    if (is.null(idx)) m6a_data else m6a_data[idx, , drop = FALSE]
  })
  data_selected_now <- reactive({
    sel <- input$tbl_rows_selected
    base <- data_filtered_now()
    if (is.null(sel) || !length(sel)) base[0, , drop = FALSE] else base[sel, , drop = FALSE]
  })
  
  output$dt_rows <- renderValueBox({
    valueBox(scales::comma(nrow(data_filtered_now())), "Rows (after filters)", icon = icon("database"), color = "light-blue")
  })
  output$dt_cols <- renderValueBox({
    valueBox(ncol(m6a_data), "Columns", icon = icon("columns"), color = "teal")
  })
  output$dt_posrate <- renderValueBox({
    d <- data_filtered_now()
    rate <- if (nrow(d)) mean(d$m6A_status == "Positive") else 0
    valueBox(scales::percent(rate, accuracy = 0.1), "Positive rate (filtered)", icon = icon("check"), color = "aqua")
  })
  observeEvent(input$dt_clear_sel, { DT::selectRows(DT::dataTableProxy("tbl"), NULL) })
  
  output$dt_dl_csv <- downloadHandler(
    filename = function() {
      scope <- if (identical(input$dt_dl_scope, "Selected rows")) "selected" else "filtered"
      paste0("m6A_", scope, ".csv")
    },
    content = function(file) {
      if (identical(input$dt_dl_scope, "Selected rows")) {
        out <- data_selected_now(); if (!nrow(out)) out <- data_filtered_now()
      } else out <- data_filtered_now()
      write.csv(out, file, row.names = FALSE)
    }
  )
  
  # ---------- Visualize ----------
  output$pie_plot <- renderPlot({
    req(input$pie_col)
    label_threshold <- 0.20
    df <- m6a_data %>%
      dplyr::count(level = .data[[input$pie_col]], name = "n") %>%
      dplyr::mutate(prop = n / sum(n)) %>%
      dplyr::arrange(dplyr::desc(prop)) %>%
      dplyr::mutate(
        lab  = ifelse(prop >= label_threshold,
                      paste0(level, "\n", scales::percent(prop, accuracy = 0.1)),
                      NA_character_),
        emph = ifelse(prop == max(prop), "max", "other")
      )
    lvls <- unique(df$level)
    pal <- if (input$pie_col == "m6A_status") {
      setNames(c("#F06A67", "#12B7B6"), c("Negative","Positive"))
    } else if (input$pie_col == "RNA_region") {
      setNames(c("#3B82F6", "#10B981", "#F59E0B", "#EF4444")[seq_along(lvls)], lvls)
    } else {
      setNames(scales::hue_pal(l = 65, c = 100)(length(lvls)), lvls)
    }
    deg <- 35; lab_y <- sum(df$n) * deg/360; lab_x <- 2.6
    ggplot2::ggplot(df, ggplot2::aes(x = 2, y = n, fill = level, alpha = emph)) +
      ggplot2::geom_col(width = 1, color = "white", linewidth = 0.6) +
      ggplot2::coord_polar(theta = "y", clip = "off") +
      ggplot2::xlim(0.5, 2.6) +
      ggplot2::geom_text(ggplot2::aes(label = lab),
                         position = ggplot2::position_stack(vjust = 0.5),
                         size = 4, lineheight = 0.9, na.rm = TRUE) +
      ggplot2::annotate("text", x = lab_x, y = lab_y,
                        label = paste0("N = ", sum(df$n), "\n", input$pie_col),
                        hjust = 0, vjust = 0, fontface = 2, size = 5, lineheight = 1.1) +
      ggplot2::scale_fill_manual(values = pal, guide = ggplot2::guide_legend(title = input$pie_col)) +
      ggplot2::scale_alpha_manual(values = c(max = 1, other = 0.55), guide = "none") +
      ggplot2::labs(title = paste("Distribution of", input$pie_col)) +
      ggplot2::theme_void(base_size = 14) +
      ggplot2::theme(
        plot.title   = ggplot2::element_text(face = "bold", size = 16, hjust = 0.5,
                                             margin = ggplot2::margin(b = 6)),
        legend.title = ggplot2::element_text(face = "bold"),
        legend.text  = ggplot2::element_text(size = 11),
        plot.margin  = ggplot2::margin(t = 8, r = 30, b = 8, l = 8)
      )
  })
  
  filtered <- reactive({
    d <- m6a_data
    if (!is.null(input$f_status) && input$f_status!="No filter") d <- d %>% filter(m6A_status==input$f_status)
    if (!is.null(input$f_region) && input$f_region!="No filter") d <- d %>% filter(RNA_region==input$f_region)
    if (!is.null(input$f_type)   && input$f_type  !="No filter") d <- d %>% filter(RNA_type  ==input$f_type)
    d
  })
  
  output$hist_plot <- renderPlot({
    d <- filtered()
    validate(need(nrow(d) > 0, "No data for current filters"))
    mu  <- mean(d$m6A_level, na.rm = TRUE)
    md  <- median(d$m6A_level, na.rm = TRUE)
    q05 <- quantile(d$m6A_level, 0.05, na.rm = TRUE)
    q95 <- quantile(d$m6A_level, 0.95, na.rm = TRUE)
    ggplot(d, aes(x = m6A_level)) +
      annotate("rect", xmin = q05, xmax = q95, ymin = 0, ymax = Inf,
               fill = "#3B82F6", alpha = 0.08) +
      geom_histogram(bins = input$bins,
                     fill = "#60A5FA", color = "white", linewidth = 0.3, alpha = 0.9) +
      geom_vline(xintercept = mu, color = "#111827", linewidth = 0.8) +
      geom_vline(xintercept = md, color = "#EF4444", linewidth = 0.8, linetype = "dashed") +
      annotate("text", x = mu, y = Inf, vjust = -0.6,
               label = sprintf("mean = %.2f", mu),
               color = "#111827", size = 3.2, fontface = "bold") +
      annotate("text", x = md, y = Inf, vjust = -2.0,
               label = sprintf("median = %.2f", md),
               color = "#EF4444", size = 3.0, fontface = "bold") +
      labs(x = "m6A_level", y = "Count") +
      scale_y_continuous(labels = scales::comma) +
      coord_cartesian(clip = "off") +
      theme_minimal(base_size = 12) +
      theme(
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.title.x = element_text(face = "bold", size = 11),
        axis.title.y = element_text(face = "bold", size = 11),
        plot.margin = ggplot2::margin(t = 4, r = 10, b = 4, l = 6)
      )
  }, res = 140)
  
  # ---------- Download ----------
  down_filtered_raw <- reactive({
    d <- m6a_data
    if (!is.null(input$down_type)   && length(input$down_type))   d <- d %>% dplyr::filter(.data$RNA_type   %in% input$down_type)
    if (!is.null(input$down_region) && length(input$down_region)) d <- d %>% dplyr::filter(.data$RNA_region %in% input$down_region)
    if (!is.null(input$down_status) && length(input$down_status)) d <- d %>% dplyr::filter(.data$m6A_status %in% input$down_status)
    if (!is.null(input$down_level) && length(input$down_level) == 2) {
      d <- d %>% dplyr::filter(.data$m6A_level >= input$down_level[1],
                               .data$m6A_level <= input$down_level[2])
    }
    d
  })
  down_subset <- reactive({
    d <- down_filtered_raw()
    cols <- input$down_cols
    if (is.null(cols) || !length(cols)) cols <- names(d)
    d[, cols, drop = FALSE]
  })
  output$vb_nrows <- renderValueBox({
    valueBox(scales::comma(nrow(down_filtered_raw())), "Rows in subset", icon = icon("database"), color = "light-blue")
  })
  output$vb_ncols <- renderValueBox({
    valueBox(scales::comma(ncol(down_subset())), "Columns selected", icon = icon("columns"), color = "teal")
  })
  output$vb_types <- renderValueBox({
    d <- down_filtered_raw(); k <- length(unique(d$RNA_type))
    valueBox(paste0(k, " RNA types"), "Distinct RNA_type", icon = icon("filter"), color = "aqua")
  })
  output$tbl_down_preview <- renderDT({
    datatable(utils::head(down_subset(), 20),
              options = list(scrollX = TRUE, pageLength = 10),
              rownames = FALSE)
  })
  output$dl_subset <- downloadHandler(
    filename = function() {
      nm   <- if (nzchar(input$down_name)) input$down_name else "m6A_subset"
      ext  <- switch(input$down_fmt, CSV = "csv", TSV = "tsv", RDS = "rds")
      gzok <- (input$down_fmt %in% c("CSV","TSV")) && identical(input$down_comp, "gzip")
      paste0(nm, ".", ext, if (gzok) ".gz" else "")
    },
    content = function(file) {
      d <- down_subset(); fmt <- input$down_fmt
      if (fmt == "RDS") {
        saveRDS(d, file = file, compress = "gzip")
      } else {
        sep <- if (fmt == "CSV") "," else "\t"
        if (identical(input$down_comp, "gzip")) {
          con <- gzfile(file, "wt"); on.exit(close(con), add = TRUE)
          write.table(d, file = con, sep = sep, row.names = FALSE, col.names = TRUE, quote = FALSE)
        } else {
          write.table(d, file = file, sep = sep, row.names = FALSE, col.names = TRUE, quote = FALSE)
        }
      }
    }
  )
  output$dl_all     <- downloadHandler("m6A_all.csv",        function(file) write.csv(m6a_data, file, row.names = FALSE))
  output$dl_q_mrna  <- downloadHandler("m6A_mRNA.csv",       function(file) write.csv(dplyr::filter(m6a_data, RNA_type == "mRNA"),       file, row.names = FALSE))
  output$dl_q_linc  <- downloadHandler("m6A_lincRNA.csv",    function(file) write.csv(dplyr::filter(m6a_data, RNA_type == "lincRNA"),    file, row.names = FALSE))
  output$dl_q_lnc   <- downloadHandler("m6A_lncRNA.csv",     function(file) write.csv(dplyr::filter(m6a_data, RNA_type == "lncRNA"),     file, row.names = FALSE))
  output$dl_q_pseudo<- downloadHandler("m6A_pseudogene.csv", function(file) write.csv(dplyr::filter(m6a_data, RNA_type == "pseudogene"), file, row.names = FALSE))
  
  # ---------- Prediction Multiple ----------
  raw_df  <- reactiveVal(NULL); pred_on <- reactiveVal(FALSE)
  observeEvent(input$file, { req(input$file); raw_df(read.csv(input$file$datapath, stringsAsFactors=FALSE)); pred_on(FALSE) })
  observeEvent(input$btn_example, { raw_df(example_data); pred_on(FALSE) })
  
  required_cols <- c("gc_content","RNA_type","RNA_region","exon_length",
                     "distance_to_junction","evolutionary_conservation","DNA_5mer")
  get_clean <- reactive({
    dat <- raw_df()
    if (is.null(dat)) return(NULL)
    # 只取需要列，缺列则在 UI 徽章里可见；Predict 时再用 req 拦截
    dat[, intersect(required_cols, colnames(dat)), drop=FALSE]
  })
  make_pred <- function(dat, thr) prediction_multiple(rf_fit, dat, positive_threshold=thr)
  
  observeEvent(input$btn_predict, {
    req(get_clean())
    # 确保所有必需列具备
    miss <- setdiff(required_cols, colnames(get_clean()))
    validate(need(!length(miss), paste("Missing required columns:", paste(miss, collapse=", "))))
    raw_df(make_pred(get_clean(), input$thresh))
    pred_on(TRUE)
  })
  observeEvent(input$thresh, {
    if (isTRUE(pred_on()) && !is.null(get_clean())) raw_df(make_pred(get_clean(), input$thresh))
  }, ignoreInit=TRUE)
  
  output$btn_template <- downloadHandler(
    filename = function() "m6A_template.csv",
    content = function(file){
      tmpl <- data.frame(
        gc_content = numeric(),
        RNA_type = character(),
        RNA_region = character(),
        exon_length = numeric(),
        distance_to_junction = numeric(),
        evolutionary_conservation = numeric(),
        DNA_5mer = character(),
        stringsAsFactors = FALSE
      )
      write.csv(tmpl, file, row.names = FALSE)
    }
  )
  
  output$col_check <- renderUI({
    dat  <- raw_df()
    cols <- if (is.null(dat)) character(0) else colnames(dat)
    tags$div(
      id = "reqcols",  # ← 用上面的 CSS 作用域
      lapply(required_cols, function(x){
        ok <- x %in% cols
        span(
          class = paste("pill", if (ok) "pill-green" else "pill-red"),
          HTML(paste0(
            if (ok) "<i class='fa fa-check'></i> " else "<i class='fa fa-times'></i> "),
            x
          )
        )
      })
    )
  })
  
  
  output$thr_badge <- renderUI({
    span(class="pill pill-blue", sprintf("Threshold = %.2f", input$thresh))
  })
  
  output$multi_summary <- renderUI({
    dat <- raw_df()
    if (is.null(dat) || !("predicted_m6A_status" %in% names(dat))) return(NULL)
    n <- nrow(dat)
    pos <- sum(dat$predicted_m6A_status == "Positive", na.rm = TRUE)
    rate <- if (n > 0) pos / n else 0
    tags$div(
      span(class="pill pill-blue",  paste0("Rows: ", scales::comma(n))), HTML("&nbsp;"),
      span(class="pill pill-green", paste0("Positive: ", pos, " (", scales::percent(rate, accuracy = 0.1), ")")), HTML("&nbsp;"),
      span(class="pill pill-gray",  sprintf("Current threshold: %.2f", input$thresh))
    )
  })
  
  output$tbl_pred <- renderDT({
    dat <- raw_df(); if (is.null(dat)) dat <- data.frame()
    dt <- datatable(dat, options = list(scrollX = TRUE, pageLength = 15), rownames = FALSE)
    if ("predicted_m6A_prob" %in% names(dat)) {
      dt <- dt %>% DT::formatRound("predicted_m6A_prob", 3) %>%
        DT::formatStyle("predicted_m6A_prob",
                        backgroundColor = DT::styleInterval(input$thresh, c("#FEE2E2", "#DCFCE7")))
    }
    if ("predicted_m6A_status" %in% names(dat)) {
      dt <- dt %>% DT::formatStyle("predicted_m6A_status",
                                   color = DT::styleEqual(c("Positive","Negative"), c("#0D9488","#DC2626")),
                                   fontWeight = "bold")
    }
    dt
  })
  output$multi_panel <- renderUI({
    dat <- raw_df()
    
    # 1) 未上传文件：显示空态占位，避免一整片空白
    if (is.null(dat)) {
      return(
        div(class="empty-state",
            tags$i(class="fa fa-upload"),
            h4("Upload a CSV to begin"),
            p("Use the ", tags$b("Template CSV"), " to avoid column issues."),
            span(class="pill pill-blue", sprintf("Threshold = %.2f", input$thresh))
        )
      )
    }
    
    # 2) 已上传但尚未 Predict：给个简短提示 + 小预览（前 8 行）
    if (!isTRUE(pred_on())) {
      return(tagList(
        div(class="btn-row",
            span(class="pill pill-blue", paste("Rows:", nrow(dat))),
            span(class="pill pill-gray", "Click 'Predict' to generate probabilities.")
        ),
        tags$hr(),
        DT::datatable(utils::head(dat, 8),
                      options = list(dom = 't', pageLength = 8, scrollX = TRUE),
                      rownames = FALSE)
      ))
    }
    
    # 3) 已预测：展示完整结果表（沿用你现有的 output$tbl_pred）
    DTOutput("tbl_pred", height = "430px")
  })
  
  
  output$btn_download <- downloadHandler(
    filename = function() "predicted.csv",
    content = function(file){
      req(get_clean())
      out <- prediction_multiple(rf_fit, get_clean(), input$thresh)
      write.csv(out, file, row.names = FALSE)
    }
  )
  
  # ---------- Prediction Single ----------
  single_pred <- reactive({
    validate(need(is_valid_5mer(input$s_5mer), "DNA_5mer must be 5 characters of A/T/C/G."))
    prediction_single(rf_fit,
                      gc_content=input$s_gc, RNA_type=input$s_RNA_type, RNA_region=input$s_RNA_region,
                      exon_length=input$s_exon, distance_to_junction=input$s_dist,
                      evolutionary_conservation=input$s_cons, DNA_5mer=input$s_5mer,
                      positive_threshold=input$s_thr)
  })
  
  output$single_status_badge <- renderUI({
    res <- single_pred()
    is_pos <- identical(res[["predicted_m6A_status"]], "Positive")
    span(class = paste("pill", if (is_pos) "pill-green" else "pill-red"),
         HTML(paste0(if (is_pos) "<i class='fa fa-thumbs-up'></i> " else "<i class='fa fa-thumbs-down'></i> "),
              "Status: ", res[["predicted_m6A_status"]]))
  })
  output$single_prob_badge <- renderUI({
    res <- single_pred()
    span(class="pill pill-blue", paste0("Prob: ", res[["predicted_m6A_prob"]]))
  })
  output$single_thr_badge <- renderUI({
    span(class="pill pill-gray", sprintf("Thr: %.2f", input$s_thr))
  })
  
  output$single_out <- renderUI({
    res <- single_pred()
    tagList(
      h5("predicted_m6A_prob: ", strong(res[["predicted_m6A_prob"]])),
      h5("predicted_m6A_status: ", strong(res[["predicted_m6A_status"]]))
    )
  })
  # --- 顶部概率条形 UI（根据 single_pred 的结果动态渲染） ---
  output$single_prob_bar <- renderUI({
    res <- single_pred()
    p <- suppressWarnings(as.numeric(res[["predicted_m6A_prob"]]))
    thr <- suppressWarnings(as.numeric(input$s_thr))
    if (is.na(p) || is.na(thr)) return(NULL)
    pct <- round(p * 100)
    cls <- if (p >= thr) "progress-bar-green" else "progress-bar-red"  # 绿/红
    div(class="progress", style="height:14px; margin-top:6px;",
        div(class=paste("progress-bar", cls), role="progressbar",
            style=paste0("width:", pct, "%;"),
            `aria-valuenow`=pct, `aria-valuemin`=0, `aria-valuemax`=100,
            paste0(pct, "%")
        )
    )
    
    
  })
  
  
  # --- Compare 多次尝试：累积表 ---
  # 保存对比条目
  single_store <- reactiveVal(data.frame())
  
  # “Add to compare”
  observeEvent(input$single_add, {
    res <- single_pred()  # 你已有的函数，返回带 prob/status 的命名向量
    row <- data.frame(
      RNA_type = input$s_RNA_type,
      RNA_region = input$s_RNA_region,
      gc_content = input$s_gc,
      exon_length = input$s_exon,
      distance_to_junction = input$s_dist,
      evolutionary_conservation = input$s_cons,
      DNA_5mer = toupper(input$s_5mer),
      thr = input$s_thr,
      predicted_m6A_prob   = as.numeric(res[["predicted_m6A_prob"]]),
      predicted_m6A_status = as.character(res[["predicted_m6A_status"]]),
      stringsAsFactors = FALSE
    )
    cur <- single_store()
    single_store(rbind(cur, row))
  })
  
  # “Clear list”
  observeEvent(input$single_clear, {
    single_store(data.frame())
  })
  
  # 右侧下方表格
  output$single_compare_tbl <- DT::renderDT({
    dat <- single_store()
    if (is.null(dat) || nrow(dat) == 0) {
      dat <- data.frame(Message = "List is empty. Click 'Add to compare'.")
    }
    tbl <- DT::datatable(dat, options = list(scrollX = TRUE, pageLength = 6), rownames = FALSE)
    if ("predicted_m6A_prob" %in% names(dat)) {
      tbl <- DT::formatRound(tbl, "predicted_m6A_prob", 3)
    }
    tbl
  })
  
  # 下载 compare.csv（注意 id 要与 UI 的 single_dl 一致）
  output$single_dl <- downloadHandler(
    filename = function() "compare.csv",
    content = function(file) {
      dat <- single_store()
      write.csv(dat, file, row.names = FALSE)
    }
  )
}

# ---- Run ----
app <- shinyApp(ui, server)
app
