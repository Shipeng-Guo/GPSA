# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

library(dplyr)
library(tidyr)
library(stringr)
library(enrichplot)
library(clusterProfiler)

### plot
geneplot <- function(dataset,gene){
  data = alltidydata[[dataset]]
  library(ggpubr)
  my_comparisons <- list(c("treat", "control"))
  ggboxplot(
    data, x = "group", y = gene,
    color = "group", palette = "jco",
    add = "jitter"
  )+
    stat_compare_means(comparisons = my_comparisons, method = "t.test")
}

trimKOSA <- function(y){
  dd = as.data.frame(y)
  out = dd %>%
    filter(p.adjust < 0.05) %>%
    ## Species
    mutate(Species = str_match(Description,".*(human|mouse|HUMAN|MOUSE).*")[,2]) %>%
    mutate(Species = toupper(Species)) %>%
    filter(Species =="HUMAN") %>%
    ## 增加一列upanddown
    mutate(upanddown = str_match(Description,".*_(.*)")[,2]) %>%
    mutate(upanddown = ifelse(upanddown=="UP",1,-1)) %>%
    ## 增加处理信息
    mutate(treat = -1) %>%
    ## NES结果
    mutate(NESdirect = sign(NES) ) %>%
    ## 得到1或者-1
    mutate(result = -treat*upanddown*NESdirect) %>%
    ## 用于排序
    mutate(NES_new = -treat*upanddown*NES) %>%
    ### 提取分组信息
    mutate(group = str_match(Description,"(.*)_(UP|DOWN)")[,2]) %>%
    ### 分组信息包含内容需要是2
    group_by(group) %>%
    mutate(number=n()) %>%
    filter(number==2) %>%
    ### 再次基础上，两个结果要统一，汇总信息
    summarise_at(c("result","NES_new"),sum) %>%
    ### 删掉非配对数据,很重要
    filter(result!=0) %>%
    ungroup() %>%
    ### 添加信息
    mutate(effect=ifelse(result>0,"Inhibited","Activated")) %>%
    mutate(Gene = str_match(group,".*_(.*)")[,2]) %>%
    mutate(Gene = toupper(Gene)) %>%
    dplyr::select(Gene,group,NES_new,effect) %>%
    ### 绝对值排序
    arrange(desc(abs(NES_new))) %>%
    rename(simScore=NES_new)
}

### plottrimdata
plotytrimdata <- function(trimdata=out,showCategory = 30){
  nrow = nrow(trimdata)
  showCategory = min(showCategory,nrow)
  trimdata = trimdata[1:showCategory,]
  ggplot(trimdata,aes(simScore,forcats::fct_reorder(group,simScore))) +
    ## 画横线
    geom_segment(aes(xend=0, yend = group)) +
    ## 画点
    geom_point(aes(color=simScore,size = abs(simScore))) +
    ## 调整颜色的区间,begin越大，整体颜色越明艳
    #scale_color_viridis_c(begin = 0.3, end = 1) +
    scale_color_gradient(low="#2ECC71", high="#7D3C98")+
    ## 调整泡泡的大小
    scale_size_continuous(range=c(2,8)) +
    theme_bw() +
    xlab("simScore") +
    ylab(NULL) +
    ggtitle("")
}

### 不美观，需要修改。
plotGSEApair <- function(term,y){
  Gene = str_match(term,".*_(.*)")[,2]
  treat = "KD"
  GSE = str_match(term,"(.*?)_.*")[,2]

  termfull1 = paste0(term,"_DOWN")
  title1 = paste(Gene,treat,GSE,"DOWN",sep = " ")
  y1 = filter(y,Description %in% termfull1)
  y1@result$Description <- ""
  g1 <- gseaplot2(y1,termfull1,title = title1,color = "red",pvalue_table = T)


  termfull2 = paste0(term,"_UP")
  title2 = paste(Gene,treat,GSE,"UP",sep = " ")
  y2 = filter(y,Description %in% termfull2)
  y2@result$Description <- ""
  g2 <- gseaplot2(y2,termfull2,title = title2,color = "red",pvalue_table = T)

  plot_grid(g1,g2,nrow = 2,labels = letters[1:2])
}


plotGSEApair2 <- function(term,y,trimdata=out,color = c("#2ECC71","#7D3C98")){
  Gene = str_match(term,".*_(.*)")[,2]
  treat = "KD"
  GSE = str_match(term,"(.*?)_.*")[,2]
  title = paste(Gene,treat,GSE,sep = " ")

  simScore = trimdata$simScore[trimdata$group==term]
  simScore = round(simScore,2)
  effect = trimdata$effect[trimdata$group==term]
  title = paste0(title,", simScore=",simScore,", ",effect)
  termfull1 = paste0(term,"_DOWN")
  termfull2 = paste0(term,"_UP")
  y1 = filter(y,Description %in% c(termfull1,termfull2))
  y1@result$Description <- str_match(y1@result$Description,".*_(.*)")[,2]
  gseaplot2(y1,c(termfull1,termfull2),title = title,pvalue_table = T,color=color)
}

### GPSA函数
### maxGSSize=500 依靠这个来筛选数据集
GPSA <- function(geneList,TERM2GENE = kosageneset,maxGSSize=500){
  kosageneset=dplyr::filter(kosageneset,rank <= maxGSSize)
  y=GSEA(geneList=geneList,TERM2GENE = kosageneset,maxGSSize = maxGSSize)
  return(y)
}

