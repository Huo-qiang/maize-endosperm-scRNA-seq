##使用国内镜像安装包

options("repos"= c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
options(BioC_mirror="http://mirrors.ustc.edu.cn/bioc/")
install.packages("Cairo")
install.packages("extrafont")
##加载包

library(readr)
library(dplyr)
## 
## Attaching package: 'dplyr'
## The following objects are masked from 'package:stats':
## 
##     filter, lag
## The following objects are masked from 'package:base':
## 
##     intersect, setdiff, setequal, union
library(ggplot2)
library(stringr)
library(magrittr)
library(purrr)
## 
## Attaching package: 'purrr'
## The following object is masked from 'package:magrittr':
## 
##     set_names
library(tidyr)
## 
## Attaching package: 'tidyr'
## The following object is masked from 'package:magrittr':
## 
##     extract
library(tibble)
library(ggpubr)
library(RColorBrewer)
library(Cairo)
library(grid)

Sys.setenv(LANGUAGE = "en") 
options(stringsAsFactors = FALSE) #禁止chr转成factor
##
##输入文件的预处理
java -mx1024M -jar stem.jar -b bulk_RNA-seq_maize_endosperm.txt output
# gene expression
exp_all <- read.delim("bulk_RNA-seq_ave678dap.txt", check.names = F,
                       colClasses=c("SPOT"="character"))
                       
# time point
tp <- colnames(exp_all)[1:ncol(exp_all)]
# profile
profiletable <- read.delim("bulk_RNA-seq_profiletable.txt", check.names = F)
genetable <- read.delim("bulk_RNA-seq_genetable.txt", check.names = F)
colnames(genetable)[4:ncol(genetable)] <- as.character(seq(1:length(tp) - 1))
##Repeat data values will be averaged with the values from the original data file using the median.
colnames(exp_all)[2:ncol(exp_all)] <- as.character(seq(1:length(tp) - 1))
# 进入STEM profile的GeneSymbol
GeneSymbolp <- base::intersect(genetable$`Gene Symbol`, exp_all$`Gene Symbol`)
length(GeneSymbolp)
#只保留进入profile的Gene Symbol
exp_all <- exp_all[exp_all$`Gene Symbol` %in% GeneSymbolp,]
exp_all <- aggregate(.~`Gene Symbol`, exp_all, median) 
dim(exp_all)
#向表达矩阵添加基因所在的profile
exp_all.pro <-  genetable %>% select(`Gene Symbol`, Profile) %>% 
  right_join(exp_all, by = "Gene Symbol") 
  profiletable %<>% separate(col = `Profile Model`, 
                          into = as.character(seq(1:length(tp) - 1)), 
                          sep = ",", convert = TRUE) %>% 
  select(Profile = `Profile ID`, `Cluster (-1 non-significant)`, `# Genes Assigned`, `p-value`, as.character(seq(1:length(tp) - 1)))
myfun <- function(df) {
  df %<>% gather(key = "x", value = "y", as.character(seq(1:length(tp) - 1))) %>% 
    mutate(x1 = as.numeric(x)) %>% select(-x) %>% rename(x = x1)
  return(df)
}

sig_num <- profiletable %>% select(Profile, `Cluster (-1 non-significant)`, `# Genes Assigned`, `p-value`)

profiletable %<>% select(- `Cluster (-1 non-significant)`, -`# Genes Assigned`, -`p-value`) %>% 
  group_by(Profile) %>% nest() 
profiletable$red <- map(profiletable$data, myfun)
profiletable %<>% ##
  select(Profile, red) %>% 
  left_join(sig_num, by = "Profile")


genetable %<>% ##
  group_by(Profile) %>% nest()
genetable$grey <- map(genetable$data, myfun)
genetable %<>% select(Profile, grey) 

exp_all.pro %<>% ##
  group_by(Profile) %>% nest()
exp_all.pro$box <- map(exp_all.pro$data, myfun)
exp_all.pro %<>%  select(Profile, box)
nest <- genetable %>% left_join(profiletable, by = "Profile") %>% 
  left_join(exp_all.pro, by = "Profile") %>%
  arrange(`Cluster (-1 non-significant)`) 
# 
Profile <- nest$Profile %>% as.character()
Profile_2 <- as.factor(Profile)
levels(Profile_2) <- Profile
nest$Profile <- Profile_2 %>% sort()
###
plot_line <- 
  ggplot(data = nest_part %>% select(Profile, grey, `# Genes Assigned`) %>% unnest()) + 
  geom_line(aes(x = x, y = y, group = `Gene Symbol`), color = "grey") + 
  geom_line(data = nest_part %>% select(Profile, red, `# Genes Assigned`) %>% unnest(), 
            aes(x = x, y = y), color = "red") + 
  facet_grid(rows = vars(Profile)) + 
  scale_x_continuous(breaks = seq(1:length(tp) - 1), 
                     labels = tp) + 

  theme(panel.background = element_rect(fill = "white", color = "black"),
        axis.title = element_blank(), 
        axis.text.x = element_text(angle = 45, hjust = .5, vjust = .5), 
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.line = element_blank(),
        panel.border = element_rect(linetype = "solid", color = "black", fill = NA)) 
##
x2 <- as.character(nest_part_new$x) %>% as.factor()
levels(x2) <- seq(1:length(tp) - 1)
nest_part_new$x2 <- x2

plot_box <- 
  nest_part_new %>% 
  ggplot(aes(x = x2, y = y, fill = x2, group = x2)) + 
  geom_boxplot(size = 0.25, 
               outlier.color = NA, 
               
               #outlier.size = 1, outlier.color = "grey", 
               show.legend = FALSE) + 
  facet_grid(rows = vars(Profile)) + 
  scale_x_discrete(breaks = seq(1:length(tp) - 1), 
                     labels = tp) + 
  scale_fill_brewer(palette = "Set2") + 
  theme(panel.background = element_rect(fill = "white", color = "black"),
        axis.title = element_blank(), 
        axis.text.x = element_text(angle = 45, hjust = .5, vjust = .5),
        axis.line.x.bottom = element_line(color = "black"),
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.line = element_blank(),
        panel.border = element_rect(linetype = "solid", color = "black", fill = NA))
#plot_box
# profile的编号和每个profile内的基因数
plot_text_1 <- 
  nest_part %>% select(Profile, `# Genes Assigned`) %>% unnest() %>% 
  add_column(., x = rep(1, nrow(.))) %>% 
  add_column(., y = rep(1, nrow(.))) %>% 
  ggplot() + 
  geom_text(mapping = aes(x = x, y = x, label = paste0("U", Profile, "\n(", `# Genes Assigned`, ")"))) + 
  facet_grid(rows = vars(Profile)) + 
  theme(panel.background = element_blank(),
        axis.title = element_blank(), 
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank()) 
#plot_text_1

# 折线图的y轴title
plot_text_2 <- nest_part %>% select(Profile, `# Genes Assigned`) %>% unnest() %>% 
  add_column(., x = rep(1, nrow(.))) %>% 
  add_column(., y = rep(1, nrow(.))) %>% 
  ggplot() + 
  geom_text(mapping = aes(x = x, y = x), 
            size = 3, 
            label = expression('Log'[2]*'FC'), angle = 90) + 
  facet_grid(rows = vars(Profile)) + 
  theme(panel.background = element_blank(),
        axis.title = element_blank(), 
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank()) 
#plot_text_2

# box plot的y轴title
plot_text_3 <- nest_part %>% select(Profile, `# Genes Assigned`) %>% unnest() %>% 
  add_column(., x = rep(1, nrow(.))) %>% 
  add_column(., y = rep(1, nrow(.))) %>% 
  ggplot() + 
  geom_text(mapping = aes(x = x, y = x), 
            size = 3,
            label = expression('Log'[2]*'(CPM)'), angle = 90) + 
  facet_grid(rows = vars(Profile)) + 
  theme(panel.background = element_blank(),
        axis.title = element_blank(), 
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank()) 
#plot_text_3
##拼图并输出
CairoPDF(file = "STEMbox.pdf", width = 6, height = 14)

grid.newpage() 
layout_1 <- grid.layout(nrow = 1, ncol = 6, 
                        widths = c(0.3, 0.2, 1, 0.2, 1, 0.5))
pushViewport(viewport(layout = layout_1)) 
print(plot_text_1, vp = viewport(layout.pos.col = 1))
print(plot_text_2, vp = viewport(layout.pos.col = 2))
print(plot_line, vp = viewport(layout.pos.col = 3))
print(plot_text_3, vp = viewport(layout.pos.col = 4))
print(plot_box, vp = viewport(layout.pos.col = 5))
dev.off()
###############################################################################################################################
