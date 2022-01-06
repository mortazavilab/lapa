library(readr)
library(dplyr)
library(ggplot2)

data <- read_csv(snakemake@input[[1]])
ymax = 42000

data$novelty <- factor(data$novelty,
  levels = rev(c("Known", "prefix ISM", "suffix ISM", 
                 "other ISM", "NIC", "NNC",
                 "Antisense", "Intergenic", "Genomic")))

data$support <- as.factor(data$support)
freqs <- data %>% count(support, novelty) %>%
  group_by(novelty) %>%
  mutate(freq = n / sum(n), total = sum(n))
freqs$novelty <- factor(freqs$novelty)

freqs$percent <- round(freqs$freq*100)
freqs[freqs$support == "no", "percent"] <- NA
freqs$tcolor_grp <- as.factor(ifelse(freqs$percent > 20, "white", "black"))

colors <- c("Known" = "#009E73","ISM" = "#0072B2", "NIC" = "#D55E00", 
            "NNC" = "#E69F00", "Antisense" = "#000000", "Intergenic" = "#CC79A7",
            "prefix ISM" = "#56B4E9", "suffix ISM" = "#698bac", "other ISM" = "#003366",
            "Genomic" = "#F0E442", "Novel" = "lightblue")

fname <- snakemake@output[[1]]
xlabel <- "Transcript category"
ylabel <- "Number of transcript models"

png(filename = fname,
    width = 2000, height = 1500, units = "px",
    bg = "white",  res = 300)

label_pad = max(freqs$total) * 0.1 + 1000
print(label_pad)
g = ggplot(freqs, aes(x = novelty, y = n, fill = novelty,
                      alpha = support)) +
  geom_bar(stat="identity", color = "black") +
  xlab(xlabel) + ylab(ylabel) +
  theme(legend.text = element_text(color="black", size = rel(1)),
        legend.title = element_text(color="black", size=rel(1))) +
  theme_bw(base_family = "Helvetica", base_size = 18) +
  scale_fill_manual("", values = colors) +
  scale_alpha_manual(values=c(0,1), name = "CAGE support") +
  theme_bw(base_family = "Helvetica", base_size = 18) +
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5),
        axis.text.x = element_text(color="black", size = rel(1.5)),
        axis.text.y = element_text(color="black", size = rel(1.5)),
        axis.title.x = element_text(color="black", size=rel(1.25)),
        axis.title.y = element_blank()) +
  coord_flip(ylim=c(0,ymax)) + guides(fill=FALSE, alpha = FALSE) +
  geom_text(aes(y = total + label_pad,
                label = paste0(percent, "%"), color = novelty),
            position = position_dodge(0.2), size = 8) +
  scale_color_manual(values = colors) +
  guides(colour=FALSE, fill=FALSE) +
  theme(legend.position=c(0.8,0.2),
        legend.title = element_blank(),
        legend.background = element_rect(fill="white", color = "black"),
        legend.key = element_rect(fill="transparent"),
        legend.text = element_text(colour = 'black', size = 16)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
print(g)
dev.off()
