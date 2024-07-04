library(tidyverse)
library(tidytidbits)
library(survivalAnalysis)
survival::lung %>%
  mutate(ecog=recode_factor(ph.ecog, `0`="0", `1`="1", `2`="2-3", `3`="2-3")) %>%
  analyse_survival(vars(time, status), by=ecog) ->
  result
kaplan_meier_plot(result,
                  break.time.by="breakByQuarterYear",
                  xlab=".OS.months",
                  legend.title="ECOG Status",
                  hazard.ratio=T,
                  risk.table=TRUE,
                  table.layout="clean",
                  ggtheme=ggplot2::theme_bw(10))
KMPlot <- kaplan_meier_plot(result,
                            break.time.by="breakByQuarterYear",
                            xlab=".OS.months",
                            legend.title="ECOG Status",
                            hazard.ratio=T,
                            risk.table=TRUE,
                            table.layout="clean",
                            ggtheme=ggplot2::theme_bw(10))
KMPlot
KMPlot
save_pdf(plot = KMPlot, folder="Plot", fileBaseName=".pdf", width=8, height=8)
save_pdf(plot = KMPlot, folder="Plot", fileBaseName="Survival.pdf", width=8, height=8)