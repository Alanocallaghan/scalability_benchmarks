library("ggplot2")
theme_set(theme_bw())
load_all("../../BASiCS")




b <- readRDS("outputs/basics_zeisel.rds")
bf <- readRDS("outputs/basics_zeisel_fixed.rds")





free <- BASiCS_showFit(b) + labs(x = bquote(log(mu[i])), y = bquote(log(delta[i])))
fixed <- BASiCS_showFit(bf, ml = ml) + labs(x = bquote(log(mu[i])), y = NULL)

ggsave(
  free + theme_bw(), 
  file = "/home/alan/Documents/github/latex-documents/figures/scalability_paper_figs/trend_free.pdf", 
  width = 3, 
  height = 3
)
ggsave(
  fixed + theme_bw(), 
  file = "/home/alan/Documents/github/latex-documents/figures/scalability_paper_figs/trend_fixed.pdf", 
  width = 3, 
  height = 3
)








de <- BASiCS_TestDE(b, bf, GroupLabel1 = "Adaptive", GroupLabel2 = "Fixed")


vars <- c("mu", "delta", "epsilon")
xl <- list(bquote(log2(mu[i]^Adapt / mu[i]^Fixed)),
    bquote(log2(delta[i]^Adapt / delta[i]^Fixed)),
    bquote(epsilon[i]^Adapt - epsilon[i]^Fixed)
)


for (i in 1:3) {
    g1 <- DEPlots_ResultDE(de@Results[[i]], Which = "GridPlot")
    ggsave(g1, 
      file = paste0(
        "/home/alan/Documents/github/latex-documents/figures/scalability_paper_figs/", 
        vars[[i]], 
        "_grid_plot.pdf"
      ),
      width = 3,
      height = 3
    )
    g2 <- DEPlots_ResultDE(de@Results[[i]], Which = "MAPlot")
    ggsave(g2, 
      file = paste0(
        "/home/alan/Documents/github/latex-documents/figures/scalability_paper_figs/", 
        vars[[i]], 
        "_ma_plot.pdf"
      ),
      width = 3,
      height = 3
    )

    g3 <- DEPlots_ResultDE(de@Results[[i]], Which = "VolcanoPlot") + labs(x = xl[[i]])

    ggsave(g3, 
      file = paste0(
        "/home/alan/Documents/github/latex-documents/figures/scalability_paper_figs/", 
        vars[[i]], 
        "_volcano_plot.pdf"
      ),
      width = 3,
      height = 3
    )
}



