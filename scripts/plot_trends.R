

get_file <- function(by, nchains, seed) {
  paste0(
    "/home/alan/Documents/benchmarkdata/", 
    by, "/", 
    sprintf("%02d", nchains),
    "_chains_seed_", 
    seed, 
    ".rds")
}
c1 <- readRDS("outputs/basics_zeisel.rds")

g1 <- plot_fits(readRDS(get_file("gene", 2, 42)), c1, alpha = 0.6)
g2 <- plot_fits(readRDS(get_file("gene", 128, 42)), c1, alpha = 0.6)

ggsave(g1 + labs(x = bquote(log(mu[i])), y = bquote(log(delta[i]))),
  file = "/home/alan/Documents/github/latex-documents/figures/scalability_paper_figs/trend_2.pdf", 
  width = 3, height = 3)

ggsave(g2 + labs(x = bquote(log(mu[i])), y = bquote(log(delta[i]))),
  file = "/home/alan/Documents/github/latex-documents/figures/scalability_paper_figs/trend_128.pdf", 
  width = 3, height = 3)

