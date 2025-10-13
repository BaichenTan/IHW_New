## ----setup, echo = FALSE------------------------------------------------------
knitr::opts_chunk$set(tidy = FALSE, cache = TRUE, autodep = TRUE, dev = c("png", "pdf"), dpi = 360)

## ----loadlibs, message = FALSE, warning = FALSE-------------------------------
library("DESeq2")
library("dplyr")
data("airway", package = "airway")
dds <- DESeqDataSet(se = airway, design = ~ cell + dex) %>% DESeq
deRes <- as.data.frame(results(dds))

## ----colnames_deRes-----------------------------------------------------------
colnames(deRes)

## ----loadihw, message = FALSE, warning = FALSE--------------------------------
library("IHW")
ihwRes <- ihw(pvalue ~ baseMean,  data = deRes, alpha = 0.1)

## ----rejections---------------------------------------------------------------
rejections(ihwRes)

## ----headadjp-----------------------------------------------------------------
head(adj_pvalues(ihwRes))
sum(adj_pvalues(ihwRes) <= 0.1, na.rm = TRUE) == rejections(ihwRes)

## ----compareBH1---------------------------------------------------------------
padjBH <- p.adjust(deRes$pvalue, method = "BH")
sum(padjBH <= 0.1, na.rm = TRUE)

## ----compareBH2, echo = FALSE-------------------------------------------------
stopifnot( sum(padjBH <= 0.1, na.rm = TRUE) < 0.9 * rejections(ihwRes) )

## ----headweights--------------------------------------------------------------
head(weights(ihwRes))

## ----nbinsnfolds, echo = FALSE------------------------------------------------
## somehow inlining this code caching does not work when code is 
nbins_ihwRes <- nbins(ihwRes)
nfolds_ihwRes <- nfolds(ihwRes)

## -----------------------------------------------------------------------------
c(nbins(ihwRes), nfolds(ihwRes))

## -----------------------------------------------------------------------------
weights(ihwRes, levels_only = TRUE)

## ----plot_ihwRes_weights, fig.width = 6, fig.height = 3.6, out.width = "600px"----
plot(ihwRes)

## ----plot_ihwRes_decisionboundary, fig.width = 6, fig.height = 3.6, out.width = "600px"----
plot(ihwRes, what = "decisionboundary") 

## ----plot_ihwRes_rawvsadj, fig.width = 6, fig.height = 3.6, warning = FALSE, out.width = "600px", message = FALSE, warning = FALSE----
library("ggplot2")
gg <- ggplot(as.data.frame(ihwRes), aes(x = pvalue, y = adj_pvalue, col = group)) + 
  geom_point(size = 0.25) + scale_colour_hue(l = 70, c = 150, drop = FALSE)
gg
gg %+% subset(as.data.frame(ihwRes), adj_pvalue <= 0.2)

## -----------------------------------------------------------------------------
ihwResDf <- as.data.frame(ihwRes)
colnames(ihwResDf)

## ----Bonferroni---------------------------------------------------------------
ihwBonferroni <- ihw(pvalue ~ baseMean, data = deRes, alpha = 0.1, adjustment_type = "bonferroni")

## ----rkpvslogp, fig.width = 6, fig.height = 3, out.width = "800px"------------
deRes <- na.omit(deRes)
deRes$geneid <- as.numeric(gsub("ENSG[+]*", "", rownames(deRes)))

# set up data frame for ggplotting
rbind(data.frame(pvalue = deRes$pvalue, covariate = rank(deRes$baseMean)/nrow(deRes), 
                 covariate_type="base mean"),
      data.frame(pvalue = deRes$pvalue, covariate = rank(deRes$geneid)/nrow(deRes), 
                 covariate_type="gene id")) %>%
ggplot(aes(x = covariate, y = -log10(pvalue))) + geom_hex(bins = 100) + 
   facet_grid( . ~ covariate_type) + ylab(expression(-log[10]~p))

## ----pvalhistogram, fig.width = 3, fig.height = 3, out.width = "350px"--------
ggplot(deRes, aes(x = pvalue)) + geom_histogram(binwidth = 0.025, boundary = 0)

## ----goodhist, fig.width = 8, fig.height = 5, out.width = "1000px"------------
deRes$baseMeanGroup <- groups_by_filter(deRes$baseMean, 8)

ggplot(deRes, aes(x=pvalue)) + 
  geom_histogram(binwidth = 0.025, boundary = 0) +
  facet_wrap( ~ baseMeanGroup, nrow = 2)

## ----goodecdf, fig.width = 5, fig.height = 3, out.width = "600px"-------------
ggplot(deRes, aes(x = pvalue, col = baseMeanGroup)) + stat_ecdf(geom = "step") 

## ----badhist, fig.width = 8, fig.height = 5, out.width = "1000px"-------------
deRes$lfcGroup <- groups_by_filter(abs(deRes$log2FoldChange),8)

ggplot(deRes, aes(x = pvalue)) + 
  geom_histogram(binwidth = 0.025, boundary = 0) +
  facet_wrap( ~ lfcGroup, nrow=2)

## ----badecdf, fig.width = 5, fig.height = 3, out.width = "600px"--------------
ggplot(deRes, aes(x = pvalue, col = lfcGroup)) + stat_ecdf(geom = "step") 

## ----sim----------------------------------------------------------------------
sim <- tibble(
  X = runif(100000, min = 0, max = 2.5),         # covariate
  H = rbinom(length(X), size = 1, prob = 0.1),   # hypothesis true or false
  Z = rnorm(length(X), mean = H * X),            # Z-score
  p = 1 - pnorm(Z))                              # pvalue

## ----simBH--------------------------------------------------------------------
sim <- mutate(sim, padj = p.adjust(p, method="BH"))
sum(sim$padj <= 0.1)

## ----simSub-------------------------------------------------------------------
reporting_threshold <- 0.1
sim <- mutate(sim, reported = (p <= reporting_threshold))
simSub <- filter(sim, reported)

## ----compareselected1, fig.width = 3, fig.height = 3, out.width = "300px"-----
simSub = mutate(simSub, padj2 = p.adjust(p, method = "BH", n = nrow(sim)))
ggplot(simSub, aes(x = padj, y = padj2)) + geom_point(cex = 0.2)

## ----justcheck----------------------------------------------------------------
stopifnot(with(simSub, max(abs(padj - padj2)) <= 0.001))

## ----m_groups-----------------------------------------------------------------
nbins <- 20
sim <- mutate(sim, group = groups_by_filter(X, nbins))
m_groups  <- table(sim$group)

## ----ihw2---------------------------------------------------------------------
ihwS <- ihw(p ~ group, alpha = 0.1, data = filter(sim, reported), m_groups = m_groups)
ihwS

## ----ihw12--------------------------------------------------------------------
ihwF <- ihw(p ~ group, alpha = 0.1, data = sim)

## ----compareihw1, fig.width = 8, fig.height = 3, out.width = "700px"----------
gridExtra::grid.arrange( 
  plot(ihwS), 
  plot(ihwF), 
  ncol = 2) 
c(rejections(ihwS),
  rejections(ihwF))

## ----compareihw2, fig.width = 3, fig.height = 3, out.width = "300px", warning = FALSE----
qplot(adj_pvalues(ihwF)[sim$reported], adj_pvalues(ihwS), cex = I(0.2), 
      xlim = c(0, 0.1), ylim = c(0, 0.1)) + coord_fixed()

