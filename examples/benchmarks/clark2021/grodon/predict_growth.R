suppressMessages({library(gRodon); library(Biostrings)})
files <- sort(list.files("/home/cumbof/muODE/data/clark2021/cds", pattern="\\.ffn$", full.names=TRUE))
cat("code\tn_genes\tn_ribo\tdoubling_h\tmumax_per_h\n")
for (f in files) {
  code <- sub("\\.ffn$","",basename(f))
  genes <- readDNAStringSet(f)
  he <- grepl("ribosomal_protein", names(genes))
  pred <- tryCatch(predictGrowth(genes, highly_expressed=he, mode="full"),
                   error=function(e){cat("#ERR",code,conditionMessage(e),"\n"); NULL})
  d <- if (is.null(pred)) NA else pred$d
  cat(sprintf("%s\t%d\t%d\t%.4f\t%.4f\n", code, length(genes), sum(he), d, log(2)/d))
}
