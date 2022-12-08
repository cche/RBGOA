
# GO_MWU uses continuous measure of significance (such as fold-change or
# -log(p-value) ) to identify GO categories that are significantly enriches
# with either up- or down-regulated genes. The advantage - no need to impose
# arbitrary significance cutoff.

# If the measure is binary (0 or 1) the script will perform a typical "GO
# enrichment" analysis based Fisher's exact test: it will show GO categories
# over-represented among the genes that have 1 as their measure.

# On the plot, different fonts are used to indicate significance and color
# indicates enrichment with either up (red) or down (blue) regulated genes.
# No colors are shown for binary measure analysis.

# The tree on the plot is hierarchical clustering of GO categories based on
# shared genes. Categories with no branch length between them are subsets of each other.

# The fraction next to GO category name indicates the fracton of "good" genes
# in it; "good" genes being the ones exceeding the arbitrary absValue cutoff
# (option in gomwuPlot). For Fisher's based test, specify absValue = 0.5.
# This value does not affect statistics and is used for plotting only.

# Stretch the plot manually to match tree to text

# Mikhail V. Matz, UT Austin, February 2015; matz@utexas.edu


# setup R error handling to go to stderr
options(show.error.messages = FALSE, error = function() {
    cat(geterrmessage(), file = stderr())
    q("no", 1, FALSE)
})

# we need that to not crash galaxy with an UTF8 error on German LC settings.
loc <- Sys.setlocale("LC_MESSAGES", "en_US.UTF-8")

library("getopt")
library("tools")
library(argparser, quietly = TRUE)
options(stringAsFactors = FALSE, useFancyQuotes = FALSE)

getScriptPath <- function() {
    cmd.args <- commandArgs()
    m <- regexpr("(?<=^--file=).+", cmd.args, perl = TRUE)
    script.dir <- dirname(regmatches(cmd.args, m))
    if (length(script.dir) == 0) stop("can't determine script dir: please call the script with Rscript")
    if (length(script.dir) > 1) stop("can't determine script dir: more than one '--file' argument detected")
    return(script.dir)
}

# Create a parser
p <- arg_parser("Perform a GO_MWU analysis")

# Add command line arguments
p <- add_argument(p, "--input", short = "-i", help = "Gene with associated value", type = "character")
p <- add_argument(p, "--goAnnotation", short = "-a", help = "goAnnotation file", type = "character")
p <- add_argument(p, "--goDatabase", short = "-g", help = "goDatabase (go.obo) file", type = "character")
p <- add_argument(p, "--goDivision", short = "-d", help = "goDivision", type = "character")
p <- add_argument(p, "--largest", short = "-l", help = "Keep GO terms that contain less than this fraction of genes",
                  type = "numeric", default = 0.1)
p <- add_argument(p, "--smallest", short = "-s", help = "Smallest number of genes in GO category to considered",
                  type = "numeric", default = 5)
p <- add_argument(p, "--absValue", short = "-b", help = "absValue", type = "numeric")
p <- add_argument(p, "--clusterheight", short = "-c", help = "Cluster height", type = "numeric", default = 0.25)
p <- add_argument(p, "--textsize", short = "-e", help = "Text size in Graph", type = "numeric", default = 1.0)
p <- add_argument(p, "--pcut", short = "-p", help = "pcut", type = "numeric", default = 1e-2)
p <- add_argument(p, "--hcut", short = "-t", help = "hcut", type = "numeric", default = 0.9)
p <- add_argument(p, "--l1", help = "First level for graph", type = "numeric", default = 0.1)
p <- add_argument(p, "--l2", help = "Second level for graph", type = "numeric", default = 0.05)
p <- add_argument(p, "--l3", help = "Third level for graph", type = "numeric", default = 0.01)
p <- add_argument(p, "--outGraph", short = "-o", help = "Graph output filename",
                  type = "numeric", default = "Rplots.png")
p <- add_argument(p, "--version", short = "-v", help = "Version", type = "logical", flag = TRUE)

# Parse the command line arguments
argv <- parse_args(p)

# if help was asked then print a friendly message
# and exit with a non-zero error code
if (argv$help) {
    print(p)
    q(status = 1)
}

if (argv$version) {
    cat("0.3.0\n")
    q(status = 1)
}

scriptdir <- getScriptPath()

### if running in Rstudio ###
# Uncomment this block, change the values and run from here.
# argv = list()
# argv$input = "heats.csv"
# argv$goAnnotation = "amil_defog_iso2go.tab"
# argv$goDatabase = "go.obo"
# argv$goDivision = "BP"

### optionally uncomment the options below that you want to change and change the default value.
# argv$clusterheight <- 0.25
# argv$textsize <- 1.2
# argv$largest <- 0.1
# argv$smallest <- 5
# argv$absValue <- -log(0.05,10)
# argv$l1 <- 0.1
# argv$l2 <- 0.05
# argv$l3 <- 0.01
# argv$pcut <- 1e-2
# argv$hcut <- 0.9


source_local <- function(fname) {
    # argv <- commandArgs(trailingOnly = FALSE)
    base_dir <- scriptdir
    source(paste(base_dir, fname, sep = "/"))
}

source_local("gomwu.functions.R")


nn <- strsplit(argv$input, "[/.]")
if (length(nn[[1]]) == 3) {
    dir <- nn[[1]][1]
    name <- nn[[1]][2]
    ext <- nn[[1]][3]
} else if (length(nn[[1]]) == 2) {
    dir <- "."
    name <- nn[[1]][1]
    ext <- nn[[1]][2]
}
# It might take a few minutes for MF and BP. Do not rerun it if you just want
# to replot the data with different cutoffs, go straight to gomwuPlot. If you
# change any of the numeric values below, delete the files that were generated
# in previos runs first.

gomwuStats(argv$input, argv$goDatabase, argv$goAnnotation, argv$goDivision, scriptdir,
    # replace with full path to perl executable if it is not in your system's PATH already
    perlPath = "perl",
    # a GO category will not be considered if it contains more than this fraction of the total number of genes
    largest = argv$largest,
    smallest = argv$smallest, # a GO category should contain at least this many genes to be considered.
    clusterCutHeight = argv$clusterheight, # threshold for merging similar (gene-sharing) terms. See README for details.
    # Alternative = "g" # by default the MWU test is two-tailed;
    # specify "g" or "l" of you want to test for "greater" or "less" instead.
    # Module = TRUE,Alternative="g" # un-remark this if you are analyzing a
    # SIGNED WGCNA module (values: 0 for not in module genes, kME for in-module genes).
    # In the call to gomwuPlot below, specify absValue=0.001 (count number of "good genes" that fall into the module)
    # Module = TRUE # un-remark this if you are analyzing an UNSIGNED WGCNA module
)
# do not continue if the printout shows that no GO terms pass 10% FDR.


# ----------- Plotting results

# change this to a pdf output
png(paste0(dir, "/", argv$outGraph), width = 6, height = 8, units = "in", res = 300)
# png(paste0(dir,"/","Rplots.png"),res=100)

results <- gomwuPlot(argv$input, argv$goAnnotation, argv$goDivision,
    # genes with the measure value exceeding this will be counted as "good genes".
    # This setting is for signed log-pvalues. Specify absValue=0.001 if you are doing
    # Fisher's exact test for standard GO enrichment or analyzing a WGCNA module (all non-zero genes = "good genes").
    absValue = argv$absValue,
    # 	absValue = 1, # un-remark this if you are using log2-fold changes
    # FDR threshold for plotting. Specify level1=1 to plot all GO categories containing genes exceeding the absValue.
    level1 = argv$l1,
    level2 = argv$l2, # FDR cutoff to print in regular (not italic) font.
    level3 = argv$l3, # FDR cutoff to print in large bold font.
    # decrease to fit more on one page, or increase (after rescaling the plot so the tree fits the text)
    # for better "word cloud" effect
    txtsize = argv$textsize,
    treeHeight = 0.5, # height of the hierarchical clustering tree
    # 	colors = c("dodgerblue2","firebrick1","skyblue2","lightcoral")
    # these are default colors, uncomment and change if needed
)
dev.off()

# text representation of results, with actual adjusted p-values
write.table(results[[1]], paste0(dir, "/", "results.tsv"), sep = "\t")


# ------- extracting representative GOs

# this module chooses GO terms that best represent *independent* groups of significant GO terms

pcut <- argv$pcut
hcut <- argv$hcut

# plotting the GO tree with the cut level (uncomment the next two lines to plot)
# plot(results[[2]],cex = 0.6)
# abline(h = hcut,col="red")

# cutting
ct <- cutree(results[[2]], h = hcut)
annots <- c()
ci <- 1
for (ci in unique(ct)) {
    message(ci)
    rn <- names(ct)[ct == ci]
    obs <- grep("obsolete", rn)
    if (length(obs) > 0) {
        rn <- rn[-obs]
    }
    if (length(rn) == 0) {
        next
    }
    rr <- results[[1]][rn, ]
    bestrr <- rr[which(rr$pval == min(rr$pval)), ]
    best <- 1
    if (nrow(bestrr) > 1) {
        nns <- sub(" .+", "", row.names(bestrr))
        fr <- c()
        for (i in 1:length(nns)) {
            fr <- c(fr, eval(parse(text = nns[i])))
        }
        best <- which(fr == max(fr))
    }
    if (bestrr$pval[best] <= pcut) {
        annots <- c(annots, sub("\\d+\\/\\d+ ", "", row.names(bestrr)[best]))
    }
}

mwus <- read.table(paste0(dir, "/", paste("MWU", argv$goDivision, name, sep = "_"), ".", ext), header = TRUE)
best_GOs <- mwus[mwus$name %in% annots, ]
write.table(best_GOs, paste0(dir, "/", "best_go.tsv"), sep = "\t", row.names = FALSE)
