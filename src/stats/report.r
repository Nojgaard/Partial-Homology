library(ggplot2)
library(latex2exp)
library(extrafont)
library(gridExtra)
library(grid)
library(extrafont)
library(tikzDevice)

options(tikzMetricPackages = c("\\usepackage[utf8]{inputenc}",
    "\\usepackage[T1]{fontenc}", "\\usetikzlibrary{calc}",
    "\\usepackage{amssymb}",
    "\\usepackage{amsfonts}"
	))

data_path <- "../../data/stats/"

reldiff_in	  <- c("25",
				   "50",
				   "100"
				   )
tbl_rand <- read.table("../../data/stats/rand_rate.dat", header = TRUE)
tbl_rand$prop <- tbl_rand$prop * 100
tbl_rand$diff <- tbl_rand$diff * 100

# settings for zoom plot
zoomtheme <- theme(legend.position="none", axis.line=element_blank(),
			# axis.text.x=element_blank(),
            # axis.text.y=element_blank(),
			axis.ticks=element_blank(),
           axis.title.x=element_blank(),
			axis.title.y=element_blank(),
            # panel.grid.major = element_blank(), 
			panel.grid.minor = element_blank(), 
            panel.background = element_rect(color='black'),
            plot.margin = unit(c(0,0,-6,-6),"mm"))
# settings for full plot
fulltheme <- theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
            panel.background = element_blank(), 
            axis.text.x=element_blank(), axis.text.y=element_blank(),
            axis.ticks=element_blank(),
            axis.title.x=element_blank(), axis.title.y=element_blank())

# zoom bounding box
xlim <- c(20,80); ylim <- c(0,12.5)


vp <- viewport(width = 0.3, height = 0.3, x = 0.45,
			   y = 0.6, just = c("right", "bottom"))


for (nleaf in reldiff_in) {
	fn <- paste("sp", nleaf, "_reldiff", sep = "")
	file_path <- paste(data_path, fn, ".dat", sep = "")
	print(file_path)

	tbl_md <- read.table(file_path, header = TRUE)
	tbl_md$prop <- tbl_md$prop * 100
	tbl_md$diff <- tbl_md$diff * 100
	tbl_md$sl <- 100 - (tbl_md$dl + tbl_md$hl)
	tbl_md$dlhl<- with(tbl_md, ifelse(sl!=34, paste0("$",sl,"/",dl,"/", hl,"$"),"$\\frac{1}{3}/\\frac{1}{3}/\\frac{1}{3}$"))


	# title <- paste("reconstructing unknown relations (",nleaf," leaves)", sep="")
	title <- paste0("$\\mathcal{I}(\\mathfrak{s}/\\mathfrak{d}/\\mathfrak{t};p;",nleaf,")$")
	title_colour <- "$\\mathfrak{s}/\\mathfrak{d}/\\mathfrak{t}$"
	print(tbl_md)
	p.zoom <- ggplot() + 
		geom_point(data=tbl_md, aes(prop, diff,colour=factor(dlhl), shape=factor(dlhl))) + 
		geom_line(data=tbl_md, aes(prop, diff,colour=factor(dlhl), shape=factor(dlhl))) + 
		geom_line(data=tbl_rand, aes(prop, diff), linetype="dashed") +
	    coord_cartesian(xlim=xlim, ylim=ylim) + zoomtheme

	p <- ggplot() + 
		geom_point(data=tbl_md, aes(prop, diff,colour=factor(dlhl), shape=factor(dlhl))) + 
		geom_line(data=tbl_md, aes(prop, diff,colour=factor(dlhl))) + 
		geom_line(data=tbl_rand, aes(prop, diff), linetype="dashed") +
		# scale_shape("", guide = FALSE) +
		# scale_color_discrete("") +
		# scale_colour_manual("", values = c("chartreuse4", "chocolate", "slateblue4", "chocolate")) +
		# scale_shape_manual("", values=c(16, 0, 3, 17)) +
	ggtitle(title) +
	labs(x="expected \\% of unknown genepairs", y="rel. diff. of original and recovered relations (in \\%)", colour=title_colour, shape=title_colour)
	  # theme(legend.text = element_text(size=30, 
                                  # vjust=1, family="Fraktur"))
	# g <- ggplotGrob(p.zoom)
# print(g)
	x_offs <- 5 ; y_offs <- 80
	# full <- p + 
	# 	annotation_custom(grob = g, xmin = min(tbl_md$prop), xmax = min(tbl_md$prop) + x_offs, ymin = max(tbl_md$diff)- y_offs, ymax = tbl_md$diff)
	# pdf(paste0(fn, ".pdf"))
	tikz(file = paste0("out/", fn, ".tex"), width = 5, height = 5)
	print(p)
	print(p.zoom, vp = vp)
	dev.off()
	# ggsave(paste0(fn, ".pdf"))

}

tbl_md <- read.table("../../data/stats/sizes_reldiff.dat", header = TRUE)
tbl_md$prop <- tbl_md$prop * 100
tbl_md$diff <- tbl_md$diff * 100
title <- paste0("$\\mathcal{I}(\\frac{1}{3}/\\frac{1}{3}/\\frac{1}{3};p;|L|)$")
p <- ggplot() + geom_point(data=tbl_md, aes(prop, diff,colour=factor(size), shape=factor(size))) + geom_line(data=tbl_md, aes(prop, diff,colour=factor(size))) + 
	geom_line(data=tbl_rand, aes(prop, diff), linetype="dashed") +
	ggtitle(title) +
	labs(x="expected \\% of unknown genepairs", y="rel. diff. of original and recovered relations (in \\%)", colour="$|L|$", shape="$|L|$")

tikz(file = paste0("out/", "sizes_reldiff", ".tex"), width = 5, height = 5)
print(p)
dev.off()
# ggsave(paste0("sizes_reldiff", ".pdf"))


# tbl_md <- read.table("../../data/stats/matrix_diff.dat", header = TRUE)
# tbl_md <- read.table("../../data/stats/sp50_reldiff.dat", header = TRUE)
# tbl_rand <- read.table("../../data/stats/rand_rate.dat", header = TRUE)
# tbl_md$dlhl<- with(tbl_md, paste0(dl,"/", hl))
# tbl_md$prop <- tbl_md$prop * 100
# tbl_md$diff <- tbl_md$diff * 100
# tbl_rand$prop <- tbl_rand$prop * 100
# tbl_rand$diff <- tbl_rand$diff * 100
# # print(tbl_md)
# print(tbl_rand)
# ggplot() + geom_point(data=tbl_md, aes(prop, diff,colour=factor(dlhl))) + geom_line(data=tbl_md, aes(prop, diff,colour=factor(dlhl))) + 
# 	geom_line(data=tbl_rand, aes(prop, diff), linetype="dashed") +
# 	ggtitle("reconstructing unknown relations (50 leaves)") +
# 	labs(x="unknown probablity (in %)", y="relative difference between true relations and predicted (in %)", colour="S/D/L dist")

# tbl_md <- read.table("../../data/stats/sizes_reldiff.dat", header = TRUE)
# tbl_rand <- read.table("../../data/stats/rand_rate.dat", header = TRUE)
# tbl_md$dlhl<- with(tbl_md, paste0(dl,"/", hl))
# tbl_md$prop <- tbl_md$prop * 100
# tbl_md$diff <- tbl_md$diff * 100
# tbl_rand$prop <- tbl_rand$prop * 100
# tbl_rand$diff <- tbl_rand$diff * 100
# # print(tbl_md)
# print(tbl_rand)
# ggplot() + geom_point(data=tbl_md, aes(prop, diff,colour=factor(size))) + geom_line(data=tbl_md, aes(prop, diff,colour=factor(size))) + 
# 	geom_line(data=tbl_rand, aes(prop, diff), linetype="dashed") +
# 	ggtitle("reconstructing unknown relations with distribution 33/33/33") +
# 	labs(x="unknown probablity (in %)", y="relative difference between true relations and predicted (in %)", colour="sizes")
# ggplot(tbl_md, aes(prop, diff,colour=factor(rate))) + geom_point() + geom_line() + geom_line(tbl_rand, aes(prop, diff)) +
# 	ggtitle("reconstructing unknown relations") +
# 	labs(x="unknown probablity (in %)", y="relative difference between true relations and predicted (in %)", colour="S/D/L dist")

options(warn=1)

# zoom bounding box
xlim <- c(60,80); ylim <- c(0,12.5)

order_in <- list(
				list("50", 10, 10),
				list("50", 45, 10),
				list("50", 10, 45),
				list("50", 33, 33)
			)

for (row in order_in) {
	nleaf <- row[[1]]
	D <- row[[2]]
	H <- row[[3]]
	fn <- paste("sp", nleaf, "_D",D,"_H",H, "_orderdiff", sep = "")
	file_path <- paste(data_path, fn, ".dat", sep = "")
	print(file_path)
	tbl_md <- read.table(file_path, header = TRUE)
	tbl_md$prop <- tbl_md$prop * 100
	tbl_md$diff <- tbl_md$diff * 100
	dist_name <- paste0((100-(D+H)), "/", D, "/", H)
	if (D == 33) {
		dist_name <- "\\frac{1}{3}/\\frac{1}{3}/\\frac{1}{3}"
	}
	# title <- paste0("recover unknown relations (distributions ",
	# 				(100-(D+H)), "/", D, "/", H, ")")
	tbl_md$order <- factor(tbl_md$order, level = c("S/D/T", "S/T/D", "D/S/T", "D/T/S", "T/S/D", "T/D/S", "RAND"))
	title <- paste0("$\\mathcal{I}(",dist_name,";p;50)$")
	p <- ggplot(tbl_md, aes(prop, diff,colour=order, shape=order)) + geom_point() + geom_line() +
		ggtitle(title) +
		scale_shape_manual(values = c(16,16,15,15,17,17,7)) +
		labs(x="expected \\% of unknown genepairs", y="rel. diff. of original and recovered relations (in \\%)", colour="Rule Order", shape="Rule Order")

	p.zoom <- ggplot(tbl_md, aes(prop, diff,colour=order, shape=order)) + geom_point() + geom_line() +
		scale_shape_manual(values = c(16,16,15,15,17,17,7)) +
		coord_cartesian(xlim=xlim, ylim=ylim) + zoomtheme
	# print(p)
	# ggsave(paste0(fn, ".pdf"))
	tikz(file = paste0("out/", fn, ".tex"), width = 5, height = 5)
	print(p)
	print(p.zoom, vp = vp)
	dev.off()
}
forbid_in <- list(
				list("50", 10, 10),
				list("50", 45, 10),
				list("50", 10, 45),
				list("50", 33, 33)
			)

for (row in forbid_in) {
	nleaf <- row[[1]]
	D <- row[[2]]
	H <- row[[3]]
	fn <- paste("sp", nleaf, "_D",D,"_H",H, "_forbiddiff", sep = "")
	file_path <- paste(data_path, fn, ".dat", sep = "")
	print(file_path)
	tbl_md <- read.table(file_path, header = TRUE)
	tbl_md$prop <- tbl_md$prop * 100
	tbl_md$diff <- tbl_md$diff * 100
	dist_name <- paste0((100-(D+H)), "/", D, "/", H)
	if (D == 33) {
		dist_name <- "\\frac{1}{3}/\\frac{1}{3}/\\frac{1}{3}"
	}
	# title <- paste0("recover unknown relations with forbidden relations (70% unknowns, ",
	# 				(100-(D+H)), "/", D, "/", H, ")")
	title <- paste0("$\\mathcal{I}(",dist_name,";70;50)$")
	tbl_md$order <- factor(tbl_md$order, level = c("S/D/T", "S/T/D", "D/S/T", "D/T/S", "T/S/D", "T/D/S", "RAND"))
	p <- ggplot(tbl_md, aes(prop, diff,colour=order, shape=order)) + geom_point() + geom_line() +
		ggtitle(title) +
		labs(x="expected \\% of forbidden genepairs", y="rel. diff. of original and recovered relations (in \\%)", colour="Rule Order", shape="Rule Order") +
		scale_shape_manual(values = c(16,16,15,15,17,17,7)) +
		ylim(0,6)
	tikz(file = paste0("out/", fn, ".tex"), width = 5, height = 5)
	print(p)
	# ggsave(paste0(fn, ".pdf"))
	dev.off()
}

# tbl_md <- read.table("../../data/stats/sp25_D10_H10_orderdiff.dat", header = TRUE)
# tbl_md$prop <- tbl_md$prop * 100
# tbl_md$diff <- tbl_md$diff * 100
# ggplot(tbl_md, aes(prop, diff,colour=factor(order))) + geom_point() + geom_line() +
# 	ggtitle("reconstructing unknown relations (distribution 45/10/45)") +
# 	labs(x="unknown probablity (in %)", y="relative difference between true relations and predicted (in %)", colour="Rule Order")

# tbl_md <- read.table("../../data/stats/sp50_D33_H33_forbiddiff.dat", header = TRUE)
# tbl_md$prop <- tbl_md$prop * 100
# tbl_md$diff <- tbl_md$diff * 100
# ggplot(tbl_md, aes(prop, diff,colour=factor(order))) + geom_point() + geom_line() +
# 	ggtitle("reconstruction with forbidden relations (70% unknowns, 33/33/33)") +
# 	labs(x="forbid relation probablity (in %)", y="relative difference between true relations and predicted (in %)", colour="Rule Order")

# tbl_ud <- read.table("../../data/stats/unknown_diff.dat", header = TRUE)
# tbl_ud$dlhl<- with(tbl_ud, paste0(dl,"_", hl))
# ggplot(tbl_ud, aes(prop, diff,colour=factor(dlhl))) + geom_point() + geom_line() +
# 	ggtitle("correctly set unknowns/total unknowns vs. unknown probability") +
# 	labs(colour="dl_hl")

# tbl_trip <- read.table("out.dat", header = TRUE)
# ggplot(tbl_trip, aes(prop, triple_diff)) + geom_point() + geom_line() +
# 	ggtitle("unknown probability vs. triple difference")

