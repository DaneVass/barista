# Connected barplot script
# Dane Vassiliadis
# 31-10-19

# Takes a dataframe of barcode counts and plots barcode composition per timepoint using stacked and connected barplots.

# inspired by Jean Fan
# https://jef.works/blog/2017/09/21/connected-barplot-for-series-data-visualization/

#dat <- as.matrix(KRAS.cpmnorm.ordered.14[c(1:100), c(1:3)])
#space = 0.5
#alpha = 0.5
#color=rainbow(nrow(dat))

# function call
connectedBarplot <- function(dat, color=rainbow(nrow(dat)), space=1, alpha=0.5, ...) {
    b <- barplot(dat, col=color, space = space)

    for (i in seq_len(ncol(dat) - 1)) {
        lines(c(b[i]+0.5, b[i+1]-0.5), c(0, 0)) ## creates bottom line of each sample

        for (j in seq_len(nrow(dat))) { ## creates barplot for each row element stacking one on top of the previous
            print(j)

            # why are separate calls to if required for 1, 2 or many rows?
            if (j == 1) {
                lines(c(b[i]+0.5, b[i+1]-0.5), c(dat[j,i], dat[j,i+1]), lwd = 1)
                polygon(c(b[i]+0.5, b[i]+0.5, b[i+1]-0.5, b[i+1]-0.5),
                        c(0, dat[j,i], dat[j,i+1], 0),
                        col=adjustcolor(color[j], alpha.f=alpha))
            }
            if (j == 2) {
                lines(c(b[i]+0.5, b[i+1]-0.5), c(colSums(dat[1:j,])[i], colSums(dat[1:j,])[i+1]), lwd = 1)
                polygon(c(b[i]+0.5, b[i]+0.5, b[i+1]-0.5, b[i+1]-0.5),
                        c(dat[1,i], colSums(dat[1:j,])[i], colSums(dat[1:j,])[i+1], dat[1,i+1]),
                        col=adjustcolor(color[j], alpha.f=alpha))
            }
            if (j > 2) {
                lines(c(b[i]+0.5, b[i+1]-0.5), c(colSums(dat[1:j,])[i], colSums(dat[1:j,])[i+1]), lwd = 0)
                polygon(c(b[i]+0.5, b[i]+0.5, b[i+1]-0.5, b[i+1]-0.5),
                        c(colSums(dat[1:(j-1),])[i], colSums(dat[1:j,])[i], colSums(dat[1:j,])[i+1], colSums(dat[1:(j-1),])[i+1]),
                        col=adjustcolor(color[j], alpha.f=alpha))
            }
        }
    }
}
