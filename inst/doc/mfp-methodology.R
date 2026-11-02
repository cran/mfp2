## ----setup, include=FALSE-------------------------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  message = FALSE,
  warning = FALSE,
  fig.width = 7,
  fig.height = 4.5,
  out.width = "100%"
)

oldwidth <- getOption("width")
options(width = 100)

has_ggplot2 <- requireNamespace("ggplot2", quietly = TRUE)
has_geepack <- requireNamespace("geepack", quietly = TRUE)
if (has_ggplot2) library(ggplot2)
has_patchwork <- requireNamespace("patchwork", quietly = TRUE)

## ----eval=has_ggplot2, echo=FALSE-----------------------------------------------------------------
#==============================================================================
# 8 FP1 FUNCTIONS 
#===============================================================================
x <- seq(0.05, 1.05, length.out = 1000)
funx <- function(x, power){
  ifelse(power == rep(0, length(x)), log(x), x^power)
}
s <- c(-2, -1, -0.5, 0, 0.5, 1, 2, 3)
# transform x using the powers in s
outx <- lapply(s, function(s) funx(x = x, power = s))

# multiply the first 3 with -1 so that the function increase rather than decrease
# due to negative powers. these are betas i,e y = beta*x^p
k <- c(-1, -1, -1, 1, 1, 1, 1, 1)
datx <- matrix(unlist(outx), ncol = length(s))
datax <- datx %*% diag(k)
# standardize the data so that y is in the same range
dataxx <- apply(datax, 2, function(x) (x - min(x)) / (max(x) - min(x)))
colnames(dataxx) <- c(paste0("x", 1:length(s)))
dataxx <- as.data.frame(dataxx)
dataxx$x <- x
width <- 0.5
fig1 <- ggplot(dataxx, aes(x)) + 
  geom_line(aes(y = x1, colour = "x1"), linewidth = width, color = "#006400") + 
  geom_line(aes(y = x2, colour = "x2"), linewidth = width, color ="#ff0000") + 
  geom_line(aes(y = x3, colour = "x3"), linewidth = width, color ="#ffd700") +
  geom_line(aes(y = x4, colour = "x4"), linewidth = width, color ="#ff00ff") +
  geom_line(aes(y = x5, colour = "x5"), linewidth = width, color ="#ffb6c1") +
  geom_line(aes(y = x6, colour = "x6"), linewidth = width, color ="#00ff00") +
  geom_line(aes(y = x7, colour = "x7"), linewidth = width, color ="#0000ff") +
  geom_line(aes(y = x8, colour = "x8"), linewidth = width, color ="#000000") +
  geom_text(aes(x = 0.175, y = 1,size = 5, label = "-2"), color = "#006400") +
  geom_text(aes(x = 0.25, y = 0.9,size = 5, label = "-1"), color = "#ff0000") +
  geom_text(aes(x = 0.33, y = 0.83,size = 5, label = "-0.5"), color = "#ffd700") +
  geom_text(aes(x = 0.4, y = 0.725,size = 5, label = "0"), color = "#ff00ff") +
  geom_text(aes(x = 0.5, y = 0.65,size = 5, label = "0.5"), color = "#ffb6c1") +
  geom_text(aes(x = 0.59, y = 0.575,size = 5, label = "1"), color = "#00ff00") +
  geom_text(aes(x = 0.7, y = 0.475,size = 5, label = "2"), color = "#0000ff") +
  geom_text(aes(x = 0.75, y = 0.4,size = 5, label = "3"), color = "#000000") +
  labs(x="x", y="Fractional polynomial f(x)") +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    axis.title=element_text(size=18),
    axis.text.x = element_blank()) 

# Define functions
f1 <- function(x) 3 - 10*x^2 + 4*x^3
f2 <- function(x) 20 - 15.4*x^2 + 4*x^3
f3 <- function(x) -20 + 6*log(x) + 6*log(x)*log(x)
f4 <- function(x) 20 + 0.3*(x^-2) - 4*(x^-1)
f5 <- function(x) -10 + 5*(x^0.5) + 14*(x^-0.5)
f6 <- function(x) 33 + 19*log(x) - 7*(x^2)
f7 <- function(x) -10 + 10*(x - 1.5) + 10*(x - 1.5)^2

# Create data frames for each function
x <- seq(0.1, 3, by = 0.01)
df1 <- data.frame(x = x, y = f1(x))
df2 <- data.frame(x = x, y = f2(x))
df3 <- data.frame(x = x, y = f3(x))
df4 <- data.frame(x = x, y = f4(x))
df5 <- data.frame(x = x, y = f5(x))
df6 <- data.frame(x = x, y = f6(x))
df7 <- data.frame(x = x, y = f7(x))

# Plot functions
fig2 <- ggplot() +
    geom_line(data = df1, aes(x = x, y = y), color = "#e41a1c") +
    geom_line(data = df2, aes(x = x, y = y), color = "#377eb8") +
    geom_line(data = df3, aes(x = x, y = y), color = "#4daf4a") +
    geom_line(data = df4, aes(x = x, y = y), color = "#984ea3") +
    geom_line(data = df5, aes(x = x, y = y), color = "#ff7f00") +
    geom_line(data = df6, aes(x = x, y = y), color = "#FF00FF") +
    geom_line(data = df7, aes(x = x, y = y), color = "#a65628") +
    scale_x_continuous(limits = c(0, 3)) +
    scale_y_continuous(expand = c(0, 0), limits = c(-25,40)) + theme_classic() +
    geom_text(aes(x = 0.75, y = 2.4,size = 5, label = "(2, 3)"), color = "#e41a1c") +
    geom_text(aes(x = 1.6, y =1.95,size = 5, label = "(2, 3)"), color = "#377eb8")+
    geom_text(aes(x = 0.7, y = -18,size = 5, label = "(0, 0)"), color = "#4daf4a") +
    geom_text(aes(x = 1.25, y = 20,size = 5, label = "(-2, -1)"), color = "#984ea3") +
    geom_text(aes(x = 1.3, y = 12,size = 5, label = "(-0.5, 0.5)"), color = "#ff7f00") +
    geom_text(aes(x = 1.25, y = 29.5,size = 5, label = "(0, 2)"), color = "#FF00FF") +
    geom_text(aes(x = 1, y = -9.5,size = 5, label = "(1, 2)"), color = "#a65628") +
    ylab(" ") +
    theme_bw() +
    theme(
      legend.position = "none",
      panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5),
      axis.title = element_text(size=18),
      axis.ticks = element_blank(),
      axis.text.y = element_blank(),
      axis.text.x = element_blank()
    )

if (has_patchwork) figure <- patchwork::wrap_plots(fig1, fig2, 
                               ncol = 2, nrow = 1, 
                               widths = 8, heights = 2)


## ----fig1, eval=has_patchwork && has_ggplot2, echo=FALSE, fig.cap="Illustration of the flexibility of the FP family. 8 FP1 functions (left panel) and a subset of FP2 functions (right panel). FPs are global functions. See section 1.3.1.4 for a proposed extension of local features", fig.height=4.5, fig.width=7, out.width="100%", fig.align="center"----
figure

## ----cleanup, include=FALSE---------------------------------------------------
options(width = oldwidth)

