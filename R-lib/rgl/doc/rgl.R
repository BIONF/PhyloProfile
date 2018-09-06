## ----setup, echo=FALSE, results="asis"-----------------------------------
source("setup.R")
knitr::opts_chunk$set(rgl.newwindow = TRUE)
set.seed(123)

## ----echo=1--------------------------------------------------------------
with(iris, plot3d(Sepal.Length, Sepal.Width, Petal.Length, 
                  type="s", col=as.numeric(Species)))
rglwidget()

## ----persp3d, webgl=TRUE, fig.height=3, fig.width=6----------------------
library(MASS)
# from the fitdistr example
set.seed(123)
x <- rgamma(100, shape = 5, rate = 0.1)
fit <- fitdistr(x, dgamma, list(shape = 1, rate = 0.1), lower = 0.001)
loglik <- function(shape, rate) sum(dgamma(x, shape=shape, rate=rate, 
                                           log=TRUE))
loglik <- Vectorize(loglik)
xlim <- fit$estimate[1]+4*fit$sd[1]*c(-1,1)
ylim <- fit$estimate[2]+4*fit$sd[2]*c(-1,1)

mfrow3d(1, 2, sharedMouse = TRUE)
persp3d(loglik, 
        xlim = xlim, ylim = ylim,
        n = 30)
zlim <- fit$loglik + c(-qchisq(0.99, 2)/2, 0)
next3d()
persp3d(loglik, 
        xlim = xlim, ylim = ylim, zlim = zlim,
        n = 30)

## ----webgl=TRUE, fig.width=3, fig.height=3-------------------------------
triangles3d(cbind(x=rnorm(9), y=rnorm(9), z=rnorm(9)), col = "green")
decorate3d()
bg3d("lightgray")
aspect3d(1,1,1)

## ----results="hide",webgl=TRUE-------------------------------------------
open3d()
cols <- rainbow(7)
layout3d(matrix(1:16, 4,4), heights=c(1,3,1,3))
text3d(0,0,0,"tetrahedron3d"); next3d()
shade3d(tetrahedron3d(col=cols[1])); next3d()
text3d(0,0,0,"cube3d"); next3d()
shade3d(cube3d(col=cols[2])); next3d()
text3d(0,0,0,"octahedron3d"); next3d()
shade3d(octahedron3d(col=cols[3])); next3d()
text3d(0,0,0,"dodecahedron3d"); next3d()
shade3d(dodecahedron3d(col=cols[4])); next3d()
text3d(0,0,0,"icosahedron3d"); next3d()
shade3d(icosahedron3d(col=cols[5])); next3d()
text3d(0,0,0,"cuboctahedron3d"); next3d()
shade3d(cuboctahedron3d(col=cols[6])); next3d()
text3d(0,0,0,"oh3d"); next3d()
shade3d(oh3d(col=cols[7]))

## ------------------------------------------------------------------------
par3d("mouseMode")

## ----echo = 2:3----------------------------------------------------------
rgl.close()
persp3d(volcano, col = "green")
rglwidget()

## ------------------------------------------------------------------------
lattice::wireframe(volcano, col = "green", 
		   screen = rglToLattice())
angles <- rglToBase()
persp(volcano, col = "green", shade = TRUE,
      theta = angles$theta, phi = angles$phi)

## ----echo=FALSE----------------------------------------------------------
setdiff(ls("package:rgl"), documentedfns)

## ----echo=FALSE, results="asis"------------------------------------------
writeIndex(cols = 5)

