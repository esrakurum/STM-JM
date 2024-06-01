

# # Load adjacency matrix 
load("49AdjMat.rdata")
dim(W) ##Check dimension of the adjacency matrix



# Number of regions obtained from the Adjacency matrix 
nregion <- nrow(W)


# Number of neighbors for each region is stored in D matrix
D <- rowSums(W) #49
D <- diag(D)

n <- nregion 

nt <- 20
gridPoints <- seq(0, 1, length.out = nt)

nknots <- 10

##Number of observations for each subject: 
nts <- table(x = dat$subj_id)

nt.subj <- c(0, as.vector(nts))

## Time-invariant data, one row per subject:
data.id <- dat[!duplicated(dat$subj_id), ]

Time.e <- data.id$event.times

## 15-point Gauss-Kronrod rule for the integrals within the survival model
wk <- c(0.0630920926299786, 0.140653259715526, 0.190350578064785, 
        0.209482141084728, 0.190350578064785, 0.140653259715526, 
        0.0630920926299786, 0.0229353220105292, 0.10479001032225, 
        0.169004726639268, 0.204432940075299, 0.204432940075299, 
        0.169004726639268, 0.10479001032225, 0.0229353220105292)
sk <- c(-0.949107912342758, -0.741531185599394, -0.405845151377397, 
        0, 0.405845151377397, 0.741531185599394, 0.949107912342758, 
        -0.991455371120813, -0.864864423359769, -0.586087235467691, 
        -0.207784955007898, 0.207784955007898, 0.586087235467691, 
        0.864864423359769, 0.991455371120813)

ordsk <- order(sk)
sk <- sk[ordsk]
wk <- wk[ordsk]

K <- length(sk)

P <- Time.e/2
st <- outer(P, sk + 1) 

############################################
# Thin plate P-spline components

## Define the matrix of fixed effects. 
V1 <- cbind(rep(1, nt), gridPoints) 
V1.h <- cbind(rep(1, length(Time.e)), Time.e) 
V1.s <- cbind(rep(1, length(c(t(st)))), c(t(st))) 

## Define knots
knots <- quantile(unique(Time.e), seq(0, 1, length = (nknots+2))[-c(1, (nknots+2))])

## The design matrix for random coefficients 
Z_K <- (abs(outer(gridPoints, knots,"-")))^3
OMEGA_all <- (abs(outer(knots, knots,"-")))^3 ## outer creates (knots-knots)^3
svd.OMEGA_all <- svd(OMEGA_all) ## we use singular-value decomposition to find omega_K^{-1/2}
sqrt.OMEGA_all <- t(svd.OMEGA_all$v %*% (t(svd.OMEGA_all$u)*sqrt(svd.OMEGA_all$d)))
Z <- t(solve(sqrt.OMEGA_all,t(Z_K)))

Z_K.h <- (abs(outer(Time.e, knots,"-")))^3 
OMEGA_all.h <- (abs(outer(knots, knots,"-")))^3
svd.OMEGA_all.h <- svd(OMEGA_all.h)
sqrt.OMEGA_all.h <- t(svd.OMEGA_all.h$v %*% (t(svd.OMEGA_all.h$u)*sqrt(svd.OMEGA_all.h$d)))
Z.h <- t(solve(sqrt.OMEGA_all.h,t(Z_K.h)))

Z_K.s <- (abs(outer(c(t(st)), knots,"-")))^3
OMEGA_all.s <- (abs(outer(knots, knots,"-")))^3
svd.OMEGA_all.s <- svd(OMEGA_all.s)
sqrt.OMEGA_all.s <- t(svd.OMEGA_all.s$v %*% (t(svd.OMEGA_all.s$u)*sqrt(svd.OMEGA_all.s$d)))
Z.s <- t(solve(sqrt.OMEGA_all.s,t(Z_K.s)))



############################################
# Data for the MCMC

Data <- list(K = K,
             event = data.id$delta_i, 
             T1 = data.id$event.times,
             y = dat$y, 
             x1 = data.id$x1,
             z1 = data.id$z1,
             zeros = rep(0, dim(data.id)[1]),
             wk = wk,
             V1 = V1,
             V1.s = V1.s,
             V1.h = V1.h,
             Z = Z,
             Z.s = Z.s,
             Z.h = Z.h,
             nknots = nknots,
             nt = nt,
             nt.subj = nt.subj,
             dim_data.id = dim(data.id)[1],
             n.regions = nregion,
             reg.id = data.id$reg_id,
             D = D, W = W, 
             mu0 = rep(0,2),
             taus.a = matrix(c(1.0E-6, 0, 0,1.0E-6), 2,2),
             identity.b = diag(1, nknots,nknots)
)

############################################

### Parameters to keep track of during the MCMC,

parms <- c("betaY0","betaY1","gammaY1", "h0","betaS1","gammaS1", "sigmasqb1", "sigmasqb2", "rho12")

