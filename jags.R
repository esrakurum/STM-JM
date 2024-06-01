model<- function ()
{
  for (i in 1:dim_data.id) {
    b[i, 1:2] ~ dmnorm(mu0[], inv.D[, ])
    
    ## Longitudinal Part
    for(j in (sum(nt.subj[1:i])+1): (sum(nt.subj[1:(i+1)])))
    { 
      logit(prob[j]) <- (inprod(a01[1:2], V1[j - sum(nt.subj[1:i]), 1:2]) + inprod(b01[1:nknots], Z[j - sum(nt.subj[1:i]),1:nknots])) + 
        x1[i]*(inprod(a11[1:2], V1[j - sum(nt.subj[1:i]), 1:2]) + inprod(b11[1:nknots], Z[j - sum(nt.subj[1:i]),1:nknots])) + 
        z1[i]*(inprod(az11[1:2], V1[j - sum(nt.subj[1:i]), 1:2]) + inprod(bz11[1:nknots], Z[j - sum(nt.subj[1:i]),1:nknots])) +
        b[i, 1]  + Q[reg.id[i], 1]
      
      y[j] ~ dbern(prob[j])
    }  
    ## Survival part
    log.h0.T[i] <- (inprod(a.h1[1:2], V1.h[i, 1:2]) + inprod(b.h1[1:nknots], Z.h[i, 1:nknots]))
    
    haz[i] <- exp(log.h0.T[i] +  x1[i]*(inprod(a12[1:2], V1.h[i, 1:2]) + inprod(b12[1:nknots], Z.h[i, 1:nknots])) 
                  + z1[i]*(inprod(az12[1:2], V1.h[i, 1:2]) + inprod(bz12[1:nknots], Z.h[i, 1:nknots]))
                  + b[i, 2] + Q[reg.id[i], 2])
    
    
    for(k in 1:K)
    {
      log.h0.s[i, k] <- inprod(a.h1[1:2], V1.s[K * (i - 1) + k, 1:2]) + inprod(b.h1[1:nknots], Z.s[K * (i - 1) + k,1:nknots])
      
      
      surv[i, k] <- exp(log.h0.s[i, k]  
                        + x1[i]*(inprod(a12[1:2], V1.s[K * (i - 1) + k, 1:2]) + inprod(b12[1:nknots], Z.s[K * (i - 1) + k, 1:nknots]))
                        + z1[i]*(inprod(az12[1:2], V1.s[K * (i - 1) + k, 1:2]) + inprod(bz12[1:nknots], Z.s[K * (i - 1) + k, 1:nknots]))
                        + b[i, 2] + Q[reg.id[i], 2])
      
    }
    log.survival[i] <- -T1[i]/2*inprod(wk[], surv[i ,])  
    phi[i] <- 100000 - (event[i] * log(haz[i])) - log.survival[i] 
    zeros[i] ~ dpois(phi[i]) 
  }
   

########
## Priors
########  
  

  for (i in 1:n.regions){
    for (j in 1:2){
      Q[i,j] <- inprod(tPsi[1:j,i], roDis[j, 1:j])
    }
  }


   for(i in 1:n.regions){ceros[i] <- 0}
  
  for (j in 1:2){
    tPsi[j, 1:n.regions] ~ dmnorm(ceros[], (D - nu[j] * W))
    nu[j] ~ dunif(0.1, 0.9)
  }
  
   roDis[1, 2] <- 0

  for (i in 1:2){
    roDis[i,i] ~ dlnorm(0.5, 0.005)
  }

    roDis[2, 1] ~ dnorm(0, 10)

   sigmasqb1 ~ dgamma (2, 0.05)
   sigmasqb2 ~ dgamma (2, 0.1)
   rho12 ~ dunif(0.1, 1)
   
   Sigma[1,1] <- sigmasqb1
   Sigma[2,2] <- sigmasqb2
   Sigma[1,2] <- sqrt(sigmasqb1)*sqrt(sigmasqb2)*rho12
   Sigma[2,1] <- Sigma[1,2]
   
   inv.D <- inverse(Sigma)
   
   b01[1:nknots] ~ dmnorm(rep(0,nknots), taub*identity.b) 
   b11[1:nknots] ~ dmnorm(rep(0,nknots), taub*identity.b) 
   b12[1:nknots] ~ dmnorm(rep(0,nknots), taub*identity.b) 
   b.h1[1:nknots]~ dmnorm(rep(0,nknots), taub*identity.b) 
   bz11[1:nknots] ~ dmnorm(rep(0,nknots), taub*identity.b) 
   bz12[1:nknots] ~ dmnorm(rep(0,nknots), taub*identity.b) 

  a01[1:2] ~ dmnorm(rep(0, 2), taus.a) 
  a11[1:2] ~ dmnorm(rep(0, 2), taus.a) 
  a12[1:2] ~ dmnorm(rep(0, 2), taus.a) 
  a.h1[1:2]~ dmnorm(rep(0, 2), taus.a) 
  az11[1:2] ~ dmnorm(rep(0, 2), taus.a) 
  az12[1:2] ~ dmnorm(rep(0, 2), taus.a) 
  
  
   
   taub ~ dgamma(1.0E-6,1.0E-6)
   sigmab <- 1/sqrt(taub)
   
   ## Obtain the regression coefficients using the fixed and random effect components in thin-plate spline expansion
   for (j in 1:nt){
     betaY0[j] <- inprod(a01[1:2], V1[j,1:2]) + inprod(b01[1:nknots], Z[j,1:nknots])
     betaY1[j] <- inprod(a11[1:2], V1[j,1:2]) + inprod(b11[1:nknots], Z[j,1:nknots])
     gammaY1[j] <- inprod(az11[1:2], V1[j,1:2]) + inprod(bz11[1:nknots], Z[j,1:nknots])
     
     h0[j] <- inprod(a.h1[1:2], V1[j, 1:2]) + inprod(b.h1[1:nknots], Z[j, 1:nknots])
     betaS1[j] <- inprod(a12[1:2], V1[j, 1:2]) + inprod(b12[1:nknots], Z[j, 1:nknots])
     gammaS1[j] <- inprod(az12[1:2], V1[j, 1:2]) + inprod(bz12[1:nknots], Z[j, 1:nknots])     
   } 
}


