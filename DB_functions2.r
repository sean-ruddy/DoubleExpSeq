DPNormC <- function( m, lambda, theta )
                {
                  mu = lambda*m
                  y = 0:ceiling(10*max(mu))
                  yly = y*log(y)
                  yly[is.na(yly)] = 0
                  lgam.y <- lgamma(y+1)
                  yly.theta <- yly*(1-theta)
                  y.theta <- theta*y
                  l.mu <- log(mu)
                  mu.theta <- mu*theta
                  num = matrix( 0 , nr=length(y) , nc=length(m) )
                  for( i in 1:length(m) )
                   {
                     num[,i] = -(mu.theta[i] + y) + yly.theta + y.theta*(l.mu[i] + 1) - lgam.y
                   }
                  normC <- colSums(exp(num))
                  return(normC)
                }

ddpois <- function( y , m, lambda , theta )
           {
             mu = lambda*m
             yly = y*log(y)
             yly[is.na(yly)] = 0
             lgam.y <- lgamma(y+1)
             yly.theta <- yly*(1-theta)
             y.theta <- theta*y
             num = matrix( 0 , nr=length(y) , nc=length(m) )
             l.mu <- log(mu)
             mu.theta <- mu*theta
             for( i in 1:length(m) )
              {
                num[,i] = -(mu.theta[i] + y) + yly.theta + y.theta*(l.mu[i] + 1) - lgam.y
              }
             num <- exp(num)
             p <- num/colSums(num)
             return(p)
           }

pdpois <- function( y , m, lambda , theta ) colSums( ddpois(0:y, m, lambda, theta) )

rdpois <- function( m, lambda , theta )
           {
             mu = lambda*m
             y = 0:ceiling(10*max(mu))
             system.time(probs <- ddpois( y , m , lambda , theta ))
             rn <- apply( probs , 2 , FUN=function(vec) sample(y,size=1,replace=TRUE,prob=vec) )
             return(rn)
           }
#rdpois <- function( m, lambda , theta ) #this function didn't perform well, but the way the code is written is important due to its quickness in the for-loop
#           {
#             unifs = runif(length(m))
#             mu = lambda*m
#             y = 0:ceiling(3*max(mu))
##             pdpois(2*max(mu),m,lambda,theta)
#             dens <- ddpois(y,m,lambda,theta)
#             P = matrix(1,nr=nrow(dens),nc=ncol(dens))
#             P[1,] = dens[1,]
#             sum = P[1,]
#             system.time(for(i in 2:nrow(dens))
#              {
#                sum = sum + dens[i,]
#                if( sum <= 1-1e-6 ) P[i,]=sum #simply putting an if statement cuts the running time in half, perhaps because the assignment is done in C and not R.
#                else break
#              })
#             options(warn=-1)
#             r = sapply( 1:length(unifs) , FUN=function(i) max(which(P[,i]<=unifs[i])))
#             r[r==-Inf] = 0
#             return(r)
#           }


DBNormCGetTerms <- function( n , cutoff = 1e6 , approx = "exact" )
                {
                   if( approx=="exact" )
                    {
                       kvals = 0:n
                       tmp1 = kvals*log(kvals)
                       tmp1[1] = 0
                       n_kvals = n-kvals
                       tmp2 = (n_kvals)*log((n_kvals))
                       tmp2[n+1] = 0
                       return(list(lgamma(n+1)-lgamma(kvals + 1)-lgamma(n_kvals + 1),(tmp1 + tmp2 - n*log(n)),kvals,(n_kvals)))
                       #lgamma(k) = lgamma(k+1) - log(k)
                    }
                   else return(NULL)
                }


#DBNormCOne <- function( n , p, phi, etap = NULL , psi = NULL , approx = c("exact","approx","1") , k = 0 , cutoff = 1e6 , log = FALSE , substitute = c("1","approx") , terms = NULL )
#                {
#                   if( any( c(length(n),length(p),length(phi)) > 1 ) ) stop( "This function requires scalar inputs for 'n', 'p' & 'phi'. For vector inputs use DBNormC()" )
#
#                   if( !is.null(etap) ) p <- exp(etap-log(1+exp(etap)))
#                   if( !is.null(psi) )
#                    {
#                      p <- exp(psi/phi-log(1+exp(psi/phi)))
#                      phi <- phi
#                    }
#                   if( n == 0 & log & k == 0 ) return( 0 )
#                    else if( n == 0 & k == 0 ) return( 1 )
#                     else if( n == 0 & k > 0 ) return( 0 )
#                      else if( n == 0 & log & k > 0 ) stop( "Must set log to FALSE if n == 0 and k > 0" )
#                   approx = match.arg( approx , c("exact","approx","1") )
#                   substitute = match.arg( substitute , c("1","approx") )
#                   if( n > cutoff & k>0 ) stop( "Must set cutoff higher if k>0" )
#                   if( k > 0 & approx!="exact" ) stop( "If k>0, 'approx' must be set to exact." )
#                   if( approx == "exact" & n > cutoff ) approx = substitute
#
#                   exFun <- function( n , p , phi , k = 0 , terms = NULL )
#                             {
#                               if( is.null( terms ) ) terms <- DBNormCGetTerms( n , cutoff = cutoff , approx = "exact" )
#                               ljexpk <- 0
#                               if( k > 0 ) ljexpk <- log(terms[[3]]^k)
#                               return( sum(exp(terms[[1]] + (1-phi)*terms[[2]] + terms[[3]]*phi*log(p) + terms[[4]]*phi*log(1-p) + ljexpk)) )
#                             }
#                   NormC <- switch( approx , "1" = phi^(-0.5) , "approx" = (( 1 + 1/(12*n) * (1-phi)/phi * ( 1 - 1/(p*(1-p))) ) * sqrt(phi))^(-1) ,
#                                             "exact" = exFun( n=n , p=p , phi=phi , k=k , terms=terms ) )
#
#                   return( exFun( n=n , p=p , phi=phi , k=k , terms=terms ))
#                   if( log ) return( log(NormC) )
#                    else return( NormC )
#                }

DBNormCOne <- function( n , p, phi, etap = NULL , psi = NULL , approx = c("exact","approx","1") , k = 0 , cutoff = 1e6 , log = FALSE , substitute = c("1","approx") , terms = NULL )
                {
                   if( is.null( terms ) ) terms <- DBNormCGetTerms( n , cutoff = cutoff , approx = "exact" )
                   ljexpk <- 0
                   if( k > 0 ) ljexpk <- log(terms[[3]]^k)
                   summ <- terms[[1]] + (1-phi)*terms[[2]] + terms[[3]]*phi*log(p) + terms[[4]]*phi*log(1-p) + ljexpk
                   tm <- max(summ)
                   return( exp(tm + log(sum(exp(summ - tm)))) ) #more stable than the below return()
#                   return( sum(exp(terms[[1]] + (1-phi)*terms[[2]] + terms[[3]]*phi*log(p) + terms[[4]]*phi*log(1-p) + ljexpk)) )
                }

DBNormCOnePois <- function( mu , theta , approx = c("exact","approx","1") )
                {
                   approx = match.arg( approx , c("exact","approx","1") )

                   exFun <- function( mu , theta )
                             {
                               y1 = 0:1e6
                               y2 = y1 - y1*log(y1)
                               y2[1] = 0
                               sum(exp( y1*theta*log(mu) + (theta-1)*y2 - (theta*mu - 0.5*log(theta)) - lgamma(y1+1) ))
                             }
                   NormC <- switch( approx , "1" = 1 , "approx" = (1 + (1-theta)/(12*mu*theta) * ( 1 + 1/(mu*theta) )) ,
                                             "exact" = exFun( mu=mu , theta=theta ) )
                   return( NormC )
                }

DBNormC <- function( n , p , phi , approx = c("exact","approx","1") , k = 0 , cutoff = 1e6 , log = FALSE , terms = NULL )
            {
              approx = match.arg( approx , c("exact","approx","1") )
              if( approx == 1 ) return(rep(phi^(-0.5),length(n)))
              if( is.null( terms ) ) terms <- lapply( n , FUN=function(off) DBNormCGetTerms( off , cutoff = cutoff , approx = approx ) )
              NormC <- sapply(1:length(n),FUN=function(i)
                      {
                        DBNormCOne( n=n[i] , p=p , phi=phi , approx=approx , k=k , cutoff=cutoff , log=log , terms=terms[[i]] )
                      })
              return(NormC)
            }


ddbinom <- function( n , p , theta , x = NULL , limit = NULL , approx="exact")
            {
              mid <- n*p
              limit.def <- min(n,ceiling(mid+sqrt(n*p*(1-p))*100))
              if( is.null(limit) ) limit = limit.def
               else limit = max(limit,limit.def)
              kvals = 0:limit
              if( is.null(x) ) x <- 0:limit
              tmp1 = kvals*log(kvals)
              tmp1[is.na(tmp1)] = 0
              tmp2 = (n-kvals)*log((n-kvals))
              tmp2[is.na(tmp2)] = 0
              num <- exp(lgamma(n+1) - lgamma(kvals + 1) - lgamma(n - kvals + 1) + (1-theta)*(tmp1 + tmp2 - n*log(n)) + kvals*theta*log(p) + (n-kvals)*theta*log(1-p) )
              probs <- switch( approx , "exact"=num/sum(num), "1"=sqrt(theta)*num)
              return(probs[x+1])
            }

DerivCpsiOne <- function( n , p , phi , approx = c("exact","approx","1") , cutoff = Inf , der=c("1","2","both") , normC = NULL , terms = NULL )
            {
              approx = match.arg( approx , c("exact","approx","1") )
              der = as.character( der )
              der = match.arg( der , c("1","2","both") )

              if( is.null(normC) ) norm0 <- DBNormCOne( n, p, phi, k=0, approx=approx, cutoff=cutoff, terms=terms)
              else norm0 <- normC

              if( der == "1" )
               {
                  if( approx == "1" ) return( 0 )
                    else if( approx == "approx" | n > cutoff ) return( -norm0^2 * 1/(12*n)*(1-phi)/phi^2 * (1-2*p)/p/(1-p) * sqrt(phi) )
                      else
                       {
                         norm1 <- DBNormCOne( n, p, phi, k=1, approx="exact", cutoff=cutoff,terms=terms )
                         return( norm1 - n*p*norm0 )
                       }
               }
                else if( der == "2" )
                 {
                    if( approx == "1" ) return( 0 )
                      else if( approx == "approx" | n > cutoff )
                       {
                         der1 = -norm0^2 * 1/(12*n)*(1-phi)/phi^2 * (1-2*p)/p/(1-p) * sqrt(phi)
                         return( 2*der1^2/norm0 + norm0^2 * 1/(12*n)*(1-phi)/phi^3 * ( 2 + (1-2*p)^2/p/(1-p) ) * sqrt(phi) )
                       }
                      else
                       {
                         norm1 <- DBNormCOne( n , p , phi , k = 1, approx="exact" , cutoff = cutoff , terms = terms)
                         norm2 <- DBNormCOne( n , p , phi , k = 2, approx="exact" , cutoff = cutoff , terms = terms )
                         return( norm2 - 2*n*p*norm1 + (n^2*p^2 - n*p/phi + n*p^2/phi)*norm0 )
                       }
                 }
                else
                 {
                    if( approx == "1" ) return( c(0,0) )
                      else if( approx == "approx" | n > cutoff )
                       {
                         der1 = -norm0^2 * 1/(12*n)*(1-phi)/phi^2 * (1-2*p)/p/(1-p) * sqrt(phi)
                         return( c( der1 , 2*der1^2/norm0 + norm0^2 * 1/(12*n)*(1-phi)/phi^3 * ( 2 + (1-2*p)^2/p/(1-p) ) * sqrt(phi) ) )
                       }
                      else
                       {
                         norm1 <- DBNormCOne( n , p , phi , k = 1, approx="exact" , cutoff = cutoff , terms = terms )
                         norm2 <- DBNormCOne( n , p , phi , k = 2, approx="exact" , cutoff = cutoff , terms = terms )
                         return( c( norm1 - n*p*norm0 , norm2 - 2*n*p*norm1 + (n^2*p^2 - n*p/phi + n*p^2/phi)*norm0 ) )
                       }
                 }
            }

DerivCetapOne <- function( n , p , phi , approx = c("exact","approx","1") , cutoff = Inf , der=c("1","2","both") , normC = NULL , terms = NULL )
            {
              approx = match.arg( approx , c("exact","approx","1") )
              der = as.character( der )
              der = match.arg( der , c("1","2","both") )

              if( is.null(normC) ) norm0 <- DBNormCOne( n, p, phi, k=0, approx=approx, cutoff=cutoff, terms=terms)
              else norm0 <- normC

              if( der == "1" )
               {
                  if( approx == "1" ) return( 0 )
                    else if( approx == "approx" | n > cutoff ) return( norm0^2 * phi^(0.5) * 1/(12*n) * (1-phi)/phi * (2*p-1)/(p*(1-p)) )
                      else
                       {
                         norm1 <- DBNormCOne( n, p, phi, k=1, approx="exact", cutoff=cutoff , terms = terms )
                         return( (norm1 - n*p*norm0)*phi )
                       }
               }
                else if( der == "2" )
                 {
                    if( approx == "1" ) return( 0 )
                      else if( approx == "approx" | n > cutoff )
                       {
                         der1 = norm0^2 * phi^(0.5) * 1/(12*n) * (1-phi)/phi * (2*p-1)/(p*(1-p))
                         return( (1-phi)/phi^(0.5) * (1/(12*n)) * ( 2*norm0*(2*p-1)*der1/(p*(1-p)) + norm0^2*(p/(1-p)+(1-p)/p)) )
                       }
                      else
                       {
                         norm1 <- DBNormCOne( n, p, phi, k=1, approx="exact", cutoff=cutoff , terms = terms)
                         norm2 <- DBNormCOne( n, p, phi, k=2, approx="exact", cutoff=cutoff , terms = terms)
                         return( phi^2*norm2 - phi^2*n*p*norm1 - n*norm0*p*(1-p)*phi - n*p*phi^2*(norm1-n*p*norm0) )
                       }
                 }
                else
                 {
                    if( approx == "1" ) return( c(0,0) )
                      else if( approx == "approx" | n > cutoff )
                       {
                         der1 = norm0^2 * phi^(0.5) * 1/(12*n) * (1-phi)/phi * (2*p-1)/(p*(1-p))
                         return( c( der1 , (1-phi)/phi^(0.5) * (1/(12*n)) * ( 2*norm0*(2*p-1)*der1/(p*(1-p)) + norm0^2*(p/(1-p)+(1-p)/p)) ) )
                       }
                      else
                       {
                         norm1 <- DBNormCOne( n, p, phi, k=1, approx="exact",cutoff=cutoff , terms = terms)
                         norm2 <- DBNormCOne( n, p, phi, k=2, approx="exact",cutoff=cutoff , terms = terms)
                         return( c( (norm1 - n*p*norm0)*phi , phi^2*norm2 - phi^2*n*p*norm1 - n*norm0*p*(1-p)*phi - n*p*phi^2*(norm1-n*p*norm0) ) )
                       }
                 }
            }

DerivCpOne <- function( n , p , phi , approx = c("exact","approx","1") , cutoff = Inf , der=c("1","2","both") , normC = NULL , terms = NULL )
            {
              approx = match.arg( approx , c("exact","approx","1") )
              der = as.character( der )
              der = match.arg( der , c("1","2","both") )

              if( der == "1" ) return( DerivCetapOne( n, p, phi, approx=approx, cutoff=cutoff, der="1" , normC = normC  , terms = terms ) / (p*(1-p)) )
                else if( der == "2" )
                 {
                   derivs <- DerivCetapOne( n, p, phi, approx=approx, cutoff=cutoff, der="both" , normC = normC , terms = terms )
                   return( (derivs[2] + (2*p-1)*derivs[1])/(p*(1-p))^2 )
                 }
                else
                 {
                   derivs <- DerivCetapOne( n, p, phi, approx=approx, cutoff=cutoff, der="both" , normC = normC , terms = terms )
                   return( c( derivs[1]/(p*(1-p)) , (derivs[2] + (2*p-1)*derivs[1])/(p*(1-p))^2 ) )
                 }
            }

DerivCphiOne <- function( n , p , phi , approx = c("exact","approx","1") , cutoff = Inf , der=c("1","2","both") , normC = NULL , terms = NULL )
              {
                approx = match.arg( approx , c("exact","approx","1") )
                der = as.character( der )
                der = match.arg( der , c("1","2","both") )

                if( approx == "approx" | (n > cutoff & approx == "exact") )
                 if( is.null(normC) ) norm0 <- DBNormCOne( n, p, phi, k=0, approx=approx, cutoff=cutoff , terms=terms )
                  else norm0 <- normC

                if( der == "1" )
                 {
                   if( approx == "1" ) return(-0.5*phi^(-1.5))
                     else if( approx == "approx" | n > cutoff ) return( norm0^2*phi^(-0.5)/2 * (1/(12*n) * (1+phi)/phi * ( 1 - 1/(p*(1-p))) - 1) )
                       else
                        {
                          singles <- apply(compositions(n, 2)[2:1,2:n], 2, DBsingle,p=p,phi=phi,k=0)
                          termn <- n*log(p)*exp(n*phi*log(p))
                          term0 <- n*log(1-p)*exp(n*phi*log(1-p))
                          k <- (1:(n-1))
                          deriv = sum( c( term0 , (n*log(n)+n*log(1-p)-n*log(n-k)+k*log(p)-k*log(1-p)-k*log(k)+k*log(n-k))*singles , termn ) )
                          return( deriv )
                        }
                 }
                else if( der == "2" )
                 {
                   if( approx == "1" ) return(0.75*phi^(-2.5))
                     else if( approx == "approx" | n > cutoff )
                      {
                         der1 = norm0^2*phi^(-0.5)/2 * (1/(12*n) * (1+phi)/phi * ( 1 - 1/(p*(1-p))) - 1)
                         return( norm0*der1*( phi^(-3/2)*(1+phi)/(12*n) * (1 - 1/(p*(1-p))) - phi^(-0.5) ) + 0.5*norm0^2*(1/(12*n)*(1-1/(p*(1-p)))*(-1.5*phi^(-5/2)*(1+phi)+phi^(-3/2))+0.5*phi^(-3/2)) )
                      }
                     else
                      {
                        singles <- apply(compositions(n, 2)[2:1,2:n], 2, DBsingle,p=p,phi=phi,k=0)
                        termn <- n^2*log(p)^2*exp(n*phi*log(p))
                        term0 <- n^2*log(1-p)^2*exp(n*phi*log(1-p))
                        k <- (1:(n-1))
                        deriv = sum( c( term0 , (n*log(n)+n*log(1-p)-n*log(n-k)+k*log(p)-k*log(1-p)-k*log(k)+k*log(n-k))^2*singles , termn ) )
                        return( deriv )
                      }
                 }
                else
                 {
                   if( approx == "1" ) return( c(-0.5*phi^(-1.5),0.75*phi^(-2.5)) )
                     else if( approx == "approx" | n > cutoff )
                      {
                         der1 = norm0^2*phi^(-0.5)/2 * (1/(12*n) * (1+phi)/phi * ( 1 - 1/(p*(1-p))) - 1)
                         der2 = norm0*der1*( phi^(-3/2)*(1+phi)/(12*n) * (1 - 1/(p*(1-p))) - phi^(-0.5) ) + 0.5*norm0^2*(1/(12*n)*(1-1/(p*(1-p)))*(-1.5*phi^(-5/2)*(1+phi)+phi^(-3/2))+0.5*phi^(-3/2))
                         return( c( der1 , der2 ) )
                      }
                     else
                      {
                        singles <- apply(compositions(n, 2)[2:1,2:n], 2, DBsingle,p=p,phi=phi,k=0)
                        termn1 <- n*log(p)*exp(n*phi*log(p))
                        term01 <- n*log(1-p)*exp(n*phi*log(1-p))
                        termn2 <- n^2*log(p)^2*exp(n*phi*log(p))
                        term02 <- n^2*log(1-p)^2*exp(n*phi*log(1-p))
                        k <- (1:(n-1))
                        der1 = sum( c( term01 , (n*log(n)+n*log(1-p)-n*log(n-k)+k*log(p)-k*log(1-p)-k*log(k)+k*log(n-k))*singles , termn1 ) )
                        der2 = sum( c( term02 , (n*log(n)+n*log(1-p)-n*log(n-k)+k*log(p)-k*log(1-p)-k*log(k)+k*log(n-k))^2*singles , termn2 ) )
                        return( c(der1,der2) )
                      }
                 }
              }

DerivC <- function( n , p , phi , approx = c("exact","approx","1") , cutoff = Inf , der=c("0","1","2","both") , wrt = c( "psi" , "etap" , "phi" , "p" ) , normC = NULL , terms=NULL )
              {
                approx = match.arg( approx , c("exact","approx","1") )
                der = match.arg( der , c("0","1","2","both") )

                if( is.null(normC) ) normC <- DBNormC( n, p, phi, k=0, approx=approx, cutoff=cutoff , terms=terms )
                if( der == "0" ) return( normC )

                wrt = match.arg( wrt , c( "psi" , "etap" , "phi" , "p" ) )
                
                FUN = get( paste( "DerivC" , wrt , "One" , sep="" ) )

                if( is.null( terms ) ) terms <- lapply( n , FUN=function(off) DBNormCGetTerms( off , cutoff = cutoff  , approx = approx) )
                deriv <- sapply(1:length(n),FUN=function(i)
                      {
                        FUN( n=n[i] , p=p , phi=phi , approx=approx , cutoff=cutoff , der=der , normC = normC[i] , terms = terms[[i]] )
                      })
                return(deriv)
              }

DBAOne <- function( n, p, phi, psi = NULL , approx = c("exact","approx","1") , cutoff = Inf , der=c("0","1","2","both") , normC = NULL , terms = NULL )
            {
              der = match.arg( der , c("0","1","2","both") )
              approx = match.arg( approx , c("exact","approx","1") )
              
              if( !is.null(psi) ) 
               {
                 p <- exp(psi/phi-log(1+exp(psi/phi)))
                 phi <- phi
               }  
              
              if( is.null(normC) ) norm0 <- DBNormCOne( n, p, phi, k=0, approx=approx, cutoff=cutoff, terms=terms)
              else norm0 <- normC

              if( der == "0" ) return(  log( norm0 ) - n*phi*log(1-p) )
                else if( der == "1" ) return( n*p + 1/norm0 * DerivCpsiOne(n,p,phi,approx=approx,cutoff=cutoff,der="1",normC=norm0,terms=terms) )
                  else if( der == "2" )
                   {
                     derivs <- DerivCpsiOne(n,p,phi,approx=approx,cutoff=cutoff,der="both",normC=norm0,terms=terms)
                     return( n*p*(1-p)/phi + derivs[2]/norm0 - (derivs[1] / norm0)^2 )
                   }
                  else
                   {
                     derivs <- DerivCpsiOne(n,p,phi,approx=approx,cutoff=cutoff,der="both",normC=norm0,terms=terms)
                     return( c( n*p + 1/norm0 * derivs[1] , n*p*(1-p)/phi + derivs[2]/norm0 - (derivs[1] / norm0)^2 ) )
                   }
            }

DBA <- function( n , p , phi , approx = c("exact","approx","1") , cutoff = Inf , der=c("0","1","2","both") , normC = NULL , terms = NULL )
        {
          approx = match.arg( approx , c("exact","approx","1") )
          der = match.arg( der , c("0","1","2","both") )

          if( is.null(normC) ) norm0 <- DBNormCOne( n, p, phi, k=0, approx=approx, cutoff=cutoff, terms=terms)
          else norm0 <- normC

          if( is.null( terms ) ) terms <- lapply( n , FUN=function(off) DBNormCGetTerms( off , cutoff = cutoff , approx = approx ) )
          A <- sapply(1:length(n),FUN=function(i)
                {
                  DBAOne(n[i],p,phi,approx=approx,cutoff=cutoff,der=der,normC=normC[i],terms=terms[[i]])
                })
          d2A <- sum(A)
          return(d2A)
        }
        
rDB <- function (size, n, p , phi , burnin = 4 * n, every = 4 * n, start = NULL)
        {
            p <- c( p , 1-p )
            w <- length(p)
            update <- function(xin) {
                jj <- sample(w, 2, replace = FALSE)
                kernel <- rep(0L, w)
                kernel[jj] <- c(-1L, +1L)
                proposed <- xin + kernel
                if (any(proposed < 0)) {
                    num <- 0
                }
                else {
                    num <- DBsingle(proposed, p , phi )
                }
                den <- DBsingle(xin, p , phi )
                if ((num == 0) & (den == 0)) {
                    print("this cannot happen")
                    alpha <- 0
                }
                else {
                    alpha <- min(1, num/den)
                }
                if (runif(1) < alpha) {
                    ans <- proposed
                }
                else {
                    ans <- xin
                }
                return(ans)
            }
            if (is.null(start)) {
                start <- tabulate(sample(w, size = n, prob = p , replace = TRUE), nbins = w)
            }
            for (i in seq_len(burnin)) {
                start <- update(start)
            }
            out <- matrix(0, size, w)
            colnames(out) <- names(p)
            for (i in seq_len(n)) {
                for (j in seq_len(every)) {
                    start <- update(start)
                }
                out[i, ] <- start
            }
            return(out)
        }

DB_loglik <- function( y, n, p, phi , groups = NULL , approx = c("exact","approx","1") , cutoff = Inf , normC = NULL , terms = NULL )
  {
    approx = match.arg( approx , c("exact","approx","1") )

    if( is.null( terms ) ) terms <- lapply( n , FUN=function(off) DBNormCGetTerms( off , cutoff = cutoff , approx = "exact" ) )
    if( is.null(normC) ) normC = DBNormC( n=n , p=p , phi=phi , approx=approx , cutoff=cutoff , terms = terms )

    single <- sum(sapply( 1:length(n) , FUN=function(i) terms[[i]][[1]][y[i]+1] + (1-phi)*terms[[i]][[2]][(y[i]+1)] + terms[[i]][[3]][y[i]+1]*phi*log(p) + terms[[i]][[4]][y[i]+1]*phi*log(1-p)))

    return( -sum(log(normC)) + single )
  }

DB_loglik_2g <- function( y, n, p, phi , groups , approx = c("exact","approx","1") , cutoff = Inf , normC = NULL , terms = NULL )
  {
    approx = match.arg( approx , c("exact","approx","1") )
    groups <- as.factor(groups)
    J = length(unique(groups))

    if( is.null( terms ) ) terms <- lapply( n , FUN=function(off) DBNormCGetTerms( off , cutoff = cutoff , approx = "exact" ) )
    yy = DBSplitIntoGroups( y , groups )[ as.character(unique(groups)) ]
    nn = DBSplitIntoGroups( n , groups )[ as.character(unique(groups)) ]
    ll <- sum(sapply( 1:J , FUN=function(ii)
               {
                 p <- p[ii]
                 y <- yy[[ii]]
                 n <- nn[[ii]]
                 terms <- terms[grep(unique(groups)[ii],groups,fixed=TRUE)]
                 normC <- DBNormC( n , p , phi , approx=approx , cutoff=cutoff ,terms=terms)
                 single <- sum(sapply( 1:length(n) , FUN=function(i) terms[[i]][[1]][y[i]+1] + (1-phi)*terms[[i]][[2]][(y[i]+1)] + terms[[i]][[3]][y[i]+1]*phi*log(p) + terms[[i]][[4]][y[i]+1]*phi*log(1-p)))
                 return( -sum(log(normC)) + single )
               } ))
    return(ll)
  }

DB_loglik_dp <- function( y, n, p, phi, approx = c("exact","approx","1") , cutoff = Inf , normC = NULL , terms = NULL )
  {
    approx = match.arg( approx , c("exact","approx","1") )
    if( is.null( terms ) ) terms <- lapply( n , FUN=function(off) DBNormCGetTerms( off , cutoff = cutoff , approx = approx ) )

    if( length(n)==1 ) n = rep(n,length(y))

    if( is.null(normC) ) normC = DBNormC( n=n , p=p , phi=phi , approx=approx , cutoff=cutoff , terms = terms )

    dc.dp <- DerivC(n=n,p=p,phi=phi,approx=approx,cutoff=cutoff,der="1",wrt="p",normC=normC,terms=terms)

    return( sum(phi*c(1/p,-1/(1-p))*c(sum(y),sum(n-y))) - sum(dc.dp/normC) )
  }

transform.p <- function( p , phi = NULL , param = c( "p" , "psi" , "etap" ) , inverse = FALSE )
                   {
                     param = match.arg( param , c( "p" , "psi" , "etap" ) )
                     if( param == "psi" & is.null(phi) ) stop( "A value for 'phi' is needed to parameterize as 'psi'" )
                     if( inverse ) p = switch( param , "p" = p , "psi" = exp(p/phi-log(1+exp(p/phi))) , etap = exp(p-log(1+exp(p))) )
                      else p = switch( param , "p" = p , "psi" = phi*(log(p)-log(1-p)) , etap = (log(p)-log(1-p)) )
                     return(p)                                     
                   }
 

DB_optimize_2g <- function (y, n, phi , groups , start = NULL, approx = c("exact","approx","1") , cutoff = Inf , terms=NULL , method = "nlm", printing = FALSE, give_fit = FALSE, ...)
  {
    approx = match.arg( approx , c("exact","approx","1") )
    if( is.null(start) )
     {
       yy = DBSplitIntoGroups( y , groups )[ as.character(unique(groups)) ]
       nn = DBSplitIntoGroups( n , groups )[ as.character(unique(groups)) ]
       start = mapply(FUN=function(y,n) sum(y)/sum(n) , yy , nn)
       if( start[1] == 1 ) start[1] = 1 - sqrt(.Machine$double.eps)
       if( start[2] == 1 ) start[2] = 1 - sqrt(.Machine$double.eps)
       if( start[1] == 0 ) start[1] = sqrt(.Machine$double.eps)
       if( start[2] == 0 ) start[2] = sqrt(.Machine$double.eps)
     }

    if( length(n)==1 ) n = rep(n,length(y))

    if( is.null(groups) ) groups = rep( 1 , length(y) )
    groups <- as.factor(groups)
    J = length(unique(groups))

    if( is.null( terms ) ) terms <- lapply( n , FUN=function(off) DBNormCGetTerms( off , cutoff = cutoff , approx = "exact" ) )

    f <- function(xin) {
        if (any(xin<0)) {
            out <- .Machine$double.xmax
        }
        else {
            out <- -DB_loglik_2g(y=y, n=n, p=xin , groups=groups , phi=phi , approx = approx , cutoff = cutoff , terms=terms )
        }
        return(out)
    }
    if (method == "Nelder") {
        out <- optim(start, f,...)
        short_out <- out$par
    }
    else if (method == "nlm") {
        out <- nlm(f, start , ...)
        short_out <- out$estimate
    }
    else {
        print("method must be Nelder or nlm")
        stop()
    }
    return(short_out)
  }

DB_optimize <- function (y, n, phi = NULL , groups = NULL , param = c( "p" , "psi" , "etap" ) , start = NULL, approx = c("exact","approx","1") , cutoff = Inf , terms=NULL , method = "nlm", printing = FALSE, give_fit = FALSE, ...)
  {
    approx = match.arg( approx , c("exact","approx","1") )
    param = match.arg( param , c( "p" , "psi" , "etap" ) )

    if( length(n)==1 ) n = rep(n,length(y))

    if( is.null(groups) ) groups = rep( 1 , length(y) )
    groups <- as.factor(groups)
    J = length(unique(groups))

    if( is.null( terms ) ) terms <- lapply( n , FUN=function(off) DBNormCGetTerms( off , cutoff = cutoff , approx = "exact" ) )

    f <- function(xin,param) {
        if( length(xin)==2 )
         {
           x <- c( transform.p( xin[1] , phi = xin[2] , param = param , inverse = TRUE ) , xin[2] )
         }
         else
          {                                      
            xin <- transform.p( xin , phi = phi , param = param , inverse = TRUE )
            x <- c(xin,phi)
          }
        if (any(x<0)) {
            out <- .Machine$double.xmax
        }
        else {
            out <- -DB_loglik(y, n, x[1], x[2] , approx = approx , cutoff = cutoff , terms=terms )
        }
        return(out)
    }
    if (method == "Nelder") {
        out <- optim(c(0.5,0.5), f, param = param,...)
        short_out <- out$par
    }
    else if (method == "Brent") {
        lower <- transform.p( sqrt(.Machine$double.eps) , phi = phi , param = param , inverse = FALSE )
        upper <- transform.p( 1-sqrt(.Machine$double.eps) , phi = phi , param = param , inverse = FALSE )
        out <- optimize(f, c(lower,upper), param = param, ... )
        short_out <- out$minimum
    }
    else if (method == "nlm") {
        out <- nlm(f, 0.5, param = param, ...)
        short_out <- out$estimate
    }
    else {
        print("method must be Nelder or nlm or Brent (single variable)")
        stop()
    }
    if( is.null(phi) ) jj <- c(transform.p(short_out[1],short_out[2],param=param,inverse=TRUE),short_out[2])
      else jj <- transform.p(short_out,phi=phi,param=param,inverse=TRUE)
    return(jj)
  }

DB_optimize_acml <- function (y, n, groups=NULL , start = NULL, approx = c( "all.exact" , "all.approx" , "all.1" , "prop.approx" , "prop.1" ) , cutoff = Inf , terms = NULL , profile = FALSE , method = "Brent", printing = FALSE, give_fit = FALSE, ...)
  {
    approx = match.arg( approx , c( "all.exact" , "all.approx" , "all.1" , "prop.approx" , "prop.1" ) )

    if( length(n)==1 ) n = rep(n,length(y))

    if( is.null(groups) ) groups = rep( 1 , length(y) )
    groups <- as.factor(groups)
    J = length(unique(groups))

    approx.2 <- "1"
    if( approx %in% c("all.exact","prop.1","prop.approx") ) approx.2 <- "exact"
    if( is.null( terms ) ) terms <- lapply( n , FUN=function(off) DBNormCGetTerms( off , cutoff = cutoff , approx = approx.2 ) )

    if (is.null(start)) {
        start <- 0.5
    }

    f <- function(xin) {
            x <- xin/(1-xin)
        if (findInterval(xin,c(0,1))!=1) {
            out <- .Machine$double.xmax
        }
        else {
            out <- -DBCondLogLik(x, y, n, , groups = groups , method="Brent" , profile=profile , approx = approx , cutoff = cutoff , reparam="none", terms=terms)
        }
        return(out)
    }
    if (method == "Brent") {
        out <- optimize(f, c(0,1) ,...)
        short_out <- out$minimum
    }
    else if (method == "nlm") {
        out <- nlm(f, start, ...)
        short_out <- out$estimate
    }
    else {
        print("method must be Brent or nlm")
        stop()
    }
      return(short_out/(1-short_out))
  }

#DBCondLogLik <- function( phi , y, n, name = NULL , method = "Brent" , approx = c( "all.exact" , "all.approx" , "all.1" , "prop.approx" , "prop.1" ) , terms=NULL , profile = FALSE , cutoff = 500 , reparam = c( "none" , "0-1" , "log" ) )
#                  {
#                    approx = match.arg( approx , c( "all.exact" , "all.approx" , "all.1" , "prop.approx" , "prop.1" ) )
#                    reparam = match.arg( reparam , c( "none" , "0-1" , "log" ) )
#
#                    approx.2 <- "1"
#                    if( approx %in% c("all.exact","prop.1","prop.approx") ) approx.2 <- "exact"
#                    if( is.null( terms ) ) terms <- lapply( n , FUN=function(off) DBNormCGetTerms( off , cutoff = cutoff , approx = approx.2 ) )
#
#                    phi = switch( reparam , "none" = phi , "0-1" = phi/(1-phi) , "log" = exp(phi) )
#                    approx1 <- switch( approx , all.exact = "exact" , all.approx = "approx" , all.1 = "1" , prop.approx = "approx" , prop.1 = "1" )
#                    approx2 <- switch( approx , all.exact = "exact" , all.approx = "approx" , all.1 = "1" , prop.approx = "exact" , prop.1 = "exact" )
#
#                    if( approx1 == "1" )
#                     {
#                       p.mle <- sum(y)/sum(n)
#                       if( p.mle == 1 | p.mle == 0 ) p.mle = abs(p.mle - sqrt(.Machine$double.eps))
#                     }
#                    else p.mle <- DB_optimize(y, n , phi = phi , method = "Brent" , approx = approx1 , cutoff = cutoff , terms=terms)
#
#                    normC <- DBNormC( n , p.mle , phi , approx=approx2 , cutoff=cutoff ,terms=terms)
#
#                    if( profile ) d2A <- 1
#                     else d2A <- DBA( n, p.mle, phi, approx = approx2 , cutoff = cutoff , der="2", normC = normC ,terms=terms)
#
#                    if( d2A == 0 | is.na(d2A) )
#                     {
#                       if( is.null(name) ) name <- 'not given'
#                       warning( paste( "For exon '" , name , "' d2A evaluated to" , d2A , ". Possible Numerical Instability. Setting approx = '1'." , sep = "" ) )
#                       d2A <- DBA( n, p.mle, phi, approx = "1" , cutoff = cutoff , der="2" , normC = normC , terms=terms )
#                     }
#
#                    cll <- DB_loglik(y, n, p.mle, phi, approx = approx2 , cutoff = cutoff , normC = normC , terms=terms ) + 0.5*log(d2A)
#
#                    return( cll )
#                  }

DBCondLogLik <- function( phi , y, n, name = NULL , groups = NULL , method = "Brent" , approx = c( "all.exact" , "all.approx" , "all.1" , "prop.approx" , "prop.1" ) , terms=NULL , profile = FALSE , cutoff = Inf , reparam = c( "none" , "0-1" , "log" ) )
                  {
                    approx = match.arg( approx , c( "all.exact" , "all.approx" , "all.1" , "prop.approx" , "prop.1" ) )
                    reparam = match.arg( reparam , c( "none" , "0-1" , "log" ) )

                    if( is.null(groups) ) groups = rep( 1 , length(y) )
                    groups <- as.factor(groups)
                    J = length(unique(groups))

                    phi = switch( reparam , "none" = phi , "0-1" = phi/(1-phi) , "log" = exp(phi) )
                    approx1 <- switch( approx , all.exact = "exact" , all.approx = "approx" , all.1 = "1" , prop.approx = "approx" , prop.1 = "1" )
                    approx2 <- switch( approx , all.exact = "exact" , all.approx = "approx" , all.1 = "1" , prop.approx = "exact" , prop.1 = "exact" )

                    if( is.null( terms ) ) terms <- lapply( n , FUN=function(off) DBNormCGetTerms( off , cutoff = cutoff , approx = approx2 ) )

                    yy = DBSplitIntoGroups( y , groups )[ as.character(unique(groups)) ]
                    nn = DBSplitIntoGroups( n , groups )[ as.character(unique(groups)) ]
                    cll <- sum(sapply( 1:J , FUN=function(ii)
                               {
                                 y <- yy[[ii]]
                                 n <- nn[[ii]]
                                 terms <- terms[grep(unique(groups)[ii],groups,fixed=TRUE)]
                                 if( approx1 == "1" )
                                  {
                                    p.mle <- sum(y)/sum(n)
                                    if( p.mle == 1 | p.mle == 0 ) p.mle = abs(p.mle - sqrt(.Machine$double.eps))
                                  }
                                 else p.mle <- DB_optimize(y, n , phi = phi , groups=NULL , method = "Brent" , approx = approx1 , cutoff = cutoff , terms=terms)
                                 normC <- DBNormC( n , p.mle , phi , approx=approx2 , cutoff=cutoff ,terms=terms)
                                 if( profile ) d2A <- 1
                                  else d2A <- DBA( n, p.mle, phi, approx = approx2 , cutoff = cutoff , der="2", normC = normC ,terms=terms)
                                 if( d2A == 0 | is.na(d2A) )
                                   {
                                     if( is.null(name) ) name <- 'not given'
                                     warning( paste( "For exon '" , name , "' d2A evaluated to" , d2A , ". Possible Numerical Instability. Setting approx = '1'." , sep = "" ) )
                                     d2A <- DBA( n, p.mle, phi, approx = "1" , cutoff = cutoff , der="2" , normC = normC , terms=terms )
                                   }
                                 cll <- DB_loglik(y, n, p.mle, phi, approx = approx2 , cutoff = cutoff , normC = normC , terms=terms ) + 0.5*log(d2A)
                                 return(cll)
                               } ))
                     return(cll)
                  }
              

DBEstimateJointDisp <- function (y,n, groups=NULL, names = NULL  )
 {
       tagnames <- rownames(y)
       if( is.null(groups) ) groups = rep( 1 , ncol(y) )
       groups <- as.factor(groups)
       J = length(unique(groups))
       G = nrow(y)
       Mj <- matrix( table(groups)[unique(as.character(groups))] , nr = G , nc = J )
       if( any(n==0) )
        {
          wh.0 <- which(apply(n,1,FUN=function(vec)any(vec==0)))
          Mj[wh.0,] = t(apply(n[wh.0,,drop=FALSE],1,FUN=function(vec) table(groups[which(vec!=0)])[unique(as.character(groups))]))
        }
       term1 <- y*log(y)
       term1[is.na(term1)] <- 0
       term2 <- (n-y)*log(n-y)
       term2[is.na(term2)] <- 0
       term3 <- n*log(n)
       term3[is.na(term3)] <- 0
       termsum1 <- rowSums( term1 + term2 - term3 )

       yy = DBSplitIntoGroups( y , groups )[ as.character(unique(groups)) ]
       nn = DBSplitIntoGroups( n , groups )[ as.character(unique(groups)) ]

       phat <- mapply( FUN = function( cts , offs ) rowSums(cts)/rowSums(offs) , yy , nn , SIMPLIFY = TRUE )
       phat <- matrix( phat[,as.character(groups)] , nr = G , nc = length(groups) )

       termsum2 <- rowSums(y*log( ifelse(phat==0,1,phat) ) + (n-y)*log( ifelse(phat==1,1,1-phat) ))

       denom = termsum1 - termsum2
       disp = 0.5*(rowSums(Mj))/denom
       disp[ disp > 500 ] <- 500
       disp[ disp < sqrt(.Machine$double.eps) ] <- sqrt(.Machine$double.eps)
       names(disp) <- NULL
       return(disp)
 }

DBEstimateACMLDisp <- function (y,n, groups, names = NULL , approx = c( "all.exact" , "all.approx" , "all.1" , "prop.approx" , "prop.1" ) , ncore = 1 , cutoff = Inf )
 { 
    approx = match.arg( approx , c( "all.exact" , "all.approx" , "all.1" , "prop.approx" , "prop.1" ) )
    if( is.null(groups) ) groups <- rep(1,ncol(y))
       groups <- as.factor(groups)
    
    if( approx == "all.1" )
     {
       tagnames <- rownames(y)
       J = length(unique(groups))
       G = nrow(y)
       Mj <- matrix( table(groups)[unique(as.character(groups))] , nr = G , nc = J )          
       if( any(n==0) ) 
        {
          wh.0 <- which(apply(n,1,FUN=function(vec)any(vec==0)))
          Mj[wh.0,] = t(apply(n[wh.0,,drop=FALSE],1,FUN=function(vec) table(groups[which(vec!=0)])[unique(as.character(groups))]))
        }  
       term1 <- y*log(y)
       term1[is.na(term1)] <- 0 
       term2 <- (n-y)*log(n-y)
       term2[is.na(term2)] <- 0 
       term3 <- n*log(n)
       term3[is.na(term3)] <- 0 
       termsum1 <- rowSums( term1 + term2 - term3 )

       yy = DBSplitIntoGroups( y , groups )[ as.character(unique(groups)) ]
       nn = DBSplitIntoGroups( n , groups )[ as.character(unique(groups)) ]       

       phat <- mapply( FUN = function( cts , offs ) rowSums(cts)/rowSums(offs) , yy , nn , SIMPLIFY = TRUE )
       phat <- matrix( phat[,as.character(groups)] , nr = G , nc = length(groups) )
       
       termsum2 <- rowSums(y*log( ifelse(phat==0,1,phat) ) + (n-y)*log( ifelse(phat==1,1,1-phat) ))
       
       denom = termsum1 - termsum2
       acml.dispersions = 0.5*(rowSums(Mj)-J)/denom
#       acml.dispersions[ acml.dispersions > 500 ] <- 500
#       acml.dispersions[ acml.dispersions < 1e-4 ] <- 1e-4
       acml.dispersions[ acml.dispersions > 1e16 ] <- 1e16
       acml.dispersions[ acml.dispersions < 1e-16 ] <- 1e-16
     }
    else
     {  
       if( ncore == 1 ) 
        acml.delta <- lapply( 1:nrow(y) , FUN = function( i ) optimize(DBCombCondLogLik,interval=c(1e-04, 1000/(1000 + 1)),y=y[i,,drop=FALSE],n=n[i,,drop=FALSE],groups=groups,tol=1e-4,maximum=TRUE,approx=approx,ncore=ncore,cutoff=cutoff,reparam="0-1",names=names,doSum=FALSE)$maximum )
         else
          acml.delta <- mclapply( 1:nrow(y) , FUN = function( i ) optimize(DBCombCondLogLik,interval=c(1e-04, 1000/(1000 + 1)),y=y[i,,drop=FALSE],n=n[i,,drop=FALSE],groups=groups,tol=1e-4,maximum=TRUE,approx=approx,ncore=ncore,cutoff=cutoff,reparam="0-1",names=names,doSum=FALSE)$maximum , mc.preschedule = TRUE , mc.cores = nCore )
  
       acml.dispersions <- unlist(acml.delta)/(1-unlist(acml.delta))
       acml.dispersions[ acml.dispersions > 500 ] <- 500
       acml.dispersions[ acml.dispersions < sqrt(.Machine$double.eps) ] <- sqrt(.Machine$double.eps)
     }  
    names(acml.dispersions) <- NULL
    return( acml.dispersions )
 }

DBEstimateEBDisp <- function (y,n,groups=NULL,nsmp=NULL,method=c("MLE","MM","P1","P3","WL","USER"),delta=NULL,median=FALSE,return.delta=FALSE,return.pq=FALSE,pq=NULL,dispersions=NULL,type.exact=c("full","limited","reverselimited"),subset=NULL,splitparams=FALSE)
 {
   if( is.null(subset) ) subset = rownames(y)
   if( is.null(groups) ) groups = rep( 1 , ncol(y) )
   groups <- as.factor(groups)
   method = match.arg( method , c("MLE","MM","P1","P3","WL","USER") )
   if( !is.null(dispersions) & length(type.exact)==2 ) warning( "Exact method is being implemented with the default for type.exact='full'" )
   type.exact = match.arg( type.exact , c("full","limited","reverselimited") )
   J = length(unique(groups))

   if( is.null(nsmp) )
    {
      nsmp <- rep(ncol(y),nrow(y))
      names(nsmp) <- rownames(y)
    }

   if( is.null(delta) ) delta = 20/(nsmp-J)

   if( ( is.null(pq) | length(pq)!=2 ) & method == "USER" ) stop("method 'USER' requires user to input a vector of 2 numeric values for argument 'pq'")
   if( is.numeric(pq) & length(pq)==2 ) method = "USER"

   CGMMEst <- function( S , alpha )
           {
             m1.S <- mean(S)
             m2.S <- mean(S^2)
             qmm <- m2.S*m1.S / ( alpha*m2.S - m1.S^2*(1+alpha) )
             bmm <- 1 + qmm*alpha / m1.S
             return(c(bmm,qmm))
           }
   transform1 <- function(x) x/(1-x)
   disp.EB <- vector( length = nrow(y) )
   S <- DBS( y , n , groups )
   if( (type.exact == "full" | type.exact == "reverselimited") & !is.null(dispersions) ) S = ncol(y)/(dispersions*2)
   S.sub <- S[subset]
   nsmp.sub <- nsmp[subset]
      names(S.sub) <- as.character(nsmp.sub)
      names(S) <- as.character(nsmp)
   N = length(S)
   a = (nsmp-J)/2
   a.sub <- a[subset]
   m.S = tapply( S.sub , names(S.sub) , mean )
   m.S = m.S[match( names(S) , names(m.S) )]
   m.S[is.na(m.S)] <- max(m.S,na.rm=TRUE)
   if( splitparams )
    {
      wh.notbdry <- which(nsmp==ncol(y))
      wh.bdry <- which(nsmp<ncol(y))
      wh.notbdry.sub <- which(nsmp.sub==ncol(y))
      wh.bdry.sub <- which(nsmp.sub<ncol(y))
      m.S.bdry <- m.S[wh.bdry]
      m.S.notbdry <- m.S[wh.notbdry]
      S.bdry <- S[wh.bdry]
      S.notbdry <- S[wh.notbdry]
      a.bdry <- a.sub[wh.bdry]
      a.notbdry <- a.sub[wh.notbdry]
      a.sub.bdry <- a.sub[wh.bdry.sub]
      a.sub.notbdry <- a.sub[wh.notbdry.sub]
      S.sub.bdry <- S.sub[wh.bdry.sub]
      S.sub.notbdry <- S.sub[wh.notbdry.sub]
      pars.bdry = switch( method , "MLE" = nlm( nll.gbp2 , p = c( 1 , 1 ) , a = a.sub.bdry , S = S.sub.bdry , hessian = FALSE , iterlim = 1000 )$estimate ,
                             "MM" = CGMMEst( S.sub.bdry , (ncol(y)-J)/2 ) ,
                             "P1" = transform1(optim( par = 0.5 , fn = nll.gbp.delta2 , gr = gr.nlgbp.delta , a = a.sub.bdry , S = S.sub.bdry , method = "Brent" , lower = 0 , upper = 1 )$par) * cbind(a.bdry,m.S.bdry) ,
                             "P3" = optim( par = c( CGMMEst( S.sub.bdry , (ncol(y)-J)/2 ) , 1 ) , fn = nll.gbp.3 , gr = gr.nlgbp.3 , a = (ncol(y)-J)/2 , S = S.sub.bdry )$par[1:2] ,
                             "WL" =  delta * cbind(a.sub.bdry,m.S.bdry) ,
                             "USER" = pq )
      pars.notbdry = switch( method , "MLE" = nlm( nll.gbp2 , p = c( 1 , 1 ) , a = a.sub.notbdry , S = S.sub.notbdry , hessian = FALSE , iterlim = 1000 )$estimate ,
                             "MM" = CGMMEst( S.sub.notbdry , (ncol(y)-J)/2 ) ,
                             "P1" = transform1(optim( par = 0.5 , fn = nll.gbp.delta2 , gr = gr.nlgbp.delta , a = a.sub.notbdry , S = S.sub.notbdry , method = "Brent" , lower = 0 , upper = 1 )$par) * cbind(a.notbdry,m.S.notbdry) ,
                             "P3" = optim( par = c( CGMMEst( S.sub.notbdry , (ncol(y)-J)/2 ) , 1 ) , fn = nll.gbp.3 , gr = gr.nlgbp.3 , a = (ncol(y)-J)/2 , S = S.sub.notbdry )$par[1:2] ,
                             "WL" =  delta * cbind(a.sub.notbdry,m.S.notbdry) ,
                             "USER" = pq )
       if( !median )
         if(is.matrix(pars.bdry)) disp.EB.bdry <- (pars.bdry[,1] + a.bdry)/(pars.bdry[,2]+S.bdry)
         else disp.EB.bdry <- (pars.bdry[1] + a.bdry)/(pars.bdry[2]+S.bdry)
       else disp.EB.bdry <- qgamma(0.5,shape=(pars.bdry[1]+a.bdry),rate=(pars.bdry[2]+S.bdry))
       if( !is.null(dispersions) & type.exact == "limited" )
         if(is.matrix(pars.bdry)) disp.EB.bdry <- (pars.bdry[,1] + a.bdry)/(pars.bdry[,2]+ncol(y)/(dispersions*2))
          else disp.EB.bdry <- (pars.bdry[1] + a.bdry)/(pars.bdry[2]+ncol(y)/(dispersions*2))
       names(disp.EB.bdry)<-rownames(y)[wh.bdry]

       if( !median )
         if(is.matrix(pars.notbdry)) disp.EB.notbdry <- (pars.notbdry[,1] + a.notbdry)/(pars.notbdry[,2]+S.notbdry)
         else disp.EB.notbdry <- (pars.notbdry[1] + a.notbdry)/(pars.notbdry[2]+S.notbdry)
       else disp.EB.notbdry <- qgamma(0.5,shape=(pars.notbdry[1]+a.notbdry),rate=(pars.notbdry[2]+S.notbdry))
       if( !is.null(dispersions) & type.exact == "limited" )
         if(is.matrix(pars.notbdry)) disp.EB.notbdry <- (pars.notbdry[,1] + a.notbdry)/(pars.notbdry[,2]+ncol(y)/(dispersions*2))
          else disp.EB.notbdry <- (pars.notbdry[1] + a.notbdry)/(pars.notbdry[2]+ncol(y)/(dispersions*2))
       names(disp.EB.notbdry) <- rownames(y)[wh.notbdry]

       disp.EB <- c(disp.EB.bdry,disp.EB.notbdry)
       disp.EB <- disp.EB[rownames(y)]
   if( return.delta ) return( c("bdry"=pars.bdry[1]/sample(a.sub.bdry,1),"notbdry"=pars.notbdry[1]/sample(a.sub.notbdry,1)) )
   if( return.pq ) return( matrix(c(pars.bdry,pars.notbdry),nc=2,nr=2,dimnames=list(c("alpha","beta"),c("bdry","notbdry"))) )
    }
   else
    {
   pars = switch( method , "MLE" = nlm( nll.gbp2 , p = c( 1 , 1 ) , a = a.sub , S = S.sub , hessian = FALSE , iterlim = 1000 )$estimate ,
                           "MM" = CGMMEst( S.sub , (ncol(y)-J)/2 ) ,
                           "P1" = transform1(optim( par = 0.5 , fn = nll.gbp.delta2 , gr = gr.nlgbp.delta , a = a.sub , S = S.sub , method = "Brent" , lower = 0 , upper = 1 )$par) * cbind(a,m.S) ,
                           "P3" = optim( par = c( CGMMEst( S , (ncol(y)-J)/2 ) , 1 ) , fn = nll.gbp.3 , gr = gr.nlgbp.3 , a = (ncol(y)-J)/2 , S = S.sub )$par[1:2] ,
                           "WL" =  delta * cbind(a.sub,m.S) ,
                           "USER" = pq )
   if( type.exact == "reverselimited" ) S <- DBS( y , n , groups )
   if( !median )
     if(is.matrix(pars)) disp.EB <- (pars[,1] + a)/(pars[,2]+S)
     else disp.EB <- (pars[1] + a)/(pars[2]+S)
   else disp.EB <- qgamma(0.5,shape=(pars[1]+a),rate=(pars[2]+S))
#   S = n/2*phi, phi=2*S/n, theta=n/(2*S), 2*theta/n = 1/S , S=n/(2*theta)=ncol(y)/(2*dispersions)
   if( !is.null(dispersions) & type.exact == "limited" )
     if(is.matrix(pars)) disp.EB <- (pars[,1] + a)/(pars[,2]+ncol(y)/(dispersions*2))
      else disp.EB <- (pars[1] + a)/(pars[2]+ncol(y)/(dispersions*2))
   disp.EB[disp.EB<=sqrt(.Machine$double.eps)] <- sqrt(.Machine$double.eps)
   if( return.delta ) return( pars[1]/sample(a,1) )
   if( return.pq ) return( pars )
    }
   return(disp.EB)
 }

#WEBSeqFilter <- function( y , n , groups )
#                 {
#                    if( any(is.na(y)) )
#                     {
#                       y[is.na(y)] <- 0
#                       warning( "NA values detected in count matrix 'y', replaced with 0" )
#                     }
#                    if( any(is.na(n)) )
#                     {
#                       n[is.na(n)] <- 0
#                       warning( "NA values detected in offset matrix 'n', replaced with 0" )
#                     }
#
#                    cols1 <- groups == unique(groups)[1]
#                    cols2 <- groups == unique(groups)[2]
#                    wh.kp <- rowSums(y)!=rowSums(n) & rowSums(y)!=0 & rowSums(n[,cols1,drop=FALSE])!=0 & rowSums(n[,cols2,drop=FALSE])!=0
#                    y <- y[wh.kp,]
#                    n <- n[wh.kp,]
#                    return( list("Counts"=y,"Offsets"=n) )
#                 }

DBGLM1 <- function ( y , n , groups , type=c("common","tagwise","EB","EBMM","EB3","WEB","ACML","joint","MAX","WEBOver","QWEB","USER") , fitDNull=TRUE , median = FALSE , dispersion = NULL , nulldispersion = NULL , prior.df = 20 , maxDisp = NULL , sigtest = c("LRT", "LRT2" , "WALD","WALDt","F","ALL") , combotest = FALSE , comboval = 0.01 , underdispersion = TRUE , ttest=TRUE , trend=FALSE , bc = TRUE , return.disp = FALSE , delta=NULL , pq=NULL)
 { 
   logitF <- function( x , inverse = FALSE ,replacement=sqrt(.Machine$double.eps))
    {
      if( !inverse & (any(x<0,na.rm=TRUE) | any(x>1,na.rm=TRUE)) ) stop( "x must be between 0 & 1 inclusive" )
      x[x==0] <- replacement
      x[x==1] <- 1-replacement
      val = switch( as.character(inverse) , "FALSE" = log(x) - log(1-x) , "TRUE" = exp( x - log(1+exp(x)) ) )
      return(val)
    }  
   if( comboval >0.5 ) comboval = 1-comboval
   tagnames <- rownames(y) 
   type = match.arg( type , c("common","tagwise","EB","EBMM","EB3", "WEB" ,"ACML","joint","MAX","WEBOver","QWEB","USER") )
   sigtest = match.arg( sigtest , c("LRT","LRT2","WALD","WALDt","F","ALL") )

   groups <- as.factor(groups)
   J = length(unique(groups))
   G = nrow(y)
   Mj <- matrix( table(groups)[unique(as.character(groups))] , nr = G , nc = J )
   if( any(n==0) )
    {
      wh.0 <- which(apply(n,1,FUN=function(vec)any(vec==0)))
      Mj[wh.0,] = t(apply(n[wh.0,,drop=FALSE],1,FUN=function(vec) table(groups[which(vec!=0)])[unique(as.character(groups))]))
    }
   rownames(Mj) <- rownames(y)
   term1 <- y*log(y)
   term1[is.na(term1)] <- 0
   term2 <- (n-y)*log(n-y)
   term2[is.na(term2)] <- 0
   term3 <- n*log(n)
   term3[is.na(term3)] <- 0
   termsum1 <- rowSums( term1 + term2 - term3 )

   if( !bc ) nsmp <- rowSums(Mj)
   else
    {
      prop <- y/n
      nsmp <- apply(prop,1,FUN=function(vec) sum(vec!=1 & vec!=0,na.rm=TRUE))
      nsmp[nsmp<=J] <- J+1
    }

   if( type == "USER" & is.null(pq) ) pq = DBEstimateEBDisp(y=y,n=n,groups=groups,method="MLE",median=median,nsmp=NULL,return.pq=TRUE)

   if( is.null( dispersion ) & is.null(maxDisp) )
      D = switch( type , "common" = DBEstimateCommonDisp(y=y,n=n,groups=groups,names=names,approx="all.1"),
                         "tagwise" = DBEstimateEBDisp(y=y,n=n,groups=groups,method="WL",median=median,nsmp=nsmp,delta=delta) ,
                         "EB" = DBEstimateEBDisp(y=y,n=n,groups=groups,method="MLE",median=median,nsmp=nsmp) ,
                         "USER" = DBEstimateEBDisp(y=y,n=n,groups=groups,method="MLE",median=median,nsmp=nsmp,pq=pq) ,
                         "EBMM" = DBEstimateEBDisp(y=y,n=n,groups=groups,method="MM",median=median,nsmp=nsmp) ,
                         "EB3" = DBEstimateEBDisp(y=y,n=n,groups=groups,method="P3",median=median,nsmp=nsmp) ,
                         "WEB" = DBEstimateEBDisp(y=y,n=n,groups=groups,method="P1",median=median,nsmp=nsmp) ,
                         "QWEB" = DBEstimateEBDisp(y=y,n=n,groups=NULL,method="P1",median=median,nsmp=nsmp) ,
                         "ACML" = DBEstimateACMLDisp(y=y,n=n,groups=groups,names=names,approx="all.1") ,
                         "joint" = DBEstimateJointDisp(y=y,n=n,groups=groups,names=names) ,
                         "MAX" = pmin(DBEstimateACMLDisp(y=y,n=n,groups=groups,names=names,approx="all.1"),
                                      DBEstimateEBDisp(y=y,n=n,groups=groups,method="P1",median=median) ) )
    else D = dispersion                              

  if( fitDNull & is.null(maxDisp) & is.null(nulldispersion) ) DNull = switch( type , "common" = DBEstimateCommonDisp(y=y,n=n,groups=NULL,names=names,approx="all.1"),
                     "tagwise" = DBEstimateTagwiseDisp(y=y,n=n,groups=NULL,names=names,approx="all.1",prior.df=prior.df) ,
                     "EB" = DBEstimateEBDisp(y=y,n=n,groups=NULL,method="MLE",median=median,nsmp=nsmp) ,
                         "USER" = DBEstimateEBDisp(y=y,n=n,groups=groups,method="MLE",median=median,nsmp=nsmp,pq=pq) ,
                     "EBMM" = DBEstimateEBDisp(y=y,n=n,groups=NULL,method="MM",median=median,nsmp=nsmp) ,
                     "EB3" = DBEstimateEBDisp(y=y,n=n,groups=NULL,method="P3",median=median,nsmp=nsmp) ,
                     "WEB" = DBEstimateEBDisp(y=y,n=n,groups=NULL,method="P1",median=median,nsmp=nsmp) ,
                     "QWEB" = DBEstimateEBDisp(y=y,n=n,groups=NULL,method="P1",median=median,nsmp=nsmp) ,
                     "ACML" = DBEstimateACMLDisp(y=y,n=n,groups=NULL,names=names,approx="all.1") ,
                     "joint" = DBEstimateJointDisp(y=y,n=n,groups=NULL,names=names) ,
                     "MAX" = pmin(DBEstimateACMLDisp(y=y,n=n,groups=groups,names=names,approx="all.1"),
                                  DBEstimateEBDisp(y=y,n=n,groups=groups,method="P1",median=median) ) )
    else if( !is.null(nulldispersion) ) DNull = nulldispersion 
      else DNull = D
   D[D==500] <- 500 - sqrt(.Machine$double.eps)
   DNull[DNull==500] <- 500 - sqrt(.Machine$double.eps)
   D[D==0] <- sqrt(.Machine$double.eps)
   DNull[DNull==0] <- sqrt(.Machine$double.eps)

   if( length(D) == nrow(y) ) names(D) <- tagnames
   if( length(DNull) == nrow(y) ) names(DNull) <- tagnames
   
   if( !underdispersion ) D[D>1] <- 1
   if( !underdispersion ) DNull[DNull>1] <- 1

   if( return.disp )
    {
      x <- cbind( D , DNull )
      colnames(x) <- c("D","DNull")
      return(x)
    }
   yy = DBSplitIntoGroups( y , groups )[unique(as.character(groups))]
   nn = DBSplitIntoGroups( n , groups )[unique(as.character(groups))]
   ysums <- sapply( yy , rowSums )
   nsums <- sapply( nn , rowSums )
   phat <- mapply( FUN = function( cts , offs ) rowSums(cts)/rowSums(offs) , yy , nn , SIMPLIFY = TRUE )
   phat <- matrix( phat[,as.character(groups)] , nr = G , nc = length(groups) )
   lphat = log(phat) - log(1-phat)       
   nysums <- nsums - ysums
   nysums[nysums==0] <- sqrt(.Machine$double.eps)
   ysums[ysums==0] <- sqrt(.Machine$double.eps)
   B = log(ysums) - log(nysums)
   B[,2] <- B[,2] - B[,1]
   if( type != "joint" ) 
    {
      V = 1/D*exp( log(nsums) - log(ysums) - log(nysums) ) 
      V[,2] = rowSums(V)
    }
   else 
    {
      phat <- mapply( FUN = function( cts , offs ) rowSums(cts)/rowSums(offs) , yy , nn , SIMPLIFY = TRUE )
      phat <- matrix( phat[,as.character(groups)] , nr = G , nc = length(groups) )       
      termsum2 <- rowSums(y*log( ifelse(phat==0,1,phat) ) + (n-y)*log( ifelse(phat==1,1,1-phat) ))
      denom = termsum1 - termsum2 
      
      p1 = logitF( rowSums(B) , inverse = TRUE )
      p0 = logitF( B[,1] , inverse = TRUE )
      dldb1 = D*(ysums[,2] - nsums[,2]*p1)
      dldb0 = dldb1 + D*(ysums[,1] - nsums[,1]*p0)
      d2ldb12 = -D*p1*(1-p1)*nsums[,2]
      d2ldb02 = d2ldb12 - D*p0*(1-p0)*nsums[,1]
      d2ldt2 = -0.5*rowSums(Mj)*D/500*(500-D)/500 + 500*(-denom)*D/500*(500-D)/500*((500-D)/500-D/500)
      d2ldb0db1 = d2ldb12
      d2ldb1dt = -denom*(500-D)/500*dldb1     
      d2ldb0dt = -denom*(500-D)/500*dldb0
      
      det1 = d2ldb12 * d2ldt2 - d2ldb1dt * d2ldb1dt
      det5 = d2ldb02 * d2ldt2 - d2ldb0dt * d2ldb0dt
      det9 = d2ldb02 * d2ldb12 - d2ldb0db1 * d2ldb0db1
      det0 = d2ldb02*det1 - d2ldb12*det5 + d2ldt2*det9
      V = -cbind( 1/det0 * det1 , 1/det0 * det5 , 1/det0 * det9 )
    }
   if( type == "joint" ) 
    {
      B = cbind( B , log( D ) - log( 500 - D ) )
      colnames(B) <- c( "Intercept" , "Group" , "Intercept:LogitDisp:Max500"  )
      colnames(V) <- c( "Intercept" , "Group" , "Intercept:LogitDisp:Max500"  )
    }  
   else 
    {
      colnames(B) <- c( "Intercept" , "Group" )
      colnames(V) <- c( "Intercept" , "Group" )
    }
    
   Z <- B/sqrt(V)
   F <- matrix(c(rep(logitF(B[,1],inverse=TRUE),table(as.character(groups))[unique(as.character(groups))][1]),rep(logitF(B[,1]+B[,2],inverse=TRUE),table(as.character(groups))[unique(as.character(groups))][2])),nr=nrow(y),nc=ncol(y))
   rownames(F)<-rownames(y)
   R <- (y/n - F)
   BNull = log(rowSums(y)) - log(pmax(sqrt(.Machine$double.eps),rowSums(n)-rowSums(y)))
   llAlt <- 0.5*log(D)*rowSums(Mj) - D*termsum1 + D*B[,1]*rowSums(y)+D*B[,2]*ysums[,2]-D*log(1+exp(B[,1]))*nsums[,1]-D*log(1+exp(B[,1]+B[,2]))*nsums[,2]
   llNull <- 0.5*log(DNull)*rowSums(Mj) - DNull*termsum1 + DNull*rowSums(y)*BNull - DNull*log(1+exp(BNull))*rowSums(n)
   llNull2 <- 0.5*log(D)*rowSums(Mj) - D*termsum1 + D*rowSums(y)*BNull - D*log(1+exp(BNull))*rowSums(n)
   LRT = -2*(llNull - llAlt)
   LRT2 = -2*(llNull2 - llAlt)

#   LRT[LRT<=0] = LRT2[LRT<=0] #doesn't have an effect since if LRT is negative then likely LRT2 is near 0 so pvalue is never that small to begin with.

   z <- y/n
   z[z==0] <- sqrt(.Machine$double.eps)
   z[z==1] <- 1-sqrt(.Machine$double.eps)
   phat[phat==0] <- sqrt(.Machine$double.eps)
   phat[phat==1] <- 1 - sqrt(.Machine$double.eps)
   gammahat <- logitF(phat)
   DEV.F <- rowSums(2*n*( z*(logitF(z)-gammahat) + log(1-z) - log(1-phat)),na.rm=TRUE)
   DEV.F.scaled <- DEV.F * D #times since D here is inverse of the standard dispersion parameterization

#DEV.F.scaled = DEV.F * 1/(Dev.F/n)=n
#DEV.R.scaled = DEV.R * 1/(DEV.F/n) = n*DEV.R/DEV.F
#DEV.R.scaled - DEV.F.scaled = n*DEV.R/DEV.F - n = n(DEV.R/DEV.F-1) = LRT2 = 2*n*(S.R/S.F - 0.5)
#DEV.R.scaled - DEV.F.scaled = DEV.R*DNULL - DEV.F*D = DEV.R*n/DEV.R - DEV.F*n/DEV.F = n - n = 0.

   gammahatNull <- BNull
   phatNull <- matrix( logitF(BNull,inverse=TRUE) , nr = G , nc = length(groups) )
   phatNull[phatNull==0] <- sqrt(.Machine$double.eps)
   phatNull[phatNull==1] <- 1 - sqrt(.Machine$double.eps)
   DEV.R <- rowSums(2*n*( z*(logitF(z)-gammahatNull) + log(1-z) - log(1-phatNull) ),na.rm=TRUE)
   DEV.R.scaled <- DEV.R * D #times since D here is inverse of the standard dispersion parameterization


   #If you do joint MLE estimates then 1/D = 1/n * DEV.F & 1/DNull = 1/n * DEV.R.
   #If you do ACML MLE estimates then 1/D = 1/(n-K) * DEV.F & 1/DNull = 1/(n-K) * DEV.R.
   FSTAT = (DEV.R.scaled - DEV.F.scaled) / (J-1)

      P <- switch( sigtest , "LRT" = pchisq(LRT,df=J-1,lower.tail=FALSE) , "LRT2" = pchisq(LRT2,df=J-1,lower.tail=FALSE) ,
                             "WALD" = 2*pnorm(abs(Z[,2]),lower.tail=FALSE) , "WALDt" = 2*pt(abs(Z[,2]),df=rowSums(Mj)-J,lower.tail=FALSE) ,
                             "F" = pf(FSTAT,df1=J-1,df2=rowSums(Mj)-J,lower.tail=FALSE) ,
                             "ALL" = cbind( "LRT" = pchisq(LRT,df=J-1,lower.tail=FALSE) ,
                                            "LRT2" = pchisq(LRT2,df=J-1,lower.tail=FALSE) ,
                                            "WALD" = 2*pnorm(abs(Z[,2]),lower.tail=FALSE) ,
                                            "WALDt" = 2*pt(abs(Z[,2]),df=rowSums(Mj)-J,lower.tail=FALSE) ,
                                            "F" = pf(FSTAT,df1=J-1,df2=rowSums(Mj)-J,lower.tail=FALSE) ) )
   
     if( combotest )
      {
        if( !ttest ) pwald <- 2*pnorm(abs(Z[,2]),lower.tail=FALSE)
         else pwald <- 2*pt(abs(Z[,2]),df=rowSums(Mj)-J,lower.tail=FALSE) #WALDt
        P <- pchisq(LRT,df=J-1,lower.tail=FALSE)
        P = ifelse( wh.boundary , pwald , P )
      }
     if( trend )
      {
        boundary <- rowMeans(y/n,na.rm=TRUE)
        wh.boundary <- boundary<=comboval | boundary>=1-comboval
        y.boundary <- y[wh.boundary,]
        n.boundary <- n[wh.boundary,]
        boundary.pvals <- DBGLM1( y.boundary , n.boundary , groups , type="WEB" , underdispersion=FALSE )$lrtpvals
      }
     if( trend ) P[wh.boundary] <- boundary.pvals

#   sum(pF1<=0.05)
#   nms1 <- names(pF1)[pF1<=0.05]
#   table(sapply(strsplit(nms1,"_",fixed=TRUE),"[",1))

#   llAlt <- 0.5*log(D)*rowSums(Mj) - D*termsum1 + D*B[,1]*rowSums(y)+D*B[,2]*ysums[,2]-D*log(1+exp(B[,1]))*nsums[,1]-D*log(1+exp(B[,1]+B[,2]))*nsums[,2]
#   llNull <- 0.5*log(D)*rowSums(Mj) - D*termsum1 + D*rowSums(y)*BNull - D*log(1+exp(BNull))*rowSums(n)
#   LRT = -2*(llNull - llAlt)
#   pLRT = p.adjust(pchisq(LRT,df=1,lower.tail=FALSE),method="BH")
#   sum(pLRT<=0.05)
#   nmsLRT <- names(pLRT)[pLRT<=0.05]
#   table(sapply(strsplit(nmsLRT,"_",fixed=TRUE),"[",1))
#
#
#   z <- y/n
#   z[z==0] <- sqrt(.Machine$double.eps)
#   z[z==1] <- 1-sqrt(.Machine$double.eps)
#   phat[phat==0] <- sqrt(.Machine$double.eps)
#   phat[phat==1] <- 1 - sqrt(.Machine$double.eps)
#   gammahat <- logit(phat)
#   lambdahat <- gammahat * D
#   liks.F <- lambdahat*y + D*n*(-log(1-z)-z*logitF(z))-D*n*log(1+exp(lambdahat/D)) + 0.5*log(D)
#
#   gammahatNull <- BNull
#   lambdahatNull <- gammahatNull*DNull
#   liks.R <- lambdahatNull*y + DNull*n*(-log(1-z)-z*logitF(z))-DNull*n*log(1+exp(lambdahatNull/DNull)) + 0.5*log(DNull)
#
#   wN2 = apply(liks.F - liks.R,2,var)
#
#   llAlt <- 0.5*log(D)*rowSums(Mj) - D*termsum1 + D*B[,1]*rowSums(y)+D*B[,2]*ysums[,2]-D*log(1+exp(B[,1]))*nsums[,1]-D*log(1+exp(B[,1]+B[,2]))*nsums[,2]
#   llNull <- 0.5*log(DNull)*rowSums(Mj) - DNull*termsum1 + DNull*rowSums(y)*BNull - DNull*log(1+exp(BNull))*rowSums(n)
#   test = ( llAlt - llNull - (2-1)/2 * log(ncol(y)) ) / (sqrt(ncol(y)*wN2))
#   ptest <- p.adjust(2*pnorm(abs(test),lower.tail=FALSE),method="BH")
#   sum(ptest<=0.05)
#   nms <- names(ptest)[ptest<=0.05]
#   table(sapply(strsplit(nms,"_",fixed=TRUE),"[",1))

   if( type == "joint" ) 
    {
      B = B[ , c(1,3,2) ]
      V = V[ , c(1,3,2) ]
    }  
   if( type == "QWEB" )
    {
      des <- model.matrix(~groups)
      yy <- split( y , 1:nrow(y))
      mm <- split( m , 1:nrow(m))
      lrtpvals <- mapply( FUN=function(cts,offs,d) summary( glm(cbind(cts,offs-cts)~des-1,family="quasibinomial") , dispersion=d ) , yy , mm , D )
      mclist <- mapply( list , yy , mm , D , SIMPLIFY=FALSE )
      lrtpvals <- mapply( FUN=function(cts,offs,d) summary( glm(cbind(cts,offs-cts)~des-1,family="quasibinomial") , dispersion=d ) , yy[1:1000] , mm[1:1000] , D[1:1000] )
      lrtpvals <- mclapply( mclist , FUN=function(data) summary( glm(cbind(data[[1]],data[[2]]-data[[1]])~des-1,family="quasibinomial") , dispersion=1/data[[3]] )$coeff[2,4] , mc.preschedule=TRUE , mc.cores = 15 )
      lrtpvals2 <- mclapply( mclist , FUN=function(data){ fit <- glm(cbind(data[[1]],data[[2]]-data[[1]])~des-1,family="quasibinomial"); anova( fit , dispersion=1/data[[3]] ,test="Rao")[6][2,1] }, mc.preschedule=TRUE , mc.cores = 15 )
      lrtpvals2 <- mclapply( mclist , FUN=function(data){ fit <- glm(cbind(data[[1]],data[[2]]-data[[1]])~des-1,family="quasibinomial"); anova( fit , dispersion=summary(fit)$disp ,test="Rao")[6][2,1] }, mc.preschedule=TRUE , mc.cores = 15 )
      lrtpvals2 <- mclapply( mclist , FUN=function(data){ fit <- glm(cbind(data[[1]],data[[2]]-data[[1]])~des-1,family="quasibinomial"); anova( fit , dispersion=summary(fit)$disp ,test="LRT")[5][2,1] }, mc.preschedule=TRUE , mc.cores = 15 )
      lrtpvals2 <- mclapply( mclist , FUN=function(data){ fit <- glm(cbind(data[[1]],data[[2]]-data[[1]])~des-1,family="quasibinomial"); anova( fit , dispersion=1/data[[3]],test="Chisq")[5][2,1] }, mc.preschedule=TRUE , mc.cores = 15 )
    }
   out <- list( B , D , DNull , V , F , R , P , Z , LRT , FSTAT , LRT2 )
   names( out ) <- c( "betas" , "disps" , "nulldisps" , "vars" , "fits" , "resids" , "lrtpvals" , "zvals" , "lrtstat" , "fstat" , "lrt2stat" )
   return(out)
 }

    


DBEstimateCommonDisp <- function (y,n, groups, names = NULL , tol = 1e-04, approx = c( "all.exact" , "all.approx" , "all.1" , "prop.approx" , "prop.1" ) , ncore = 1 , cutoff = Inf , reparam = c("0-1","log") )
 {                                                                                              
    approx = match.arg( approx , c( "all.exact" , "all.approx" , "all.1" , "prop.approx" , "prop.1" ) )
    reparam = match.arg( reparam , c( "0-1" , "log" ) )
    
     if( is.null(groups) ) groups = rep( 1 , ncol(y) )
     groups <- as.factor(groups)

    if( approx == "all.1" )
     {
       tagnames <- rownames(y)
       J = length(unique(as.character(groups)))
       G = nrow(y)
       Mj <- matrix( table(groups)[unique(as.character(groups))] , nr = G , nc = J )          
       if( any(n==0) ) 
        {
          wh.0 <- which(apply(n,1,FUN=function(vec)any(vec==0)))
          Mj[wh.0,] = t(apply(n[wh.0,,drop=FALSE],1,FUN=function(vec) table(groups[which(vec!=0)])[unique(as.character(groups))]))
        }  
   prop <- y/n
   nsmp <- apply(prop,1,FUN=function(vec) sum(vec!=1 & vec!=0,na.rm=TRUE))
   nsmp[nsmp<=J] <- J+1

       term1 <- y*log(y)
       term1[is.na(term1)] <- 0 
       term2 <- (n-y)*log(n-y)
       term2[is.na(term2)] <- 0 
       term3 <- n*log(n)
       term3[is.na(term3)] <- 0                                           
       termsum1 <- sum( term1 + (term2 - term3) )

       yy = DBSplitIntoGroups( y , groups )[ as.character(unique(groups)) ]
       nn = DBSplitIntoGroups( n , groups )[ as.character(unique(groups)) ]       

       phat <- mapply( FUN = function( cts , offs ) rowSums(cts)/rowSums(offs) , yy , nn , SIMPLIFY = TRUE )
       phat <- matrix( phat[,as.character(groups)] , nr = G , nc = length(groups) )
       
       termsum2 <- sum(y*log( ifelse(phat==0,1,phat) ) + (n-y)*log( ifelse(phat==1,1,1-phat) ))
       
       denom = termsum1 - termsum2

       common.disp = 0.5*(sum(nsmp)-G*J)/denom
        names(common.disp) <- NULL
     }
    else
     {  
       delta <- switch( reparam , "0-1" = optimize(DBCombCondLogLik,interval=c(1e-04, 100/(100 + 1)),y=y,n=n,groups=groups,tol=tol,maximum=TRUE,approx=approx,ncore=ncore,cutoff=cutoff,reparam=reparam,names=names)$maximum ,
                                  "log" = optimize(DBCombCondLogLik,interval=c(log(1e-04),log(100)),y=y,n=n,groups=groups,tol=tol,maximum=TRUE,approx=approx,ncore=ncore,cutoff=cutoff,reparam=reparam,names=names)$maximum)
       common.disp = switch( reparam , "0-1" = delta/(1 - delta) , "log" = exp(delta) )
        names(common.disp) <- NULL
     }  
    return( common.disp )
 }

DBSplitIntoGroups <- function (y, groups)
  {
    if( is.null(groups) ) groups <- rep(1,ncol(y))
    groups <- as.factor( groups )
    if( is.null(nrow(y)) ) nrows = 1
     else nrows <- nrow(y)
    yy <- lapply(split(t(y), groups), FUN = function(u) matrix(u, nrow = nrows, byrow = TRUE))[ unique(as.character(groups)) ]
    for (i in 1:length(unique(groups))) rownames(yy[[i]]) <- rownames(y)
    return( yy )
  }

DBCombCondLogLik <- function (y,n,groups=NULL, names=NULL , delta, approx = c( "all.exact" , "all.approx" , "all.1" , "prop.approx" , "prop.1" ) , ncore = 1 , doSum = TRUE , cutoff = Inf , reparam = c("0-1","log") )
  {
    approx = match.arg( approx , c( "all.exact" , "all.approx" , "all.1" , "prop.approx" , "prop.1" ) )
    reparam = match.arg( reparam , c( "0-1" , "log" ) )

    if( is.null(groups) ) groups <- rep(1,ncol(y))
      groups <- as.factor(groups)
      tagnames <- rownames(y)
      yy = DBSplitIntoGroups( y , groups )
      nn = DBSplitIntoGroups( n , groups )
      Y = mapply( FUN = function( cts , offs )
                         {
                           list = lapply( 1:nrow(cts) , FUN = function( ii ) cbind( S=cts[ii,],F=offs[ii,]-cts[ii,] ) )
                           names( list ) <- tagnames
                           return( list )
                         } , yy , nn , SIMPLIFY = FALSE )

    phis = switch( reparam , "0-1" = delta/(1-delta) , "log" = exp(delta) )

    m0 = matrix(0,nrow=ifelse(doSum,1,nrow(y)),ncol=length(phis))
    for( j in 1:length(phis) )
     {
        L0 = 0
        phi = phis[j]
        for( i in 1:length(Y) ) #number of groups
         {
            YY = Y[[i]]
            if( ncore == 1 )
             {
               l0 <- lapply( YY , FUN = function( mtx ) DBCondLogLik( phi=phi , y=mtx[,1] , n=rowSums(mtx) , approx = approx , cutoff = cutoff , reparam="none" ) )
             }
            else
             {                                 
               l0 <- mclapply( YY , FUN = function( mtx ) DBCondLogLik( phi=phi , y=mtx[,1] , n=rowSums(mtx) , approx = approx , cutoff = cutoff , reparam="none" ) , mc.preschedule = TRUE , mc.cores = ncore )
             }
            L0 <- do.call( "rbind" , l0 ) + L0
         }
        if( any( is.na( L0 ) ) ) L0 <- apply( L0 , 2 , FUN = function( col ){ col[ is.na( col ) ] <- mean( col , na.rm = TRUE ) ; return( col ) } )
        m0[,j] <- matrix(if(doSum) colSums(L0) else L0, nrow = ifelse(doSum,1,nrow(y)), ncol = 1)
     }
    return( m0 )
  }


DBEstimateTagwiseDisp <-  function (y,n,groups=NULL , names=NULL , common.dispersion=NULL ,  tol = 1e-04 , approx = c( "all.exact" , "all.approx" , "all.1" , "prop.approx" , "prop.1" ) , ncore = 1 , prior.df = 20 , grid.length = 11, grid.range = c(-6, 6), cutoff = Inf )
  {
    approx = match.arg( approx , c( "all.exact" , "all.approx" , "all.1" , "prop.approx" , "prop.1" ) )

    if( is.null(groups) ) groups <- rep(1,ncol(y))
      groups <- as.factor(groups)

    J = length(unique(groups))
    delta = prior.df / (ncol(y) - J)
    
    ntags = nrow(y)

    if( approx == "all.1" )
     {
       tagnames <- rownames(y)
       G = nrow(y)
       Mj <- matrix( table(groups)[unique(as.character(groups))] , nr = G , nc = J )          
       if( any(n==0) ) 
        {
          wh.0 <- which(apply(n,1,FUN=function(vec)any(vec==0)))
          Mj[wh.0,] = t(apply(n[wh.0,,drop=FALSE],1,FUN=function(vec) table(groups[which(vec!=0)])[unique(as.character(groups))]))
        }  
       
       term1 <- y*log(y)
       term1[is.na(term1)] <- 0 
       term2 <- (n-y)*log(n-y)
       term2[is.na(term2)] <- 0 
       term3 <- n*log(n)
       term3[is.na(term3)] <- 0 

       termsum1 <- rowSums( term1 + term2 - term3 )

       yy = DBSplitIntoGroups( y , groups )[ as.character(unique(groups)) ]
       nn = DBSplitIntoGroups( n , groups )[ as.character(unique(groups)) ]       

       phat <- mapply( FUN = function( cts , offs ) rowSums(cts)/rowSums(offs) , yy , nn , SIMPLIFY = TRUE )
       phat <- matrix( phat[,as.character(groups)] , nr = G , nc = length(groups) )
       
       termsum2 <- rowSums(y*log( ifelse(phat==0,1,phat) ) + (n-y)*log( ifelse(phat==1,1,1-phat) ))
       
       denom = termsum1 - termsum2

       tagwise.dispersion = 0.5*(rowSums(Mj)-J + (sum(Mj)/G-J)*delta) / (denom + delta*mean(denom))
       tagwise.dispersion[ tagwise.dispersion > 500 ] <- 500
       tagwise.dispersion[ tagwise.dispersion < sqrt(.Machine$double.eps) ] <- sqrt(.Machine$double.eps)
     }
    else
     {
       if( is.null(common.dispersion) ) stop( "Must give common dispersion value" )
       spline.pts <- seq(from = grid.range[1], to = grid.range[2], length = grid.length)
       spline.disp <- common.dispersion * 2^spline.pts #the common dispersion is a good starting place to draw a neighborhood in which the max of the WL would be found. It doesn't need to be super accurate
       grid.vals <- spline.disp/(1 + spline.disp)
     
       l0 <- DBCombCondLogLik(y=y,n=n,groups=groups,delta=grid.vals,ncore=ncore,approx=approx,doSum=FALSE,cutoff=cutoff,reparam="0-1")
  
       m0 <- matrix(colSums(l0), ntags, grid.length, byrow = TRUE)
         
       l0a <- l0 + delta/ntags * m0
  
       d <- rep(0, ntags)
       for (j in 1:ntags) d[j] <- maximizeInterpolant(spline.pts,l0a[j, ])
       tagwise.dispersion <- common.dispersion * 2^d
       tagwise.dispersion[ tagwise.dispersion > 500 ] <- 500
       tagwise.dispersion[ tagwise.dispersion < sqrt(.Machine$double.eps) ] <- sqrt(.Machine$double.eps)
     }  
    names(tagwise.dispersion) <- NULL
    return( tagwise.dispersion )
  }

rho <- function( vec )
        {
          out = -log(1-vec) - vec*(log(vec)-log(1-vec)) 
          out[is.na(out)] = 0
          return(out)
        }

#DBZ <- function( Y , m , groups=NULL )
#           {
#              if( is.null(groups) ) groups <- rep(1,ncol(Y))
#                groups <- as.factor(groups)
#             if( !is.null(dim(Y)) ) if( dim(Y)[2] != 1 & dim(Y)[1] != 1) stop( "Count data 'Y' should be a vector or 1-column (or 1-row) matrix" )
#               else Y<-as.vector(Y)
#             if( !is.null(dim(m)) ) if( dim(m)[2] != 1 & dim(m)[1] != 1) stop( "Count data 'm' should be a vector or 1-column (or 1-row) matrix" ) else m<-as.vector(m)
#             if( any(Y<0) ) stop( "Count data has negative values" )
#             if( !all(findInterval(Y,c(-1e-8,1+1e-8))==1) ) Y <- Y/m
#
#             if( any(m==0) ) 
#              {
#                Y<-Y[-which(m==0)]
#                m<-m[-which(m==0)]
#              }  
#
#
#             M = sum(m)
#             stats <- rho(1/M*sum(m*Y)) - 1/M*sum(m*rho(Y)) 
#             Z <- M * stats
#
#             return( Z )
#           }

DBS <- function( y , n , groups=NULL )
           {
             if( is.null(groups) ) groups <- rep(1,ncol(y))
               groups <- as.factor(groups)
  
             if( any(y<0) ) stop( "Count data has negative values" )

             z <- y/n
             z[is.na(z)]<-0
             zz = DBSplitIntoGroups( z , groups )[ as.character(unique(groups)) ]
             nn = DBSplitIntoGroups( n , groups )[ as.character(unique(groups)) ]       
             
             N <- rowSums(n) 
             U <- (1/N) * rowSums(n*rho(z))
             Nk = sapply( nn , rowSums )
             Tk <- sapply( mapply( "*" , zz , nn , SIMPLIFY = FALSE ), rowSums ) / Nk
             R <- (1/N) * rowSums(Nk*rho(Tk))
             
             S <- N*(R - U)

             S = pmax(S,sqrt(.Machine$double.eps))
             
             return( S )
           }


nll.gbp <- function( pars , a , S ) 
           {
                 b = pars[1]
                 q = pars[2]
                 N = length(S)
                 S1 = log( S )
                 S2 = log( S + q )
                 lB = lbeta( a , b )
                 llval = sum((a-1)*S1) - sum((a+b)*S2) + N*b*log(q) - sum(lB)
                 return( -llval )
           }
nll.gbp2 <- function( pars , a , S )
           {
             nllval <- nll.gbp( pars , a , S )
             attr(nllval, "gradient") <- gr.nlgbp( pars = pars , a = a, S = S )
             attr(nllval, "hessian") <- Jcb( pars = pars , a = a, S = S )
             return( nllval )
           }

gr.nlgbp <- function( pars , a , S )
           {
             require( numDeriv , quietly = TRUE )
             b = pars[1]
             q = pars[2]
             N = length(S)
             S2 = log( S + q )
             S3 = 1/(S+q)
             dg = digamma( b ) - digamma( a + b )

             derivs = c( -sum(S2) + N*log(q) - sum(dg) , sum(-(a+b)*S3) + N*b/q )

             gr = c( derivs[1] , derivs[2] )
             return( -gr )
           }

Jcb <- function( pars , a , S )
        {
          require( numDeriv , quietly = TRUE )
          b = pars[1]
          q = pars[2]
          N = length(S)
          S2 = log( S + q )
          S3 = 1/(S+q)
          S4 = 1/(S+q)^2
          tg = trigamma(b) - trigamma(a+b)

          mtx = matrix( c( -sum(tg) , -sum(S3) + N/q , -sum(S3) + N/q , sum((a+b)*S4) - N*b/q^2 ) , 2 , 2 )
          return( -mtx )
        }

nll.gbp.3 <- function( pars , a , S )
           {
             b = pars[1]
             q = pars[2]
             p = pars[3]
             N = length(S)
             tmp1 = log(S)
             S1 = sum( tmp1 )
             S2 = sum( log(1 + exp( p*(log(q)-tmp1) )) + p*tmp1 )
             lB = lbeta( a , b )
             llval = (a*p-1)*S1 - (a+b)*S2 + N*log(p) + N*p*b*log(q) - N*lB
             return( -llval )
           }
nll.gbp2.3 <- function( pars , a , S )
           {
             nllval <- nll.gbp.3( pars , a , S )
             attr(nllval, "gradient") <- gr.nlgbp.3( pars = pars , a = a, S = S )
             attr(nllval, "hessian") <- Jcb.3( pars = pars , a = a, S = S )
             return( nllval )
           }

gr.nlgbp.3 <- function( pars , a , S )
           {
             require( numDeriv , quietly = TRUE )
             b = pars[1]
             q = pars[2]
             p = pars[3]
             N = length(S)
             tmp1 = log(S)
             tmp2 = 1 + exp( p*(log(q)-tmp1) )
             tmp2b = log(tmp2)
             S1 = sum( tmp1 )
             S2 = sum( tmp2b )
             S3 = sum( exp( p*(log(q)-tmp1) - tmp2b ) )
             S4 = sum( tmp1/tmp2 )
             dg = digamma( b ) - digamma( a + b )

             derivs = c( -S2 - p*S1 + N*p*log(q) - N*dg , -((a+b)*p/q)*S3 + N*b*p/q , N/p + N*b*log(q) + a*S1 - (a+b)*log(q)*S3 - (a+b)*S4 )

             gr = c( derivs[1] , derivs[2] , derivs[3] )
             return( -gr )
           }

Jcb.3 <- function( pars , a , S )
        {
          require( numDeriv , quietly = TRUE )
          b = pars[1]
          q = pars[2]
          p = pars[3]
          N = length(S)

          tmp1 = log(S)
          tmp2 = 1 + exp( p*(log(q)-tmp1) )
          tmp2b = log(tmp2)
          S1 = sum( tmp1 )
          S2 = sum( tmp2b )
          S3 = sum( exp( p*(log(q)-tmp1) - tmp2b ) )
          S4 = sum( exp( log(tmp1) - tmp2b ) )
          tmp5 = exp( 2*tmp2b )
          S5 = sum( exp( 2*(p*(log(q)-tmp1) - tmp2b) ) )
          tmp6 = exp( p*(log(q)-tmp1) - 2*tmp2b )*(log(q)-tmp1)
          S6 = sum(tmp6)
          S7 = sum( tmp6*(log(q)-tmp1) )
          tg = trigamma(b) - trigamma(a+b)

          mtx = matrix( c( -N*tg , -(p/q)*S3 + N*p/q , N*log(q) - log(q)*S3 - S4 ,
                           -(p/q)*S3 + N*p/q , -((a+b)*p*(p-1)/q^2)*S3 + ((a+b)*(p/q)^2)*S5 - N*b*p/q^2 , N*b/q - ((a+b)/q)*(S3 + p*S6) ,
                           N*log(q) - log(q)*S3 - S4 , N*b/q - ((a+b)/q)*(S3 + p*S6) , -N/p^2 - (a+b)*S7 ) , 3 , 3 )
          return( -mtx )
        }

nll.gbp.delta <- function( delta , a , S )
           {
             b = delta*a
             m.S = tapply( S , names(S) , mean )
             m.S = m.S[match( names(S) , names(m.S) )]
             q = delta*m.S
             N = length(S)
             S1 = log( S )
             S2 = log( S + q )
             lB = lbeta( a , b )
             llval = sum((a-1)*S1 - (a+b)*S2 + b*log(q) - lB)
             return( -llval )
           }

nll.gbp.delta2 <- function( gamma , a , S )
           {
             delta = gamma / (1 - gamma)
             b = delta*a
             m.S = tapply( S , names(S) , mean )
             m.S = m.S[match( names(S) , names(m.S) )]
             q = delta*m.S
             N = length(S)
             S1 = log( S )
             S2 = log( S + q )
             lB = lbeta( a , b )
             llval = sum((a-1)*S1 - (a+b)*S2 + b*log(q) - lB)
             return( -llval )
           }

gr.nlgbp.delta <- function( delta , a , S )

           {
             b = delta*a
             m.S = tapply( S , names(S) , mean )
             m.S = m.S[match( names(S) , names(m.S) )]
             q = delta*m.S
             N = length(S)
             S2 = log( S + q )
             S3 = 1/(S+q)
             dg = digamma( b ) - digamma( a + b )

             gr = sum(-S2 - a*(1+delta)*m.S*S3 + a*(1+log(q)) - a*dg)

             return( -gr )
           }

mydexpbinomial <- function (lmean = "logit", ldispersion = "logit", idispersion = 0.25,
    zero = 2 , approx = c("1","approx","exact") )
  {

    lmean <- as.list(substitute(lmean))
    emean <- link2list(lmean)
    lmean <- attr(emean, "function.name")

    ldisp <- as.list(substitute(ldispersion))
    edisp <- link2list(ldisp)
    ldisp <- attr(edisp, "function.name")

    idisp <- idispersion

    approx = match.arg( approx , c("1","approx","exact") )

    if (!is.Numeric(idispersion, positive = TRUE)) stop("bad input for 'idispersion'")


    new( "vglmff" ,
         blurb = c( "Double Exponential Binomial distribution\n\n",
                    "Link:     ", namesof("mean", lmean, earg = emean), ", ",
                    namesof("dispersion", ldisp, earg = edisp), "\n", "Mean:     ",
                    "mean\n") ,
         constraints = eval(substitute(expression({
                                                   constraints = VGAM:::cm.zero.vgam(constraints, x, .zero, M)
                                                   }) ,
                                       list(.zero = zero))) ,

         initialize = eval(substitute(expression({
                                                   if (!all(w == 1)) extra$orig.w = w
                                                   if (ncol(cbind(w)) != 1) stop("'weights' must be a vector or a one-column matrix")
                                                    NCOL = function(x) if (is.array(x) && length(dim(x)) >
                                                        1 || is.data.frame(x)) ncol(x) else as.integer(1)
                                                    if (NCOL(y) == 1) {
                                                        if (is.factor(y)) y = (y != levels(y)[1])
                                                        nvec = rep(1, n)
                                                        y[w == 0] <- 0
                                                        if (!all(y == 0 || y == 1)) stop("response values 'y' must be 0 or 1")
                                                        init.mu = mustart = (0.5 + w * y)/(1 + w)
                                                        no.successes = y
                                                        if (min(y) < 0) stop("Negative data not allowed!")
                                                        if (any(abs(no.successes - round(no.successes)) >
                                                            1e-08)) stop("Number of successes must be integer-valued")
                                                    } else if (NCOL(y) == 2) {
                                                        if (min(y) < 0) stop("Negative data not allowed!")
                                                        if (any(abs(y - round(y)) > 1e-08)) stop("Count data must be integer-valued")
                                                        y = round(y)
                                                        nvec = y[, 1] + y[, 2]
                                                        y = ifelse(nvec > 0, y[, 1]/nvec, 0)
                                                        w = w * nvec
                                                        init.mu = mustart = (0.5 + nvec * y)/(1 + nvec)
                                                    } else stop("for the dexpbinomial family, response 'y' must be a ",
                                                        "vector of 0 and 1's\n", "or a factor (first level = fail, ",
                                                        "other levels = success),\n", "or a 2-column matrix where col 1 is the no. of ",
                                                        "successes and col 2 is the no. of failures")
                                                    dn2 = if (is.matrix(y)) dimnames(y)[[2]] else NULL
                                                    dn2 = if (length(dn2)) {
                                                        paste("E[", dn2, "]", sep = "")
                                                    } else {
                                                        "mu"
                                                    }
                                                        predictors.names <- c(namesof(dn2, .lmean, earg = .emean,
                                                            short = TRUE), namesof("dispersion", .ldisp, earg = .edisp,
                                                            short = TRUE))
                                                        tmp2 <- rep(.idisp, len = n)
                                                        if (!length(etastart)) etastart = cbind(theta2eta(init.mu,
                                                            .lmean, earg = .emean), theta2eta(tmp2, .ldisp, earg = .edisp))
                                                    }),
                                      list(.lmean = lmean, .emean = emean, .ldisp = ldisp,
                                           .edisp = edisp, .idisp = idisp))) ,

         linkinv = eval(substitute(function(eta, extra = NULL) {
                                     eta2theta(eta[, 1], link = .lmean, earg = .emean)
                                      }, list(.lmean = lmean, .emean = emean, .ldisp = ldisp, .edisp = edisp))),
         last = eval(substitute(expression({
                                              misc$expected <- TRUE
                                              misc$link <- c(mean = .lmean, dispersion = .ldisp)
                                              misc$earg <- list(mean = .emean, dispersion = .edisp)
                                              misc$approx <- approx
                                              misc$nvec <- nvec
                                              misc$cutoff <- cutoff
                                          }),
                                list(.lmean = lmean, .emean = emean, .ldisp = ldisp, .edisp = edisp))),


         loglikelihood = eval(substitute(function(mu, Disper=NULL , y, w, nvec, cutoff = Inf, residuals = FALSE, eta, extra = NULL , approx = c("1","approx","exact") ) {
                                            prob = eta2theta(eta[, 1], link = .lmean, earg = .emean)

                                            Disper = eta2theta(eta[, 2], link = .ldisp, earg = .edisp)

                                            approx = match.arg( approx , c("1","approx","exact") )

                                            if (residuals) stop("loglikelihood residuals ", "not implemented yet") else {
                                                temp1 = y * log(ifelse(y > 0, y, 1))
                                                temp2 = (1 - y) * log1p(ifelse(y < 1, -y, 0))

                                            normc = DBNormC( nvec , mu , Disper , approx = approx , cutoff = cutoff )

                                             sum( - log(normc) + w * (y * Disper * log(prob) + (1 - y) * Disper * log1p(-prob) + temp1 * (1 - Disper) + temp2 * (1 - Disper)))
                                            }
                                          },
                              list(.lmean = lmean, .emean = emean, .ldisp = ldisp, .edisp = edisp))),
         vfamily = "dexpbinomial",

         deriv = eval(substitute(expression({
                                              prob = eta2theta(eta[, 1], link = .lmean, earg = .emean)
                                              Disper = eta2theta(eta[, 2], link = .ldisp, earg = .edisp)

                                              temp1 = y * log(ifelse(y > 0, y, 1))
                                              temp2 = (1 - y) * log1p(ifelse(y < 1, -y, 0))
                                              temp3 = prob * (1 - prob)
                                              temp3 = pmax(temp3, .Machine$double.eps * 10000)

                                              normc = DBNormC( nvec , mu , Disper , approx = approx , cutoff = cutoff )
                                              dxnormcdpx = DerivC( nvec , mu , Disper , approx = approx , cutoff = cutoff , der = "both" , wrt = "p" , normC = normc )
                                              dxnormcdphix = DerivC( nvec , mu , Disper , approx = approx , cutoff = cutoff , der = "both" , wrt = "phi" , normC = normc )
                                              
                                              dnormcdp <- dxnormcdpx[1,]
                                              d2normcdp2 <- dxnormcdpx[2,]

                                              dnormcdphi <- dxnormcdphix[1,]
                                              d2normcdphi2 <- dxnormcdphix[2,]

                                              dl.dprob = w * Disper * (y - prob)/temp3 - 1/normc * dnormcdp
                                              dl.dDisper = - 1/normc * dnormcdphi + w * (y * log(prob) + (1 - y) * log1p(-prob) - temp1 - temp2)

                                              dprob.deta = dtheta.deta(theta = prob, .lmean, earg = .emean)
                                              dDisper.deta = dtheta.deta(theta = Disper, .ldisp, earg = .edisp)

                                              cbind(dl.dprob * dprob.deta, dl.dDisper * dDisper.deta)
                                          }),
                                  list(.lmean = lmean, .emean = emean, .ldisp = ldisp, .edisp = edisp))),

         weight = eval(substitute(expression({
                                              wz = matrix(as.numeric(NA), nrow = n, ncol = 2)
                                              wz[, iam(1, 1, M)] = ( w * (Disper/temp3) + (normc*d2normcdp2-dnormcdp^2)/normc^2 )* dprob.deta^2
                                              wz[, iam(2, 2, M)] = ( (normc*d2normcdphi2-dnormcdphi^2)/normc^2 )* dDisper.deta^2
                                              wz
                                             }),
                                  list(.lmean = lmean, .emean = emean, .ldisp = ldisp, .edisp = edisp))))
  }


mydexpbinomial.profile <- function (lmean = "logit")
  {

    zero=NULL
    lmean <- as.list(substitute(lmean))
    emean <- link2list(lmean)
    lmean <- attr(emean, "function.name")

    new( "vglmff" ,
         blurb = c( "Double Exponential Binomial distribution\n\n",
                    "Link:     ", namesof("mean", lmean, earg = emean), ", ") ,

         constraints = eval(substitute(expression({ #what to do with zero and constraints. Where and how is it used? I think it checks M = length(eta) which = 1 here and zero must be <=1 or NULL. It seems it should be NULL.
                                                   constraints = VGAM:::cm.zero.vgam(constraints, x, .zero, M)
                                                   }) ,
                                       list(.zero = zero))) ,

         initialize = eval(substitute(expression({
                                                   if (!all(w == 1)) extra$orig.w = w
                                                   if (ncol(cbind(w)) != 1) stop("'weights' must be a vector or a one-column matrix")
                                                    NCOL = function(x) if (is.array(x) && length(dim(x)) >
                                                        1 || is.data.frame(x)) ncol(x) else as.integer(1)
                                                    if (NCOL(y) == 1) {
                                                        if (is.factor(y)) y = (y != levels(y)[1])
                                                        nvec = rep(1, n)
                                                        y[w == 0] <- 0
                                                        if (!all(y == 0 || y == 1)) stop("response values 'y' must be 0 or 1")
                                                        init.mu = mustart = (0.5 + w * y)/(1 + w)
                                                        no.successes = y
                                                        if (min(y) < 0) stop("Negative data not allowed!")
                                                        if (any(abs(no.successes - round(no.successes)) >
                                                            1e-08)) stop("Number of successes must be integer-valued")
                                                    } else if (NCOL(y) == 2) {
                                                        if (min(y) < 0) stop("Negative data not allowed!")
                                                        if (any(abs(y - round(y)) > 1e-08)) stop("Count data must be integer-valued")
                                                        y = round(y)
                                                        nvec = y[, 1] + y[, 2]
                                                        y = ifelse(nvec > 0, y[, 1]/nvec, 0)
                                                        w = w * nvec
                                                        init.mu = mustart = (0.5 + nvec * y)/(1 + nvec)
                                                    } else stop("for the dexpbinomial family, response 'y' must be a ",
                                                        "vector of 0 and 1's\n", "or a factor (first level = fail, ",
                                                        "other levels = success),\n", "or a 2-column matrix where col 1 is the no. of ",
                                                        "successes and col 2 is the no. of failures")
                                                    dn2 = if (is.matrix(y)) dimnames(y)[[2]] else NULL
                                                    dn2 = if (length(dn2)) {
                                                        paste("E[", dn2, "]", sep = "")
                                                    } else {
                                                        "mu"
                                                    }
                                                        predictors.names <- c(namesof(dn2, .lmean, earg = .emean,
                                                            short = TRUE))
                                                        if (!length(etastart)) etastart = theta2eta(init.mu,
                                                            .lmean, earg = .emean)
                                                    }),
                                      list(.lmean = lmean, .emean = emean))) ,

         linkinv = eval(substitute(function(eta, extra = NULL) {
                                     eta2theta(eta, link = .lmean, earg = .emean)
                                      }, list(.lmean = lmean, .emean = emean))),
         last = eval(substitute(expression({
                                              misc$expected <- TRUE
                                              misc$link <- c(mean = .lmean)
                                              misc$earg <- list(mean = .emean)
                                              misc$approx <- approx
                                              misc$nvec <- nvec
                                              misc$cutoff <- cutoff
                                              misc$Disper <- Disper
                                          }),
                                list(.lmean = lmean, .emean = emean))),


         loglikelihood = eval(substitute(function(mu, Disper , y, w, nvec, cutoff = Inf, residuals = FALSE, eta, extra = NULL , approx = c("1","approx","exact") ) {
                                            prob = eta2theta(eta, link = .lmean, earg = .emean)

                                            approx = match.arg( approx , c("1","approx","exact") )

                                            if (residuals) stop("loglikelihood residuals ", "not implemented yet") else {
                                                temp1 = y * log(ifelse(y > 0, y, 1))
                                                temp2 = (1 - y) * log1p(ifelse(y < 1, -y, 0))

                                            normc = DBNormC( nvec , mu , Disper , approx = approx , cutoff = cutoff )

                                             sum( -log(normc) + w * (y * Disper * log(prob) + (1 - y) * Disper * log1p(-prob) + temp1 * (1 - Disper) + temp2 * (1 - Disper)))
                                            }
                                          },
                              list(.lmean = lmean, .emean = emean ))),
         vfamily = "dexpbinomial",

         deriv = eval(substitute(expression({
                                              prob = eta2theta(eta, link = .lmean, earg = .emean)

                                              temp1 = y * log(ifelse(y > 0, y, 1))
                                              temp2 = (1 - y) * log1p(ifelse(y < 1, -y, 0))
                                              temp3 = prob * (1 - prob)
                                              temp3 = pmax(temp3, .Machine$double.eps * 10000)

                                              normc = DBNormC( nvec , mu , Disper , approx = approx , cutoff = cutoff )
                                              dxnormcdpx = DerivC( nvec , mu , Disper , approx = approx , cutoff = cutoff , der = "both" , wrt = "p" , normC = normc )

                                              dnormcdp <- dxnormcdpx[1,]
                                              d2normcdp2 <- dxnormcdpx[2,]

                                              dl.dprob = w * Disper * (y - prob)/temp3 - 1/normc * dnormcdp
                                              dprob.deta = dtheta.deta(theta = prob, .lmean, earg = .emean)

                                              cbind(dl.dprob * dprob.deta)
                                          }),
                                  list(.lmean = lmean, .emean = emean))),

         weight = eval(substitute(expression({
                                              wz = matrix(as.numeric(NA), nrow = n, ncol = 1)
                                              wz[, iam(1, 1, M)] = ( w * (Disper/temp3) + (normc*d2normcdp2-dnormcdp^2)/normc^2 )* dprob.deta^2 
                                              wz
                                             }),
                                  list(.lmean = lmean, .emean = emean))))
  }


myvglm <- function (formula, family, data = list(), nvec, cutoff=500 , Disper=NULL, approx = c("1","approx","exact") , weights = NULL, subset = NULL,
    na.action = na.fail, etastart = NULL, mustart = NULL, coefstart = NULL,
    control = vglm.control(...), offset = NULL, method = "vglm.fit",
    model = FALSE, x.arg = TRUE, y.arg = TRUE, contrasts = NULL,
    constraints = NULL, extra = list(), form2 = NULL, qr.arg = TRUE,
    smart = TRUE, ...)
{
    approx = match.arg( approx , c("1","approx","exact") )

    dataname <- as.character(substitute(data))
    function.name <- "vglm"
    ocall <- match.call()
    if (smart)
        setup.smart("write")
    if (missing(data))
        data <- environment(formula)
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "weights", "na.action",
        "etastart", "mustart", "offset"), names(mf), 0)
    mf <- mf[c(1, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    switch(method, model.frame = return(mf), myvglm.fit = 1, stop("invalid 'method': ",
        method))
    mt <- attr(mf, "terms")
    xlev = .getXlevels(mt, mf)
    y <- model.response(mf, "any")
    x <- if (!is.empty.model(mt))
        model.matrix(mt, mf, contrasts)
    else matrix(, NROW(y), 0)
    attr(x, "assign") = VGAM:::attrassigndefault(x, mt)
    if (!is.null(form2)) {
        if (!is.null(subset))
            stop("argument 'subset' cannot be used when ", "argument 'form2' is used")
        retlist = shadowvglm(formula = form2, family = family,
            data = data, na.action = na.action, control = vglm.control(...),
            method = method, model = model, x.arg = x.arg, y.arg = y.arg,
            contrasts = contrasts, constraints = constraints,
            extra = extra, qr.arg = qr.arg)
        Ym2 <- retlist$Ym2
        Xm2 <- retlist$Xm2
        if (length(Ym2)) {
            if (nrow(as.matrix(Ym2)) != nrow(as.matrix(y)))
                stop("number of rows of 'y' and 'Ym2' are unequal")
        }
        if (length(Xm2)) {
            if (nrow(as.matrix(Xm2)) != nrow(as.matrix(x)))
                stop("number of rows of 'y' and 'Ym2' are unequal")
        }
    }
    else {
        Xm2 = Ym2 = NULL
    }
    offset <- model.offset(mf)
    if (is.null(offset))
        offset <- 0
    w <- model.weights(mf)
    if (!length(w)) {
        w <- rep(1, nrow(mf))
    }
    else if (ncol(as.matrix(w)) == 1 && any(w < 0))
        stop("negative weights not allowed")
    if (is.character(family))
        family <- get(family)
    if (is.function(family))
        family <- family()
    if (!inherits(family, "vglmff")) {
        stop("'family = ", family, "' is not a VGAM family function")
    }
    eval(vcontrol.expression)
    if (length(slot(family, "first")))
        eval(slot(family, "first"))
    vglm.fitter <- get(method)
    fit <- vglm.fitter(x = x, y = y, Disper=Disper, w = w, nvec = nvec , cutoff = cutoff , offset = offset, approx = approx ,
        Xm2 = Xm2, Ym2 = Ym2, etastart = etastart, mustart = mustart,
        coefstart = coefstart, family = family, control = control,
        constraints = constraints, criterion = control$criterion,
        extra = extra, qr.arg = qr.arg, Terms = mt, function.name = function.name,
        ...)
    fit$misc$dataname <- dataname
    if (smart) {
        fit$smart.prediction <- get.smart.prediction()
        wrapup.smart()
    }
    answer <- new(Class = "vglm", assign = attr(x, "assign"),
        call = ocall, coefficients = fit$coefficients, constraints = fit$constraints,
        criterion = fit$crit.list, df.residual = fit$df.residual,
        df.total = fit$df.total, dispersion = 1, effects = fit$effects,
        family = fit$family, misc = fit$misc, model = if (model)
            mf
        else data.frame(), R = fit$R, rank = fit$rank, residuals = as.matrix(fit$residuals),
        rss = fit$rss, smart.prediction = as.list(fit$smart.prediction),
        terms = list(terms = mt) )
    if (!smart)
        answer@smart.prediction <- list(smart.arg = FALSE)
    if (qr.arg) {
        class(fit$qr) = "list"
        slot(answer, "qr") = fit$qr
    }
    if (length(attr(x, "contrasts")))
        slot(answer, "contrasts") = attr(x, "contrasts")
    if (length(fit$fitted.values))
        slot(answer, "fitted.values") = as.matrix(fit$fitted.values)
    slot(answer, "na.action") = if (length(aaa <- attr(mf, "na.action")))
        list(aaa)
    else list()
    if (length(offset))
        slot(answer, "offset") = as.matrix(offset)
    if (length(fit$weights))
        slot(answer, "weights") = as.matrix(fit$weights)
    if (x.arg)
        slot(answer, "x") = fit$x
    if (x.arg && length(Xm2))
        slot(answer, "Xm2") = Xm2
    if (y.arg && length(Ym2))
        slot(answer, "Ym2") = as.matrix(Ym2)
    if (!is.null(form2))
        slot(answer, "callXm2") = retlist$call
    answer@misc$formula = formula
    answer@misc$form2 = form2
    if (length(xlev))
        slot(answer, "xlevels") = xlev
    if (y.arg)
        slot(answer, "y") = as.matrix(fit$y)
    slot(answer, "control") = fit$control
    slot(answer, "extra") = if (length(fit$extra)) {
        if (is.list(fit$extra))
            fit$extra
        else {
            warning("'extra' is not a list, therefore placing ",
                "'extra' into a list")
            list(fit$extra)
        }
    }
    else list()
      slot(answer, "iter") = fit$iter
    slot(answer, "post") = fit$post
    fit$predictors = as.matrix(fit$predictors)
    if (length(fit$misc$predictors.names) == ncol(fit$predictors))
        dimnames(fit$predictors) = list(dimnames(fit$predictors)[[1]],
            fit$misc$predictors.names)
    slot(answer, "predictors") = fit$predictors
    if (length(fit$prior.weights))
        slot(answer, "prior.weights") = as.matrix(fit$prior.weights)
    answer
}


  myvglm.fit <- function (x, y, Disper=NULL , w = rep(1, length(x[, 1])), nvec , cutoff = Inf , approx = c("1","approx","exact"), X_vlm_arg = NULL,
      Xm2 = NULL, Ym2 = NULL, etastart = NULL, mustart = NULL,
      coefstart = NULL, offset = 0, family, control = vglm.control(),
      criterion = "coefficients", qr.arg = FALSE, constraints = NULL,
      extra = NULL, Terms = Terms, function.name = "vglm", ...)
  {
      approx = match.arg( approx , c("1","approx","exact") )
      specialCM = NULL
      post = list()
      check.rank <- TRUE
      nonparametric <- FALSE
      epsilon <- control$epsilon
      maxit <- control$maxit
      save.weight <- control$save.weight
      trace <- control$trace
      orig.stepsize <- control$stepsize
      minimize.criterion <- control$min.criterion
      n <- dim(x)[1]
      new.s.call <- expression({
          if (c.list$one.more) {
              fv <- c.list$fit
              new.coeffs <- c.list$coeff
              if (length(slot(family, "middle"))) eval(slot(family,
                  "middle"))
              eta <- fv + offset
              mu <- slot(family, "linkinv")(eta, extra)
              if (length(slot(family, "middle2"))) eval(slot(family,
                  "middle2"))
              old.crit <- new.crit
              new.crit <- switch(criterion, coefficients = new.coeffs,
                  tfun(mu = mu, Disper=Disper ,y = y, w = w, nvec = nvec , cutoff = cutoff , res = FALSE, eta = eta,
                    extra, approx = approx ))
              if (trace && orig.stepsize == 1) {
                  cat("VGLM    linear loop ", iter, ": ", criterion,
                    "= ")
                  UUUU = switch(criterion, coefficients = format(new.crit,
                    dig = round(2 - log10(epsilon))), format(round(new.crit,
                    4)))
                  switch(criterion, coefficients = {
                    if (length(new.crit) > 2) cat("\n")
                    cat(UUUU, fill = TRUE, sep = ", ")
                  }, cat(UUUU, fill = TRUE, sep = ", "))
              }
              {
                  take.half.step = (control$half.stepsizing && length(old.coeffs)) &&
                    ((orig.stepsize != 1) || (criterion != "coefficients" && (if (minimize.criterion) new.crit > old.crit else new.crit < old.crit)))
                  if (!is.logical(take.half.step)) take.half.step = TRUE
                  if (take.half.step) {
                    stepsize <- 2 * min(orig.stepsize, 2 * stepsize)
                    new.coeffs.save <- new.coeffs
                    if (trace) cat("Taking a modified step")
                    repeat {
                      if (trace) {
                        cat(".")
                        flush.console()
                      }
                      stepsize <- stepsize/2
                      if (too.small <- stepsize < 0.001) break
                      new.coeffs <- (1 - stepsize) * old.coeffs +
                        stepsize * new.coeffs.save
                      if (length(slot(family, "middle"))) eval(slot(family,
                        "middle"))
                      fv <- X_vlm_save %*% new.coeffs
                      if (M > 1) fv <- matrix(fv, n, M, byrow = TRUE)
                      eta <- fv + offset
                      mu <- slot(family, "linkinv")(eta, extra)
                      if (length(slot(family, "middle2"))) eval(slot(family,
                        "middle2"))
                      new.crit <- switch(criterion, coefficients = new.coeffs,
                        tfun(mu = mu, Disper = Disper, y = y, w = w, nvec = nvec, cutoff = cutoff, res = FALSE,
                          eta = eta, extra,approx=approx))
                      
                      if ((criterion == "coefficients") || (minimize.criterion &&
                        new.crit < old.crit) || (!minimize.criterion &&
                        new.crit > old.crit)) break
                    }
                    if (trace) cat("\n")
                    if (too.small) {
                      warning("iterations terminated because ",
                        "half-step sizes are very small")
                      one.more <- FALSE
                    } else {
                      if (trace) {
                        cat("VGLM    linear loop ", iter, ": ",
                          criterion, "= ")
                        UUUU = switch(criterion, coefficients = format(new.crit,
                          dig = round(2 - log10(epsilon))), format(round(new.crit,
                          4)))
                        switch(criterion, coefficients = {
                          if (length(new.crit) > 2) cat("\n")
                          cat(UUUU, fill = TRUE, sep = ", ")
                        }, cat(UUUU, fill = TRUE, sep = ", "))
                      }
                      one.more <- eval(control$convergence)
                    }
                  } else {
                    one.more <- eval(control$convergence)
                  }
              }
              flush.console()
              if (!is.logical(one.more)) one.more = FALSE
              if (one.more) {
                  iter <- iter + 1
                  deriv.mu <- eval(slot(family, "deriv"))
                  wz <- eval(slot(family, "weight"))
                  if (control$checkwz) wz = VGAM:::checkwz(wz, M = M,
                    trace = trace, wzepsilon = control$wzepsilon)
                  U <- VGAM:::vchol(wz, M = M, n = n, silent = !trace)
                  tvfor <- VGAM:::vforsub(U, as.matrix(deriv.mu), M = M,
                    n = n)
                  z <- eta + VGAM:::vbacksub(U, tvfor, M = M, n = n) -
                    offset
                  c.list$z <- z
                  c.list$U <- U
                  if (copy_X_vlm) c.list$X_vlm <- X_vlm_save
              }
              c.list$one.more <- one.more
              c.list$coeff = runif(length(new.coeffs))
              old.coeffs <- new.coeffs
          }
          c.list
      })
      copy_X_vlm <- FALSE
      stepsize <- orig.stepsize
      old.coeffs <- coefstart
      intercept.only <- ncol(x) == 1 && dimnames(x)[[2]] == "(Intercept)"
      y.names <- predictors.names <- NULL
      n.save <- n
      if (length(slot(family, "initialize")))
          eval(slot(family, "initialize"))
      if (length(etastart)) {
          eta <- etastart
          mu <- if (length(mustart))
              mustart
          else if (length(body(slot(family, "linkinv"))))
              slot(family, "linkinv")(eta, extra)
          else warning("argument 'etastart' assigned a value ",
              "but there is no 'linkinv' slot to use it")
      }
      if (length(mustart)) {
          mu <- mustart
          if (length(body(slot(family, "linkfun")))) {
              eta <- slot(family, "linkfun")(mu, extra)
          }
          else {
              warning("argument 'mustart' assigned a value ", "but there is no 'link' slot to use it")
          }
      }
      M <- if (is.matrix(eta))
          ncol(eta)
      else 1
      if (length(slot(family, "constraints")))
          eval(slot(family, "constraints"))
      Blist <- VGAM:::process.constraints(constraints, x, M, specialCM = specialCM)
      ncolBlist <- unlist(lapply(Blist, ncol))
      dimB <- sum(ncolBlist)
      X_vlm_save = if (length(X_vlm_arg))
          X_vlm_arg
      else lm2vlm.model.matrix(x, Blist, xij = control$xij, Xm2 = Xm2)
      if (length(coefstart)) {
          eta <- if (ncol(X_vlm_save) > 1)
              X_vlm_save %*% coefstart + offset
          else X_vlm_save * coefstart + offset
          eta <- if (M > 1)
              matrix(eta, ncol = M, byrow = TRUE)
          else c(eta)
          mu <- slot(family, "linkinv")(eta, extra)
      }
      if (criterion != "coefficients") {
          tfun <- slot(family, criterion)
      }
      iter <- 1
      new.crit <- switch(criterion, coefficients = 1, tfun(mu = mu,
          Disper = Disper , y = y, w = w, nvec = nvec , cutoff = cutoff , res = FALSE, eta = eta, extra,approx=approx))
      old.crit <- if (minimize.criterion)
          10 * new.crit + 10
      else -10 * new.crit - 10
      deriv.mu <- eval(slot(family, "deriv"))
      wz <- eval(slot(family, "weight"))
      if (control$checkwz)
          wz = VGAM:::checkwz(wz, M = M, trace = trace, wzepsilon = control$wzepsilon)
      U <- VGAM:::vchol(wz, M = M, n = n, silent = !trace)
      tvfor <- VGAM:::vforsub(U, as.matrix(deriv.mu), M = M, n = n)
      z <- eta + VGAM:::vbacksub(U, tvfor, M = M, n = n) - offset
      c.list <- list(z = as.double(z), fit = as.double(t(eta)),
          one.more = TRUE, coeff = as.double(rep(1, ncol(X_vlm_save))),
          U = as.double(U), copy_X_vlm = copy_X_vlm, X_vlm = if (copy_X_vlm) as.double(X_vlm_save) else double(3))
      dX_vlm <- as.integer(dim(X_vlm_save))
      nrow_X_vlm <- dX_vlm[[1]]
      ncol_X_vlm <- dX_vlm[[2]]
      if (nrow_X_vlm < ncol_X_vlm)
          stop(ncol_X_vlm, "parameters but only ", nrow_X_vlm,
              " observations")
      bf.call <- expression(VGAM:::vlm.wfit(xmat = X_vlm_save, z, Blist = NULL,
          U = U, matrix.out = FALSE, is.vlmX = TRUE, qr = qr.arg,
          xij = NULL))
      while (c.list$one.more) {
          tfit <- eval(bf.call)
          c.list$coeff <- tfit$coefficients
          tfit$predictors <- tfit$fitted.values
          c.list$fit <- tfit$fitted.values
          c.list <- eval(new.s.call)
          NULL
      }
      if (maxit > 1 && iter >= maxit && !control$nowarning)
          warning("convergence not obtained in ", maxit, " iterations")
      dnrow_X_vlm <- labels(X_vlm_save)
      xnrow_X_vlm <- dnrow_X_vlm[[2]]
      ynrow_X_vlm <- dnrow_X_vlm[[1]]
      if (length(slot(family, "fini")))
          eval(slot(family, "fini"))
      if (M > 1)
          tfit$predictors <- matrix(tfit$predictors, n, M)
      coefs <- tfit$coefficients
      asgn <- attr(X_vlm_save, "assign")
      names(coefs) <- xnrow_X_vlm
      rank <- tfit$rank
      cnames <- xnrow_X_vlm
      if (check.rank && rank < ncol_X_vlm)
          stop("vglm only handles full-rank models (currently)")
      R <- tfit$qr$qr[1:ncol_X_vlm, 1:ncol_X_vlm, drop = FALSE]
      R[lower.tri(R)] <- 0
      attributes(R) <- list(dim = c(ncol_X_vlm, ncol_X_vlm), dimnames = list(cnames,
          cnames), rank = rank)
      effects <- tfit$effects
      neff <- rep("", nrow_X_vlm)
      neff[seq(ncol_X_vlm)] <- cnames
      names(effects) <- neff
      dim(tfit$predictors) <- c(n, M)
      dn <- labels(x)
      yn <- dn[[1]]
      xn <- dn[[2]]
      residuals <- z - tfit$predictors
      if (M == 1) {
          tfit$predictors <- as.vector(tfit$predictors)
          residuals <- as.vector(residuals)
          names(residuals) <- names(tfit$predictors) <- yn
      }
      else {
          dimnames(residuals) <- dimnames(tfit$predictors) <- list(yn,
              predictors.names)
      }
      if (is.matrix(mu)) {
          if (length(dimnames(y)[[2]])) {
              y.names <- dimnames(y)[[2]]
          }
          if (length(dimnames(mu)[[2]])) {
              y.names <- dimnames(mu)[[2]]
          }
          dimnames(mu) <- list(yn, y.names)
      }
      else {
          names(mu) <- names(fv)
      }
      df.residual <- nrow_X_vlm - rank
      fit <- list(assign = asgn, coefficients = coefs, constraints = Blist,
          df.residual = df.residual, df.total = n * M, effects = effects,
          fitted.values = mu, offset = offset, rank = rank, residuals = residuals,
          R = R, terms = Terms )
      if (qr.arg) {
          fit$qr <- tfit$qr
          dimnames(fit$qr$qr) <- dnrow_X_vlm
      }
      if (M == 1) {
          wz <- as.vector(wz)
      }
      fit$weights <- if (save.weight)
          wz
      else NULL
      misc <- list(colnames.x = xn, colnames.X_vlm = xnrow_X_vlm,
          criterion = criterion, function.name = function.name,
          intercept.only = intercept.only, predictors.names = predictors.names,
          M = M, n = n, nonparametric = nonparametric, nrow_X_vlm = nrow_X_vlm,
          orig.assign = attr(x, "assign"), p = ncol(x), ncol_X_vlm = ncol_X_vlm,
          ynames = dimnames(y)[[2]])
      crit.list <- list()
      if (criterion != "coefficients")
          crit.list[[criterion]] <- fit[[criterion]] <- new.crit
      for (ii in names(VGAM:::.min.criterion.VGAM)) {
          if (ii != criterion && any(slotNames(family) == ii) &&
              length(body(slot(family, ii)))) {
              fit[[ii]] <- crit.list[[ii]] <- (slot(family, ii))(mu = mu,
                  Disper = Disper , y = y, w = w, nvec = nvec, cutoff = cutoff, res = FALSE, eta = eta, extra)
          }
      }
      if (w[1] != 1 || any(w != w[1]))
          fit$prior.weights <- w
      if (length(slot(family, "last")))
          eval(slot(family, "last"))
      structure(c(fit, list(predictors = tfit$predictors, contrasts = attr(x,
          "contrasts"), control = control, crit.list = crit.list,
          extra = extra, family = family, iter = iter, misc = misc,
          post = post, rss = tfit$rss, x = x, y = y)), vclass = slot(family,
          "vfamily"))
  }



 vcontrol.expression <- expression({
    control <- control
    mylist <- family@vfamily
    for (i in length(mylist):1) {
        for (ii in 1:2) {
            temp <- paste(if (ii == 1)
                ""
            else paste(function.name, ".", sep = ""), mylist[i],
                ".control", sep = "")
            tempexists = if (is.R())
                exists(temp, envir = VGAM:::VGAMenv)
            else exists(temp, inherit = TRUE)
            if (tempexists) {
                temp <- get(temp)
                temp <- temp(...)
                for (k in names(temp)) control[[k]] <- temp[[k]]
            }
        }
    }
    orig.criterion = control$criterion
    if (control$criterion != "coefficients") {
        try.crit = c(names(VGAM:::.min.criterion.VGAM), "coefficients")
        for (i in try.crit) {
            if (any(slotNames(family) == i) && ((is.R() && length(body(slot(family,
                i)))) || ((!is.R() && length(slot(family, i)) >
                1)))) {
                control$criterion <- i
                break
            }
            else control$criterion <- "coefficients"
        }
    }
    control$min.criterion <- control$min.criterion[control$criterion]
    for (ii in 1:2) {
        temp <- paste(if (ii == 1)
            ""
        else paste(function.name, ".", sep = ""), family@vfamily[1],
            ".", control$criterion, ".control", sep = "")
        if (exists(temp, inherit = T)) {
            temp <- get(temp)
            temp <- temp(...)
            for (k in names(temp)) control[[k]] <- temp[[k]]
        }
    }
})



db.free.estimate <- function(data.matrix , X, dispersions = NULL , nulldispersions = NULL , offsets = NULL , MC = FALSE , nCore = NULL , approx = c("1","approx","exact") , null.model = TRUE , limit=500 , cutoff = Inf , betas.true = NULL )
                        {
                           require( MM )
                           if( !is.list(X) ) X = list(X)

                           approx = match.arg( approx , c("1","approx","exact") )
                           if( !is.null(betas.true) ) betas <- as.vector( betas.true )

                           if( !is.list( data.matrix ) ) stop("Data set must be a list")

                            if( MC )
                             {
                               require( parallel )
                               if( is.null( nCore ) ) stop( "must specify number of cores (nCore) when multicoring" )

                             if( is.null(betas.true) )
                              {
                                if( is.null( dispersions ) ) combList <- mapply( list , X , data.matrix , offsets , SIMPLIFY = FALSE )
                                 else combList <- mapply( list , X , data.matrix , offsets , split( dispersions , 1:length(dispersions) ) , split( nulldispersions , 1:length(nulldispersions) ) , SIMPLIFY = FALSE )
                              }   
                             else
                              {
                                if( is.null( dispersions ) ) combList <- mapply( list , X , data.matrix , offsets , split(betas,1:length(betas)) , SIMPLIFY = FALSE )
                                 else combList <- mapply( list , X , data.matrix , offsets , split( dispersions , 1:length(dispersions) ) , split( nulldispersions , 1:length(nulldispersions) ) , split(betas,1:length(betas)) , SIMPLIFY = FALSE )
                              } 
                             
                             names( combList ) <- names( offsets )

                               dbFit = mclapply( combList ,
                                                FUN = function( data )
                                                       {
                                                       #problem data is rows 1, 6 and 1074. For 6, you get Inf likelihood but for the null test it's fine. row 1074 is opposite of 6
                                                          des<-data[[1]]
                                                          CT<-data[[2]]
                                                          OFF<-data[[3]]
                                                          if( !is.null(dispersions) )
                                                           {
                                                             disp <- data[[4]]
                                                             nulldisp <- data[[5]]
                                                             if( !is.null(betas.true) ) beta <- data[[6]]
                                                              else beta = NULL
                                                           }  
                                                          else
                                                           {
                                                             disp = NULL
                                                             nulldisp=NULL
                                                             if( !is.null(betas.true) ) beta <- data[[4]]
                                                              else beta = NULL 
                                                           }  

                                                          if( is.null(disp) ) fam = get("mydexpbinomial")
                                                           else fam = get("mydexpbinomial.profile")

                                                          options(warn=-1)
                                                          assign( "limit2" , limit , envir=.GlobalEnv )
                                                          optim=FALSE
                                                          optim.null = FALSE
                                                          if( is.null(disp) )
                                                           {
                                                             fit = try(myvglm(cbind(CT,OFF-CT) ~ des[,2] , fam = mydexpbinomial(ldisp=elogit(min=0,max=limit2), idispersion=.25 , zero=2) , nvec = OFF , Disper=disp , approx=approx , method = "myvglm.fit" , cutoff = cutoff ),TRUE)
                                                             if( null.model & class(fit)!="try-error" ) fit.null = try(myvglm(cbind(CT,OFF-CT) ~ 1 , fam = mydexpbinomial(ldisp=elogit(min=0,max=limit2), idispersion=.25 , zero=2) , nvec = OFF , Disper=nulldisp , approx=approx , method = "myvglm.fit" , cutoff = cutoff ),TRUE)
                                                             if( !is.null(betas.true) )
                                                              {
                                                                fit.true = myvglm(cbind(CT,OFF-CT) ~ 1 , fam = mydexpbinomial(ldisp=elogit(min=0,max=limit2), idispersion=.25 , zero=2) , offset =  des[,2]*beta ,  nvec = OFF , Disper=disp , approx=approx , method = "myvglm.fit" , cutoff = cutoff )
                                                                pval.true <- pchisq(2*(fit@criterion$loglikelihood - fit.true@criterion$loglikelihood),df=1,lower.tail=FALSE)
                                                                min.conf <- 1-pval.true
                                                              }  
                                                             else min.conf = NA
                                                           }
                                                           else
                                                             {
                                                               fit = try(myvglm(cbind(CT,OFF-CT) ~ des[,2] , fam = "mydexpbinomial.profile" , nvec = OFF , Disper=disp , approx=approx , method = "myvglm.fit" , cutoff = cutoff ),TRUE)
                                                               if( class(fit) != "try-error" ) if( fit@criterion$loglikelihood == -Inf | fit@criterion$loglikelihood == Inf ) optim = TRUE
                                                               if( class(fit) == "try-error" | optim )
                                                                {
                                                                  optim = TRUE
                                                                  nlogl <- function( pars , n , y , disp , groups )
                                                                     {
                                                                       b0 <- pars[1]
                                                                       b1 <- pars[2]
                                                                       p0 <- exp(b0)/(1+exp(b0))
                                                                       p1 <- exp(b0+b1)/(1+exp(b0+b1))
                                                                       uniq <- unique(groups)
                                                                       lnormC <- c( DBNormC(n=n[groups%in%uniq[1]],p=p0,phi=disp,log=TRUE) , DBNormC(n=n[groups%in%uniq[2]],p=p1,phi=disp,log=TRUE) )
                                                                       out <- disp*sum( y*b0 + y*b1*(groups%in%uniq[2]) - n*log(1+exp(b0 + b1*(groups%in%uniq[2]))) ) - sum( lnormC )
                                                                       return( -out )
                                                                     }
                                                                   fit = try(optim( c(0,0) , fn = nlogl , n = OFF , y = CT , disp = disp , groups = des[,2] ),TRUE)
                                                                }
                                                               if( null.model & !optim ) 
                                                                { 
                                                                  fit.null = try(myvglm(cbind(CT,OFF-CT) ~ 1 , fam = "mydexpbinomial.profile" , nvec = OFF , Disper=nulldisp , approx=approx , method = "myvglm.fit" , cutoff = cutoff ),TRUE)
                                                                  if( class(fit.null) == "try-error" ) optim.null = TRUE
                                                                   else fit.null = fit.null@criterion$loglikelihood
                                                                }  
                                                               if( null.model & (optim.null | optim) )
                                                                {
                                                                  nlogl <- function( par , n , y , disp )
                                                                     {
                                                                       b0 <- par
                                                                       p0 <- exp(b0)/(1+exp(b0))
                                                                       lnormC <- DBNormC(n=n,p=p0,phi=nulldisp,log=TRUE)
                                                                       out <- disp*sum( y*b0 - n*log(1+exp(b0))) - sum( lnormC )
                                                                       return( -out )
                                                                     }
                                                                   fit.null = try(optim( 0 , fn = nlogl , n = OFF , y = CT , disp = nulldisp , lower=-20, upper=20, method = "Brent"),TRUE)
                                                                   fit.null = fit.null$value
                                                                }
                                                               if( !is.null(betas.true) )
                                                                {
                                                                  fit.true = myvglm(cbind(CT,OFF-CT) ~ 1 , fam = "mydexpbinomial.profile" , offset = des[,2]*beta, nvec = OFF , Disper=disp , approx=approx , method = "myvglm.fit" , cutoff = cutoff )
                                                                  pval.true <- pchisq(2*(fit@criterion$loglikelihood - fit.true@criterion$loglikelihood),df=1,lower.tail=FALSE)
                                                                  min.conf <- 1-pval.true
                                                                }  
                                                               else min.conf = NA 
                                                             }

                                                          options(warn=0)
                                                          if( !optim ) 
                                                           {
                                                             if( null.model ) pval <- pchisq(2*(fit@criterion$loglikelihood - fit.null@criterion$loglikelihood),df=1,lower.tail=FALSE)
                                                              else pval = NA
                                                             fitted <- fitted.values( fit )
                                                             resid <- resid( fit )
                                                             coef <- coefficients( fit )
                                                             if( is.null(disp) ) disp <- elogit( coef["(Intercept):2"] , inverse = TRUE , min = 0 , max = limit2 )
                                                             vars <- try(diag( vcov(fit) ),TRUE)
                                                             if( class(vars) == "try-error" ) vars = NA
                                                             if( !is.na(vars) ) zval <- coef/sqrt(vars)
                                                              else zval <- NA
                                                           }
                                                          else
                                                           {   
                                                             if( null.model ) pval <- pchisq(2*(fit$value - fit.null),df=1,lower.tail=FALSE)
                                                              else pval = NA                         
                                                             coef = fit$par                                     
                                                             fitted = NA; resid=NA; vars=NA; zval=NA;
                                                           }
                                                          return( list( zval = zval , fitted = fitted , vars = vars , resid = resid , coef = coef , disp = disp , pval = pval , conf = min.conf ) )
                                                       } , mc.cores = nCore , mc.preschedule = TRUE )
                             }
                            else
                             {
                               message( "Not coded yet. Set MC to TRUE" )
                             }

                            names( dbFit ) <- names( data.matrix )
                            db.fitted <- sapply( dbFit , "[[" , "fitted" , simplify = FALSE )
                              names( db.fitted ) <- names( dbFit )
                            db.resid <- sapply( dbFit , "[[" , "resid" , simplify = FALSE )
                              names( db.resid ) <- names( dbFit )
                            db.coef <- t(sapply( dbFit , "[[" , "coef" , simplify = TRUE ))
                              if( NROW( db.coef ) == 1 & length(data.matrix) > 1 ) db.coef = t(db.coef)
                              rownames( db.coef ) <- names( dbFit )
                            db.disp <- sapply( dbFit ,  "[[" , "disp" , simplify = TRUE )
                              names( db.disp ) <- names( dbFit )
                            db.pval <- sapply( dbFit ,  "[[" , "pval" , simplify = TRUE )
                              names( db.pval ) <- names( dbFit )
                            db.var <- t(sapply( dbFit , "[[" , "vars" , simplify = TRUE ))
                              if( NROW( db.var ) == 1 & length(data.matrix) > 1 ) db.var = t(db.var)
                              rownames( db.var ) <- names( dbFit )
                            db.zval <- t(sapply( dbFit , "[[" , "zval" , simplify = TRUE ))
                              if( NROW( db.zval ) == 1 & length(data.matrix) > 1 ) db.zval = t(db.zval)
                              rownames( db.zval ) <- names( dbFit )
                            db.conf <- sapply( dbFit ,  "[[" , "conf" , simplify = TRUE )
                              names( db.conf ) <- names( dbFit )

                            out <- list( db.coef , db.disp , db.var , db.fitted , db.resid , db.pval , db.zval , db.conf )
                            names( out ) <- c( "betas" , "disps" , "vars" , "fits" , "resids" , "lrtpvals" , "zvals" , "conf" )
                            return(out)
                        }

db.constrained.estimate <- function(data.matrix, type = c( "matrix" , "list") , X , Beta.qb , disp.qb , offsets=NULL , MC=FALSE , nCore = NULL , weights = NULL , trunc = FALSE , approx = c( "1" , "approx" , "exact" ) , cutoff = Inf , betas.true = NULL )
                      {
                        require(MM)
                        approx = match.arg( approx , c("1","approx","exact") )
                        require( locfit )
                        if( length( type ) == 2 ) type <- "matrix"

                        if( !is.list(X) ) X = list(X)
                        if( !is.null(betas.true) ) betas <- as.vector( betas.true )

                        if (class(data.matrix) != type & type == "matrix" ) {
                            data.matrix = as.matrix(data.matrix)
                        }
                        if (class(data.matrix) != type & type == "list" ) {
                            stop( "argument is not list" )
                        }

                        if( type == "list")
                         {
                           Beta.qb.list <- lapply( apply( Beta.qb , 1 , list ) , "[[" , 1 )
                           p = dim(X[[1]])[2]
                           n1 = sapply( X , FUN = function( Y ){ length( which( Y[ , 2 ] == 0 ) ) } , simplify = TRUE )
                           n2 = sapply( X , FUN = function( Y ){ length( which( Y[ , 2 ] == 1 ) ) } , simplify = TRUE )
                           n = mapply( sum , n1 , n2 , SIMPLIFY = FALSE )
                           Y = data.matrix
                           N = offsets
                           m = length(Y) # no. of genes
                           x = mapply( FUN = function( des , b ){ apply(des %*% as.matrix(b), 2, mean) } , X , Beta.qb.list , SIMPLIFY = TRUE )
                              names( x ) <- names(Beta.qb.list)
                           prob = 1/(1+exp(-x))
                           xx = x[ prob >= 0.001 & prob <= 0.999 & !is.na(prob) ]
                            y = disp.qb[ prob >= 0.001 & prob <= 0.999 & !is.na(prob)]
                            yy = y[ y!=0 ]
                            xxx = xx[ y!=0 ]
                            yyy = log( yy )
                            if( trunc )
                             {
                               xxx = xxx[ yy >= 1 ]
                               yyy = yyy[ yy >= 1 ]
                             }
                            if( !is.null(weights) ) weights = weights[ names(yyy) ]
                           if( is.null( weights ) )
                            {
                              disp.fit.loc <- locfit( yyy ~ xxx )
                              dispD.fit.loc <- locfit( yyy ~ xxx , deriv = 1 )
                            }
                           else
                            {
                              dispD.fit.loc <- locfit( yyy ~ xxx , weights = weights )
                              dispD.fit.loc <- locfit( yyy ~ xxx , weights = weights , deriv = 1 )
                            }
                            
                           dispersionsLog.free <- vector(length=length(x))
                           dispersionsLog.free[!is.na(x)] <- predict( disp.fit.loc , newdata = x[!is.na(x)] )
                           dispersionsLog.free[is.na(x)] <- NA

                           Beta.model = matrix(NA, p, m)
                           Bvar.model = matrix(NA, p, m)
                            if( MC )
                             {
                               require( parallel )
                               if( is.null( nCore ) ) stop( "must specify number of cores (nCore) when multicoring" )

                               if( is.null(betas.true) ) combList <- mapply( list , X , Beta.qb.list , Y , N , SIMPLIFY = FALSE )
                                else combList <- mapply( list , X , Beta.qb.list , Y , N , split(betas,1:length(betas)) , SIMPLIFY = FALSE )

                               temp = mclapply( combList ,
                                                FUN = function( data )
                                                       {
                                                         des<-data[[1]]
                                                         b<-data[[2]]
                                                         CT<-data[[3]]
                                                         OFF<-data[[4]]
                                                         if( !is.null(betas.true) ) beta <- data[[5]]
                                                          else beta = NULL 
                                                         fit <- try(optim(par = b , fn = db.ProfileLik ,
                                                            X = des, counts = CT , disp.fit = disp.fit.loc , dispD.fit = dispD.fit.loc , offsets = OFF , approx = approx , cutoff = cutoff ,
                                                            hessian = TRUE, method="Nelder-Mead"), TRUE)
                                                         fit.null <- try(optim(par = 0.5 , fn = db.ProfileLik ,
                                                            X = des[,1,drop=FALSE], counts = CT , disp.fit = disp.fit.loc , dispD.fit = dispD.fit.loc , offsets = OFF , approx = approx , cutoff = cutoff , int.only = TRUE ,
                                                            hessian = FALSE, method="Brent" , lower = 0 , upper = 1 ), TRUE)
                                                         if( class(fit)!="try-error" & class(fit.null)!="try-error" )
                                                          {
                                                            lrtpval <- pchisq(2*abs(fit$value-fit.null$value),df=ncol(des)-1,lower.tail=FALSE)
                                                             if( !is.null(betas.true) )
                                                              {
                                                                fit.true = try(optim(par = 0.5 , fn = db.ProfileLik , model.offsets = des[,2]*beta ,
                                                                X = des[,1,drop=FALSE], counts = CT , disp.fit = disp.fit.loc , dispD.fit = dispD.fit.loc , offsets = OFF , approx = approx , cutoff = cutoff , int.only = TRUE ,
                                                                hessian = FALSE, method="Brent" , lower = 0 , upper = 1 ), TRUE)
                                                                pval.true <- pchisq(2*abs(fit$value - fit.true$value),df=ncol(des)-1,lower.tail=FALSE)
                                                                min.conf <- 1-pval.true
                                                              } 
                                                             else min.conf = NA
                                                          } 
                                                          else 
                                                           {
                                                             lrtpval = NA  
                                                             min.conf = NA
                                                           }  
                                                         return( list( fit , lrtpval , min.conf ) )
                                                       } , mc.cores = nCore , mc.preschedule = TRUE )
                             }
                            else
                             {
                               message( "Not coded yet. Set MC to TRUE" )
                             }

                            lrtpvals <- as.vector(sapply( temp , "[[" , 2 ))
                            conf <- as.vector(sapply( temp , "[[" , 3 ))
                            temp <- lapply( temp , "[[" , 1 )
                            Beta.model <- t(sapply( temp , FUN = function(x)
                                                               {
                                                                 if(class(x)!="try-error") return( x$par[1:p] )
                                                                 else{ return( rep( NA, p ) ) }
                                                               } , simplify = TRUE ))
                           Beta.model.list <- lapply( apply( as.matrix(Beta.model) , 1 , list ) , "[[" , 1 )

                        fitted <- mapply( FUN = function( des , b ){ apply(des %*% as.matrix(b), 2, mean) } , X , Beta.model.list , SIMPLIFY = TRUE )

                           dispersionsLog.cstr <- vector(length=length(fitted))
                           if( any(is.na(fitted)) )
                            {
                              wh.na <- which(is.na(fitted)==TRUE)
                              wh.not.na <- which(!is.na(fitted)==TRUE)
                              fitted.notNA <- fitted[ wh.not.na ]
                              disps.cstr.tmp <- predict( disp.fit.loc , newdata = fitted.notNA )
                              dispersionsLog.cstr[wh.not.na] <- disps.cstr.tmp
                              dispersionsLog.cstr[wh.na] <- NA
                            } else
                               {
                                 dispersionsLog.cstr <- predict( disp.fit.loc , newdata = fitted )
                                }
                           names(dispersionsLog.cstr)<-names(fitted)

                            var <- sapply( temp , FUN = function(x)
                                                               {
                                                                 return(try(diag(solve(x$hessian)), TRUE))
                                                               } , simplify = FALSE )
                            Bvar.model <- t(sapply( var , FUN = function(x)
                                                               {
                                                                 if(class(x)!="try-error") return( x[1:p] )
                                                                 else{ return( rep( NA, p ) ) }
                                                               } , simplify = TRUE ))
                            convergence.constrained <- sapply( temp , FUN = function(x)
                                           {
                                             if(class(x)!="try-error") return( x$convergence )
                                             else{ return( rep( NA, 1 ) ) }
                                           } , simplify = TRUE )

                            k = mapply( FUN = function( x , y ){ x - 2 - apply( as.matrix(y == 0) , 2 , sum ) } , n , Y , SIMPLIFY = TRUE )
                            k[k <= 0] = 1
                            options(warn = -1)
                            pt.model = 2 * pt(-abs(as.numeric(Beta.model[,2])/sqrt(as.numeric(Bvar.model[,2]))), df = k)
                              names(pt.model)<-rownames(Beta.model)
                            options(warn = 0)
                            out = list(Beta.model, Bvar.model, pt.model , fitted , lrtpvals , conf , dispersionsLog.free , dispersionsLog.cstr , convergence.constrained )
                            names(out) = c("betas", "vars", "pvals","fits","lrtpvals","conf","dispsFit.log.free","dispsFit.log.cstr","convergence")
                            return(out)
                         }
                        else stop( "This part of code not yet written" )
                      }

db.ProfileLik <- function( Beta.qb , counts , offsets , X , disp.fit , dispD.fit = NULL , approx = c("1","approx","exact") , cutoff = Inf , int.only = FALSE , model.offsets = NULL )
             {
               approx = match.arg( approx , c("1","approx","exact") )
               dimX <- dim(as.matrix(X))
               X <- matrix(as.numeric(as.character(X)),nr=dimX[1],ncol=dimX[2])
               #if( dimX[2] == 1 ) X <- cbind(1,X) #commented out 5/24/2013. Don't believe it's needed but you never know so commenting it out.

               if( int.only ) Beta.qb <- Beta.qb/(1-Beta.qb)
               if( class(Beta.qb) != "matrix" ) Beta.qb <- as.matrix(Beta.qb)
               
               if( is.null(model.offsets) ) model.offsets = 0                 
               mu = exp( X%*%Beta.qb + model.offsets ) / (1 + exp( X%*%Beta.qb + model.offsets ) )
               newx = mean(X%*%Beta.qb + model.offsets)
               disp = exp( predict( disp.fit , newdata = newx ) )
  
               db.lik = vector(length=length(counts))
               for( i in 1:length(counts) )
                {
                  paras = paras( c(mu[i],disp) )
                  db.lik[i] = DB_loglik( counts[i] , offsets[i] , mu[i] , disp , approx = approx , cutoff = cutoff )
                }
               return( -sum(db.lik) )
             }



