#' Weak-IV Robust Anderson Rubin Test
# fork of `ivmodel::AR.test` that accepts FELM objects and fits models with FEs
#' @param data dataframe
#' @param Y outcome (string)
#' @param D treatment (string)
#' @param Z instrument (string)
#' @param X control variables (character vector)
#' @param FE fixed effects (character vector)
#' @param prec precision of CI in string
#' @importFrom lfe felm
#' @export
AR_test = function(data, Y, D, Z, X, FE, prec = 3){
  # fit IV model in FELM to get correct DoFs
  mod = felm(formula_lfe(Y = Y, W = D, X = X, Z = Z, D = FE), data = data)
  # residualise
  Ytil = partialer(Y, X = X, FE = FE, data = data, weights = weights)
  Dtil = partialer(D, X = X, FE = FE, data = data, weights = weights)
  Ztil = partialer(Z, X = X, FE = FE, data = data, weights = weights)
  alpha = 0.05; n= length(Ztil);  k=mod$p;  l=length(Z)
  ZtilQR = qr(Ztil)

  # compute F
  tmp = Ytil - 0 * Dtil # test null of \beta = 0
  Fstat         = c(sum(qr.fitted(ZtilQR, tmp)^2))/
                  c(sum(tmp^2)- sum(qr.fitted(ZtilQR, tmp)^2))*(n-k-l)/l
  p.value=1-pf(Fstat, df1=l, df2=n-k-l)

  # CI construction ingredients
  cval          = qf(1-alpha, df1=l, df2=n-k-l)*l/(n-k-l)
  coef.beta0sq  = cval*sum(Dtil^2)-(cval+1)*sum(qr.fitted(ZtilQR, Dtil)^2)
  coef.beta0    = -2*cval*sum(Dtil*Ytil)+2*(cval+1)*sum(Dtil*qr.fitted(ZtilQR, Ytil))
  coef.constant = cval*sum(Ytil^2)-(cval+1)*sum(qr.fitted(ZtilQR, Ytil)^2)
  Delta         = coef.beta0^2-4*coef.constant*coef.beta0sq
  ci = matrix(NA, ncol=2); colnames(ci)<-c("lower", "upper")

  # CI construction gymnastics
  if(coef.beta0sq==0){
      if(coef.beta0>0){
        info=c("[",round(-coef.constant/coef.beta0, prec),",Infinity)")
        ci[1,]=c(-coef.constant/coef.beta0, Inf)
      }
      if(coef.beta0<0){
        info=c("(-Infinity,",round(-coef.constant/coef.beta0, prec) ,"]");
        ci[1,]=c(-Inf, -coef.constant/coef.beta0)
      }
      if(coef.beta0==0){
        if(coef.constant>=0){
          info="Whole Real Line"
          ci[1,]=c(-Inf, Inf)
        }
        if(coef.constant<0){
          info="Empty Set"
        }
      }
  } else if(coef.beta0sq!=0){
    if(Delta<=0){
      if(coef.beta0sq>0){
        info="Whole Real Line"
        ci[1,]=c(-Inf, Inf)
      }
      if(coef.beta0sq<0){
        info="Empty Set"
      }
    }
    if(Delta>0){
      # Roots of quadratic equation
      root1=(-coef.beta0+sqrt(Delta))/(2*coef.beta0sq)
      root2=(-coef.beta0-sqrt(Delta))/(2*coef.beta0sq)
      upper.root=max(root1,root2)
      lower.root=min(root1,root2)
      if(coef.beta0sq<0){
        info=paste("[",  round(lower.root, prec),", ", round(upper.root, prec),"]",sep="")
        ci[1, ]=c(lower.root, upper.root)
      }
      if(coef.beta0sq>0){
        info= paste("(-Infinity,",round(lower.root, prec),"] union [",round(upper.root, prec),",Infinity)")
        ci[1, ]=c(-Inf, lower.root)
        ci<-rbind(ci, c(upper.root, Inf))
      }
    }
  }
  return(list(Fstat=Fstat, df=c(l, n-k-l), p.value=p.value,
              ci.info=info, ci=ci))

}
