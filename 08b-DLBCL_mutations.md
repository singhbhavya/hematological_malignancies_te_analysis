DLBCL mutations - example analysis
================
Matthew Bendall
30 January, 2025

- [Statistical Testing (example)](#statistical-testing-example)
  - [Setup](#setup)
  - [Load and wrangle data](#load-and-wrangle-data)
  - [Chi square](#chi-square)
  - [Logistic regression](#logistic-regression)
    - [Cluster as predictor (variant ~
      cluster)](#cluster-as-predictor-variant--cluster)
    - [Variant as predictor (cluster ~
      variant)](#variant-as-predictor-cluster--variant)

<style type="text/css">
.scroll-150 {
  max-height: 150px;
  overflow-y: auto;
  background-color: inherit;
}
</style>

# Statistical Testing (example)

## Setup

## Load and wrangle data

``` r
load("r_outputs/05o-DLBCL_pca_ccp_clusters_metadata.Rdata")

dlbcl.tbl <- DLBCL_metadata %>% 
  rownames_to_column() %>%
  as_tibble %>%
  mutate(
    subtype.herv = factor(clust.retro.k7, levels = str_c('HC', 1:7))
  )

mut.simple.df <- readRDS('r_outputs/08-DLBCL_mutations/mut.simple.df.Rds')
fus.df <- readRDS('r_outputs/08-DLBCL_mutations/fus.df.Rds')

# factorize
mut.simple.df %<>% mutate(across(matches("."), ~ factor(.x, levels=c('.', 'x'))))
fus.df %<>% mutate(across(matches("."), ~ factor(.x, levels=c('.', 'x'))))

mutfus.tbl <- (
    dlbcl.tbl %>%
      select(case, subtype.herv) %>%
      left_join(
          mut.simple.df %>% rownames_to_column('case') %>% as_tibble(),
          by="case"
      ) %>%
      left_join(
          fus.df %>% rownames_to_column('case') %>% as_tibble(),
          by="case"
      )
)
```

``` r
mutfus.tbl
```

    ## # A tibble: 529 × 3,255
    ##    case   subtype.herv KMT2D PIM1  MYD88 TP53  HLA.B BTG2  CREBBP TNFAIP3 TMSB4X
    ##    <chr>  <fct>        <fct> <fct> <fct> <fct> <fct> <fct> <fct>  <fct>   <fct> 
    ##  1 DLBCL… HC1          .     .     .     x     .     .     .      .       .     
    ##  2 DLBCL… HC2          .     x     .     .     x     x     .      .       .     
    ##  3 DLBCL… HC3          x     .     .     .     .     .     x      .       .     
    ##  4 DLBCL… HC1          .     x     .     x     .     .     x      .       .     
    ##  5 DLBCL… HC4          .     .     x     .     .     .     .      .       .     
    ##  6 DLBCL… HC1          .     .     .     .     .     .     .      .       .     
    ##  7 DLBCL… HC5          x     .     .     x     .     .     .      .       .     
    ##  8 DLBCL… HC4          .     .     .     x     x     .     .      .       .     
    ##  9 DLBCL… HC1          .     .     .     .     .     .     .      .       .     
    ## 10 DLBCL… HC3          x     .     .     .     x     .     .      .       .     
    ## # ℹ 519 more rows
    ## # ℹ 3,244 more variables: OSBPL10 <fct>, HIST1H1E <fct>, B2M <fct>,
    ## #   CARD11 <fct>, DTX1 <fct>, TNFRSF14 <fct>, HLA.A <fct>, BTG1 <fct>,
    ## #   CD79B <fct>, SOCS1 <fct>, TBL1XR1 <fct>, DUSP2 <fct>, BCL6 <fct>,
    ## #   CCND3 <fct>, IRF8 <fct>, PRDM1 <fct>, KLHL6 <fct>, MEF2B <fct>, SGK1 <fct>,
    ## #   SPEN <fct>, BCL2 <fct>, CD58 <fct>, STAT3 <fct>, MPEG1 <fct>, SETD1B <fct>,
    ## #   FAS <fct>, HIST1H1C <fct>, IRF2BP2 <fct>, ARID1A <fct>, BCL10 <fct>, …

## Chi square

###### Is there a significant difference in the distribution of a mutation across the HERV subtypes?

(Could wrap this in an lapply)

``` r
chisq.myfunc <- function(g) {
      
    cat(paste0("## Gene: ", g, "\n"))
    obs.tbl <- as.matrix(table(mutfus.tbl$subtype.herv, mutfus.tbl[[g]]))
    
    cat("\n## Observed")
    print(obs.tbl)
    
    res <- chisq.test(obs.tbl)
    print(res)
    
    inttext <- paste0(
      "\n## Interpretation:\n",
      "## There is ", 
      ifelse(res$p.value<0.05, "", "not "),
      "a signficant difference in the distribution of ",
      g,
      " mutations\n## across HERV subtypes.\n"
    )
    
    cat(inttext)
    
    # cat("\n## Observed from test")
    # res$observed
    
    cat("\n## Expected")
    print(res$expected)
    
    cat("\n## Residuals")
    print(res$residuals)
    data.frame(res$residuals) %>% 
      mutate(Var2 = ifelse(Var2=='x', 'VAR', 'WT')) %>% 
      pivot_wider(names_from = Var2, values_from=Freq) %>%
      mutate(msg = ifelse(WT<0 & VAR>0, "more mutations", "fewer mutations")) %>%
      mutate(msg = paste0(Var1, ": ", msg, " than expected")) %>%
      pull(msg) %>%
      str_c(collapse='\n##     ') %>% cat('\n##     ', ., '\n', sep = '')

    contrib <- (res$observed - res$expected)^2 / res$expected
    pct.contrib <- (contrib / res$statistic) * 100
    
    cat("\n## Percent contribution")
    print(pct.contrib)
    
    if(res$p.value<0.05) {
      inttext2 <- paste0("\n## Interpretation:")
      pc <- data.frame(pct.contrib) %>% arrange(desc(Freq))
      pc <- {if(max(pc$Freq)>10) pc %>% filter(Freq>10) else pc %>% slice_max(Freq, n=3)}
      msgs <- sapply(1:nrow(pc), function(i) paste0(
        'For ', pc$Var1[i], ", ",
        ifelse(pc$Var2[i]=="x", "a mutation in ", "wildtype "), g, 
        " contributes ", 
        sprintf('%.1f', pc$Freq[i]), 
        "% to the overall signficant test."
      ))
      inttext2 <- paste0(inttext2, str_c(str_c('\n##     ', msgs), collapse=""))
      cat(inttext2)
    }
}
```

``` r
chisq.myfunc('MYD88')
```

    ## ## Gene: MYD88
    ## 
    ## ## Observed     
    ##        .  x
    ##   HC1 64 19
    ##   HC2 25 47
    ##   HC3 58 18
    ##   HC4 89 10
    ##   HC5 48  9
    ##   HC6 34  3
    ##   HC7 70 27
    ## 
    ##  Pearson's Chi-squared test
    ## 
    ## data:  obs.tbl
    ## X-squared = 81.694, df = 6, p-value = 1.596e-15
    ## 
    ## 
    ## ## Interpretation:
    ## ## There is a signficant difference in the distribution of MYD88 mutations
    ## ## across HERV subtypes.
    ## 
    ## ## Expected     
    ##              .         x
    ##   HC1 61.81190 21.188100
    ##   HC2 53.61996 18.380038
    ##   HC3 56.59885 19.401152
    ##   HC4 73.72745 25.272553
    ##   HC5 42.44914 14.550864
    ##   HC6 27.55470  9.445298
    ##   HC7 72.23800 24.761996
    ## 
    ## ## Residuals     
    ##                .          x
    ##   HC1  0.2783115 -0.4753583
    ##   HC2 -3.9084611  6.6756851
    ##   HC3  0.1862436 -0.3181056
    ##   HC4  1.7786768 -3.0379952
    ##   HC5  0.8519736 -1.4551782
    ##   HC6  1.2278494 -2.0971773
    ##   HC7 -0.2633164  0.4497467
    ## 
    ## ##     HC1: fewer mutations than expected
    ## ##     HC2: more mutations than expected
    ## ##     HC3: fewer mutations than expected
    ## ##     HC4: fewer mutations than expected
    ## ##     HC5: fewer mutations than expected
    ## ##     HC6: fewer mutations than expected
    ## ##     HC7: more mutations than expected
    ## 
    ## ## Percent contribution     
    ##                 .           x
    ##   HC1  0.09481387  0.27659986
    ##   HC2 18.69912617 54.55083425
    ##   HC3  0.04245926  0.12386610
    ##   HC4  3.87261024 11.29753965
    ##   HC5  0.88850932  2.59204224
    ##   HC6  1.84543999  5.38368960
    ##   HC7  0.08487223  0.24759719
    ## 
    ## ## Interpretation:
    ## ##     For HC2, a mutation in MYD88 contributes 54.6% to the overall signficant test.
    ## ##     For HC2, wildtype MYD88 contributes 18.7% to the overall signficant test.
    ## ##     For HC4, a mutation in MYD88 contributes 11.3% to the overall signficant test.

``` r
chisq.myfunc('fusion.BCL6')
```

    ## ## Gene: fusion.BCL6
    ## 
    ## ## Observed     
    ##        .  x
    ##   HC1 72 11
    ##   HC2 58 14
    ##   HC3 43 33
    ##   HC4 91  8
    ##   HC5 51  6
    ##   HC6 33  4
    ##   HC7 74 23
    ## 
    ##  Pearson's Chi-squared test
    ## 
    ## data:  obs.tbl
    ## X-squared = 44.579, df = 6, p-value = 5.675e-08
    ## 
    ## 
    ## ## Interpretation:
    ## ## There is a signficant difference in the distribution of fusion.BCL6 mutations
    ## ## across HERV subtypes.
    ## 
    ## ## Expected     
    ##              .        x
    ##   HC1 67.22841 15.77159
    ##   HC2 58.31862 13.68138
    ##   HC3 61.55854 14.44146
    ##   HC4 80.18810 18.81190
    ##   HC5 46.16891 10.83109
    ##   HC6 29.96929  7.03071
    ##   HC7 78.56814 18.43186
    ## 
    ## ## Residuals     
    ##                 .           x
    ##   HC1  0.58195177 -1.20150510
    ##   HC2 -0.04172216  0.08614010
    ##   HC3 -2.36537323  4.88357999
    ##   HC4  1.20738859 -2.49279001
    ##   HC5  0.71100180 -1.46794346
    ##   HC6  0.55361287 -1.14299626
    ##   HC7 -0.51536628  1.06403186
    ## 
    ## ##     HC1: fewer mutations than expected
    ## ##     HC2: more mutations than expected
    ## ##     HC3: more mutations than expected
    ## ##     HC4: fewer mutations than expected
    ## ##     HC5: fewer mutations than expected
    ## ##     HC6: fewer mutations than expected
    ## ##     HC7: more mutations than expected
    ## 
    ## ## Percent contribution     
    ##                 .           x
    ##   HC1  0.75970862  3.23835390
    ##   HC2  0.00390487  0.01664500
    ##   HC3 12.55082932 53.49949467
    ##   HC4  3.27014646 13.93941216
    ##   HC5  1.13400370  4.83383394
    ##   HC6  0.68752015  2.93064145
    ##   HC7  0.59580628  2.53969949
    ## 
    ## ## Interpretation:
    ## ##     For HC3, a mutation in fusion.BCL6 contributes 53.5% to the overall signficant test.
    ## ##     For HC4, a mutation in fusion.BCL6 contributes 13.9% to the overall signficant test.
    ## ##     For HC3, wildtype fusion.BCL6 contributes 12.6% to the overall signficant test.

## Logistic regression

### Cluster as predictor (variant ~ cluster)

Cluster as the predictor of variant (i.e. variant is the dependent
variable, it depends on the cluster assignment).

``` r
xtabs(~subtype.herv + MYD88, data=mutfus.tbl)
```

    ##             MYD88
    ## subtype.herv  .  x
    ##          HC1 64 19
    ##          HC2 25 47
    ##          HC3 58 18
    ##          HC4 89 10
    ##          HC5 48  9
    ##          HC6 34  3
    ##          HC7 70 27

``` r
model.logit <- glm(MYD88 ~ subtype.herv, data=mutfus.tbl, family="binomial")
summary(model.logit)
```

    ## 
    ## Call:
    ## glm(formula = MYD88 ~ subtype.herv, family = "binomial", data = mutfus.tbl)
    ## 
    ## Coefficients:
    ##                 Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)     -1.21444    0.26126  -4.648 3.34e-06 ***
    ## subtype.hervHC2  1.84572    0.35991   5.128 2.92e-07 ***
    ## subtype.hervHC3  0.04437    0.37557   0.118   0.9060    
    ## subtype.hervHC4 -0.97161    0.42365  -2.293   0.0218 *  
    ## subtype.hervHC5 -0.45953    0.44744  -1.027   0.3044    
    ## subtype.hervHC6 -1.21330    0.65622  -1.849   0.0645 .  
    ## subtype.hervHC7  0.26179    0.34580   0.757   0.4490    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for binomial family taken to be 1)
    ## 
    ##     Null deviance: 591.92  on 520  degrees of freedom
    ## Residual deviance: 515.57  on 514  degrees of freedom
    ##   (8 observations deleted due to missingness)
    ## AIC: 529.57
    ## 
    ## Number of Fisher Scoring iterations: 4

### Variant as predictor (cluster ~ variant)

Variant as the predictor of cluster (i.e. cluster is the dependent
variable, it depends on the variant). Allows testing of multiple
variants.

``` r
xtabs(~MYD88 + subtype.herv, data=mutfus.tbl)
```

    ##      subtype.herv
    ## MYD88 HC1 HC2 HC3 HC4 HC5 HC6 HC7
    ##     .  64  25  58  89  48  34  70
    ##     x  19  47  18  10   9   3  27

``` r
model.logit <- glm(subtype.herv ~ MYD88, data=mutfus.tbl, family="binomial")
summary(model.logit)
```

    ## 
    ## Call:
    ## glm(formula = subtype.herv ~ MYD88, family = "binomial", data = mutfus.tbl)
    ## 
    ## Coefficients:
    ##             Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)   1.6219     0.1368   11.86   <2e-16 ***
    ## MYD88x        0.1699     0.2830    0.60    0.548    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for binomial family taken to be 1)
    ## 
    ##     Null deviance: 456.94  on 520  degrees of freedom
    ## Residual deviance: 456.57  on 519  degrees of freedom
    ##   (8 observations deleted due to missingness)
    ## AIC: 460.57
    ## 
    ## Number of Fisher Scoring iterations: 4

``` r
model.logit <- glm(subtype.herv ~ MYD88*CD79B, data=mutfus.tbl, family="binomial")
summary(model.logit)
```

    ## 
    ## Call:
    ## glm(formula = subtype.herv ~ MYD88 * CD79B, family = "binomial", 
    ##     data = mutfus.tbl)
    ## 
    ## Coefficients:
    ##               Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)    1.85238    0.15525  11.932  < 2e-16 ***
    ## MYD88x        -0.03531    0.33693  -0.105   0.9165    
    ## CD79Bx        -1.73460    0.37704  -4.601 4.21e-06 ***
    ## MYD88x:CD79Bx  1.65212    0.65396   2.526   0.0115 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for binomial family taken to be 1)
    ## 
    ##     Null deviance: 456.94  on 520  degrees of freedom
    ## Residual deviance: 437.08  on 517  degrees of freedom
    ##   (8 observations deleted due to missingness)
    ## AIC: 445.08
    ## 
    ## Number of Fisher Scoring iterations: 4
