library(BCEA)
#load_all("..")
#load_all("voi")
#load_all("../../../BCEA")

if (0) {  ## TEST DISABLED WHILE ldr has been taken off CRAN

chemo_bcea <- bcea(e=chemo_cea$e, c=chemo_cea$c, wtp=30000)

test_that("Agrees with BCEA: single parameter, default GAM",{
    expect_equal(
        BCEA::evppi(parameter=c("pi1"), he=chemo_bcea, input=chemo_pars)$evppi[1],
        voi::evppi(outputs=chemo_cea, inputs=chemo_pars, poi=c("pi1"))[chemo_cea$k == 30000]
    )
})

test_that("Agrees with BCEA: single parameter, default GAM",{
    expect_equal(
        BCEA::evppi(parameter=c("pi1","rho"), he=chemo_bcea, input=chemo_pars)$evppi[1],
        voi::evppi(outputs=chemo_cea, inputs=chemo_pars, poi=c("pi1","rho"))[chemo_cea$k == 30000]
    )
})

chemo_bcea100 <- chemo_bcea
chemo_bcea100$n.sim <- 100
chemo_pars100 <- chemo_pars[1:100,]

## This is slow, should be in extra slow tests file

test_that("Agrees with BCEA: single parameter, INLA",{
    pars <- c("pi1","rho","gamma","gamma2")
    expect_equal(
        BCEA::evppi(parameter=pars, he=chemo_bcea100, input=chemo_pars, method="INLA")$evppi[1]
       ,
        voi::evppi(outputs=chemo_cea, inputs=chemo_pars, poi=pars, method="inla", nsim=100)[chemo_cea$k == 30000]
    )
})


    }
