library(wildboottestjlr)
library(JuliaConnectoR)
library(pracma)

N <- 1000
seed <- 879345
data = wildboottestjlr:::create_data(N = 1000,
                                                  N_G1 = 40,
                                                  icc1 = 0.5,
                                                  N_G2 = 20,
                                                  icc2 = 0.2,
                                                  numb_fe1 = 10,
                                                  numb_fe2 = 10,
                                                  seed = 90864369,
                                                  weights = 1:N / N)

resp = data$proposition_vote
predexog = cbind(1, as.matrix(data[,c("treatment","log_income")]))
clustid = data$group_id1
R <- matrix(c(0,1,0), nrow=1)
r <- c(0)

stopJulia()
startJuliaServer()
juliaEval(r'(push!(LOAD_PATH, raw"D:\OneDrive\Documents\Macros\WildBootTests.jl"))')
juliaEval('using WildBootTests')

tic(); juliaLet('wildboottest(R, r; resp, predexog, clustid)', R=R, r=r, resp=resp, predexog=predexog, clustid=clustid); toc()
tic(); juliaLet('wildboottest(R, r; resp, predexog, clustid)', R=R, r=r, resp=resp, predexog=predexog, clustid=clustid); toc()
