## ----
## Set parameters
## ----
n_sims=10

# MultiMAP
print("lot 1")
system(paste("source multimap3.7/bin/activate",
             "python 4a_run_multiMAP.py",
             "../2210_0simu/simulations/rna/",
             "../2210_0simu/simulations/met/",
             "4a_multiMAP/",
             1,
             n_sims))
print("lot 2")
system(paste("source multimap3.7/bin/activate",
             "python 4a_run_multiMAP.py",
             "../2210_0simu/simulations/rna/",
             "../2210_0simu/simulations/met/",
             "4a_multiMAP/",
             2,
             n_sims))
print("lot 3")
system(paste("source multimap3.7/bin/activate",
             "python 4a_run_multiMAP.py",
             "../2210_0simu/simulations/rna/",
             "../2210_0simu/simulations/met/",
             "4a_multiMAP/",
             3,
             n_sims))
