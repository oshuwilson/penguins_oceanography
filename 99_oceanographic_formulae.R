#-----------------------------------
#Oceanographic Formula Combinations
#-----------------------------------
vars <- c("depth", "slope", "dshelf", "s(sst, k = 5, bs = 'ts')", "mld", 
          "s(sal, k = 5, bs = 'ts')", "s(ssh, k = 5, bs = 'ts')", "curr", 
          "front_freq", "s(eddies, k = 5, bs = 'ts')")

formula <- ~ depth + slope + dshelf + mld + curr + front_freq +
  s(sst, k=5, bs = "ts") + s(sal, k=5, bs = "ts") + s(ssh, k=5, bs = "ts") + s(eddies, k=5, bs = "ts")

library(combinat)

combinations <- combn(vars, 2)

# Create GAM formulae
gam_formulae <- apply(combinations, 2, function(vars) {
  paste("~", vars[1], "+", vars[2], sep = " ")
})

# Save formulae
saveRDS(gam_formulae, "code/out/oceanographic_formulae.rds")
