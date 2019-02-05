# Run this on n59
# ssh bioinfo.sb-roscoff.fr
# qlogin -q formation.q@@bignode.formation
# export PATH=/opt/6.x/R-3.5.1/bin/:$PATH

library(dyno)
library(dplyr)

data("fibroblast_reprogramming_treutlein")

dim(fibroblast_reprogramming_treutlein$counts)

# Prepare the data -------------------------------------------------------------

# some methods require pure counts, others require normalized expression
# we there have to provide both
dataset <- wrap_expression(
  counts = fibroblast_reprogramming_treutlein$counts,
  expression = fibroblast_reprogramming_treutlein$expression
)

# Choosing what method to ruun -------------------------------------------------
# Go to http://35.232.29.209/ or guidelines.dynverse.org
# When running this on your laptop, you can do `guidelines_shiny()`

# Running methods --------------------------------------------------------------

# infer a trajectory using default parameters
model <- infer_trajectory(dataset, "slingshot")

# the trajectory is contained in the following two variables
model$milestone_network
model$progressions

# infer a trajectory with prior information
model <- infer_trajectory(dataset, "paga")

dataset <- dataset %>% add_prior_information(start_id = dataset$counts[,"Vim"] %>% which.max() %>% names())

model <- infer_trajectory(dataset, "paga")

# check out parameters of a method
?ti_mst

# infer a trajectory with some different parameters
model <- infer_trajectory(dataset, "mst", ndim = 4)

# Save models ------------------------------------------------------------------
saveRDS(model, "model.rds")
saveRDS(dataset, "dataset.rds")
