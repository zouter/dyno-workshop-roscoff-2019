# Run this on n59
# > ssh bioinfo.sb-roscoff.fr
# > qlogin -q formation.q@@bignode.formation # login onto n59
# > export PATH=/opt/6.x/R-3.5.1/bin/:$PATH # add R libraries to PATH

library(dyno)
library(dplyr)

data("fibroblast_reprogramming_treutlein")

# we'll use an example dataset which is quite small
# most methods will probably take much longer on your dataset
dim(fibroblast_reprogramming_treutlein$counts)

# Prepare the data -------------------------------------------------------------

# some methods require pure counts, others require normalized expression
# we therefore have to provide both
dataset <- wrap_expression(
  counts = fibroblast_reprogramming_treutlein$counts,
  expression = fibroblast_reprogramming_treutlein$expression
)

# To select genes or not to select genes? That is the questions
# Pro:
# - Some methods don't scale well with a lot of genes
# - Some methods include their own feature selection
# - Noise can make it really difficult to infer a trajectory
# Contra:
# - Every step adds yet another bias

# As Antonio said:
# => Try first without filtering, if that doesn't work try with filtering

# Choosing what method to run -------------------------------------------------
# Go to http://guidelines.dynverse.org
# When running this on your laptop, you can do `guidelines_shiny()`

# you can find a list of all available methods using
dynmethods::methods %>% select(id, name, source)

# Running methods --------------------------------------------------------------

# infer a trajectory using default parameters
# this will download the container of the method first
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
