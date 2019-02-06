library(dyno)
library(dplyr)
library(ggplot2)

dataset <- readRDS("dataset.rds")
model <- readRDS("model.rds")

# Plotting ---------------------------------------------------------------------

# in dimensionality reduction
# will use by default the dimensionality reduction of the TI method
plot_dimred(model)

# trajectories can be projected onto any dimensionality reduction
dimred = dyndimred::dimred_pca(dataset$expression)
plot_dimred(model, dimred = dimred)

dimred = dyndimred::dimred_tsne(dataset$expression)
plot_dimred(model, dimred = dimred)

# it's possible that you'll have to install umap manually
# devtools::install_github("jlmelville/uwot")
dimred = uwot::umap(dataset$expression)
rownames(dimred) = rownames(dataset$expression)
plot_dimred(model, dimred = dimred)

dimred = uwot::umap(dataset$expression, n_neighbors = 300)
rownames(dimred) = rownames(dataset$expression)
plot_dimred(model, dimred = dimred)

# color cells according to pseudotime
# this uses the root of the trajectory if provided, see later
plot_dimred(model, color_cells = "pseudotime")

# plot a grouping/clustering
plot_dimred(model, grouping = dyno::fibroblast_reprogramming_treutlein$grouping)

# plot the grouping as background density
plot_dimred(
  model,
  color_cells = "grouping",
  color_density = "grouping",
  grouping = dyno::fibroblast_reprogramming_treutlein$grouping,
  label_milestones = FALSE
) + theme(legend.position = "none")

# plot the expression of a feature
plot_dimred(
  model,
  color_cells = "feature",
  feature_oi = "Vim",
  expression_source = dataset
)

# plotting the trajectory topology itself
# useful if the topology is too complex to be plotted inside a 2D dimensionality reduction
plot_graph(model)

# plot as a dendrogram
# useful for trees
plot_dendro(model)

# plot a heatmap (in some rstudio version, pressing zoom is necessary)
# see later for how these features are chosen
plot_heatmap(model, expression_source = dataset)

plot_heatmap(model, expression_source = dataset, features_oi = 500)

# Annotation -------------------------------------------------------------------
# you can annotate the trajectory model based on your biological knowledge

# rooting if you know the root milestone
plot_dimred(model, label_milestones = TRUE)
model <- model %>% add_root(root_milestone_id = "1")
plot_dimred(model, label_milestones = TRUE)

# rooting if you know a marker gene of the progenitor cells
model <- model %>% add_root_using_expression("Vim", expression_source = dataset)
plot_dimred(model, label_milestones = TRUE)

# labelling milestones
model <- model %>% label_milestones(c("3" = "start", "4" = "end1", "5" = "end2"))
plot_dimred(model, label_milestones = TRUE)

# labelling milestones using expression
model <- label_milestones_markers(
  model,
  markers = list(
    MEF = c("Vim"),
    Myocyte = c("Myl1"),
    Neuron = c("Stmn3")
  ),
  expression_source = dataset
)
plot_dimred(model)

# Trajectory differential expression -------------------------------------------
# dyno includes a small toolkit for doing trajectory differential expression
# it looks at how well the gene expression can predict certain aspects of the trajectory
# it well then rank these genes acording to their predictive performance (= feature importance)

# if you just want to know genes that change the strongest in the trajectory
# -> overall feature importance
overall_feature_importance <- calculate_overall_feature_importance(model, expression_source = dataset)

plot_heatmap(
  model,
  expression_source = dataset,
  features_oi = overall_feature_importance %>% top_n(40, importance) %>% pull(feature_id)
)

# find features that are predictive for a particular branch, either up or downregulated
branch_feature_importance <- calculate_branch_feature_importance(model, expression_source=dataset$expression)

neuron_features <- branch_feature_importance %>%
  filter(to == names(model$milestone_labelling)[which(model$milestone_labelling =="Neuron")]) %>%
  top_n(50, importance) %>%
  pull(feature_id)

plot_heatmap(
  model,
  expression_source = dataset$expression,
  features_oi = neuron_features
)

# features that are predictive for the branching point
# i.e. they change it the branching point
branching_milestone <- model$milestone_network %>% group_by(from) %>% filter(n() > 1) %>% pull(from) %>% first()

branch_feature_importance <- calculate_branching_point_feature_importance(model, expression_source=dataset$expression, milestones_oi = branching_milestone)

branching_point_features <- branch_feature_importance %>% top_n(20, importance) %>% pull(feature_id)

plot_heatmap(
  model,
  expression_source = dataset$expression,
  features_oi = branching_point_features
)

# for giggles, let's plot them in the dimensionality reduction
space <- dyndimred::dimred_mds(dataset$expression)
purrr::map(branching_point_features[1:9], function(feature_oi) {
  plot_dimred(model, dimred = space, expression_source = dataset$expression, feature_oi = feature_oi, label_milestones = FALSE) +
    ggplot2::theme(legend.position = "none") +
    ggplot2::ggtitle(feature_oi)
}) %>% patchwork::wrap_plots()


# Some final tips when inferring trajectories ----------------------------------

# Start slowly: select a subset of cells with a linear trajectory and try to interpret that
# Explore dimensionality reductions, try to see the trajectory yourself
# Do a small in silico validation by looking at the most relevant genes
# A trajectory is still a model - experimental validation is necessary
