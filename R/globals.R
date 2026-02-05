utils::globalVariables(c("lag_end", "lag_start", "value"))

.glmmtmb_family <- c("nbinom2", "nbinom1", "nbinom12", "compois", "truncated_compois",
                     "genpois", "truncated_genpois", "truncated_poisson", "truncated_nbinom2",
                     "truncated_nbinom1", "beta_family", "betabinomial", "tweedie", "skewnormal", "lognormal",
                     "ziGamma", "t_family", "ordbeta", "bell")

.glm_family <- c("gaussian", "binomial", "Gamma",
                   "inverse.gaussian", "poisson", "quasi", "quasibinomial", "quasipoisson")
