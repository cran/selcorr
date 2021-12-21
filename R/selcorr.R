#' Post-Selection Inference for Generalized Linear Models
#'
#' \code{selcorr} returns (unconditional) post-selection confidence intervals and p-values
#' for the coefficients of (generalized) linear models.
#'
#' When \code{boot.repl = 0}, an approximate asymptotic distribution of the test statistic
#' is used to calculate p-values and calibrate the profile-likelihood confidence intervals.
#' This approach is faster, but p-values and confidence intervals can be more precisely calibrated
#' by parametrically bootstrapping the test statistic (with \code{boot.repl} the number of replicates).
#' Parallel computing can be used to speed up the bootstrapping: see Examples.
#'
#' @param object an object representing a model of an appropriate class.
#' This is used as the initial model in a (bidirectional) stepwise model selection.
#' @param fixed.vars the names of all independent variables that must be included
#' in the selected model. The default is none.
#' @param further.vars the names of all independent variables that can be included
#' in the selected model, but are not part of \code{object}. The default is none.
#' @param boot.repl a number or list of bootstrap replicates. The default is no bootstrapping.
#' See Details and Examples for clarification.
#' @param k the multiple of the number of degrees of freedom used as penalty in the model selection.
#' The default \code{k = 2} corresponds to the AIC.
#' @param conf.level the level of the confidence intervals.
#' @param quiet if \code{TRUE}, then \code{selcorr} does not generate an output.
#'
#' @return the selected model is returned, without correction for model-selection,
#' but with up to two additional components. There is an \code{output} component corresponding to
#' the post-selection inference, which is also printed unless \code{quiet = TRUE}.
#' When \code{boot.repl} is not \code{0}, there is also a \code{boot.repl} component
#' corresponding to the bootstrap replicates.
#'
#' @examples
#' ## linear regression:
#' selcorr(lm(Fertility ~ ., swiss))
#'
#' ## logistic regression:
#' swiss.lr = within(swiss, Fertility <- (Fertility > 70))
#' selcorr(glm(Fertility ~ ., binomial, swiss.lr))
#'
#' ## parallel bootstrapping:
#' \dontrun{
#' library(future.apply)
#' plan(multisession)
#' boot.repl = future_replicate(8, selcorr(lm(Fertility ~ ., swiss), boot.repl = 1000,
#'                                         quiet = TRUE)$boot.repl, simplify = FALSE)
#' plan(sequential)
#' selcorr(lm(Fertility ~ ., swiss), boot.repl = do.call("rbind", boot.repl))}
#' @export
#' @importFrom methods is
#' @importFrom stats anova approx confint drop1 extractAIC formula pnorm
#' qnorm setNames simulate step update
selcorr = function(object, fixed.vars = NULL, further.vars = NULL,
                   boot.repl = 0, k = 2, conf.level = 0.95, quiet = FALSE) {
  suppressWarnings(suppressMessages({
    pos.part = \(x) ifelse(x > 0, x, 0)
    potential.vars = union(all.vars(formula(object)[[3]]), further.vars)
    if (!all(apply(update(object, paste(c("~ .", further.vars),
                                        collapse = " + "))$model[, potential.vars],
                   2, is.numeric))) stop("All independent variables should be numeric.")
    modaic = step(object, list(lower = paste(c("~ 1", fixed.vars), collapse = " + "),
                               upper = paste(c("~ .", further.vars), collapse = " + ")), trace = 0)
    output = summary(modaic)$coef
    colnames(output)[2:4] = c("CI lower", "CI upper", "p")
    if ("(Intercept)" %in% rownames(output)) {
      output["(Intercept)", 2:3] = confint(modaic, "(Intercept)", conf.level)
      if (is(modaic, "glm")) {
        output["(Intercept)", 4] = anova(modaic, update(modaic, ~ . -1),
                                         test = "Chisq")[2, "Pr(>Chi)"]}}
    modaic.wos = lapply(setdiff(rownames(output), "(Intercept)"), \(iv) {
      step(update(modaic, paste("~ . -", iv)),
           list(lower = paste(c("~ 1", setdiff(fixed.vars, iv)), collapse = " + "),
                upper = paste(c("~ .", setdiff(potential.vars, iv)), collapse = " + ")), trace = 0)})
    if (length(boot.repl) == 1 && boot.repl == 0) {
      output[setdiff(rownames(output), "(Intercept)"), 4] = sapply(modaic.wos, \(modaic.wo) {
        2 * pnorm(-sqrt(pos.part(extractAIC(modaic.wo, k = k)[2] - extractAIC(modaic, k = k)[2]) +
                          k))})
    } else {
      if (length(boot.repl) == 1) {
        modaic.wos.vars = lapply(modaic.wos, \(modaic.wo) sort(all.vars(formula(modaic.wo)[[3]])))
        modaic.sims = lapply(unique(modaic.wos.vars), \(modaic.wo.vars) {
          modaic.wo = modaic.wos[[which.max(sapply(modaic.wos.vars, identical, modaic.wo.vars))]]
          lapply(simulate(modaic.wo, ceiling(abs(boot.repl))), \(sim) {
            datasim = cbind(setNames(as.data.frame(sim), all.vars(formula(object)[[2]])),
                            update(object, paste(c("~ .", further.vars),
                                                 collapse = " + "))$model[potential.vars])
            modaic.wo.up = update(modaic.wo, paste(c("~ .", fixed.vars), collapse = " + "),
                                  data = datasim)
            attr(modaic.wo.up$terms, ".Environment") = environment()
            step(modaic.wo.up, list(lower = paste(c("~ 1", fixed.vars), collapse = " + "),
                                    upper = paste(c("~ .", potential.vars), collapse = " + ")),
                 trace = 0)})})
        simdiffs = sapply(setdiff(rownames(output), "(Intercept)"), \(iv) {
          sapply(modaic.sims[[which.max(sapply(unique(modaic.wos.vars), identical,
                                               sort(all.vars(formula(modaic.wos[[which.max(
                                                 iv == setdiff(rownames(output),
                                                               "(Intercept)"))]])[[3]]))))]],
                 \(modaic.sim) {
                   datasim = cbind(modaic.sim$model[all.vars(formula(object)[[2]])],
                                   update(object, paste(c("~ .", further.vars),
                                                        collapse = " + "))$model[potential.vars])
                   attr(modaic.sim$terms, ".Environment") = environment()
                   modaic.sim.wo = step(update(modaic.sim, paste("~ . -", iv)),
                                        list(lower = paste(c("~ 1", setdiff(fixed.vars, iv)),
                                                           collapse = " + "),
                                             upper = paste(c("~ .", setdiff(potential.vars, iv)),
                                                           collapse = " + ")), trace = 0)
                   ifelse(iv %in% all.vars(formula(modaic.sim)[[3]]),
                          pos.part(extractAIC(modaic.sim.wo, k = k)[2] -
                                     extractAIC(modaic.sim, k = k)[2]), 0)})})
      } else simdiffs = boot.repl
      output[setdiff(rownames(output), "(Intercept)"), 4] =
        sapply(setdiff(rownames(output), "(Intercept)"), \(iv) {
          approx(c(0, 2 * pnorm(-sqrt(sort(simdiffs[, iv], decreasing = TRUE) + k))),
                 0:length(simdiffs[, iv]) / length(simdiffs[, iv]),
                 2 * pnorm(-sqrt(pos.part(extractAIC(modaic.wos[[
                   which.max(iv == setdiff(rownames(output), "(Intercept)"))]], k = k)[2] -
                     extractAIC(modaic, k = k)[2]) + k)), ties = "min", rule = 2)$y})
      modaic$boot.repl = simdiffs}
    output[setdiff(rownames(output), "(Intercept)"), 2:3] =
      t(sapply(setdiff(rownames(output), "(Intercept)"), \(iv) {
        q = max(1e-300, 2 * pnorm(-sqrt(pos.part(diff(drop1(modaic, iv, k = k)$AIC)) + k) *
                                    qnorm((1 - conf.level) / 2) / qnorm(output[iv, 4] / 2)))
        r = 1e-15
        if (q > r) ci = confint(modaic, iv, 1 - q) else ci = NA
        while(any(is.na(ci)) & r < 1) {
          ci = output[iv, 1] + (confint(modaic, iv, 1 - r) - output[iv, 1]) *
            qnorm(q / 2) / qnorm(r / 2)
          r = r * 10}
        ci}))}))
  modaic$output = output
  if (!quiet) print(output)
  invisible(modaic)}
