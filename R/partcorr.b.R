# This file is a generated template, your changes will not be overwritten

# `self$data` contains the data
# `self$options` contains the options
# ..corrvars - variables to calculate correlations for
# ..ctrlvars - variables to control for
# ..shwSig  - show significance level for correlations (p-values)
# ..flgSig  - flag significant correlations
# ..sidSig  - one- or two-tailed significance calculations
# `self$results` contains the results object (to populate)

partcorrClass <- if (requireNamespace('jmvcore')) R6::R6Class(
    "partcorrClass",
    inherit = partcorrBase,
    private = list(
        # ====================================================
        .init = function() {
            # get variables
            matrix <- self$results$get('matrix')
            varCrr <- self$options$get('corrvars')
            varCtl <- self$options$get('ctrlvars')

            # set title according to whether the procedure is controlling for variables or not
            matrix$setTitle(ifelse(length(varCtl) > 0, 'Partial Correlation Matrix', 'Correlation Matrix'))

            # initialize the results table (add the required number of columns)
            for (i in seq_along(varCrr)) {
                matrix$addColumn(name=paste0(varCrr[[i]], '[r]'),  title=varCrr[[i]], type='number', format='zto')
                matrix$addColumn(name=paste0(varCrr[[i]], '[rp]'), title=varCrr[[i]], type='number', format='zto,pvalue', visible='(shwSig)')
            }

            # initialize the results table (empty cells above and put "-" in the main diagonal)
            for (i in seq_along(varCrr)) {
                values <- list()
                for (j in seq(i, length(varCrr))) {
                    values[[paste0(varCrr[[j]], '[r]')]]  <- ''
                    values[[paste0(varCrr[[j]], '[rp]')]] <- ''
                }
                values[[paste0(varCrr[[i]], '[r]')]]  <- '\u2014'
                values[[paste0(varCrr[[i]], '[rp]')]] <- '\u2014'
                matrix$setRow(rowKey=varCrr[[i]], values)
            }
            # initialize the results table (assign the notes underneath)
            matrix$setNote('ctlNte', ifelse(length(varCtl) > 0, paste0('Controlling for ', paste(varCtl, collapse=", ")), 'Not controlling for any variables, i.e. the table shows correlations'))
            matrix$setNote('sigNte', paste0(ifelse(self$options$get('sidSig') == 'onetailed', 'One-tailed significance', 'Two-tailed significance'), ifelse(self$options$get('flgSig'), ': * p < .05, ** p < .01, *** p < .001', '')))
        },
        # ====================================================
        .run = function() {
            # get variables
            matrix <- self$results$get('matrix')
            varCrr <- self$options$get('corrvars')
            lngCrr <- length(varCrr)
            varCtl <- self$options$get('ctrlvars')
            lngCtl <- length(varCtl)

            # there are at least two variables required to calculate a (partial) correlation
            # if the variable to control for is empty, calculate a standard correlation
            if (lngCrr > 1) {
               m  <- cor(do.call("cbind", lapply(self$data[, c(varCrr, varCtl)], jmvcore::toNumeric)), use='pairwise', method='pearson')
               X  <- m[varCrr, varCrr]
               if (lngCtl > 0) {
     	           Y  <-       m[varCrr, varCtl]
                   pi <- solve(m[varCtl, varCtl])
                   Rp <- cov2cor(X - Y %*% pi %*% t(Y))
               } else {
                   Rp <- X
               }

               df <- dim(self$data)[1] - lngCtl
               Rt <- (Rp * sqrt(df - 2)) / sqrt(1 - Rp ^ 2)
               if (self$options$sidSig == 'onetailed') { nt = 1 } else { nt = 2 }
               Pp <- -nt *  expm1(pt(abs(Rt), (df - 2), log.p = TRUE))

               # populate results
               for (i in 2:lngCrr) {
                    for (j in seq_len(i - 1)) {
                        values <- list()
                        values[[paste0(varCrr[[j]], '[r]')]]  <- Rp[i, j]
                        values[[paste0(varCrr[[j]], '[rp]')]] <- Pp[i, j]
                        matrix$setRow(rowNo=i, values)
                        if (self$options$get('flgSig')) {
                            if      (Pp[i, j] < .001)
                                matrix$addSymbol(rowNo=i, paste0(varCrr[[j]], '[r]'), '***')
                            else if (Pp[i, j] < .01)
                                matrix$addSymbol(rowNo=i, paste0(varCrr[[j]], '[r]'), '**')
                            else if (Pp[i, j] < .05)
                                matrix$addSymbol(rowNo=i, paste0(varCrr[[j]], '[r]'), '*')
                        }
                    }
                }
            }
        }
    )
)
