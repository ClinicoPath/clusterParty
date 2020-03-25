proximityClass <- if (requireNamespace("jmvcore")) R6::R6Class(
    "proximityClass",
    inherit = proximityBase,
    private = list(
        # ====================================================
        .init = function() {
            self$results$txtPfm$setVisible(FALSE)
            self$results$mtxSbj$setVisible(FALSE)
            self$results$mtxVar$setVisible(FALSE)
        },
        # ====================================================
        .run = function() {
            # get variables
            pxmVar = self$options$get("vars")

            if (length(pxmVar) > 1) {
                # get variables
                pxmLbl = self$options$get("label")
                mtxSbj = self$results$get("mtxSbj")
                mtxVar = self$results$get("mtxVar")

                # extract level name, method cetgory (simil. / dissim.), method name and description
                lvlNme = tolower(substr(self$options$get("lvlMsr"), 4, 6))
                mthNme = self$options$get(paste0(lvlNme, substr(self$options$get("disSim"), 4, 6)))
                mthDsc = paste0(nmePxm(mthNme), " (SPSS: ", spsPxm(mthNme), ")")

                # prepare data matrix (check for missing values, etc.)
                numSbj = dim(self$data)[1]
                dtaMtx = self$data[, pxmVar]
                sbjLbl = as.list(sprintf(paste0("C %", sprintf("%d", ceiling(log10(numSbj))), "d"), seq(1:numSbj)))
                if (length(pxmLbl) == 1) {
                    rawLbl = self$data[, pxmLbl]
                    if (length(levels(rawLbl)) == numSbj)
                        sbjLbl = as.list(as.character(rawLbl))
                }
                blnVld = !apply(is.na(dtaMtx), 1, any)
                dtaMtx = dtaMtx[blnVld, ]
                sbjLbl = sbjLbl[blnVld]
                numIvN = length(blnVld[!blnVld])
                # for binary data, ensure that these only contain the two categories defined by binAbs and binPrs
                numIvB = 0
                mthExp = c(0, 0)
                if      (lvlNme == "bin") {
                    binAbs = as.numeric(self$options$get("binAbs"))
                    binPrs = as.numeric(self$options$get("binPrs"))
                    # determine matrix positions for present and absent and replace them with 0 (for absent) and 1 (for present)
                    blnAbs = (dtaMtx == binAbs)
                    blnPrs = (dtaMtx == binPrs)
                    dtaMtx[blnAbs] = 0
                    dtaMtx[blnPrs] = 1
                    # determine lines that contain cell that are neither coded as present or absent and remove those
                    blnBnC = apply(blnAbs | blnPrs, 1, all)
                    dtaMtx = dtaMtx[blnBnC, ]
                    sbjLbl = sbjLbl[blnBnC]
                    numIvB = length(blnBnC[!blnBnC])
                    # adjust methods description (footnote in results table): replace placeholder for present and absent with their codes
                    mthDsc = gsub("_BA_", sprintf('%d', binAbs),    gsub("_BP_", sprintf('%d', binPrs),    mthDsc))
                }
                else if (lvlNme == "int") {
                    # assign method exponent and adjust methods description (footnote in results table): replace placeholder for method exponent with the values used 
                    mthExp = c(self$options$get("intPwr"), self$options$get("intRot"))
                    mthDsc = gsub("_IP_", sprintf('%d', mthExp[1]), gsub("_IR_", sprintf('%d', mthExp[2]), mthDsc))
                }

                # calculate proximity measures
                # dtaMtx is the data matrix with subjects in rows and variables in columns;
                #     to calculate similarities between subjects / cases the matrix needs to be transposed 
                if (self$options$get("btwDir") == "btwSbj")
                    dtaMtx = t(dtaMtx)
                # xfmRwM and xfmRwD indicate whether a standardization should be applied to the data values for either cases or variables
                # before computing proximities (NB: not applicable to binary data)
                # xfmRwM is a character, indication the transformation method: xfmNon - none, xfmZsc - z-standardization, 
                #     xfmRNP - range -1 to 1, xfmRZP - range 0 to 1, xfmMag - max. magnitude of 1, xfmAvr - mean of 1, xfmStd - std. dev. of 1
                xfmRwM = self$options$get("xfmMth")
                # xfmRwD indicates whether the transformation 
                xfmRwD = (substr(self$options$get("btwDir"), 4, 6) == substr(self$options$get("xfmDir"), 4, 6))
                # mthNme is composed of variable type ([1:3]: int, cnt, bin) and method ([4:6], e.g., Euc → Euclidian)
                # mthExp is only required for Minkowski and Custom within integer dissimilarities
                # xfmRes is a vector of booleans, indicating which transformation should be applied after calculating the proximity measures
                xfmRes = c(self$options$get("xfmAbs"), self$options$get("xfmInv"), self$options$get("xfmRsc"))
                resMtx = algPxm(dtaMtx, xfmRwM, xfmRwD, mthNme, mthExp, xfmRes)

                # prepare notes for the table footer (about included participants and used method) and assign results
                mthNte = paste0(gsub("clcSim", "Similarities ", gsub("clcDis", "Dissimilarities ", self$options$get("disSim"))), " between ", 
                                gsub("btwSbj", "cases",         gsub("btwVar", "variables",        self$options$get("btwDir"))), "; Method: ", mthDsc)
                incNte = paste0("Total number of cases: ", sprintf("%d", numSbj), "; excluded because of invalid data (NA): ", sprintf("%d", numIvN), 
                                ifelse(lvlNme == "bin", paste0(" and of invalid categories: ", sprintf("%d", numIvB)), ""), "; Included cases: ",    sprintf("%d", numSbj - numIvN - numIvB))

                if      (self$options$get("btwDir") == "btwVar") {
                    # WARNING: if there is a too high number of participants or variables
#                   if length(pxmVar) > 50 {
#                       # get an Preformatted 
#                       stop("your error message here")
#                   }

                    # initialize the results table (add the required number of columns)
                    for (i in seq_along(pxmVar)) {
                        mtxVar$addColumn(name=pxmVar[[i]], title=pxmVar[[i]], type="number", format="zto")
                    }

                    # populate results
                    for (i in seq_along(pxmVar)) {
                        values <- vector("list", length(pxmVar))
                        names(values) = pxmVar
                        values[] = ""
                        for (j in seq_len(i)) {
#                       for (j in seq_along(pxmVar)) {
#                           values[[pxmVar[[j]]]] = ifelse(i > j, resMtx[i, j], ifelse(i < j, "", "\u2014"))
#                           values[[pxmVar[[j]]]] = ifelse(i >= j, resMtx[i, j], "")
                            values[[pxmVar[[j]]]] = resMtx[i, j]
                        }
                        mtxVar$setRow(rowNo=i, values)
                    }
                    # add notes and make the table visible 
                    mtxVar$setNote("mthNte", mthNte)
                    mtxVar$setNote("incNte", incNte)
                    mtxVar$setVisible(TRUE)
                }
                else if (self$options$get("btwDir") == "btwSbj") {
                    # initialize the results table (add the required number of rows and columns)
                    for (i in seq_along(sbjLbl)) { mtxSbj$addColumn(name=sbjLbl[[i]], title=sbjLbl[[i]], type="number", format="zto") }
                    for (i in seq_along(sbjLbl)) { mtxSbj$addRow(rowKey=sbjLbl[[i]], list('.name' = sbjLbl[[i]])) }

                    # populate results
                    for (i in seq_along(sbjLbl)) {
                        values <- vector("list", length(sbjLbl))
                        names(values) = sbjLbl
                        values[] = ""
                        for (j in seq_len(i)) {
#                       for (j in seq_along(sbjLbl)) {
#                           values[[sbjLbl[[j]]]] = ifelse(i > j,  resMtx[i, j], ifelse(i < j, "", "\u2014"))
#                           values[[sbjLbl[[j]]]] = ifelse(i >= j, resMtx[i, j], "")
                            values[[sbjLbl[[j]]]] = resMtx[i, j]
                        }
                        mtxVar$setRow(rowNo=i, values)
                    }
                    # add notes and make the table visible 
                    mtxSbj$setNote("mthNte", mthNte)
                    mtxSbj$setNote("incNte", incNte)
                    mtxSbj$setVisible(TRUE)
                }
            }
        }
    )
)

# =====================================================================================================================
# Algorithms to calculate proximities (incl. transformations before and after this calculation)
# =====================================================================================================================
algPxm <- function(dtaMtx, xfmRwM, xfmRwD, mthNme, mthExp, xfmRes) {
    # =================================================================================================================
    # transforming data before calculating the proximity measures
    # =================================================================================================================
    # standardize data values for either cases or variables before computing proximities (NB: not applicable to binary data)
    # standardization methods are z scores, range –1 to 1, range 0 to 1, max. magnitude of 1, mean of 1, and std. dev. of 1
    if (mthNme[1:3] != 'bin' && xfmRwM != 'xfmNon') {
        xfmDim = ifelse(xfmRwD, 1, 2)
        # z-scores - Z - correct
        if      (xfmRwM == "xfmZsc")
            dtaMtx = sweep(sweep(dtaMtx, xfmDim, as.vector(apply(dtaMtx, xfmDim, mean)),  "-"), xfmDim, as.vector(apply(dtaMtx, xfmDim, sd)),          "/")
        # range -1 to 1 - RANGE - over subjects correct, over variables possibly wrong in SPSS
        else if (xfmRwM == "xfmRNP")
            dtaMtx = sweep(dtaMtx, xfmDim, as.vector(diff(apply(dtaMtx, xfmDim, range))), "/")
        # range 0 to 1 - RESCALE - over subjects correct, over variables possibly wrong in SPSS
        else if (xfmRwM == "xfmRZP")
            dtaMtx = sweep(sweep(dtaMtx, xfmDim, as.vector(apply(dtaMtx, xfmDim, min)),   "-"), xfmDim, as.vector(diff(apply(dtaMtx, xfmDim, range))), "/")
        # maximum magnitude of 1 - MAX - correct
        else if (xfmRwM == "xfmMag")
            dtaMtx = sweep(dtaMtx, xfmDim, as.vector(apply(dtaMtx, xfmDim, max)),  "/")
        # mean of 1 - MEAN - correct
        else if (xfmRwM == "xfmAvr")
            dtaMtx = sweep(dtaMtx, xfmDim, as.vector(apply(dtaMtx, xfmDim, mean)), "/")
        # standard deviation of 1 - SD - correct
        else if (xfmRwM == "xfmStd")
            dtaMtx = sweep(dtaMtx, xfmDim, as.vector(apply(dtaMtx, xfmDim, sd)),   "/")
    }

    # =================================================================================================================
    # integer measures that implemented using R functions (cor, dist)
    # =================================================================================================================
    # integer - similarity - Pearson correlation - CORRELATION
    if      (mthNme == "intCrr") {
        resMtx = cor(dtaMtx)
    }
    # integer - dissimilarity - Euclidian distance, squared Euclidian distance, Chebychev, Block, Minkowski, Customized
    else if (mthNme == "intEuc" || mthNme == "intSqE" || mthNme == "intChb" || mthNme == "intBlk" || mthNme == "intMnk" || mthNme == "intCst") {
       # distance measures from the R-function dist
       dstMth = ifelse(mthNme == "intEuc" || mthNme == "intSqE", "euclidian", ifelse(mthNme == "intChb", "maximum", ifelse(mthNme == "intBlk", "manhattan", ifelse(mthNme == "intMnk" || mthNme == "intCst", "minkowski", "error"))))
       dstExp = ifelse(mthNme == "intSqE", 2, ifelse(mthNme == "intCst", (mthExp[1] / mthExp[2]), 1))
       resMtx = as.matrix(dist(t(dtaMtx), upper=T, diag=T, method=dstMth, p = mthExp[1]) ^ dstExp)
    }
    # =================================================================================================================
    # integer, count and binary measures; implemented based upon:
    # www.ibm.com/support/knowledgecenter/SSLVMB_22.0.0/com.ibm.spss.statistics.algorithms/alg_proximities.htm
    # see www.itl.nist.gov/div898/software/dataplot/refman2/auxillar/binmatch.htm for further possible measures
    # =================================================================================================================
    else {
        # create result matrix
        resMtx = array(0, c(ncol(dtaMtx), ncol(dtaMtx)))
        dimnames(resMtx) = list(colnames(dtaMtx), colnames(dtaMtx))

        # preparation: binary measures
        if (mthNme == "binShp")
            dtaMtx = scale(dtaMtx, center=T, scale=F)

        # for most calculations, the main diagonal doesn't need to be calculated
        # however, there are three measures that represent an exception from that
        # those methods are the same where the main diagonal is divided by two at
        #  the very end of this algorithm
        if (mthNme == "binAnD" || mthNme == "binDsp" || mthNme == "binRnR")
            minI = 1
        else
            minI = 2

        for (i in minI:ncol(dtaMtx)) {
            for (j in 1:(i - minI + 1)) {
                # =====================================================================================================
                # calculate match vector and t1 / t2 required for binary measures
                # =====================================================================================================
                if (mthNme == "binAnD" || mthNme == "binDsp" || mthNme == "binKc2" || mthNme == "binLmb" || mthNme == "binOch" || mthNme == "binPh4" || mthNme == "binPtD" || 
                    mthNme == "binSk4" || mthNme == "binSk5" || mthNme == "binSzD" || mthNme == "binVar" || mthNme == "binYlQ" || mthNme == "binYlY") {
                    mtcVec = c(as.vector(table(dtaMtx[, i] == 1 & dtaMtx[, j] == 1)["TRUE"]), 
                               as.vector(table(dtaMtx[, i] == 1 & dtaMtx[, j] != 1)["TRUE"]),
                               as.vector(table(dtaMtx[, i] != 1 & dtaMtx[, j] == 1)["TRUE"]),
                               as.vector(table(dtaMtx[, i] != 1 & dtaMtx[, j] != 1)["TRUE"]))
                    if (mthNme == "binAnD" || mthNme == "binLmb") {
                        t1 = max(mtcVec[1], mtcVec[2], na.rm=T) + max(mtcVec[3], mtcVec[4], na.rm=T) + max(mtcVec[1], mtcVec[3], na.rm=T) + max(mtcVec[2], mtcVec[4], na.rm=T)
                        t2 = max(sum(mtcVec[1], mtcVec[3], na.rm=T), sum(mtcVec[2], mtcVec[4], na.rm=T), na.rm=T) + max(sum(mtcVec[1], mtcVec[2], na.rm=T), sum(mtcVec[3], mtcVec[4], na.rm=T), na.rm=T)
                    }
                }

                # =====================================================================================================
                # measures for integer data
                # =====================================================================================================
                # integer - similarity - Cosine - COSINE
                if      (mthNme == "intCos") {
                    resMtx[i, j] = crossprod(dtaMtx[, i], dtaMtx[, j]) / sqrt(crossprod(dtaMtx[, i]) * crossprod(dtaMtx[, j]))
                }
#               PLACEHOLDER FOR FUTURE IMPLEMENTATIONS
#               else if (mthNme == "int") {
#               }    
                # =====================================================================================================
                # measures for count data
                # =====================================================================================================
                # dissimilarities - counts - Chi-square measure: disChi - CHISQ
                else if (mthNme == "cntChi") {
                    resMtx[i, j] = NA # to implement
                }
                # dissimilarities - counts - Phi-square measure: disPhi - PH2
                else if (mthNme == "cntPhi") {
                    resMtx[i, j] = NA # to implement
                }
#               PLACEHOLDER FOR FUTURE IMPLEMENTATIONS
#               else if (mthNme == "cnt") {
#               }
                # =====================================================================================================
                # measures for binary data
                # =====================================================================================================
                # binary - similarity - Anderberg's D: binAnD - D
                else if (mthNme == "binAnD") {
                    resMtx[i, j] = (t1 - t2) / (2 * nrow(dtaMtx))
                }
                # binary - similarity - Dice: binDic - DICE
                else if (mthNme == "binDic") {
                    resMtx[i, j] = (as.vector(table(dtaMtx[, i] & dtaMtx[, j])["TRUE"]) * 2) / (as.vector(table(dtaMtx[, i] | dtaMtx[, j])["TRUE"]) + as.vector(table(dtaMtx[, i] & dtaMtx[, j])["TRUE"]))
                }
                # binary - similarity - Dispersion: binDsp - DISPER
                else if (mthNme == "binDsp") {
                    resMtx[i, j] = ((mtcVec[1] * mtcVec[4]) - ifelse(is.na(mtcVec[2]) || is.na(mtcVec[3]), 0, (mtcVec[2] * mtcVec[3]))) / (nrow(dtaMtx) ^ 2)
                }
                # binary - dissimilarity - Euclidian distance: binEuc - BEUCLID
                else if (mthNme == "binEuc") {
                    resMtx[i, j] = sqrt(as.vector(table(dtaMtx[, i] != dtaMtx[, j])["TRUE"]))
                }
                # binary - similarity - Jaccard: binJcc - JACCARD
                else if (mthNme == "binJcc") {
                    resMtx[i, j] = as.vector(table(dtaMtx[, i] & dtaMtx[, j])["TRUE"]) / as.vector(table(dtaMtx[, i] | dtaMtx[, j])["TRUE"])
                }
                # binary - similarity - Hamann: binHmn - HAMANN
                else if (mthNme == "binHmm") {
                    resMtx[i, j] = (as.vector(table(dtaMtx[, i] == dtaMtx[, j])["TRUE"]) - as.vector(table(dtaMtx[, i] != dtaMtx[, j])["TRUE"])) / nrow(dtaMtx)
                }
                # binary - similarity - Kulczynski 1: binKc1 - K1
                else if (mthNme == "binLmb") {
                    resMtx[i, j] = (t1 - t2) / (2 * nrow(dtaMtx) - t2)
                }
                # binary - similarity - Kulczynski 2: binKc2 - K2
                else if (mthNme == "binKc1") {
                    resMtx[i, j] = as.vector(table(dtaMtx[, i] & dtaMtx[, j])["TRUE"]) / (as.vector(table(dtaMtx[, i] | dtaMtx[, j])["TRUE"]) - as.vector(table(dtaMtx[, i] & dtaMtx[, j])["TRUE"]))
                }
                # binary - similarity - Lambda: binLmb - LAMBDA
                else if (mthNme == "binKc2") {
                    resMtx[i, j] = (mtcVec[1] / (mtcVec[1] + mtcVec[2]) + mtcVec[1] / (mtcVec[1] + mtcVec[3])) / 2
                }
                # binary - dissimilarity - Lance and Williams: binLnW - BLWMN
                else if (mthNme == "binLnW") {
                    resMtx[i, j] = sum(abs(dtaMtx[, i] - dtaMtx[, j])) / sum(dtaMtx[, i] + dtaMtx[, j])
                }
                # binary - similarity - Ochiai: binOch - OCHIAI
                else if (mthNme == "binOch") {
                    resMtx[i, j] = sqrt((mtcVec[1] / (mtcVec[1] + mtcVec[2])) * (mtcVec[1] / (mtcVec[1] + mtcVec[3])))
                }
                # binary - similarity - Phi 4-point correlation: binPh4 - PHI
                else if (mthNme == "binPh4") {
                    resMtx[i, j] = ((mtcVec[1] * mtcVec[4]) - (mtcVec[2] * mtcVec[3])) / sqrt((mtcVec[1] + mtcVec[2]) * (mtcVec[1] + mtcVec[3]) * (mtcVec[2] + mtcVec[4]) * (mtcVec[3] + mtcVec[4]))
                }
                # binary - dissimilarity - Pattern difference: binPtD - PATTERN
                else if (mthNme == "binPtD") {
                    resMtx[i, j] = (mtcVec[2] * mtcVec[3]) / (nrow(dtaMtx) ^ 2)
                }
                # binary - similarity - Russel and Rao: binRnR - RR
                else if (mthNme == "binRnR") {
                    resMtx[i, j] = as.vector(table(dtaMtx[, i] & dtaMtx[, j])["TRUE"]) / nrow(dtaMtx)
                }
                # binary - similarity - Rogers and Tanimoto: binRnT - RT
                else if (mthNme == "binRnT") {
                    resMtx[i, j] = as.vector(table(dtaMtx[, i] == dtaMtx[, j])["TRUE"]) / (as.vector(table(dtaMtx[, i] == dtaMtx[, j])["TRUE"]) + 2 * as.vector(table(dtaMtx[, i] == dtaMtx[, j])["FALSE"]))
                }
                # binary - dissimilarity - Shape: binShp - BSHAPE
                else if (mthNme == "binShp") {
                    resMtx[i, j] = sum((dtaMtx[, i] - dtaMtx[, j]) ^ 2) / nrow(dtaMtx)
                }
                # binary - similarity - Sokal and Sneath 1: binSk1 - SS1
                else if (mthNme == "binSk1") {
                    resMtx[i, j] = (2 * as.vector(table(dtaMtx[, i] == dtaMtx[, j])["TRUE"])) / (2 * as.vector(table(dtaMtx[, i] == dtaMtx[, j])["TRUE"]) + as.vector(table(dtaMtx[, i] == dtaMtx[, j])["FALSE"]))
                }
                # binary - similarity - Sokal and Sneath 2: binSk2 - SS2
                else if (mthNme == "binSk2") {
                    resMtx[i, j] = as.vector(table(dtaMtx[, i] & dtaMtx[, j])["TRUE"]) / (2 * as.vector(table(dtaMtx[, i] | dtaMtx[, j])["TRUE"]) - as.vector(table(dtaMtx[, i] & dtaMtx[, j])["TRUE"]))
                }
                # binary - similarity - Sokal and Sneath 3: binSk3 - SS3
                else if (mthNme == "binSk3") {
                    resMtx[i, j] = as.vector(table(dtaMtx[, i] == dtaMtx[, j])["TRUE"]) / as.vector(table(dtaMtx[, i] == dtaMtx[, j])["FALSE"])
                }
                # binary - similarity - Sokal and Sneath 4: binSk4 - SS4
                else if (mthNme == "binSk4") {
                    resMtx[i, j] = (mtcVec[1] / (mtcVec[1] + mtcVec[2]) + mtcVec[1] / (mtcVec[1] + mtcVec[3]) + mtcVec[4] / (mtcVec[2] + mtcVec[4]) + mtcVec[4] / (mtcVec[3] + mtcVec[4])) / 4
                }
                # binary - similarity - Sokal and Sneath 5: binSk5 - SS5
                else if (mthNme == "binSk5") {
                    resMtx[i, j] = (mtcVec[1] * mtcVec[4]) / sqrt((mtcVec[1] + mtcVec[2]) * (mtcVec[1] + mtcVec[3]) * (mtcVec[2] + mtcVec[4]) * (mtcVec[3] + mtcVec[4]))
                }
                # binary - similarity - Simple matching: binSmM - SM
                else if (mthNme == "binSmM") {
                    resMtx[i, j] = as.vector(table(dtaMtx[, i] == dtaMtx[, j])["TRUE"]) / nrow(dtaMtx)
                }
                # binary - dissimilarity - Squared Euclidian distance: binSqE - BSEUCLID
                else if (mthNme == "binSqE") {
                    resMtx[i, j] = as.vector(table(dtaMtx[, i] != dtaMtx[, j])["TRUE"])
                }
                # binary - dissimilarity - Size difference: binSzD - SIZE
                else if (mthNme == "binSzD") {
                    resMtx[i, j] = (mtcVec[2] - mtcVec[3]) ^ 2 / (nrow(dtaMtx) ^ 2)
                }
                # binary - dissimilarity - Variance: binVar - VARIANCE
                else if (mthNme == "binVar") {
                    resMtx[i, j] = (mtcVec[2] + mtcVec[3]) / (4 * nrow(dtaMtx))
                }
                # binary - similarity - Yule's Q: binYlQ - Q
                else if (mthNme == "binYlQ") {
                    resMtx[i, j] = (mtcVec[1] * mtcVec[4] - mtcVec[2] * mtcVec[3]) / (mtcVec[1] * mtcVec[4] + mtcVec[2] * mtcVec[3])
                }
                # binary - similarity - Yule's Y: binYlY - Y
                else if (mthNme == "binYlY") {
                    resMtx[i, j] = (sqrt(mtcVec[1] * mtcVec[4]) - sqrt(mtcVec[2] * mtcVec[3])) / (sqrt(mtcVec[1] * mtcVec[4]) + sqrt(mtcVec[2] * mtcVec[3]))
                }
#               PLACEHOLDER FOR FUTURE IMPLEMENTATIONS
#               else if (mthNme == "bin") {
#               }
            }
        }

        # =============================================================================================================
        # transpose and add the matrix (only the bottom triangle is calculated, th upper is "mirrored")
        # =============================================================================================================
        resMtx = resMtx + t(resMtx)
        # handle the main diagonal
        if      (mthNme == "binEuc" || mthNme == "binLnW" || mthNme == "binPtD" || mthNme == "binShp" || mthNme == "binSqE" || mthNme == "binSzD" || mthNme == "binVar")
            diag(resMtx) = 0
        else if (mthNme == "intCos" || mthNme == "cntChi" || mthNme == "cntChi" || mthNme == "binDic" || mthNme == "binJcc" || mthNme == "binHmm" || mthNme == "binKc2" || 
                 mthNme == "binLmb" || mthNme == "binOch" || mthNme == "binPh4" || mthNme == "binRnT" || mthNme == "binSmM" || mthNme == "binSk1" || mthNme == "binSk2" ||
                 mthNme == "binSk4" || mthNme == "binSk5" || mthNme == "binYlY" || mthNme == "binYlQ")
            diag(resMtx) = 1
        else if (mthNme == "binKc1" || mthNme == "binSk3")
            diag(resMtx) = NA
        else if (mthNme == "binAnD" || mthNme == "binDsp" || mthNme == "binRnR")
            diag(resMtx) = diag(resMtx) / 2
    
    }
    # =================================================================================================================
    # end: own implementations
    # =================================================================================================================

    # =================================================================================================================
    # transform the values generated by the distance measure, applied after the distance measure has been computed
    # available options are [1] absolute values, [2] change sign, and [3] rescale to 0–1 range.
    # =================================================================================================================
    if (xfmRes[1])
        resMtx = abs(resMtx)
    if (xfmRes[2])
        resMtx[row(resMtx) != col(resMtx)] = -resMtx[row(resMtx) != col(resMtx)]
    if (xfmRes[3])
        resMtx[row(resMtx) != col(resMtx)] = (resMtx[row(resMtx) != col(resMtx)] - min(resMtx[row(resMtx) != col(resMtx)])) / diff(range(resMtx[row(resMtx) != col(resMtx)]))

    # =================================================================================================================
    # return matrix with results
    # =================================================================================================================
    return(resMtx)
}

# =====================================================================================================================
# Description of the proximity measure
# =====================================================================================================================
nmePxm <- function(mthNme) {
    # measures for interval data
    if      (mthNme == "intBlk")
        resNme = "Block"
    else if (mthNme == "intChb")
        resNme = "Chebychev"
    else if (mthNme == "intCst")
        resNme = "Customized, based upon Minkowski with Power: _IP_, Root: _IR_"
    else if (mthNme == "intCos")
        resNme = "Cosine"
    else if (mthNme == "intCrr")
        resNme = "Pearson correlation"
    else if (mthNme == "intEuc")
        resNme = "Euclidian distance"
    else if (mthNme == "intMnk")
        resNme = "Minkowski (Power: _IP_)"
    else if (mthNme == "intSqE")
        resNme = "Squared Euclidian distance"
    # count / frequency measures
    else if (mthNme == "cntChi")
        resNme = "Chi-squared measure"
    else if (mthNme == "cntPhi")
        resNme = "Phi-squared measure"
    # binary measures
    else if (mthNme == "binAnD")
        resNme = "Anderberg's D; present: _BP_, absent: _BA_"
    else if (mthNme == "binDic")
        resNme = "Dice; present: _BP_, absent: _BA_"
    else if (mthNme == "binDsp")
        resNme = "Dispersion; present: _BP_, absent: _BA_"
    else if (mthNme == "binEuc")
        resNme = "(Binary) Euclidian distance; present: _BP_, absent: _BA_"
    else if (mthNme == "binHmn")
        resNme = "Hamann; present: _BP_, absent: _BA_"
    else if (mthNme == "binJcc")
        resNme = "Jaccard; present: _BP_, absent: _BA_"
    else if (mthNme == "binKc1")
        resNme = "Kulczynski 1; present: _BP_, absent: _BA_"
    else if (mthNme == "binKc2")
        resNme = "Kulczynski 2; present: _BP_, absent: _BA_"
    else if (mthNme == "binLmb")
        resNme = "Lambda; present: _BP_, absent: _BA_"
    else if (mthNme == "binLnW")
        resNme = "Lance and Williams; present: _BP_, absent: _BA_"
    else if (mthNme == "binOch")
        resNme = "Ochiai; present: _BP_, absent: _BA_"
    else if (mthNme == "binPh4")
        resNme = "Phi 4-point correlation; present: _BP_, absent: _BA_"
    else if (mthNme == "binPtD")
        resNme = "Pattern difference; present: _BP_, absent: _BA_"
    else if (mthNme == "binRnR")
        resNme = "Russel and Rao; present: _BP_, absent: _BA_"
    else if (mthNme == "binRnT")
        resNme = "Rogers and Tanimoto; present: _BP_, absent: _BA_"
    else if (mthNme == "binShp")
        resNme = "Shape; present: _BP_, absent: _BA_"
    else if (mthNme == "binSk1")
        resNme = "Sokal and Sneath 1; present: _BP_, absent: _BA_"
    else if (mthNme == "binSk2")
        resNme = "Sokal and Sneath 2; present: _BP_, absent: _BA_"
    else if (mthNme == "binSk3")
        resNme = "Sokal and Sneath 3; present: _BP_, absent: _BA_"
    else if (mthNme == "binSk4")
        resNme = "Sokal and Sneath 4; present: _BP_, absent: _BA_"
    else if (mthNme == "binSk5")
        resNme = "Sokal and Sneath 5; present: _BP_, absent: _BA_"
    else if (mthNme == "binSmM")
        resNme = "Simple matching; present: _BP_, absent: _BA_"
    else if (mthNme == "binSqE")
        resNme = "(Binary) squared Euclidian distance; present: _BP_, absent: _BA_"
    else if (mthNme == "binSzD")
        resNme = "Size difference; present: _BP_, absent: _BA_"
    else if (mthNme == "binVar")
        resNme = "(Binary) variance; present: _BP_, absent: _BA_"
    else if (mthNme == "binYlQ")
        resNme = "Yule's Q; present: _BP_, absent: _BA_"
    else if (mthNme == "binYlY")
        resNme = "Yule's Y; present: _BP_, absent: _BA_"
#   else if (mthNme == "")
#       resNme = ""
    else
        resNme = ""

    return(resNme)
}

# =====================================================================================================================
# Information regarding the equivalent SPSS proximity measure
# =====================================================================================================================
spsPxm <- function(mthNme) {
    # measures for interval data
    if      (mthNme == "intBlk")
        spssMt = "BLOCK"
    else if (mthNme == "intChb")
        spssMt = "CHEBYCHEV"
    else if (mthNme == "intCst")
        spssMt = "POWER(_IP_, _IR_)"
    else if (mthNme == "intCos")
        spssMt = "COSINE"
    else if (mthNme == "intCrr")
        spssMt = "CORRELATION"
    else if (mthNme == "intEuc")
        spssMt = "EUCLID"
    else if (mthNme == "intMnk")
        spssMt = "MINKOWSKI(_IP_)"
    else if (mthNme == "intSqE")
        spssMt = "SEUCLID"
    # count / frequency measures
    else if (mthNme == "cntChi")
        spssMt = "CHISQ"
    else if (mthNme == "cntPhi")
        spssMt = "PH2"
    # binary measures
    else if (mthNme == "binAnD")
        spssMt = "D(_BP_, _BA_)"
    else if (mthNme == "binDic")
        spssMt = "DICE(_BP_, _BA_)"
    else if (mthNme == "binDsp")
        spssMt = "DISPER(_BP_, _BA_)"
    else if (mthNme == "binEuc")
        spssMt = "BEUCLID(_BP_, _BA_)"
    else if (mthNme == "binHmn")
        spssMt = "HAMANN(_BP_, _BA_)"
    else if (mthNme == "binJcc")
        spssMt = "JACCARD(_BP_, _BA_)"
    else if (mthNme == "binKc1")
        spssMt = "K1(_BP_, _BA_)"
    else if (mthNme == "binKc2")
        spssMt = "K2(_BP_, _BA_)"
    else if (mthNme == "binLmb")
        spssMt = "LAMBDA(_BP_, _BA_)"
    else if (mthNme == "binLnW")
        spssMt = "BLWMN(_BP_, _BA_)"
    else if (mthNme == "binOch")
        spssMt = "OCHIAI(_BP_, _BA_)"
    else if (mthNme == "binPh4")
        spssMt = "PHI(_BP_, _BA_)"
    else if (mthNme == "binPtD")
        spssMt = "PATTERN(_BP_, _BA_)"
    else if (mthNme == "binRnR")
        spssMt = "RR(_BP_, _BA_)"
    else if (mthNme == "binRnT")
        spssMt = "RT(_BP_, _BA_)"
    else if (mthNme == "binShp")
        spssMt = "BSHAPE(_BP_, _BA_)"
    else if (mthNme == "binSk1")
        spssMt = "SS1(_BP_, _BA_)"
    else if (mthNme == "binSk2")
        spssMt = "SS2(_BP_, _BA_)"
    else if (mthNme == "binSk3")
        spssMt = "SS3(_BP_, _BA_)"
    else if (mthNme == "binSk4")
        spssMt = "SS4(_BP_, _BA_)"
    else if (mthNme == "binSk5")
        spssMt = "SS5(_BP_, _BA_)"
    else if (mthNme == "binSmM")
        spssMt = "SM(_BP_, _BA_)"
    else if (mthNme == "binSqE")
        spssMt = "BSEUCLID(_BP_, _BA_)"
    else if (mthNme == "binSzD")
        spssMt = "SIZE(_BP_, _BA_)"
    else if (mthNme == "binVar")
        spssMt = "VARIANCE(_BP_, _BA_)"
    else if (mthNme == "binYlQ")
        spssMt = "Q(_BP_, _BA_)"
    else if (mthNme == "binYlY")
        spssMt = "Y(_BP_, _BA_)"
#   else if (mthNme == "")
#       spssMt = ""
    else
        spssMt = ""

    return(spssMt)
}
