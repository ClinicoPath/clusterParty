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
