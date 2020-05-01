#ifndef SEQR_RESULT_CREATOR_H
#define SEQR_RESULT_CREATOR_H

#include <Rcpp.h>
// [[Rcpp::depends(RcppParallel)]]
//' @importFrom  RcppParallel RcppParallelLibs
#include<RcppParallel.h>
#include "kmer_manager.h"
#include "hash/custom_hashers.h"

template<template<typename key, typename value, typename...> class kmer_dictionary_t>
class KMerMatrixCreatorWorker : public RcppParallel::Worker {
public:
    Rcpp::IntegerMatrix outputKMerCounts;

    RcppParallel::RMatrix<int> outputKMerCountsWrapper;

    KMerMatrixCreatorWorker(int nrow,
                            int ncol,
                            std::vector<KMerManager<kmer_dictionary_t>> &kmerCountsManagers,
                            kmer_dictionary_t<std::vector<int>, int> &hashIndexer,
                            Rcpp::StringVector &uniqueKMerStrings) :
            outputKMerCounts(Rcpp::IntegerMatrix(nrow, ncol)),
            outputKMerCountsWrapper(outputKMerCounts),
            kmerCountsManagers(kmerCountsManagers),
            hashIndexer(hashIndexer) {
        Rcpp::colnames(this->outputKMerCounts) = Rcpp::wrap(uniqueKMerStrings);
    }

    inline void operator()(std::size_t beginRow, std::size_t endRow) override {
        for (int r = beginRow; r < endRow; ++r) {
            for (const auto &kmerHashPair: kmerCountsManagers[r].getDictionary()) {
                int c = hashIndexer[kmerHashPair.first];
                outputKMerCountsWrapper(r, c) = kmerHashPair.second.cnt;
            }
        }
    }

private:
    std::vector<KMerManager<kmer_dictionary_t>> &kmerCountsManagers;
    kmer_dictionary_t<std::vector<int>, int> &hashIndexer;

};


#endif //SEQR_RESULT_CREATOR_H
