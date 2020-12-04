#ifndef SEQR_USER_PARAMS_H
#define SEQR_USER_PARAMS_H

#include "Rcpp.h"
#include <string>
#include <vector>

struct UserParams {

    int k;
    std::vector<int> gaps;
    bool positional;
    bool withKMerCounts;
    const std::string kMerDictionaryName;
    int batchSize;
    int hashDim;
    bool verbose;
    bool parallelMode;

    static UserParams createForContiguous(Rcpp::Environment &params) {
        UserParams res(params);
        res.gaps.resize(res.k - 1, 0);
        return std::move(res);
    }

    static UserParams createForGapped(Rcpp::Environment &params) {
        UserParams res(params);
        res.gaps = std::move(Rcpp::as<std::vector<int>>(params.get("gaps")));
        return std::move(res);
    }

private:
    UserParams(Rcpp::Environment &params) :
            k(Rcpp::as<int>(params.get("k"))),
            positional(Rcpp::as<bool>(params.get("positional"))),
            withKMerCounts(Rcpp::as<bool>(params.get("with_kmer_counts"))),
            kMerDictionaryName(Rcpp::as<std::string>(params.get("kmer_dictionary_name"))),
            batchSize(Rcpp::as<bool>(params.get("batch_size"))),
            hashDim(Rcpp::as<int>(params.get("hash_dim"))),
            verbose(Rcpp::as<bool>(params.get("verbose"))),
            parallelMode(Rcpp::as<bool>(params.get("parallel_mode"))) {
    }
};


#endif //SEQR_USER_PARAMS_H
