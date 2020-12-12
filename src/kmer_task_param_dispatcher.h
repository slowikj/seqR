#ifndef SEQR_KMER_TASK_PARAM_DISPATCHER_H
#define SEQR_KMER_TASK_PARAM_DISPATCHER_H

#include <Rcpp.h>
#include "user_params.h"
#include "count_kmers_specific/count_kmers_integer_matrix.h"
#include "count_kmers_specific/count_kmers_numeric_matrix.h"
#include "count_kmers_specific/count_kmers_string_matrix.h"
#include "count_kmers_specific/count_kmers_string_vector.h"
#include "dictionary/dictionaries.h"

template<class sequences_t,
        class alphabet_t,
        class algorithm_params_t>
inline Rcpp::List countKMers(
        sequences_t &sequences,
        alphabet_t &alphabet,
        const UserParams &userParams,
        algorithm_params_t &algorithmParams
);

template<class sequences_t,
        class alphabet_t,
        class algorithm_params_t>
inline Rcpp::List countKMersDictionaryDispatch(
        sequences_t &sequences,
        alphabet_t &alphabet,
        const UserParams &userParams,
        algorithm_params_t &algorithmParams);

template<class sequences_t,
        class alphabet_t,
        class algorithm_params_t,
        template<typename key, typename value, typename...> class kmer_dictionary_t>
inline Rcpp::List countKMersParallelModeDispatch(
        sequences_t &sequences,
        alphabet_t &alphabet,
        const UserParams &userParams,
        algorithm_params_t &algorithmParams);

// ------------------ IMPLEMENTATION ------------------

template<class sequences_t,
        class alphabet_t,
        class algorithm_params_t>
inline Rcpp::List countKMers(
        sequences_t &sequences,
        alphabet_t &alphabet,
        const UserParams &userParams,
        algorithm_params_t &algorithmParams
) {
    return countKMersDictionaryDispatch<sequences_t, alphabet_t, algorithm_params_t>(
            sequences, alphabet, userParams, algorithmParams);
}

template<class sequences_t, class alphabet_t, class algorithm_params_t>
inline Rcpp::List countKMersDictionaryDispatch(
        sequences_t &sequences,
        alphabet_t &alphabet,
        const UserParams &userParams,
        algorithm_params_t &algorithmParams) {
    if (userParams.kMerDictionaryName == dictionary::names::STL_UNORDERED_MAP) {
        return countKMersParallelModeDispatch<sequences_t, alphabet_t, algorithm_params_t, dictionary::StlUnorderedMapWrapper>(
                sequences, alphabet, userParams, algorithmParams);
    } else if (userParams.kMerDictionaryName == dictionary::names::LINEAR_LIST) {
        return countKMersParallelModeDispatch<sequences_t, alphabet_t, algorithm_params_t, dictionary::LinearListDictionary>(
                sequences, alphabet, userParams, algorithmParams);
    } else if (userParams.kMerDictionaryName == dictionary::names::STL_ORDERED_MAP) {
        return countKMersParallelModeDispatch<sequences_t, alphabet_t, algorithm_params_t, dictionary::StlOrderedMapWrapper>(
                sequences, alphabet, userParams, algorithmParams);
    } else if (userParams.kMerDictionaryName == dictionary::names::EMILIB_HASH_MAP) {
        return countKMersParallelModeDispatch<sequences_t, alphabet_t, algorithm_params_t, dictionary::EmilibHashMapWrapper>(
                sequences, alphabet, userParams, algorithmParams);
    } else if (userParams.kMerDictionaryName == dictionary::names::MARTINUS_ROBIN_HOOD_DICTIONARY) {
        return countKMersParallelModeDispatch<sequences_t, alphabet_t, algorithm_params_t, dictionary::MartinusRobinHoodDictionary>(
                sequences, alphabet, userParams, algorithmParams);
    } else {
        std::string errorMessage = "unsupported k-mer dictionary name: " + userParams.kMerDictionaryName;
        throw Rcpp::exception(errorMessage.c_str());
    }
}

template<class sequences_t, class alphabet_t, class algorithm_params_t,
        template<typename key, typename value, typename...> class kmer_dictionary_t>
inline Rcpp::List countKMersParallelModeDispatch(
        sequences_t &sequences,
        alphabet_t &alphabet,
        const UserParams &userParams,
        algorithm_params_t &algorithmParams) {
    if (userParams.parallelMode) {
        return parallelCountKMersSpecific<decltype(algorithmParams), kmer_dictionary_t>(
                sequences, alphabet, userParams, algorithmParams);
    } else {
        return sequentialCountKMersSpecific<decltype(algorithmParams), kmer_dictionary_t>(
                sequences, alphabet, userParams, algorithmParams);
    }
}

#endif //SEQR_KMER_TASK_PARAM_DISPATCHER_H
