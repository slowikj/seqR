#pragma once

#include "count_kmers_specific/count_kmers_specific.h"
#include "dictionary/dictionaries.h"
#include "kmer_manager.h"
#include "user_params.h"

template <class sequences_t,
          class alphabet_t,
          class algorithm_params_t>
inline Rcpp::List countKMers(
    sequences_t &sequences,
    alphabet_t &alphabet,
    const UserParams &userParams,
    algorithm_params_t &algorithmParams);

template <class sequences_t, class alphabet_t, class algorithm_params_t,
          template <typename key, typename value, typename...> class result_dictionary_t>
inline Rcpp::List countKMersKMerManagerDispatch(
    sequences_t &sequences,
    alphabet_t &alphabet,
    const UserParams &userParams,
    algorithm_params_t &algorithmParams);

template <class sequences_t,
          class alphabet_t,
          class algorithm_params_t,
          template <typename key, typename value, typename...> class result_dictionary_t>
inline Rcpp::List countKMersParallelModeDispatch(
    sequences_t &sequences,
    alphabet_t &alphabet,
    const UserParams &userParams,
    algorithm_params_t &algorithmParams);

// ------------------ IMPLEMENTATION ------------------

template <class sequences_t,
          class alphabet_t,
          class algorithm_params_t>
inline Rcpp::List countKMers(
    sequences_t &sequences,
    alphabet_t &alphabet,
    const UserParams &userParams,
    algorithm_params_t &algorithmParams) {
  return countKMersKMerManagerDispatch<sequences_t, alphabet_t, algorithm_params_t,
                                       dictionary::MartinusRobinHoodDictionary>(
      sequences, alphabet, userParams, algorithmParams);
}

template <class sequences_t, class alphabet_t, class algorithm_params_t,
          template <typename key, typename value, typename...> class result_dictionary_t>
inline Rcpp::List countKMersKMerManagerDispatch(
    sequences_t &sequences,
    alphabet_t &alphabet,
    const UserParams &userParams,
    algorithm_params_t &algorithmParams) {
  if (userParams.withKMerCounts) {
    return countKMersParallelModeDispatch<sequences_t, alphabet_t, algorithm_params_t,
                                          CountingKMerManager<result_dictionary_t>,
                                          result_dictionary_t>(
        sequences, alphabet, userParams, algorithmParams);
  } else {
    return countKMersParallelModeDispatch<sequences_t, alphabet_t, algorithm_params_t,
                                          PresenceKMerManager<result_dictionary_t>,
                                          result_dictionary_t>(
        sequences, alphabet, userParams, algorithmParams);
  }
}

template <class sequences_t, class alphabet_t, class algorithm_params_t,
          class kmer_manager_t,
          template <typename key, typename value, typename...> class result_dictionary_t>
inline Rcpp::List countKMersParallelModeDispatch(
    sequences_t &sequences,
    alphabet_t &alphabet,
    const UserParams &userParams,
    algorithm_params_t &algorithmParams) {
  if (userParams.parallelMode) {
    return parallelCountKMersSpecific<algorithm_params_t, kmer_manager_t, result_dictionary_t>(
        sequences, alphabet, userParams, algorithmParams);
  } else {
    return sequentialCountKMersSpecific<algorithm_params_t, kmer_manager_t, result_dictionary_t>(
        sequences, alphabet, userParams, algorithmParams);
  }
}
