#pragma once

#include <Rcpp.h>
#include "user_params.h"
#include "count_kmers_specific/count_kmers_specific.h"
#include "dictionary/dictionaries.h"

template <class sequences_t,
		  class alphabet_t,
		  class algorithm_params_t>
inline Rcpp::List countKMers(
	sequences_t &sequences,
	alphabet_t &alphabet,
	const UserParams &userParams,
	algorithm_params_t &algorithmParams);

template <class sequences_t,
		  class alphabet_t,
		  class algorithm_params_t>
inline Rcpp::List countKMersDictionaryDispatch(
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
	algorithm_params_t &algorithmParams)
{
	return countKMersDictionaryDispatch<sequences_t, alphabet_t, algorithm_params_t>(
		sequences, alphabet, userParams, algorithmParams);
}

template <class sequences_t,
		  class alphabet_t,
		  class algorithm_params_t>
inline Rcpp::List countKMersDictionaryDispatch(
	sequences_t &sequences,
	alphabet_t &alphabet,
	const UserParams &userParams,
	algorithm_params_t &algorithmParams)
{
	if (userParams.kMerDictionaryName == dictionary::names::STL_UNORDERED_MAP)
	{
		return countKMersKMerManagerDispatch<sequences_t, alphabet_t, algorithm_params_t,
											 dictionary::StlUnorderedMapWrapper>(
			sequences, alphabet, userParams, algorithmParams);
	}
	else if (userParams.kMerDictionaryName == dictionary::names::LINEAR_LIST)
	{
		return countKMersKMerManagerDispatch<sequences_t, alphabet_t, algorithm_params_t,
											 dictionary::LinearListDictionary>(
			sequences, alphabet, userParams, algorithmParams);
	}
	else if (userParams.kMerDictionaryName == dictionary::names::STL_ORDERED_MAP)
	{
		return countKMersKMerManagerDispatch<sequences_t, alphabet_t, algorithm_params_t,
											 dictionary::StlOrderedMapWrapper>(
			sequences, alphabet, userParams, algorithmParams);
	}
	else if (userParams.kMerDictionaryName == dictionary::names::EMILIB_HASH_MAP)
	{
		return countKMersKMerManagerDispatch<sequences_t, alphabet_t, algorithm_params_t,
											 dictionary::EmilibHashMapWrapper>(
			sequences, alphabet, userParams, algorithmParams);
	}
	else if (userParams.kMerDictionaryName == dictionary::names::MARTINUS_ROBIN_HOOD_DICTIONARY)
	{
		return countKMersKMerManagerDispatch<sequences_t, alphabet_t, algorithm_params_t,
											 dictionary::MartinusRobinHoodDictionary>(
			sequences, alphabet, userParams, algorithmParams);
	}
	else if (userParams.kMerDictionaryName == dictionary::names::FLAT_HASHMAP)
	{
		return countKMersKMerManagerDispatch<sequences_t, alphabet_t, algorithm_params_t,
											 dictionary::FlatHashMapWrapper>(
			sequences, alphabet, userParams, algorithmParams);
	}
	else
	{
		std::string errorMessage = "unsupported k-mer dictionary name: " + userParams.kMerDictionaryName;
		throw Rcpp::exception(errorMessage.c_str());
	}
}

template <class sequences_t, class alphabet_t, class algorithm_params_t,
		  template <typename key, typename value, typename...> class result_dictionary_t>
inline Rcpp::List countKMersKMerManagerDispatch(
	sequences_t &sequences,
	alphabet_t &alphabet,
	const UserParams &userParams,
	algorithm_params_t &algorithmParams)
{
	if (userParams.withKMerCounts)
	{
		return countKMersParallelModeDispatch<sequences_t, alphabet_t, algorithm_params_t,
			CountingKMerManager<result_dictionary_t>,
			result_dictionary_t>(
			sequences, alphabet, userParams, algorithmParams);
	}
	else
	{
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
	algorithm_params_t &algorithmParams)
{
	if (userParams.parallelMode)
	{
		return parallelCountKMersSpecific<decltype(algorithmParams), kmer_manager_t, result_dictionary_t>(
			sequences, alphabet, userParams, algorithmParams);
	}
	else
	{
		return sequentialCountKMersSpecific<decltype(algorithmParams), kmer_manager_t, result_dictionary_t>(
			sequences, alphabet, userParams, algorithmParams);
	}
}
