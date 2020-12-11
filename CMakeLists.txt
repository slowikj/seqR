cmake_minimum_required(VERSION 3.15)

project(seqR)
set(CMAKE_CXX_STANDARD 17)

set(R_LIBRARY "~/R/x86_64-pc-linux-gnu-library/3.6/")
set(R_INCLUDE "/usr/share/R/include/")
include_directories(${R_LIBRARY}/Rcpp/include)
include_directories(${R_LIBRARY}/RcppParallel/include)
include_directories(${R_INCLUDE})

add_library(seqR
        src/utils.h
        src/hash/single_hasher.h
        src/hash/polynomial_single_hasher.h
        src/hash/complex_hasher.h
        src/hash/rolling_window.h src/dictionary/unordered_map_wrapper.h src/alphabet_encoder/default_alphabet_encoder.h
        src/kmer_manager.h
        src/contiguous_kmer_counter.h
        src/gapped_kmer_counter.h
        src/kmer_strings_creator.h
        src/count_kmers_specific/rmatrix_row_getter.h
        src/proxy_count_gapped_kmers.cpp
        src/proxy_count_contiguous_kmers.cpp src/count_kmers_specific/safe_sequences_wrapper.h src/dictionary/iterator_wrapper.h src/dictionary/linear_list_dictionary.h src/dictionary/supported_dict_names.h src/kmer_task_config.h src/kmer_task_solver.h src/kmer_counting_result.h src/count_kmers_specific/count_kmers_string_matrix.h src/count_kmers_specific/count_kmers_integer_matrix.h src/count_kmers_specific/count_kmers_numeric_matrix.h src/count_kmers_specific/count_kmers_string_list.h src/hash/custom_vector_hasher.h src/common_config.h src/alphabet_encoder/identity_alphabet_encoder.h src/count_kmers_specific/encoded_sequence_row.h src/user_params.h src/kmer_task_param_dispatcher.h src/hash/globals.h inst/include/fast_modulo.h inst/include/emilib/hash_map.hpp src/hash/prefix_sequence_polynomial_hasher.h src/dictionary/ordered_map_wrapper.h src/dictionary/emilib_hash_map_wrapper.h inst/thirdparty/robin_hood_martinus.h src/dictionary/martinus_robin_hood_dictionary.h src/dictionary/stl_dictionary_wrapper.h)