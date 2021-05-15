cmake_minimum_required(VERSION 3.10)

project(seqR)
set(CMAKE_CXX_STANDARD 17)

set(R_LIBRARY "~/R/x86_64-pc-linux-gnu-library/4.0/")
set(R_INCLUDE "/usr/share/R/include/") 
include_directories(${R_LIBRARY}/Rcpp/include)
include_directories(${R_LIBRARY}/RcppParallel/include)

include_directories(${R_INCLUDE})

add_library(seqR
        src/utils.h
        src/hash/single_hasher.h
        src/hash/polynomial_single_hasher.h
        src/hash/complex_hasher.h
        src/hash/rolling_window.h
        src/kmer_manager.h
        src/contiguous_kmer_counter.h
        src/gapped_kmer_counter.h
        src/count_kmers_specific/count_kmers_string_list.h
        src/count_kmers_specific/encoded_sequence/encoded_sequence_proxy.h
        src/count_kmers_specific/encoded_sequence/raw_encoded_sequences_list.h
        src/kmer_strings_creator.h
        src/dictionary/helper/iterator_wrapper.h
        src/dictionary/supported_dict_names.h
        src/kmer_task_config.h
        src/kmer_task_solver.h
        src/kmer_counting_result.h
        src/count_kmers_specific/count_kmers_string_vector.h
        src/hash/custom_vector_hasher.h
        src/common_config.h
        src/user_params.h
        src/kmer_task_param_dispatcher.h
        src/hash/globals.h
        src/hash/prefix_sequence_polynomial_hasher.h
        inst/thirdparty/dictionaries/robin_hood_martinus.h
        src/dictionary/martinus_robin_hood_dictionary.h
        src/dictionary/dictionaries.h
        src/proxy_count_gapped_kmers.cpp
        src/proxy_count_contiguous_kmers.cpp
        src/proxy_merge_kmer_results.cpp)