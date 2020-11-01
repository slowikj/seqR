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
        src/hash/complex_hasher.h src/hash/custom_hashers.h
        src/hash/rolling_window.h src/dictionary/unordered_map_wrapper.h src/alphabet_encoder.h
        src/proxy_hashing.cpp
        src/input_to_string_item_converter.h
        src/kmer_manager.h
        src/kmer_counting_common_algorithm.h
        src/kmer_counter.h
        src/gapped_kmer_counter.h
        src/kmer_strings_creator.h
        src/sequence_getter.h
        src/proxy_alphabet_encoder.cpp
        src/proxy_count_gapped_kmers.cpp
        src/proxy_count_kmers.cpp src/rcpp_to_cpp_converters.h src/safe_sequences_wrapper.h src/dictionary/iterator.h src/dictionary/linear_list_dictionary.h src/dictionary/supported_dict_names.h src/kmer_task_config.h src/kmer_task_solver.h src/kmer_counting_result.h src/kmer_position_info.h src/kmer_task_solver_type_specific.h)