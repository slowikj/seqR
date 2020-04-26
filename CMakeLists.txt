cmake_minimum_required(VERSION 3.15)

project(seqR)
set(CMAKE_CXX_STANDARD 17)

set(R_LIBRARY "~/R/x86_64-pc-linux-gnu-library/3.6/")
set(R_INCLUDE "/usr/share/R/include/")
include_directories(${R_LIBRARY}/Rcpp/include)
include_directories(${R_LIBRARY}/RcppParallel/include)
include_directories(${R_INCLUDE})

add_library(seqR
        src/app_conf.h
        src/utils.h
        src/utils.cpp
        src/hash/single_hasher.h
        src/hash/polynomial_single_hasher.h
        src/hash/complex_hasher.h src/hash/custom_hashers.h
        src/hash/rolling_window.h src/dictionary.h src/alphabet_encoder.h
        src/hashing_proxy.cpp
        src/input_to_string_item_converter.h
        src/input_to_string_item_converter.cpp
        src/kmer_counts_manager.h
        src/kmer_counting_common.h
        src/kmer_counter.h
        src/gapped_kmer_counter.h
        src/kmer_hash_indexer.h
        src/kmer_hash_indexer.cpp
        src/kmer_strings_creator.h
        src/sequence_getter.h
        src/tidysq_encoded_sequence.h
        src/tidysq_encoded_sequence.cpp
        src/gapped_kmer_counter.cpp
        src/alphabet_encoder_proxy.cpp
        src/count_gapped_kmers_proxy.cpp
        src/count_kmers_proxy.cpp src/rcpp_to_cpp_converters.h src/safe_sequences_wrapper.h)