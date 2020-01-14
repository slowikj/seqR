
// [[Rcpp::plugins("c++1z")]]
#include<Rcpp.h>
#include<string>
#include<unordered_map>
#include<vector>
#include<functional>
#include<tuple>
#include<algorithm>

// ------------------------ SEQUENCE INTEGER ENCODING ------------------------

typedef short ITEM_ENCODING_TYPE;

template<class ALPHABET_VECTOR_TYPE, class CPP_ITEM_TYPE, class RCPP_ITEM_TYPE>
std::tuple<std::unordered_map<CPP_ITEM_TYPE, ITEM_ENCODING_TYPE>*,
           std::vector<std::string>*>
enumerate_alphabet(ALPHABET_VECTOR_TYPE& alphabet,
                   std::function<CPP_ITEM_TYPE(const RCPP_ITEM_TYPE&)>& rcpp2cpp_converter,
                   std::function<std::string(const CPP_ITEM_TYPE&)>& cpp2string_converter) {
  auto val2num_encoder = new std::unordered_map<CPP_ITEM_TYPE, ITEM_ENCODING_TYPE>();
  auto num2str_decoder = new std::vector<std::string>();
  for(const auto& alphabet_item: alphabet) {
    CPP_ITEM_TYPE cpp_alfabet_item = rcpp2cpp_converter(alphabet_item);
    if(val2num_encoder->find(cpp_alfabet_item) == val2num_encoder->end()) {
      (*num2str_decoder).push_back(cpp2string_converter(cpp_alfabet_item));
      (*val2num_encoder)[cpp_alfabet_item] = (*num2str_decoder).size();
    }
  }
  return {val2num_encoder, num2str_decoder};
}

const ITEM_ENCODING_TYPE NOT_ALLOWED_CHARACTER_CODE = -1;

template<class VECTOR_TYPE, class CPP_ITEM_TYPE, class RCPP_ITEM_TYPE>
std::vector<ITEM_ENCODING_TYPE> enumerate_sequence(VECTOR_TYPE& seq,
                                                   std::unordered_map<CPP_ITEM_TYPE, ITEM_ENCODING_TYPE>* val2num_encoder,
                                                   std::function<CPP_ITEM_TYPE(const RCPP_ITEM_TYPE&)>& rcpp2cpp_converter) {
  std::vector<ITEM_ENCODING_TYPE> res;
  res.reserve(seq.size());
  for(const auto& seq_item: seq) {
    CPP_ITEM_TYPE item_cpp = rcpp2cpp_converter(seq_item);
    ITEM_ENCODING_TYPE item_encoded = val2num_encoder -> find(item_cpp) == val2num_encoder -> end()
      ? NOT_ALLOWED_CHARACTER_CODE
      : (*val2num_encoder)[item_cpp];
    res.push_back(item_encoded);
  }
  return res;
}

template<class VECTOR_TYPE, class CPP_ITEM_TYPE, class RCPP_ITEM_TYPE>
std::tuple<std::vector<ITEM_ENCODING_TYPE>,
           std::unordered_map<CPP_ITEM_TYPE, ITEM_ENCODING_TYPE>*,
           std::vector<std::string>*>
enumerate_sequence_nonnull(VECTOR_TYPE& sequence,
                           VECTOR_TYPE& alphabet,
                           std::function<CPP_ITEM_TYPE(const RCPP_ITEM_TYPE&)> rcpp2cpp_converter,
                           std::function<std::string(const CPP_ITEM_TYPE&)> cpp2string_converter) {
  auto [val2num_encoder, num2str_decoder] = enumerate_alphabet<VECTOR_TYPE, CPP_ITEM_TYPE, RCPP_ITEM_TYPE>(
    alphabet, rcpp2cpp_converter, cpp2string_converter);
  auto res = enumerate_sequence<VECTOR_TYPE, CPP_ITEM_TYPE, RCPP_ITEM_TYPE>(sequence, val2num_encoder, rcpp2cpp_converter);
  return { res, val2num_encoder, num2str_decoder };
}

template<class NON_NULL_TYPE>
NON_NULL_TYPE get_value_or_construct_default(Rcpp::Nullable<NON_NULL_TYPE>& nullableValue) {
  return nullableValue.isNull() ? NON_NULL_TYPE() : NON_NULL_TYPE(nullableValue);
}

template<class VECTOR_TYPE, class CPP_ITEM_TYPE, class RCPP_ITEM_TYPE>
std::tuple<std::vector<ITEM_ENCODING_TYPE>,
           std::unordered_map<CPP_ITEM_TYPE, ITEM_ENCODING_TYPE>*,
           std::vector<std::string>*>
enumerate_sequence(Rcpp::Nullable<VECTOR_TYPE>& sequence,
                   Rcpp::Nullable<VECTOR_TYPE>& alphabet,
                   std::function<CPP_ITEM_TYPE(const RCPP_ITEM_TYPE&)> rcpp2cpp_converter,
                   std::function<std::string(const CPP_ITEM_TYPE&)>& cpp2string_converter) {
  VECTOR_TYPE nonNullSequence = get_value_or_construct_default(sequence);
  VECTOR_TYPE nonNullAlphabet = get_value_or_construct_default(alphabet);
  return enumerate_sequence_nonnull<VECTOR_TYPE, CPP_ITEM_TYPE, RCPP_ITEM_TYPE>(
      nonNullSequence,
      nonNullAlphabet,
      rcpp2cpp_converter,
      cpp2string_converter
  );
}

// ------------------------ SEQUENCE INTEGER ENCODING - R WRAPPERS ------------------------

template<class VECTOR_TYPE, class CPP_ITEM_TYPE, class RCPP_ITEM_TYPE>
Rcpp::IntegerVector enumerate_sequence_and_wrap(Rcpp::Nullable<VECTOR_TYPE> sequence,
                                                Rcpp::Nullable<VECTOR_TYPE> alphabet,
                                                std::function<CPP_ITEM_TYPE(const RCPP_ITEM_TYPE&)> rcpp2cpp_converter,
                                                std::function<std::string(const CPP_ITEM_TYPE&)> cpp2string_converter) {
  auto [res, val2num_encoder, num2str_decoder] = enumerate_sequence<VECTOR_TYPE, CPP_ITEM_TYPE, RCPP_ITEM_TYPE>(
    sequence,
    alphabet,
    rcpp2cpp_converter,
    cpp2string_converter
  );
  delete val2num_encoder;
  delete num2str_decoder;
  return Rcpp::wrap(res);
}

//' @export
// [[Rcpp::export]]
Rcpp::IntegerVector enumerate_string_sequence(Rcpp::Nullable<Rcpp::StringVector> sequence,
                                              Rcpp::Nullable<Rcpp::StringVector> alphabet) {
  return enumerate_sequence_and_wrap<Rcpp::StringVector, std::string, Rcpp::String::StringProxy>(
      sequence,
      alphabet,
      [](const Rcpp::String::StringProxy& s) -> std::string {return Rcpp::as<std::string>(s);},
      [](const std::string& s) -> std::string { return s; });
}

//' @export
// [[Rcpp::export]]
Rcpp::IntegerVector enumerate_integer_sequence(Rcpp::Nullable<Rcpp::IntegerVector> sequence,
                                               Rcpp::Nullable<Rcpp::IntegerVector> alphabet) {
  return enumerate_sequence_and_wrap<Rcpp::IntegerVector, int, int>(
      sequence,
      alphabet,
      [](const int& x) -> int { return x; },
      [](const int& x) -> std::string { return std::to_string(x); }
  );
}

//' @export
// [[Rcpp::export]]
Rcpp::IntegerVector enumerate_numeric_sequence(Rcpp::Nullable<Rcpp::NumericVector> sequence,
                                               Rcpp::Nullable<Rcpp::NumericVector> alphabet) {
  return enumerate_sequence_and_wrap<Rcpp::NumericVector, double, double>(
    sequence,
    alphabet,
    [](const double& x) -> double { return x; },
    [](const double& x) -> std::string { return std::to_string(x); }
  );
}

// ------------------------ COMPUTATION OF ALLOWED SEQUENCE RANGES ------------------------

std::vector<int> get_not_allowed_sequence_positions(const std::vector<ITEM_ENCODING_TYPE>& encoded_sequence) {
  std::vector<int> res;
  res.push_back(-1); // left sentinel
  for(size_t seq_i = 0; seq_i < encoded_sequence.size(); ++seq_i) {
    if(encoded_sequence[seq_i] == NOT_ALLOWED_CHARACTER_CODE) {
      res.push_back(seq_i);
    }
  }
  res.push_back(encoded_sequence.size()); // right sentinel
  return res;
}

// ------------------------ COMPUTATION OF ALLOWED SEQUENCE RANGES - R WRAPPER ------------------------

//' @export
// [[Rcpp::export]]
Rcpp::IntegerVector get_not_allowed_sequence_positions(Rcpp::IntegerVector encoded_sequence) {
  std::vector<ITEM_ENCODING_TYPE> cpp_encoded_sequence = Rcpp::as<std::vector<ITEM_ENCODING_TYPE>>(encoded_sequence);
  auto not_allowed_sequence_positions = get_not_allowed_sequence_positions(cpp_encoded_sequence);
  return static_cast<Rcpp::IntegerVector>(Rcpp::wrap(not_allowed_sequence_positions)) + 1;
}

// ------------------------ K-MERS COMPUTATION (COMPACT SUBWORDS - ROLLING WINDOW HASHES) --------------

// ------------------------ K-MERS COMPUTATION FOR ONE SEQUENCE ----------------------------------------

struct KMerHashInfo {
  int seq_start_pos;
  int cnt;
  
  KMerHashInfo() {}
  
  KMerHashInfo(int seq_start_pos, int cnt): seq_start_pos(seq_start_pos), cnt(cnt) {}
  
  KMerHashInfo(const KMerHashInfo& second) {
    this -> seq_start_pos = second.seq_start_pos;
    this -> cnt = second.cnt;
  }
  
  ~KMerHashInfo() {}
};

int compute_hash(int previous_hash, int current_item, int P, int M) {
  return (static_cast<long long>(previous_hash) * P + current_item) % M;
}

void add_hash(std::unordered_map<int, KMerHashInfo>& kmer_counter,
              int hash,
              int begin,
              bool positional_kmer,
              int P,
              int M) {
  if(positional_kmer) {
    hash = compute_hash(hash, (begin + 1), P, M);
  }
  if(kmer_counter.find(hash) == kmer_counter.end()) {
    kmer_counter[hash] = KMerHashInfo(begin, 0);
  }
  ++kmer_counter[hash].cnt;
}

void update_kmers_for_subsequence(std::unordered_map<int, KMerHashInfo>& kmer_counter,
                                  int k,
                                  std::vector<ITEM_ENCODING_TYPE>& encoded_sequence,
                                  int sequence_begin,
                                  int sequence_end,
                                  bool positional_kmer,
                                  int P,
                                  int P_K_1,
                                  int M) {
  int hash = std::accumulate(encoded_sequence.begin() + sequence_begin,
                             encoded_sequence.begin() + sequence_begin + k,
                             0,
                             [&P, &M](int prev_hash, int current_item) -> int { return compute_hash(prev_hash, current_item, P, M); });
  add_hash(kmer_counter, hash, sequence_begin, positional_kmer, P, M);
  for(int current_sequence_begin = sequence_begin + 1; current_sequence_begin <= sequence_end; ++current_sequence_begin) {
    hash = (hash - encoded_sequence[current_sequence_begin - 1] * P_K_1 + M) % M; // remove a previous item
    hash = compute_hash(hash, encoded_sequence[current_sequence_begin + k - 1], P, M); // add a current item
    add_hash(kmer_counter, hash, current_sequence_begin, positional_kmer, P, M);
  }
}   

int compute_power_fast(int base, int power, int modulo) {
  long long res = 1;
  long long current_base_power = base;
  while(power > 0) {
    if(power & 1) {
      res = (res * current_base_power) % modulo;
    }
    power >>= 1;
    current_base_power = (current_base_power * current_base_power) % modulo;
  }
  return static_cast<int>(res);
}

std::unordered_map<int, KMerHashInfo> count_kmers(std::vector<ITEM_ENCODING_TYPE>& encoded_sequence,
                                                  int k,
                                                  bool positional_kmer,
                                                  int P,
                                                  int P_K_1,
                                                  int M) {
  std::unordered_map<int, KMerHashInfo> res;
  auto not_allowed_positions = get_not_allowed_sequence_positions(encoded_sequence);
  
  for(int not_allowed_position_ind = 0; not_allowed_position_ind < not_allowed_positions.size() - 1; ++not_allowed_position_ind) {
    if(not_allowed_positions[not_allowed_position_ind + 1] - not_allowed_positions[not_allowed_position_ind] > k) {
      int sequence_begin = not_allowed_positions[not_allowed_position_ind] + 1;
      int sequence_end = not_allowed_positions[not_allowed_position_ind + 1] - k;
      update_kmers_for_subsequence(res, k, encoded_sequence, sequence_begin, sequence_end, positional_kmer, P, P_K_1, M);
    }
  }
  return res;
}

//' @export
// [[Rcpp::export]]
Rcpp::DataFrame count_kmers_hashed(Rcpp::IntegerVector encoded_sequence,
                                   int k,
                                   bool positional_kmer,
                                   int P,
                                   int P_K_1,
                                   int M) {
  auto cpp_encoded_sequence = Rcpp::as<std::vector<ITEM_ENCODING_TYPE>>(encoded_sequence);
  auto kmers_counter_hash_map = count_kmers(cpp_encoded_sequence, k, positional_kmer, P, P_K_1, M);
  
  std::map<int, int> kmers_counter_map;
  for(const auto& kmer_info_entry: kmers_counter_hash_map) {
    kmers_counter_map[kmer_info_entry.second.seq_start_pos] = kmer_info_entry.second.cnt;
  }
  
  std::vector<short> positions;
  std::vector<int> counts;
  for(const auto& map_entry: kmers_counter_map) {
    positions.push_back(map_entry.first + 1);
    counts.push_back(map_entry.second);
  } 
  
  return Rcpp::DataFrame::create(
    Rcpp::Named("position") = Rcpp::wrap(positions),
    Rcpp::Named("cnt") = Rcpp::wrap(counts)
  );
}

// ------------------------ K-MERS (BOTH GAPPED AND NORMAL) DECODING ----------------------------------------

typedef std::function<std::string(std::string&, int)> KMER_STRING_DECORATOR; // (current_string, kmer_begin_pos) -> decorated_string

typedef std::function<int(int,int)> NEXT_KMER_POSITION_GENERATOR; // (current_pos, current_position) -> next position

int get_total_size_of_kmer(std::vector<ITEM_ENCODING_TYPE>& encoded_sequence,
                           int begin_position,
                           int k,
                           std::vector<std::string>& num2str_decoder,
                           NEXT_KMER_POSITION_GENERATOR& next_sequence_position) {
  int res = 0;
  int current_position = begin_position;
  for(int i = 0; i < k; ++i) {
    int item_size = num2str_decoder[encoded_sequence[current_position]].size();
    res += i == 0 ? item_size : item_size + 1;
    current_position = next_sequence_position(current_position, i);
  }
  return res;
}

int get_total_size_of_kmer(std::vector<ITEM_ENCODING_TYPE>& encoded_kmer,
                           std::vector<std::string>& num2str_decoder) {
  return std::accumulate(encoded_kmer.begin() + 1,
                         encoded_kmer.end(),
                         num2str_decoder[encoded_kmer[0]].size(),
                         [&num2str_decoder](int prev_sum, int code) -> int { return prev_sum + num2str_decoder[code].size() + 1; });
}

const std::string KMER_ITEM_SEPARATOR = ".";

const std::string KMER_POSITION_SEPARATOR = "_";

const int NO_POSITION = -1;

std::string decode_kmer(std::vector<ITEM_ENCODING_TYPE>& encoded_sequence,
                        int begin_position,
                        int k,
                        std::vector<std::string>& decoder,
                        NEXT_KMER_POSITION_GENERATOR& next_sequence_position,
                        KMER_STRING_DECORATOR& kmer_string_decorator) {
  int kmer_size = get_total_size_of_kmer(encoded_sequence, begin_position, k, decoder, next_sequence_position);
  std::string res;
  res.reserve(kmer_size);
  int current_position = begin_position;
  for(int i = 0; i < k; ++i) {
    res += decoder[encoded_sequence[current_position]] + (i == k - 1 ? "" : KMER_ITEM_SEPARATOR);
    current_position = next_sequence_position(current_position, i);
  }
  return kmer_string_decorator(res, begin_position);
}

std::string decode_kmer(std::vector<ITEM_ENCODING_TYPE>& encoded_kmer,
                        std::vector<std::string>& num2str_decoder) {
  std::string res;
  res.reserve(get_total_size_of_kmer(encoded_kmer, num2str_decoder));
  res += num2str_decoder[encoded_kmer[0]];
  for(int ek_i = 1; ek_i < encoded_kmer.size(); ++ek_i) {
    res += KMER_ITEM_SEPARATOR + num2str_decoder[encoded_kmer[ek_i]];
  }
  return res;
}

std::string decorated_with_position(std::string& kmer, int position, bool positional_kmer) {
  return positional_kmer ?
    std::to_string(position + 1) + KMER_POSITION_SEPARATOR + kmer :
    kmer;
}

int generate_next_kmer_position(int cur_pos, int d_i, Rcpp::IntegerVector& d) {
  return cur_pos + d[d_i] + 1;
}

//' @export
// [[Rcpp::export]]
Rcpp::StringVector decode_kmer(Rcpp::IntegerVector encoded_sequence,
                               Rcpp::IntegerVector d,
                               int begin_position,
                               Rcpp::StringVector decoder,
                               bool positional_kmer) {
  std::vector<std::string> num2str_decoder = Rcpp::as<std::vector<std::string>>(decoder);
  std::vector<ITEM_ENCODING_TYPE> encoded_sequence_cpp;
  for(const auto& elem: encoded_sequence) {
    encoded_sequence_cpp.push_back(elem - 1);
  }
  NEXT_KMER_POSITION_GENERATOR next_pos = [&d](int cur_pos, int i) -> int {
    return generate_next_kmer_position(cur_pos, i, d);  
  };
  KMER_STRING_DECORATOR decorator = [&positional_kmer](std::string& cur_str, int begin_pos) -> std::string {
   return decorated_with_position(cur_str, begin_pos, positional_kmer);
  };
  
  return Rcpp::wrap(  
    decode_kmer(encoded_sequence_cpp, begin_position - 1, d.size() + 1, num2str_decoder, next_pos, decorator)
  );
}

// ------------------------ ALL K-MERS GENERATION ----------------------------------------

void add_encoded_kmer(std::vector<ITEM_ENCODING_TYPE>& encoded_kmer,
                      std::vector<std::string>& num2str_decoder,
                      int hash,
                      std::unordered_map<int, KMerHashInfo>& found_kmer_hashes,
                      std::vector<std::string>& used_kmers, // with count > 0
                      std::vector<std::string>& unused_kmers) { // with count == 0
  std::string str_kmer = decode_kmer(encoded_kmer, num2str_decoder);
  if(found_kmer_hashes.find(hash) == found_kmer_hashes.end()) {
    unused_kmers.push_back(std::move(str_kmer));
  } else {
    used_kmers.push_back(std::move(str_kmer));
  }
}

void generate_all_kmers(size_t k,
                        std::vector<ITEM_ENCODING_TYPE>& current_kmer,
                        int current_hash,
                        int P,
                        int M,
                        std::vector<std::string>& used_kmers,
                        std::vector<std::string>& unused_kmers,
                        std::unordered_map<int, KMerHashInfo>& found_kmer_hashes,
                        std::vector<std::string>& num2str_decoder) {
  if(current_kmer.size() == k) {
    add_encoded_kmer(current_kmer, num2str_decoder, current_hash, found_kmer_hashes, used_kmers, unused_kmers);
  } else {
    for(int code = 0; code < num2str_decoder.size(); ++code) {
      int next_hash = compute_hash(current_hash, code, P, M);
      current_kmer.push_back(code);
      generate_all_kmers(k, current_kmer, next_hash, P, M, used_kmers, unused_kmers, found_kmer_hashes, num2str_decoder);
      current_kmer.pop_back();
    }
  }
}

//' @export
// [[Rcpp::export]]
std::vector<std::string> generate_all_kmers(int k,
                                      int P,
                                      int M,
                                      Rcpp::StringVector& decoder) {
  std::vector<std::string> used_kmers, unused_kmers;
  std::unordered_map<int, KMerHashInfo> found_kmer_hashes;
  std::vector<std::string> num2str_decoder = Rcpp::as<std::vector<std::string>>(decoder);
  std::vector<ITEM_ENCODING_TYPE> current_kmer;
  generate_all_kmers(k, current_kmer, 0, P, M, used_kmers, unused_kmers, found_kmer_hashes, num2str_decoder);
  return unused_kmers;
}
