
// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>
// [[Rcpp::plugins("cpp11")]]
#include<Rcpp.h>
#include<string>
#include<unordered_map>
#include<vector>
#include<iostream>
#include<algorithm>
#include<functional>
#include<mutex>


// a helper typedef used for adding presentation details for k-mer
// the first parameter denotes input k-mer string
// the second parameter denotes the position of k-mer in the given sequence
typedef std::function<std::string(std::string&, int)> KMER_DECORATOR;

// a prime number used as a base for hashing function
const int HASH_CONST = 233;

// a modulo constant used for hashing function
const int MOD = 1e9 + 33;

const std::string KMER_ITEM_SEPARATOR = ".";
  
//' @name get_window_length
//' @title Get k-mer window length
//' 
//' @description Compute a k-mer window length. The window length is the total size 
//' used by the k-mer - the number of elements and the size of gaps.
//' 
//' @param d  \code{integer} vector with distances between consequent elements
//' @return \code{integer} representing the total window length
//' @export
// [[Rcpp::export]]
int get_window_length(const Rcpp::IntegerVector& d) {
  int res = 1;
  for(const int& item: d) {
    res += item + 1;
  }
  return res;
}

//' @name get_hash
//' @title Get hash of k-mer that is in a given sequence
//' 
//' @param s  \code{integer} vector representing input sequence
//' @param d  \code{integer} vector representing the gaps between consecutive elements of k-mer
//' @param begin_index  \code{integer} value representing the begin index of the k-mer
//' @param pos  \code{logical} value representing whether k-mer is positional
//' 
//' @return \code{integer} value representing the result of a hashing function
//' @export
// [[Rcpp::export]]
int get_hash(const std::vector<int>& s,
             const Rcpp::IntegerVector& d,
             int begin_index,
             bool pos) {
  long long lres = ((pos ? begin_index : 0) * HASH_CONST + s[begin_index]) % MOD;
  int current_letter_index = begin_index;
  for(int i = 0; i < d.size(); ++i) {
    current_letter_index += d[i] + 1;
    lres = ((lres * HASH_CONST) + s[current_letter_index]) % MOD;
  }
  return (int)lres;
}

//' @name get_hash_for_word
//' @title Get hash of a sequence
//' 
//' @param kmer  \code{integer} vector representing a word to be hashed
//' @return \code{integer} value representing the result of a hashing function
//' @export
// [[Rcpp::export]]
int get_hash_for_word(const std::vector<int>& kmer) {
  long long hash = 0;
  for(int i = 0; i < kmer.size(); ++i) {
    hash = ((hash * HASH_CONST) + kmer[i]) % MOD;
  }
  return (int)hash;
}

//' @name get_total_size_of_kmer
//' @title Get the total size of k-mer
//' 
//' @description Computes the number of characters of the result k-mer
//' taking into account the base alphabet. 
//' 
//' @param s  \code{integer} vector of encoded elements of a sequence (see Details)
//' @param d  \code{integer} vector which denotes the gaps in k-mer
//' @param begin_index  \code{integer} representing the begin index (in \code{s}) of the k-mer
//' @param num2str  a \code{hash map} representing the encoding between number and \code{string} representation of each alphabet's item
//' @return \code{int} denoting the total size (number of characters) of k-mer
//' @details Each element of a sequence is previously encoded to an integer in order to make hashing computation
//' more convenient
int get_total_size_of_kmer(const std::vector<int>& s,
                           const Rcpp::IntegerVector& d,
                           int begin_index,
                           std::unordered_map<int, std::string>& num2str) {
  int res = num2str[s[begin_index]].size();
  int current_index = begin_index;
  for(int d_i = 0; d_i < d.size(); ++d_i) {
    current_index += d[d_i] + 1;
    res += num2str[s[current_index]].size() + 1; // + 1 because of the dot separator
  }
  return res;
}

//' @name get_total_size_of_kmer
//' @title Get the total size (number of characters) of a k-mer
//' 
//' @description The number of characters of the result k-mer (after decoding from \code{integer} to \code{string})
//' 
//' @param kmer  \code{integer} vector representing the encoded k-mer (\link{get_total_size_of_kmer})
//' @param num2str  \code{hash map} representing the encoding between the integer and string
//' @return the number of characters in the result \code{string} that is the result of decoding each \code{integer} from \code{kmer}
int get_total_size_of_kmer(const std::vector<int>& kmer,
                           std::unordered_map<int, std::string>& num2str) {
  int total_size = 0;
  for(const int& item: kmer) {
    total_size += num2str[item].size() + 1;
  }
  return total_size - 1;
}

//' @name create_kmer
//' @title Create k-mer
//' 
//' @description Creates k-mer (of type \code{string}) from encoded (\code{integer}) vector
//' based on encoding described in \code{num2str} and kmer_decorator
//' 
//' @param s  \code{integer} vector representing an encoded sequence
//' @param d  \code{integer} vector representing the gaps in k-mer
//' @param begin_index  \code{integer} representing the start of k-mer (in \code{s})
//' @param num2str  \code{hash map} representing encoding of sequence items between \code{integer} and \code{string}
//' @param kmer_decorator  a \code{function} that can add extra characters to k-mer (for example position information)
//' 
//' @return a \code{string} representing a result k-mer (that is used for presentation)
std::string create_kmer(const std::vector<int>& s,
                        const Rcpp::IntegerVector& d,
                        int begin_index,
                        std::unordered_map<int, std::string>& num2str,
                        KMER_DECORATOR kmer_decorator) {
  int total_kmer_size = get_total_size_of_kmer(s, d, begin_index, num2str);
  std::string res;
  res.reserve(total_kmer_size);
  res += num2str[s[begin_index]];
  int current_index = begin_index;
  for(int d_i = 0; d_i < d.size(); ++d_i) {
    current_index += d[d_i] + 1;
    res += KMER_ITEM_SEPARATOR + num2str[s[current_index]];
  }
  return kmer_decorator(res, begin_index);
}

//' @name create_kmer
//' @title Create k-mer
//' 
//' @description Creates k-mer (of type \code{string}) from encoded (\code{integer}) vector
//' based on encoding described in \code{num2str} and kmer_decorator
//' 
//' @param kmer  \code{integer} vector representing an encoded sequence
//' @param num2str  \code{hash map} representing encoding of sequence items between \code{integer} and \code{string}
//' @param kmer_decorator  a \code{function} that can add extra characters to k-mer (for example position information)
//' 
//' @return a \code{string} representing a result k-mer (that is used for presentation)
std::string create_kmer(const std::vector<int>& kmer,
                        std::unordered_map<int, std::string>& num2str,
                        KMER_DECORATOR kmer_decorator) {
  int total_kmer_size = get_total_size_of_kmer(kmer, num2str);
  std::string res;
  res.reserve(total_kmer_size);
  for(const int& elem: kmer) {
    if(res.size() == 0) {
      res += num2str[elem];
    } else {
      res += KMER_ITEM_SEPARATOR + num2str[elem];
    }
  }
  return kmer_decorator(res, -1);
}

//' @name update_kmers
//' @title Update k-mers in a \code{hash map}
//' 
//' @param kmers  a \code{hash map} reference representing the found k-mers (see details)
//' @param d  \code{integer} vector representing gaps in k-mer
//' @param s  \code{integer} vector representing an encoded sequence
//' @param kmer_hash \code{integer} representing computed hash of k-mer
//' @param kmer_begin_index \code{integer} representing the begin index of k-mer in \code{s}
//' @param num2str  \code{hash map} representing encoding of sequence items between \code{integer} and \code{string}
//' @param kmer_decorator  \code{function} that can add extra characters to \code{string} k-mer (for example position information)
//' 
//' @details k-mers \code{hashmap} contains key-value pairs: key is an \code{integer} representing a hash of k-mer,
//' whereas the value represents a pair: (k-mer \code{string} for presentation, number of k-mer occurrences)
void update_kmers(std::unordered_map<int, std::pair<std::string, int>>& kmers,
                  const Rcpp::IntegerVector& d,
                  const std::vector<int>& s,
                  int kmer_hash,
                  int kmer_begin_index,
                  std::unordered_map<int, std::string>& num2str,
                  KMER_DECORATOR kmer_decorator) {
  auto map_elem = kmers.find(kmer_hash);
  if(map_elem == kmers.end()) {
    kmers[kmer_hash] = std::make_pair(create_kmer(s, d, kmer_begin_index, num2str, kmer_decorator), 0);
  }
  ++kmers[kmer_hash].second;
}

//' @name add_kmer_if_not_exists
//' @title Add k-mer to a \code{hash map} if it does not exist
//' 
//' @param kmers  \code{hash map} containing k-mers (see \link{update_kmers})
//' @param kmer  encoded \code{integer} vector representing a k-mer
//' @param num2str  \code{hash map} representing encoding between \code{integer} and \code{string}
//' @param kmer_decorator  \code{function} that can add extra characters to the \code{string} representation of k-mer
void add_kmer_if_not_exists(std::unordered_map<int, std::pair<std::string, int>>& kmers,
                            std::vector<int> kmer,
                            std::unordered_map<int, std::string>& num2str,
                            KMER_DECORATOR kmer_decorator) {
  int hash = get_hash_for_word(kmer);
  auto map_entry = kmers.find(hash);
  if(map_entry == kmers.end()) {
    kmers[hash] = std::make_pair(create_kmer(kmer, num2str, kmer_decorator), 0);
  }
}

//' @name update_kmers_with_alphabet
//' @title Update k-mers with alphabet
//' 
//' @description Generates and add k-mers (based on the given alphabet) that do not exist in the \code{hash map}.
//' 
//' @param kmers  \code{hash map} reference that stores k-mers (see \link{update_kmers})
//' @param alphabet  \code{integer} vector representing encoded alphabet
//' @param currentKmer \code{integer} vector representing the part of currently generated k-mer
//' @param k  \code{integer} representing the number of k-mer items
//' @param num2str  \code{hash map} representing the sequence encoding between \code{integer} and \code{string}
//' @param kmer_decorator  \code{function} that can add extra characters for \code{string} representation of k-mer (for presentation reasons)
void update_kmers_with_alphabet(std::unordered_map<int, std::pair<std::string, int>>& kmers,
                                const std::vector<int>& alphabet,
                                std::vector<int>& currentKmer,
                                int k,
                                std::unordered_map<int, std::string>& num2str,
                                KMER_DECORATOR kmer_decorator) {
  if(currentKmer.size() == k) {
    add_kmer_if_not_exists(kmers, currentKmer, num2str, kmer_decorator);
  } else {
    for(char letter: alphabet) {
      currentKmer.push_back(letter);
      update_kmers_with_alphabet(kmers, alphabet, currentKmer, k, num2str, kmer_decorator);
      currentKmer.pop_back();
    }
  }
}

//' @name is_kmer_allowed
//' @title Is k-mer allowed
//' 
//' @description Checks whether all elements of the given k-mer are contained in the alphabet set
//' @param s  \code{integer} vector of encoded sequence characters
//' @param d  \code{integer} vector representing gaps between elements of k-mer
//' @param begin_index  \code{integer} representing the start index of k-mer in \code{s}
//' @param is_item_allowed \code{hash map} that answers the question whether the element is in the alphabet
//' @return \code{logical} value denoting whether k-mer is valid (contains valid characters that are in the alphabet)
bool is_kmer_allowed(const std::vector<int>& s,
                     const Rcpp::IntegerVector& d,
                     int begin_index,
                     std::unordered_map<int, bool>& is_item_allowed) {
  int current_index = begin_index;
  int i = 0;
  do {
    if(!is_item_allowed[s[current_index]]) {
      return false;
    }
    current_index += d[i] + 1;
    ++i;
  } while (i <= d.length());
  return true;
}

//' @name count_kmers_helper
//' @title Count k-mers 
//' @description Counts the occurrences of k-mers (the size of k-mer should be larger than one)
//' 
//' @param s  \code{integer} vector representing encoded input sequence
//' @param d  \code{integer} vector representing gaps in k-mer
//' @param alphabet  \code{integer} vector representing encoded alphabet
//' @param num2str  \code{hash map} representing the sequence elements encoding between \code{integer} and \code{string}
//' @param kmer_decorator \code{function} that can add extra characters in order to enhance the presentation of \code{string} k-mer
//' @param pos  \code{logical} value representing whether to count positional k-mers
//' 
//' @return \code{hash map} whose key is a \code{string} presentation of k-mer and value is the number of its occurrences
std::unordered_map<std::string, int> count_kmers_helper(const std::vector<int>& s,
                                                   const Rcpp::IntegerVector& d,
                                                   const std::vector<int>& alphabet,
                                                   std::unordered_map<int, std::string> num2str,
                                                   KMER_DECORATOR kmer_decorator,
                                                   bool pos) {
  std::unordered_map<int, bool> isItemAllowed;
  for(const int& c: alphabet) {
    isItemAllowed[c] = true;
  }
  
  std::unordered_map<int, std::pair<std::string, int>> kmers; // hash -> (kmer, count)
  
  // fill map with found (allowed) k-mers
  int window_length = get_window_length(d);
  int last_window_index = s.size() - window_length;
  for(int kmer_begin_index = 0; kmer_begin_index <= last_window_index; ++kmer_begin_index) {
    if(is_kmer_allowed(s, d, kmer_begin_index, isItemAllowed)) {
      update_kmers(kmers, d, s, get_hash(s, d, kmer_begin_index, pos), kmer_begin_index, num2str, kmer_decorator);
    }
  }
  
  if(!pos) {
    // fill map with kmers that can be created with the given alphabet
    std::vector<int> currentKmer;
    currentKmer.clear();
    update_kmers_with_alphabet(kmers, alphabet, currentKmer, (int)(d.size() + 1), num2str, kmer_decorator);
  }
  
  std::unordered_map<std::string, int> res;
  for(const auto& entry: kmers) {
    res[entry.second.first] = entry.second.second;
  }
  return res;
}

//' @name fill_items_coding_maps
//' @title Prepare encoding and decoding \code{hash maps} for sequence items
//' 
//' @details Enumerates each sequence item in order to convert \code{non-integer} values to \code{integer} ones
//' \code{B} is the template Rcpp input type
//' \code{S} is the template c++ input type of one element
//' @param elems  the input elements of a sequence
//' @param val2num  the reference to \code{hash map} representing the encoding to \code{integer} value
//' @param num2str the reference to \code{hash map} representing the (reversed to \code{val2num}) encoding to \code{string} value
//' @param lowest_not_used_num  the reference to \code{integer} value denoting current counter used to encode elements
//' @param val2str_converter  \code{function} that take a sequence item and returns its string representation that is used for presentation
template <class B, class S>
void fill_items_coding_maps(B& elems,
                            std::unordered_map<S, int>& val2num,
                            std::unordered_map<int, std::string>& num2str,
                            int& lowest_not_used_num,
                            std::function<std::string(S)> val2str_converter) {
  for(int i = 0; i < elems.size(); ++i) {
    S s_item = static_cast<S>(elems[i]);
    if(val2num.find(s_item) == val2num.end()) {
      val2num[s_item] = lowest_not_used_num;
      num2str[lowest_not_used_num] = val2str_converter(s_item);
      ++lowest_not_used_num;
    }
  }
}

//' @name fill_encoded_int_vector
//' @title Encode sequence vector (replace items to numbers)
//' 
//' @details \code{SEQ_TYPE} - the type of a sequence of the input sequence (Rcpp)
//' \code{ELEM_TYPE} - the type of an item of the input sequence (c++)
//' 
//' @param str_v  the input (Rcpp) sequence
//' @param res  the result (encoded) vector
//' @param val2int  the encoder
template <class SEQ_TYPE, class ELEM_TYPE>
void fill_encoded_int_vector(SEQ_TYPE str_v,
                             std::vector<int>& res,
                             std::unordered_map<ELEM_TYPE, int>& val2int) {
  for(int i = 0; i < str_v.size(); ++i) {
    res[i] = val2int[static_cast<ELEM_TYPE>(str_v[i])];
  }
}

//' @name get_kmers
//' @title Get k-mers
//' 
//' @description Counts the occurrences of k-mers (the size of k-mer should be larger than one)
//' \code{B} - the (Rcpp) type of an input sequence
//' \code{S} - the (c++) type of an item of the sequence
//' 
//' @param s  input sequence
//' @param d  \code{integer} vector representing gaps in k-mer
//' @param alphabet  Rcpp sequence representing alphabet
//' @param val2str_converter  \code{function} representing the conversion to string representation of an item
//' @param kmer_decorator  \code{function} that can add extra characters in order to enhance the presentation of k-mer.
//' @return \code{hash map} containing string representations of k-mers with their occurrence counts 
template <class B, class S>
std::unordered_map<std::string, int> get_kmers(B& s,
                                               Rcpp::IntegerVector& d,
                                               B& alphabet,
                                               std::function<std::string(S)> val2str_converter,
                                               KMER_DECORATOR kmer_decorator,
                                               bool pos) {
  // create S -> int (and vice versa) coding maps in order to deal with numbers
  int lowest_not_used_num = 1;
  std::unordered_map<S, int> val2num;
  std::unordered_map<int, std::string> num2str;
  fill_items_coding_maps(s, val2num, num2str, lowest_not_used_num, val2str_converter);
  fill_items_coding_maps(alphabet, val2num, num2str, lowest_not_used_num, val2str_converter);
  
  std::vector<int> int_s(s.size());
  fill_encoded_int_vector<B, S>(s, int_s, val2num);
  
  std::vector<int> int_alphabet(alphabet.size());
  fill_encoded_int_vector<B, S>(alphabet, int_alphabet, val2num);
  
  return count_kmers_helper(int_s, d, int_alphabet, num2str, kmer_decorator, pos);
}

// a helper function used for determining how to decorate k-mer
KMER_DECORATOR get_kmer_decorator(bool pos) {
  return pos ?
    [](std::string& s, int p) { return std::to_string(p) + "_" + s; } :
    [](std::string& s, int p) { return s; };
}

// a helper function for extracting a single logical value
bool is_first_true(Rcpp::LogicalVector& v) {
  return static_cast<bool>(v[0]);
}

//' @name count_kmers_str
//' @title Count k-mers for string sequences (the size of k-mer should be larger than one)
//' 
//' @param s  a \code{string} vector representing an input sequence
//' @param d  an \code{integer} vector representing gaps between consecutive elements of k-mer
//' @param alphabet a \code{string} vector representing valid elements of k-mer
//' @param pos a \code{logical} value that denotes whether positional k-mers should be generated
//' @return a named vector with counts of k-mers
//' 
//' @details K-mers that contain elements from \code{alphabet} but do not exist in the input sequence are also generated.
//' 
//' @examples
//' count_kmers_str(
//' c("a", "b", "c", "d", "x", "y", "z", "z", "a", "a"),
//' d=c(0,0),
//' c("a", "b", "c", "z"),
//' pos=FALSE)
//' @export
// [[Rcpp::export]]
std::unordered_map<std::string, int> count_kmers_str(Rcpp::StringVector& s,
                                                     Rcpp::IntegerVector& d,
                                                     Rcpp::StringVector& alphabet,
                                                     Rcpp::LogicalVector& pos) {
  bool positional = is_first_true(pos);
  return get_kmers<Rcpp::StringVector, std::string>(s, d, alphabet,
                                             [](std::string s) { return s; },
                                             get_kmer_decorator(positional),
                                             positional);
}

//' @name count_kmer_num
//' @title Count k-mers for numeric sequences (the size of k-mer should be larger than one)
//' 
//' 
//' @param s  a \code{numeric} vector representing an input sequence
//' @param d  an \code{integer} vector representing gaps between consecutive elements of k-mer
//' @param alphabet a \code{numeric} vector representing valid elements of k-mer
//' @param pos a \code{logical} value that denotes whether positional k-mers should be generated
//' @return a named vector with counts of k-mers
//' 
//' @details K-mers that contain elements from \code{alphabet} but do not exist in the input sequence are also generated.
//' 
//' @examples
//' count_kmers_str(c(1,2,3,5,3,7),
//' d=c(0,0),
//' c(1, 2, 3, 4),
//' pos=FALSE)
//' @export
// [[Rcpp::export]]
std::unordered_map<std::string, int> count_kmer_num(Rcpp::NumericVector& s,
                                                    Rcpp::IntegerVector& d,
                                                    Rcpp::NumericVector& alphabet,
                                                    Rcpp::LogicalVector& pos) {
  bool positional = is_first_true(pos);
  return get_kmers<Rcpp::NumericVector, double>(s, d, alphabet,
                                         [](double d) { return std::to_string(d); },
                                         get_kmer_decorator(positional),
                                         positional);
}

struct MapReduceWorker: public RcppParallel::Worker {
  
  std::unordered_map<std::string, int> output;
  
  std::function<std::unordered_map<std::string, int>(int)> proc;
  
  std::mutex mr_mutex;
  
  MapReduceWorker(std::unordered_map<std::string, int> output,
            std::function<std::unordered_map<std::string, int>(int)> proc)
    : output(output), proc(proc) {
  }
  
  void operator()(std::size_t begin, std::size_t end) {
    std::unordered_map<std::string, int> tmp_res;
    for(int r_ind = begin; r_ind < end; ++r_ind) {
      for(auto& pair: proc(r_ind)) {
        tmp_res[pair.first] += pair.second;
      }
    }
    
    mr_mutex.lock();
    for(auto& pair: tmp_res) {
      output[pair.first] += pair.second; 
    }
    mr_mutex.unlock();
  }
};


//' @name count_kmers_larger_than_one
//' @title Count k-mers that contains more than one item
//' 
//' @param m  \code{character} matrix - each row represents one sequence
//' @param d  an \code{integer} vector representing gaps between consecutive elements of k-mer
//' @param alphabet a \code{numeric} vector representing valid elements of k-mer
//' @param pos a \code{logical} value that denotes whether positional k-mers should be generated
//' @return a named vector with counts of k-mers
//' @example count_kmers_larger_than_one(
//' matrix(data=c("a", "b", "c", "b", "c", "a"), nrow=2),
//' c(0),
//' c("a", "b", "c"),
//' FALSE)
//' @importFrom  RcppParallel RcppParallelLibs
//' @export
// [[Rcpp::export]]
std::unordered_map<std::string, int> count_kmers_larger_than_one(Rcpp::StringMatrix& m,
                                                 Rcpp::IntegerVector& d,
                                                 Rcpp::StringVector& alphabet,
                                                 Rcpp::LogicalVector& pos) {
  std::function<std::unordered_map<std::string, int>(int)> proc = [&m, &d, &alphabet, &pos](int r_ind) {
    Rcpp::StringVector v = m(r_ind, Rcpp::_);
    return count_kmers_str(v, d, alphabet, pos);
  };
  
  std::unordered_map<std::string, int> res;
  MapReduceWorker worker(res, proc);
  RcppParallel::parallelFor(0, m.rows(), worker);
  
  return worker.output;
}

//' @name count_unigrams
//' @title Count unigrams
//' @param m  \code{string} matrix that contains one sequence in each row
//' @param alphabet  \code{string} vector that contains valid elements to construct unigrams
//' @param pos  \code{logical} vector denoting whether to count positional k-mers
//' @return named \code{integer} vector with unigrams' counts 
//' @export
// [[Rcpp::export]]
std::unordered_map<std::string, int> count_unigrams(Rcpp::StringMatrix& m,
                                                    Rcpp::StringVector& alphabet,
                                                    Rcpp::LogicalVector& pos) {
  std::unordered_map<std::string, int> is_allowed;
  for(auto& a: alphabet) {
    is_allowed[static_cast<std::string>(a)] = true;
  }
  
  bool is_positional = is_first_true(pos);
  KMER_DECORATOR decorator = get_kmer_decorator(is_positional);
  
  std::function<std::unordered_map<std::string, int>(int)> proc = [&decorator, &m, &is_allowed](int row_index) {
    std::unordered_map<std::string, int> kmers;
    Rcpp::StringMatrix::Row row = m(row_index, Rcpp::_);
    for(int c=0; c < row.size(); ++c) {
      std::string str = static_cast<std::string>(row[c]);
      if(is_allowed[str]) {
        kmers[decorator(str, c)]++;
      }
    }
    return kmers;
  };
  
  std::unordered_map<std::string, int> res;
  MapReduceWorker worker(res, proc);
  RcppParallel::parallelFor(0, m.nrow(), worker);
  
  if(!is_positional) {
    for(auto& a: alphabet) {
      std::string str = static_cast<std::string>(a);
      if(worker.output.find(str) == worker.output.end()) {
        worker.output[str] = 0;
      }
    }
  }
  return worker.output;
}


