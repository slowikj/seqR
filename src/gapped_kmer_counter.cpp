#include "gapped_kmer_counter.h"

std::vector<std::pair<int, int>> getContiguousIntervals(const Rcpp::IntegerVector& gaps) {
  std::vector<std::pair<int,int>> res;
  int currentKMerIndex = 0;
  for(int gapIndex = 0; gapIndex < gaps.size(); ++gapIndex) {
    int beginGapIndex = gapIndex;
    while(gapIndex < gaps.size() && gaps[gapIndex] == 0) {
      ++gapIndex;
    }
    res.push_back(std::make_pair(currentKMerIndex,
                                 currentKMerIndex + (gapIndex - beginGapIndex)));
    
    if(gapIndex < gaps.size()) {
      currentKMerIndex += gaps[gapIndex] + (gapIndex - beginGapIndex) + 1;
    }
    
    if(gapIndex == gaps.size() - 1) {
      res.push_back(std::make_pair(currentKMerIndex, currentKMerIndex));
    }
  }
  return res;
}
