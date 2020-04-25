// [[Rcpp::plugins("cpp17")]]
#include "gapped_kmer_counter.h"
#include <Rcpp.h>

std::vector<std::pair<int, int>> getContiguousIntervals(const Rcpp::IntegerVector &gaps) {
    std::vector<std::pair<int, int>> res;
    int currentKMerIndex = 0;
    for (int gapIndex = 0; gapIndex < gaps.size(); ++gapIndex) {
        int beginGapIndex = gapIndex;
        while (gapIndex < gaps.size() && gaps[gapIndex] == 0) {
            ++gapIndex;
        }
        res.push_back(std::make_pair(currentKMerIndex,
                                     currentKMerIndex + (gapIndex - beginGapIndex)));

        if (gapIndex < gaps.size()) {
            currentKMerIndex += gaps[gapIndex] + (gapIndex - beginGapIndex) + 1;
        }

        if (gapIndex == gaps.size() - 1) {
            res.push_back(std::make_pair(currentKMerIndex, currentKMerIndex));
        }
    }
    return res;
}

int getIntervalLength(const std::pair<int, int> &interval) {
    return interval.second - interval.first + 1;
}

bool isGappedKMerAllowed(int seqBegin,
                         const std::vector<std::pair<int, int>> &contiguousKMerIntervals,
                         const std::vector<int> &notAllowedItemsPrefixCount) {
    return std::all_of(
            std::begin(contiguousKMerIntervals),
            std::end(contiguousKMerIntervals),
            [&notAllowedItemsPrefixCount, &seqBegin](const std::pair<int, int> &interval) -> bool {
                return (notAllowedItemsPrefixCount[seqBegin + interval.second]
                        - (seqBegin + interval.first == 0 ?
                           0 :
                           notAllowedItemsPrefixCount[seqBegin + interval.first - 1])) == 0;
            }
    );
}

std::size_t getTotalKMerSize(const Rcpp::IntegerVector &gaps) {
    return Rcpp::sum(gaps + 1) + 1;
}
