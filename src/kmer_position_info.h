#ifndef SEQR_KMER_POSITION_INFO_H
#define SEQR_KMER_POSITION_INFO_H

class KMerPositionInfo {
public:
    int seqNum;

    int position;

    KMerPositionInfo(int seqNum, int position) :
            seqNum(seqNum), position(position) {
    }

    KMerPositionInfo(const KMerPositionInfo &) = default;

    KMerPositionInfo &operator=(const KMerPositionInfo &) = default;

};

#endif //SEQR_KMER_POSITION_INFO_H
