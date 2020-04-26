#ifndef SEQR_SAFE_SEQUENCES_WRAPPER_H
#define SEQR_SAFE_SEQUENCES_WRAPPER_H

#include "rcpp_to_cpp_converters.h"

template<class elem_t>
class SafeSequencesWrapper {
public:

    class Row {
    public:
        explicit Row(const std::vector<elem_t> &row) :
                row(row) {
        }

        Row() = delete;

        Row(const Row &) = delete;

        Row(Row &&) noexcept = default;

        Row &operator=(const Row &) = delete;

        Row &operator=(Row &&) = delete;

        const elem_t &operator[](std::size_t index) const {
            return row[index];
        }

        std::size_t size() const {
            return row.size();
        }

    private:
        const std::vector<elem_t> &row;
    };

    inline Row row(std::size_t index) const {
        return std::move(Row(sequences_[index]));
    }

    inline std::size_t size() const {
        return sequences_.size();
    }

protected:
    std::vector<std::vector<elem_t>> sequences_;

    SafeSequencesWrapper() = default;
};

template<class cell_t>
class SafeMatrixSequenceWrapper : public SafeSequencesWrapper<cell_t> {
public:

    template<class matrix_t>
    explicit SafeMatrixSequenceWrapper(const matrix_t &matrix) {
        initSequences(matrix);
    }

    SafeMatrixSequenceWrapper(SafeMatrixSequenceWrapper<cell_t> &&) noexcept = default;

    SafeMatrixSequenceWrapper(const SafeMatrixSequenceWrapper<cell_t> &) = delete;

    SafeMatrixSequenceWrapper &operator=(SafeMatrixSequenceWrapper<cell_t> &&) = delete;

    SafeMatrixSequenceWrapper &operator=(const SafeMatrixSequenceWrapper<cell_t> &) = delete;

    SafeMatrixSequenceWrapper() = delete;

private:

    template<class matrix_t>
    inline void initSequences(const matrix_t &matrix) {
        this->sequences_.resize(matrix.nrow());
        for (int row = 0; row < matrix.nrow(); ++row) {
            setSequenceRow(matrix, row);
        }
    }

    template<class matrix_t>
    inline void setSequenceRow(const matrix_t &matrix, int row) {
        this->sequences_[row].resize(matrix.ncol());
        for (int column = 0; column < matrix.ncol(); ++column) {
            auto matrixCell = matrix(row, column);
            this->sequences_[row][column] = convert<decltype(matrixCell), cell_t>(matrix(row, column));
        }
    }
};

class SafeTidysqSequencesWrapper : public SafeSequencesWrapper<unsigned char> {
public:
    explicit SafeTidysqSequencesWrapper(const Rcpp::List &unpackedSequences) {
        initSequences(unpackedSequences);
    }

    SafeTidysqSequencesWrapper(SafeTidysqSequencesWrapper &&) = default;

    SafeTidysqSequencesWrapper &operator=(SafeTidysqSequencesWrapper &&) = delete;

    SafeTidysqSequencesWrapper(const SafeTidysqSequencesWrapper &) = delete;

    SafeTidysqSequencesWrapper &operator=(const SafeTidysqSequencesWrapper &) = delete;

    SafeTidysqSequencesWrapper() = delete;

private:
    inline void initSequences(const Rcpp::List &unpackedSequences) {
        this->sequences_.resize(unpackedSequences.size());
        for (int seqNum = 0; seqNum < unpackedSequences.size(); ++seqNum) {
            setSequenceRow(unpackedSequences[seqNum], seqNum);
        }
    }

    inline void setSequenceRow(const Rcpp::RawVector &row, int seqNum) {
        this->sequences_[seqNum].resize(row.size());
        for (int elemIndex = 0; elemIndex < row.size(); ++elemIndex) {
            this->sequences_[seqNum][elemIndex] = row[elemIndex];
        }
    }
};


#endif //SEQR_SAFE_SEQUENCES_WRAPPER_H
