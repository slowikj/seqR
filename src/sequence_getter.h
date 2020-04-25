#ifndef SEQUENCE_GETTER_H
#define SEQUENCE_GETTER_H

#include <Rcpp.h>
#include <functional>
#include <memory>
#include <vector>

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

        Row(Row &&) = default;

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

template<class source_t, class target_t>
inline target_t convert(const source_t& source) {
    return Rcpp::as<target_t>(source);
}

template<>
inline int convert(const int& source) {
    return source;
}

template<>
inline double convert(const double& source) {
    return source;
}

template<class cell_t>
class SafeMatrixSequenceWrapper : public SafeSequencesWrapper<cell_t> {
public:

    template <class matrix_t>
    explicit SafeMatrixSequenceWrapper(const matrix_t &matrix) {
        initSequences(matrix);
    }

    SafeMatrixSequenceWrapper(SafeMatrixSequenceWrapper<cell_t> &&) = default;

    SafeMatrixSequenceWrapper(const SafeMatrixSequenceWrapper<cell_t> &) = delete;

    SafeMatrixSequenceWrapper &operator=(SafeMatrixSequenceWrapper<cell_t> &&) = delete;

    SafeMatrixSequenceWrapper &operator=(const SafeMatrixSequenceWrapper<cell_t> &) = delete;

    SafeMatrixSequenceWrapper() = delete;

private:

    template <class matrix_t>
    inline void initSequences(const matrix_t &matrix) {
        this->sequences_.resize(matrix.nrow());
        for (int row = 0; row < matrix.nrow(); ++row) {
            setSequenceRow(matrix, row);
        }
    }

    template <class matrix_t>
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
    explicit SafeTidysqSequencesWrapper(const std::vector<Rcpp::RawVector> &unpackedSequences) {
        initSequences(unpackedSequences);
    }

    SafeTidysqSequencesWrapper(SafeTidysqSequencesWrapper &&) = default;

    SafeTidysqSequencesWrapper &operator=(SafeTidysqSequencesWrapper &&) = delete;

    SafeTidysqSequencesWrapper(const SafeTidysqSequencesWrapper &) = delete;

    SafeTidysqSequencesWrapper &operator=(const SafeTidysqSequencesWrapper &) = delete;

    SafeTidysqSequencesWrapper() = delete;

private:
    inline void initSequences(const std::vector<Rcpp::RawVector> &unpackedSequences) {
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

template<class input_vector_t>
using SequenceGetter_t = std::function<input_vector_t(int)>;

SequenceGetter_t<SafeTidysqSequencesWrapper::Row> getTidysqRowGetter(SafeTidysqSequencesWrapper& safeWrapper);

template<class elem_t>
SequenceGetter_t<typename SafeSequencesWrapper<elem_t>::Row> getSafeMatrixRowGetter(SafeMatrixSequenceWrapper<elem_t>& sequenceWrapper) {
    return [&sequenceWrapper](int rowNum) -> typename SafeSequencesWrapper<elem_t>::Row {
        return std::move(sequenceWrapper.row(rowNum));
    };
}

#endif
