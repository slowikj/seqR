#ifndef SEQR_SAFE_SEQUENCES_WRAPPER_H
#define SEQR_SAFE_SEQUENCES_WRAPPER_H

#include "rcpp_to_cpp_converters.h"

template<class sequence_t, class elem_t>
class BaseSequencesWrapper {
public:

    class Row {
    public:
        explicit Row(const sequence_t &row) :
                row(row) {
        }

        Row() = delete;

        Row(const Row &) = delete;

        Row(Row &&) noexcept = default;

        Row &operator=(const Row &) = delete;

        Row &operator=(Row &&) noexcept = delete;

        const elem_t &operator[](std::size_t index) const {
            return row[index];
        }

        std::size_t size() const {
            return row.size();
        }

    private:
        const sequence_t &row;
    };

    inline Row row(std::size_t index) const {
        return std::move(Row(sequences_[index]));
    }

    inline std::size_t size() const {
        return sequences_.size();
    }

protected:
    std::vector<sequence_t> sequences_;

    BaseSequencesWrapper() = default;
};

template<class elem_t>
class SafeSequencesVectorWrapper : public BaseSequencesWrapper<std::vector<elem_t>, elem_t> {
public:
    explicit SafeSequencesVectorWrapper(std::vector<std::vector<elem_t>> &&sequences) {
        this->sequences_ = std::move(sequences);
    }

    SafeSequencesVectorWrapper(SafeSequencesVectorWrapper &&) noexcept = default;

    SafeSequencesVectorWrapper(const SafeSequencesVectorWrapper &) = delete;

    SafeSequencesVectorWrapper &operator=(SafeSequencesVectorWrapper &&) noexcept = default;

    SafeSequencesVectorWrapper &operator=(const SafeSequencesVectorWrapper &) = delete;

    SafeSequencesVectorWrapper() = delete;
};

template<class cell_t>
class SafeSequencesMatrixWrapper : public BaseSequencesWrapper<std::vector<cell_t>, cell_t> {
public:

    template<class matrix_t>
    explicit
    SafeSequencesMatrixWrapper(const matrix_t &matrix)
            : SafeSequencesMatrixWrapper(matrix, 0, matrix.nrow()) {}

    template<class matrix_t>
    SafeSequencesMatrixWrapper(const matrix_t &matrix, int beginRow, int endRow) {
        initSequences(matrix, beginRow, endRow);
    }

    SafeSequencesMatrixWrapper(SafeSequencesMatrixWrapper &&) noexcept = default;

    SafeSequencesMatrixWrapper(const SafeSequencesMatrixWrapper &) = delete;

    SafeSequencesMatrixWrapper &operator=(SafeSequencesMatrixWrapper &&) noexcept = default;

    SafeSequencesMatrixWrapper &operator=(const SafeSequencesMatrixWrapper &) = delete;

    SafeSequencesMatrixWrapper() = delete;

private:

    template<class matrix_t>
    inline void initSequences(const matrix_t &matrix, std::size_t beginRow, std::size_t endRow) {
        this->sequences_.resize(endRow - beginRow);
        for (std::size_t row = beginRow; row < endRow; ++row) {
            setSequenceRow(matrix, beginRow, row);
        }
    }

    template<class matrix_t>
    inline void setSequenceRow(const matrix_t &matrix, std::size_t rowOffset, std::size_t row) {
        size_t sequenceNum = row - rowOffset;
        this->sequences_[sequenceNum].resize(matrix.ncol());
        for (int column = 0; column < matrix.ncol(); ++column) {
            auto matrixCell = matrix(row, column);
            this->sequences_[sequenceNum][column] = std::move(
                    convert<decltype(matrixCell), cell_t>(matrix(row, column)));
        }
    }
};

class SafeSequencesStringListWrapper : public BaseSequencesWrapper<std::string, char> {
public:

    explicit SafeSequencesStringListWrapper(const Rcpp::List &inputVector)
            : SafeSequencesStringListWrapper(inputVector, 0, inputVector.size()) {}

    SafeSequencesStringListWrapper(const Rcpp::List &inputVector, std::size_t begin, std::size_t end) {
        this->initSequences(inputVector, begin, end);
    }

private:

    void initSequences(const Rcpp::List &inputList, int begin, int end) {
        this->sequences_.resize(end - begin);
        for (int i = begin; i < end; ++i) {
            Rcpp::StringVector listElem = inputList[i];
            this->sequences_[i - begin] = std::move(Rcpp::as<std::string>(listElem[0]));
        }
    }
};

#endif //SEQR_SAFE_SEQUENCES_WRAPPER_H