#pragma once

#include <vector>

template <class init_elem_t_, class encoded_elem_t_>
class _EncodedSequence
{
public:
    using init_elem_t = init_elem_t_;
    using encoded_elem_t = encoded_elem_t_;

    const static std::size INVALID_ELEM = EncodedSequences<init_elem_t, encoded_elem_t>::INVALID_ELEM;

    _EncodedSequence(
        std::size_t sequenceNum,
        const EncodedSequences<init_elem_t, encoded_elem_t> &encodedSequences)
        : _sequenceNum(sequenceNum),
          _encodedSequences(encodedSequences)
    {
    }

    _EncodedSequence() = delete;

    _EncodedSequence(const _EncodedSequence &) = default;

    _EncodedSequence &operator=(const _EncodedSequence &) = default;

    _EncodedSequence(_EncodedSequence &&) = default;

    _EncodedSequence &operator=(_EncodedSequence &&) = default;

    inline encoded_elem_t operator[](std::size_t index) const
    {
        return _encodedSequences.getElem(_sequenceNum, index);
    }

    inline init_elem_t decode(std::size_t index) const
    {
        return _encodedSequences.decode(_sequenceNum, index);
    }

    inline std::size_t size() const
    {
        return _encodedSequences.getSequenceSize(_sequenceNum);
    }

    inline bool isAllowed(std::size_t index) const
    {
        return _encodedSequences.isAllowed(_sequenceNum, index);
    }

private:
    std::size_t _sequenceNum;
    const EncodedSequences<init_elem_t, encoded_elem_t> &_encodedSequences;
};

template <class init_elem_t_, class encoded_elem_t_>
class EncodedSequences
{
public:
    using init_elem_t = init_elem_t_;
    using encoded_elem_t = encoded_elem_t_;
    using Entry = EncodedSequence<init_elem_t, encoded_elem_t>;

    const static std::size_t INVALID_ELEM = 1;

    EncodedSequences(
        std::vector<encoded_elem_t> &&items,
        std::vector<std::size_t> &&sequenceStarts,
        std::vector<init_elem_t> &&elemDecoder)
        : _items(std::move(items)),
          _sequenceStarts(std::move(sequenceStarts)),
          _elemDecoder(std::move(elemDecoder))
    {
    }

    EncodedSequences() = delete;

    EncodedSequences(const EncodedSequences &) = delete;

    EncodedSequences(EncodedSequences &&) = default;

    EncodedSequences &operator=(const EncodedSequences &other) = default;

    EncodedSequences &operator=(EncodedSequences &&other) = default;

    inline Entry operator[](std::size_t sequenceNum) const
    {
        return Entry(
            sequenceNum,
            *this);
    }

    inline encoded_elem_t getElem(std::size_t sequenceNum,
                                  std::size_t offset) const
    {
        return _items[_getRawIndex(sequenceNum, offset)];
    }

    inline init_elem_t decode(std::size_t sequenceNum,
                              std::size_t offset) const
    {
        return _elemDecoder[getElem(sequenceNum, offset)];
    }

    inline std::size_t gSequenceSize(std::size_t sequenceNum) const
    {
        return _items[sequenceNum + 1] - _items[sequenceNum];
    }

    inline std::size_t size() const
    {
        return _sequenceStarts.size() - 1;
    }

    inline bool isAllowed(std::size_t sequenceNum,
                          std::size_t offset) const
    {
        return getElem(sequenceNum, offset) != INVALID_ELEM;
    }

private:
    std::vector<encoded_elem_t> _items;
    std::vector<std::size_t> _sequenceStarts;
    std::vector<init_elem_t> _elemDecoder;

    inline std::size_t _getRawIndex(std::size_t sequenceNum,
                                    std::size_t offset) const
    {
        return _sequenceStarts[sequenceNum] + offset;
    }
};
