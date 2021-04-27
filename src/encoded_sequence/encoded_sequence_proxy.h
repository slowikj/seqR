#pragma once

#include "encoded_sequence_proxy.h"

template <class sequences_list_t>
class EncodedSequenceProxy
{
public:
    using init_elem_t = typename sequences_list_t::init_elem_t;
    using encoded_elem_t = typename sequences_list_t::encoded_elem_t;

    EncodedSequenceProxy(
        std::size_t sequenceNum,
        const sequences_list_t &encodedSequencesList)
        : _sequenceNum(sequenceNum),
          _encodedSequencesList(encodedSequencesList)
    {
    }

    EncodedSequenceProxy() = delete;

    EncodedSequenceProxy(const EncodedSequenceProxy &) = default;

    EncodedSequenceProxy &operator=(const EncodedSequenceProxy &) = default;

    EncodedSequenceProxy(EncodedSequenceProxy &&) = default;

    EncodedSequenceProxy &operator=(EncodedSequenceProxy &&) = default;

    inline encoded_elem_t operator[](std::size_t index) const
    {
        return _encodedSequencesList.getElem(_sequenceNum, index);
    }

    inline init_elem_t decode(std::size_t index) const
    {
        return _encodedSequencesList.decode(_sequenceNum, index);
    }

    inline std::size_t size() const
    {
        return _encodedSequencesList.getSequenceSize(_sequenceNum);
    }

    inline bool isAllowed(std::size_t index) const
    {
        return _encodedSequencesList.isAllowed(_sequenceNum, index);
    }

private:
    std::size_t _sequenceNum;
    const RawEncodedSequencesList<init_elem_t, encoded_elem_t> &_encodedSequencesList;
};