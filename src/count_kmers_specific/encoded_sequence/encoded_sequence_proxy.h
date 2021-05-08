#pragma once

template <class sequences_list_t>
class EncodedSequenceProxy
{
public:
    using init_elem_t = typename sequences_list_t::init_elem_t;
    using encoded_elem_t = typename sequences_list_t::encoded_elem_t;

    EncodedSequenceProxy(
        std::size_t sequenceIndex,
        const sequences_list_t &encodedSequencesList)
        : _sequenceIndex(sequenceIndex),
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
        return _encodedSequencesList.getElem(_sequenceIndex, index);
    }

    inline init_elem_t decode(std::size_t index) const
    {
        return _encodedSequencesList.decode(_sequenceIndex, index);
    }

    inline std::size_t size() const
    {
        return _encodedSequencesList.getSequenceSize(_sequenceIndex);
    }

    inline bool isAllowed(std::size_t index) const
    {
        return _encodedSequencesList.isAllowed(_sequenceIndex, index);
    }

    inline bool areAllElementsAllowed() const
    {
        return _encodedSequencesList.areAllElementsAllowed();
    }

private:
    std::size_t _sequenceIndex;
    const sequences_list_t &_encodedSequencesList;
};