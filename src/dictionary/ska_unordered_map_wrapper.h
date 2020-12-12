#ifndef SEQR_SKA_UNORDERED_MAP_H
#define SEQR_SKA_UNORDERED_MAP_H

#include "helper/stl_like_dictionary_wrapper.h"
#include "../../inst/thirdparty/dictionaries/ska/unordered_map.hpp"

namespace dictionary {

    template<class K, class V, class Hash>
    class SkaUnorderedMapWrapper : public StlLikeDictionaryWrapper<K, V, ska::unordered_map, Hash> {
    };

}

#endif //SEQR_SKA_UNORDERED_MAP_H
