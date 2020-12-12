#ifndef SEQR_SKA_FLAT_HASH_MAP_WRAPPER_H
#define SEQR_SKA_FLAT_HASH_MAP_WRAPPER_H

#include "../../inst/thirdparty/dictionaries/ska/flat_hash_map.hpp"
#include "helper/stl_like_dictionary_wrapper.h"

namespace dictionary {

    template<class K, class V, class Hash>
    class SkaFlatHashMapWrapper : public StlLikeDictionaryWrapper<K, V, ska::flat_hash_map> {
    };

}

#endif //SEQR_SKA_FLAT_HASH_MAP_WRAPPER_H
