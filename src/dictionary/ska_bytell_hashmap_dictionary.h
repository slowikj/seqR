#ifndef SEQR_SKA_BYTELL_HASHMAP_H
#define SEQR_SKA_BYTELL_HASHMAP_H

#include "../../inst/thirdparty/dictionaries/ska/bytell_hashmap.hpp"
#include "./helper/stl_like_dictionary_wrapper.h"

namespace dictionary {

    template<class K, class V, class Hash>
    class SkaBytellHashmapWrapper : public StlLikeDictionaryWrapper<K, V, ska::bytell_hash_map, Hash> {
    };

}

#endif //SEQR_SKA_BYTELL_HASHMAP_H
