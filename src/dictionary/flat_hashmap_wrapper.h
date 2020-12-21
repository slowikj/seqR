#ifndef SEQR_PARALLEL_HASHMAP_WRAPPER_H
#define SEQR_PARALLEL_HASHMAP_WRAPPER_H

#include "../../inst/thirdparty/dictionaries/parallel_hashmap/phmap.h"

namespace dictionary {

    template<class K, class V, class Hash>
    class FlatHashMapWrapper : public StlLikeDictionaryWrapper<K, V, phmap::flat_hash_map, Hash> {
    };

}
#endif //SEQR_PARALLEL_HASHMAP_WRAPPER_H
