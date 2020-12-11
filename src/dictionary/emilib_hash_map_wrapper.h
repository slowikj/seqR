#ifndef SEQR_EMILIB_HASH_MAP_WRAPPER_H
#define SEQR_EMILIB_HASH_MAP_WRAPPER_H

#include <functional>
#include "../../inst/thirdparty/emilib/hash_map.hpp"
#include "helper/stl_like_dictionary_wrapper.h"

namespace dictionary {

    namespace internal {
        template<class K, class V, class H>
        using emilib_unordered_map = robin_hood::unordered_map<K, V, H>;
    }

    template<class K, class V, class Hash=std::hash<K>>
    class EmilibHashMapWrapper : public StlLikeDictionaryWrapper<K, V, internal::emilib_unordered_map, Hash> {
    };

}

#endif //SEQR_EMILIB_HASH_MAP_WRAPPER_H
