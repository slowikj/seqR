#ifndef SEQR_EMILIB_HASH_MAP_WRAPPER_H
#define SEQR_EMILIB_HASH_MAP_WRAPPER_H

#include <functional>
#include "../../inst/thirdparty/emilib/hash_map.hpp"

namespace dictionary {

    template<class K, class V, class Hash=std::hash<K>>
    class EmilibHashMapWrapper : public StlDictionaryWrapper<K, V, emilib::HashMap, Hash> {
    };

}

#endif //SEQR_EMILIB_HASH_MAP_WRAPPER_H
