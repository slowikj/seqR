#ifndef SEQR_ORDERED_MAP_WRAPPER_H
#define SEQR_ORDERED_MAP_WRAPPER_H

#include <map>
#include "./stl_dictionary_wrapper.h"

namespace dictionary {

    template<class K, class V, class...>
    class OrderedMapWrapper : public StlDictionaryWrapper<K, V, std::map> {
    };
}

#endif //SEQR_ORDERED_MAP_WRAPPER_H
