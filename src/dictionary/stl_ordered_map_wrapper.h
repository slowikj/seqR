#ifndef SEQR_STL_ORDERED_MAP_WRAPPER_H
#define SEQR_STL_ORDERED_MAP_WRAPPER_H

#include <map>
#include "helper/stl_like_dictionary_wrapper.h"

namespace dictionary {

    template<class K, class V, class...>
    class StlOrderedMapWrapper : public StlLikeDictionaryWrapper<K, V, std::map> {
    };
}

#endif //SEQR_STL_ORDERED_MAP_WRAPPER_H
