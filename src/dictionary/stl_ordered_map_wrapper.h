#pragma once

#include <map>
#include "helper/stl_like_dictionary_wrapper.h"

namespace dictionary {

    template<class K, class V, class...>
    class StlOrderedMapWrapper : public StlLikeDictionaryWrapper<K, V, std::map> {
    };
}

