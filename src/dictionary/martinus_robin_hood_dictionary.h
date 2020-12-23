#pragma once

#include "../../inst/thirdparty/dictionaries/robin_hood_martinus.h"
#include "helper/stl_like_dictionary_wrapper.h"

namespace dictionary {

    namespace internal {
        template<class K, class V, class H>
        using martinus_unordered_map = robin_hood::unordered_map<K, V, H>;
    }

    template<class K, class V, class Hash=std::hash<K>>
    class MartinusRobinHoodDictionary : public StlLikeDictionaryWrapper<K, V, internal::martinus_unordered_map, Hash> {
    };

}
