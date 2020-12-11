#ifndef SEQR_MARTINUS_ROBIN_HOOD_DICTIONARY_H
#define SEQR_MARTINUS_ROBIN_HOOD_DICTIONARY_H

#include "../../inst/thirdparty/robin_hood_martinus.h"
#include "helper/stl_dictionary_wrapper.h"

namespace dictionary {

    namespace internal {
        template<class K, class V, class H>
        using martinus_unordered_map = robin_hood::unordered_map<K, V, H>;
    }

    template<class K, class V, class Hash=std::hash<K>>
    class MartinusRobinHoodDictionary : public StlDictionaryWrapper<K, V, internal::martinus_unordered_map, Hash> {
    };

}

#endif //SEQR_MARTINUS_ROBIN_HOOD_DICTIONARY_H
