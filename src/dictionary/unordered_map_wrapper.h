#ifndef CUSTOM_DICTIONARY_H
#define CUSTOM_DICTIONARY_H

#include <vector>
#include <memory>
#include <iterator>
#include <unordered_map>
#include <functional>
#include <algorithm>
#include "helper/stl_dictionary_wrapper.h"

namespace dictionary {

    template<class K, class V, class Hash=std::hash<K>>
    class UnorderedMapWrapper : public StlDictionaryWrapper<K, V, std::unordered_map, Hash> {
    };
}

#endif
