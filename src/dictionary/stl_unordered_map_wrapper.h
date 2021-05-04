#pragma once

#include <vector>
#include <memory>
#include <iterator>
#include <unordered_map>
#include <functional>
#include <algorithm>
#include "helper/stl_like_dictionary_wrapper.h"

namespace dictionary
{

    template <class K, class V, class Hash = std::hash<K>>
    class StlUnorderedMapWrapper : public StlLikeDictionaryWrapper<K, V, std::unordered_map, Hash>
    {
    };
}
