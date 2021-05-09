#pragma once

#include <algorithm>
#include <functional>
#include <iterator>
#include <memory>
#include <unordered_map>
#include <vector>

#include "helper/stl_like_dictionary_wrapper.h"

namespace dictionary {

template <class K, class V, class Hash = std::hash<K>>
class StlUnorderedMapWrapper : public StlLikeDictionaryWrapper<K, V, std::unordered_map, Hash> {
};
}  // namespace dictionary
