#pragma once

#include "hash/polynomial_single_hasher.h"
#include "hash/globals.h"
#include "hash/prefix_sequence_polynomial_hasher.h"
#include <vector>
#include <memory>
#include <utility>
#include <algorithm>

namespace gappedKMers
{

	template <class encoded_sequence_t,
			  class kmer_manager_t>
	inline kmer_manager_t count(
		const encoded_sequence_t &sequence,
		const std::vector<int> &gaps,
		std::size_t totalKMerSize,
		bool isPositionalKMer,
		bool withKMerCounts,
		const std::vector<hashing::PolynomialSingleHasherConfig> &hasherConfigs);

	template <class vector_t>
	inline std::vector<std::pair<int, int>> getContiguousIntervals(const vector_t &gaps);

	template <class encoded_sequence_t>
	class SequenceWrapper;

	// ------------------ IMPLEMENTATION ------------------

	template <class encoded_sequence_t,
			  class kmer_manager_t>
	inline kmer_manager_t count(
		const encoded_sequence_t &sequence,
		const std::vector<int> &gaps,
		std::size_t totalKMerSize,
		bool isPositionalKMer,
		const std::vector<hashing::PolynomialSingleHasherConfig> &hasherConfigs)
	{
		SequenceWrapper<encoded_sequence_t> sequenceWrapper(sequence, hasherConfigs);
		std::vector<std::pair<int, int>> contiguousIntervals = getContiguousIntervals(gaps);

		kmer_manager_t kMerManager;
		int lastSequenceIndex = static_cast<int>(sequence.size()) - totalKMerSize + 1;
		for (int seqInd = 0; seqInd < lastSequenceIndex; ++seqInd)
		{
			if (sequenceWrapper.isGappedKMerAllowed(seqInd, contiguousIntervals))
			{
				auto hash = sequenceWrapper.getGappedKMerHash(seqInd, contiguousIntervals, isPositionalKMer);
				kMerManager.add(std::move(hash), seqInd);
			}
		}

		return kMerManager;
	}

	template <class vector_t>
	inline std::vector<std::pair<int, int>> getContiguousIntervals(
		const vector_t &gaps)
	{
		std::vector<std::pair<int, int>> res;
		int currentKMerIndex = 0;
		for (int gapIndex = 0; gapIndex < gaps.size(); ++gapIndex)
		{
			int beginGapIndex = gapIndex;
			while (gapIndex < gaps.size() && gaps[gapIndex] == 0)
			{
				++gapIndex;
			}
			res.emplace_back(currentKMerIndex,
							 currentKMerIndex + (gapIndex - beginGapIndex));

			if (gapIndex < gaps.size())
			{
				currentKMerIndex += gaps[gapIndex] + (gapIndex - beginGapIndex) + 1;
			}

			if (gapIndex == gaps.size() - 1)
			{
				res.emplace_back(currentKMerIndex, currentKMerIndex);
			}
		}
		return res;
	}

	template <class encoded_sequence_t>
	class SequenceWrapper
	{
	public:
		SequenceWrapper(const encoded_sequence_t &sequence,
						const std::vector<hashing::PolynomialSingleHasherConfig> &hasherConfigs)
			: _sequence(sequence),
			  _hasherConfigs(hasherConfigs),
			  _sequenceHasher(sequence, hasherConfigs)
		{
			if (!sequence.areAllElementsAllowed())
			{
				_optionalNotAllowedPrefixCount = _prepareNotAllowedItemsPrefixCount(sequence);
			}
		}

		inline bool isGappedKMerAllowed(
			int seqBegin,
			const std::vector<std::pair<int, int>> &contiguousKMerIntervals) const
		{
			return _sequence.areAllElementsAllowed() ||
				std::all_of(
				   std::begin(contiguousKMerIntervals),
				   std::end(contiguousKMerIntervals),
				   [this, &seqBegin]
				   (const std::pair<int, int> &interval) -> bool {
					   int notAllowedItems = (_optionalNotAllowedPrefixCount[seqBegin + interval.second] -
											  ((seqBegin + interval.first) == 0
												   ? 0
												   : _optionalNotAllowedPrefixCount[seqBegin + interval.first - 1]));
					   return notAllowedItems == 0;
				   });
		}

		inline typename hashing::PrefixSequencePolynomialHasher<encoded_sequence_t>::hash_t
		getGappedKMerHash(
			int beginPosition,
			const std::vector<std::pair<int, int>> &contiguousKMerIntervals,
			bool isPositionalKMer) const
		{
			auto res = _sequenceHasher.getHashForSeveralIntervals(beginPosition, contiguousKMerIntervals);
			if (isPositionalKMer)
			{
				res.push_back(beginPosition);
			}
			return res;
		}

	private:
		const encoded_sequence_t &_sequence;
		const std::vector<hashing::PolynomialSingleHasherConfig> &_hasherConfigs;
		hashing::PrefixSequencePolynomialHasher<encoded_sequence_t> _sequenceHasher;
		std::vector<int> _optionalNotAllowedPrefixCount;

		inline std::vector<int> _prepareNotAllowedItemsPrefixCount(
			const encoded_sequence_t &sequence) const
		{
			std::vector<int> res;
			res.reserve(sequence.size());
			for (int i = 0; i < sequence.size(); ++i)
			{
				bool isNotPresent = !sequence.isAllowed(i);
				res[i] = (i == 0) ? isNotPresent : res[i - 1] + isNotPresent;
			}
			return res;
		}
	};

}
