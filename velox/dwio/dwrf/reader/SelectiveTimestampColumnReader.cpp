/*
 * Copyright (c) Facebook, Inc. and its affiliates.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include "velox/dwio/dwrf/reader/SelectiveTimestampColumnReader.h"
#include "velox/dwio/common/BufferUtil.h"
#include "velox/dwio/dwrf/common/DecoderUtil.h"

namespace facebook::velox::dwrf {

using namespace dwio::common;

SelectiveTimestampColumnReader::SelectiveTimestampColumnReader(
    const std::shared_ptr<const TypeWithId>& nodeType,
    DwrfParams& params,
    common::ScanSpec& scanSpec)
    : SelectiveColumnReader(nodeType->type, params, scanSpec, nodeType) {
  EncodingKey encodingKey{fileType_->id, params.flatMapContext().sequence};
  auto& stripe = params.stripeStreams();
  version_ = convertRleVersion(stripe.getEncoding(encodingKey).kind());
  auto data = encodingKey.forKind(proto::Stream_Kind_DATA);
  bool vints = stripe.getUseVInts(data);
  seconds_ = createRleDecoder</*isSigned*/ true>(
      stripe.getStream(data, params.streamLabels().label(), true),
      version_,
      memoryPool_,
      vints,
      LONG_BYTE_SIZE);
  auto nanoData = encodingKey.forKind(proto::Stream_Kind_NANO_DATA);
  bool nanoVInts = stripe.getUseVInts(nanoData);
  nano_ = createRleDecoder</*isSigned*/ false>(
      stripe.getStream(nanoData, params.streamLabels().label(), true),
      version_,
      memoryPool_,
      nanoVInts,
      LONG_BYTE_SIZE);
}

uint64_t SelectiveTimestampColumnReader::skip(uint64_t numValues) {
  numValues = SelectiveColumnReader::skip(numValues);
  seconds_->skip(numValues);
  nano_->skip(numValues);
  return numValues;
}

void SelectiveTimestampColumnReader::seekToRowGroup(uint32_t index) {
  SelectiveColumnReader::seekToRowGroup(index);
  auto positionsProvider = formatData_->seekToRowGroup(index);
  seconds_->seekToRowGroup(positionsProvider);
  nano_->seekToRowGroup(positionsProvider);
  // Check that all the provided positions have been consumed.
  VELOX_CHECK(!positionsProvider.hasNext());
}

template <bool dense>
void SelectiveTimestampColumnReader::readHelper(RowSet rows) {
  ExtractToReader extractValues(this);
  common::AlwaysTrue filter;
  DirectRleColumnVisitor<
      int64_t,
      common::AlwaysTrue,
      decltype(extractValues),
      dense>
      visitor(filter, this, rows, extractValues);

  if (version_ == velox::dwrf::RleVersion_1) {
    decodeWithVisitor<velox::dwrf::RleDecoderV1<true>>(seconds_.get(), visitor);
  } else {
    decodeWithVisitor<velox::dwrf::RleDecoderV2<true>>(seconds_.get(), visitor);
  }

  // Save the seconds into their own buffer before reading nanos into
  // 'values_'
  dwio::common::ensureCapacity<uint64_t>(
      secondsValues_, numValues_, &memoryPool_);
  secondsValues_->setSize(numValues_ * sizeof(int64_t));
  memcpy(
      secondsValues_->asMutable<char>(),
      rawValues_,
      numValues_ * sizeof(int64_t));

  // We read the nanos into 'values_' starting at index 0.
  numValues_ = 0;
  if (version_ == velox::dwrf::RleVersion_1) {
    decodeWithVisitor<velox::dwrf::RleDecoderV1<false>>(nano_.get(), visitor);
  } else {
    decodeWithVisitor<velox::dwrf::RleDecoderV2<false>>(nano_.get(), visitor);
  }
}

void SelectiveTimestampColumnReader::read(
    vector_size_t offset,
    RowSet rows,
    const uint64_t* incomingNulls) {
  prepareRead<int64_t>(offset, rows, incomingNulls);
  VELOX_CHECK(!scanSpec_->filter());
  bool isDense = rows.back() == rows.size() - 1;
  if (isDense) {
    readHelper<true>(rows);
  } else {
    readHelper<false>(rows);
  }
  readOffset_ += rows.back() + 1;
}

namespace {
void fillTimestamps(
    Timestamp* timestamps,
    const uint64_t* nullsPtr,
    const int64_t* secondsPtr,
    const uint64_t* nanosPtr,
    vector_size_t numValues) {
  for (vector_size_t i = 0; i < numValues; i++) {
    if (!nullsPtr || !bits::isBitNull(nullsPtr, i)) {
      auto nanos = nanosPtr[i];
      uint64_t zeros = nanos & 0x7;
      nanos >>= 3;
      if (zeros != 0) {
        for (uint64_t j = 0; j <= zeros; ++j) {
          nanos *= 10;
        }
      }
      auto seconds = secondsPtr[i] + EPOCH_OFFSET;
      if (seconds < 0 && nanos != 0) {
        seconds -= 1;
      }
      timestamps[i] = Timestamp(seconds, nanos);
    }
  }
}

} // namespace

void SelectiveTimestampColumnReader::getValues(RowSet rows, VectorPtr* result) {
  // We merge the seconds and nanos into 'values_'
  auto tsValues = AlignedBuffer::allocate<Timestamp>(numValues_, &memoryPool_);
  auto rawTs = tsValues->asMutable<Timestamp>();
  auto secondsData = secondsValues_->as<int64_t>();
  auto nanosData = values_->as<uint64_t>();
  auto rawNulls = nullsInReadRange_
      ? (returnReaderNulls_ ? nullsInReadRange_->as<uint64_t>()
                            : rawResultNulls_)
      : nullptr;
  fillTimestamps(rawTs, rawNulls, secondsData, nanosData, numValues_);
  values_ = tsValues;
  rawValues_ = values_->asMutable<char>();
  getFlatValues<Timestamp, Timestamp>(rows, result, fileType_->type, true);
}

} // namespace facebook::velox::dwrf