/***************************************************************************
 * A simple implementation of MSD (Most Significant Digit) and LSD (Least
 * Significant Digit) for the base types: signed and unsigned integers (that are
 * 8, 16, 32, and 64 bits long), floats, and doubles.
 */

#ifndef RADIX_SORT_H
#define RADIX_SORT_H

#include <algorithm>  // For swap
#include <cstdint>    // For system independent base types
#include <cstring>    // For memset
#include <vector>

/**
 * Sorts a given array, of length N, using the MSD version of Radix sort. As of
 * right now, this function only works on base 2. This version is in-place,
 * using O(1) extra space.
 */
template <typename ValueType>
void RadixSortMSD(ValueType* array, size_t N);

/**
 * Sorts a given array, of length N, of some known base type, using the LSD
 * version of Radix sort. This version uses O(Base*N) extra space, though it
 * could be much faster, especially for larger bases.
 */
template <size_t Base, typename ValueType>
void RadixSortLSD(ValueType* array, size_t N);

/* * * * * Implementation Below This Point * * * * */
namespace details {

/*
 * For all the basic types on which we perform radix sort, we will begin by
 * converting the type to some "radixable" form. From there, we perform radix
 * sort and finally convert back to the original type. We first define the
 * radixable types as unsigned ints of the proper length.
 */
using Radixable8 = uint8_t;
using Radixable16 = uint16_t;
using Radixable32 = uint32_t;
using Radixable64 = uint64_t;

template <size_t Bytes>
struct Radixable {};

template <>
struct Radixable<1> {
  using type = Radixable8;
};
template <>
struct Radixable<2> {
  using type = Radixable16;
};
template <>
struct Radixable<4> {
  using type = Radixable32;
};
template <>
struct Radixable<8> {
  using type = Radixable64;
};

/*
 * Here we define a set of transformation rules to convert from and to the
 * "radixable" form.
 */
template <class ValueType, class RadixableType>
struct RadixTransform {
  RadixableType TransformFrom(ValueType);
  ValueType TransformTo(RadixableType);
};

/*
 * The first group of transformation rules are by far the easiest and most
 * straight forward. Unsigned numbers are already in radixable form.
 */
template <>
struct RadixTransform<uint8_t, Radixable8> {
  Radixable8 TransformFrom(uint8_t value) { return value; }
  uint8_t TransformTo(Radixable8 rdxble) { return rdxble; }
};

template <>
struct RadixTransform<uint16_t, Radixable16> {
  Radixable16 TransformFrom(uint16_t value) { return value; }
  uint16_t TransformTo(Radixable16 rdxble) { return rdxble; }
};

template <>
struct RadixTransform<uint32_t, Radixable32> {
  Radixable32 TransformFrom(uint32_t value) { return value; }
  uint32_t TransformTo(Radixable32 rdxble) { return rdxble; }
};

template <>
struct RadixTransform<uint64_t, Radixable64> {
  Radixable64 TransformFrom(uint64_t value) { return value; }
  uint64_t TransformTo(Radixable64 rdxble) { return rdxble; }
};

/*
 * Now we get a little more complicated with signed numbers. Lets begin by
 * looking at what happens if we just reinterpret the numbers as being unsigned:
 *
 * Before: 33 36 27 -35 43 -15  36  42  -1 -29
 * After:  27 33 36  36 42  43 -35 -29 -15  -1
 *
 * The result is ascending positive numbers, followed be ascending negative
 * numbers. Positive numbers ascending makes sense. In the two's complement
 * form, negative numbers have the most significant bit set, which, in the
 * unsigned interpretation makes them valued more highly.
 *
 * By flipping the high bit during the radix sort, positive numbers attain the
 * high bit and both positives and negatives sort in ascending order.
 */
template <>
struct RadixTransform<int8_t, Radixable8> {
  Radixable8 TransformFrom(int8_t value) {
    Radixable8 rdxble = *(Radixable8*)&value;
    rdxble ^= (1 << 7);
    return rdxble;
  }

  int8_t TransformTo(Radixable8 rdxble) {
    int8_t value = *(int8_t*)&rdxble;
    value ^= (1 << 7);
    return value;
  }
};

template <>
struct RadixTransform<int16_t, Radixable16> {
  Radixable16 TransformFrom(int16_t value) {
    Radixable16 rdxble = *(Radixable16*)&value;
    rdxble ^= (1 << 15);
    return rdxble;
  }

  int16_t TransformTo(Radixable16 rdxble) {
    int16_t value = *(int16_t*)&rdxble;
    value ^= (1 << 15);
    return value;
  }
};

template <>
struct RadixTransform<int32_t, Radixable32> {
  Radixable32 TransformFrom(int32_t value) {
    Radixable32 rdxble = *(Radixable32*)&value;
    rdxble ^= (1 << 31);
    return rdxble;
  }

  int32_t TransformTo(Radixable32 rdxble) {
    int32_t value = *(int32_t*)&rdxble;
    value ^= (1 << 31);
    return value;
  }
};

template <>
struct RadixTransform<int64_t, Radixable64> {
  Radixable64 TransformFrom(int64_t value) {
    Radixable64 rdxble = *(Radixable64*)&value;
    rdxble ^= (((Radixable64)1) << 63);
    return rdxble;
  }

  int64_t TransformTo(Radixable64 rdxble) {
    int64_t value = *(int64_t*)&rdxble;
    value ^= (((uint64_t)1) << 63);
    return value;
  }
};

/*
 * Now for the most interesting group of basic types: floating point. For
 * simplicity, I'm assuming IEEE 754 format, which is the most common, though
 * not standard. It has several advantages that will make it easier to
 * transform. First off the exponent bits are greater than the mantissa, and
 * greater exponent guarentees greater value. I'm going to ignore special values
 * such as infinity and NaN for simplicity. Lets begin by looking at what
 * happens without any bit manipulation:
 *
 * Before: 33 36 27 -35 43 -15  36  42  -1 -29
 * After:  27 33 36  36 42  43  -1 -15 -29 -35
 *
 * This is similar to signed integers in that we have positive numbers followed
 * by negative numbers. Flipping the high (sign) bit will be a key to swapping
 * their regions. Positive numbers are still ascending (which is nice), but
 * negative numbers are now descending (which is less nice). This is due to the
 * fact that IEEE 754 essentially stores the absolute value and a sign bit,
 * rather than two's complement. What would be really awesome is if we could
 * invoke equivalent of two's complement negative operation on the negative
 * values (so that sorting them gives positive values in ascending order), then
 * flip the high bit on those numbers (making them negative in ascending order).
 * Lets see what that gives us:
 *
 * Before: 33 36 27 -35 43 -15  36  42  -1 -29
 * After:  27 33 36  36 42  43 -35 -29 -15  -1
 *
 * This is what signed numbers initially looked like. Flipping the high bit, we
 * can radix floating point values.
 */
template <>
struct RadixTransform<float, Radixable32> {
  Radixable32 TransformFrom(float value) {
    Radixable32 rdxble = *(Radixable32*)&value;
    // Two's complement negative numbers, but keep them negative.
    if ((rdxble >> 31) == 1) {
      rdxble *= -1;
      rdxble ^= (1 << 31);
    }
    // Flip the high bit on all numbers to swap the positive and negative
    // regions.
    rdxble ^= (1 << 31);
    return rdxble;
  }

  float TransformTo(Radixable32 rdxble) {
    rdxble ^= (1 << 31);
    if ((rdxble >> 31) == 1) {
      rdxble ^= (1 << 31);
      rdxble *= -1;
    }
    return *(float*)&rdxble;
  }
};

template <>
struct RadixTransform<double, Radixable64> {
  Radixable64 TransformFrom(double value) {
    Radixable64 rdxble = *(Radixable64*)&value;
    // Two's complement negative numbers, but keep them negative.
    if ((rdxble >> 63) == 1) {
      rdxble *= -1;
      rdxble ^= (((Radixable64)1) << 63);
    }
    // Flip the high bit on all numbers to swap the positive and negative
    // regions.
    rdxble ^= (((Radixable64)1) << 63);
    return rdxble;
  }

  double TransformTo(Radixable64 rdxble) {
    rdxble ^= (((Radixable64)1) << 63);
    if ((rdxble >> 63) == 1) {
      rdxble ^= (((Radixable64)1) << 63);
      rdxble *= -1;
    }
    return *(double*)&rdxble;
  }
};

/**
 * Here is the actual implementation of MSD Radix Sort. This algorithm can be
 * thought of as a form of quick sort. We use each bit, starting with the high
 * bit, as a pivot. We move all values with a 0 at the specified but to the left
 * half and all values with a 1 to the right. We then recursively call this
 * algorithm on both halves with the next significant bit as the pivot.
 */
template <typename RadixableType>
void RadixSortMSD(RadixableType* data, size_t N, size_t bit) {
  /*
   * We technically only need to stop if N <= 1, but certain optimizations can
   * be made here to greatly speed up execution. I kept the is_sorted check,
   * which provides early termination. Another optimization that can greatly
   * reduce practical runtime is to run insertion sort when N <= 32 (or so). For
   * small arrays, divide-and-conquer adds too much overhead.
   */
  if (N <= 1 || std::is_sorted(data, data + N)) {
    return;
  }

  /*
   * We create a pointer to the right and left edges of the data. While the left
   * pointer is pointing to a number with a 0 in the current bit, shift the left
   * pointer right. We shift the right pointer left in a similar fashion as long
   * it points to a number with a 1 in the current bit. When both pointers are
   * pointing to numbers with the wrong bit, those numbers are swapped. Then
   * back to shifting left and right pointers. When the pointers overlap, we
   * stop.
   */
  int left = 0;
  int right = N;
  while (left < right) {
    size_t left_bit = ((data[left] >> bit) & 1);

    if (left_bit == 0) {
      ++left;
      continue;
    }

    size_t right_bit = ((data[right - 1] >> bit) & 1);
    if (right_bit == 1) {
      --right;
      continue;
    }

    std::swap(data[left], data[right - 1]);
    ++left;
    --right;
  }

  /*
   * As long as more bits exist, call both sides of the pivot recursively.
   */
  if (bit > 0) {
    RadixSortMSD(data, left, bit - 1);
    RadixSortMSD(data + left, N - left, bit - 1);
  }
}

/**
 * Here is the actual implementation of LSD Radix Sort. This algorithm creates a
 * queue for each possible digit in the given base. Numbers are placed into
 * queues (called buckets) based on their least significant digit, then removed
 * from the queues so that all numbers are sorted by that digit. When we move to
 * the next digit, numbers placed into the same buckets will be ordered (within
 * the bucket) based on all less-significant digits. This is continued until the
 * numbers are sorted. If the base were to be resticted to powers of 2 then the
 * number of iterations could be easily computed without a call to the log
 * function.
 *
 * We note that this function dynamically allocates O(Base*N) additional space
 * for the buckets.
 */
template <size_t Base, typename RadixableType>
void RadixSortLSD(RadixableType* data, size_t N, size_t) {
  /*
   * We create Base buckets of length N. In the worst case, if all N numbers
   * have the same digit in some position, a bucket would need to be able to
   * hold all numbers.
   *
   * For efficiency and simplicity, we never clean the buckets. Instead we
   * create a set of Base pointers which keep track of how much of the bucket we
   * are using. For example, if pointers[2] = 5, then 5 numbers so far have a 2
   * digit in the current slot. If we encounter another 2 digit, we place the
   * number into buckets[2][pointers[2]], then increment pointers[2]. Pointers
   * are cleaned at every digit.
   */
  std::vector<RadixableType[Base]> buckets(N);
  size_t pointers[Base];

  for (size_t div = 1; !std::is_sorted(data, data + N); div *= Base) {
    std::memset(pointers, 0, sizeof(pointers));

    for (size_t i = 0; i < N; ++i) {
      size_t bkt = (data[i] / div) % Base;
      buckets[pointers[bkt]++][bkt] = data[i];
    }

    int index = 0;
    for (size_t b = 0; b < Base; ++b) {
      for (size_t p = 0; p < pointers[b]; ++p) {
        data[index++] = buckets[p][b];
      }
    }
  }
}

/**
 * This function takes in the original data, transforms it to some radixable
 * form, invokes the selected radix function, then transforms the data back to
 * the original form (now sorted). To allow passing function pointers,
 * RadixSortMSD and RadixSortLSD needed to have the same signature, so
 * RadixSortLSD was given a dummy size_t argument. Since the function pointer
 * needs to know what RadixableType type is in advance, RadixableType needs to
 * be defined prior to this function execution.
 */
template <typename RadixableType, typename ValueType>
void TransformAndInvoke(ValueType* array, size_t N,
                        void (*RadixFunc)(RadixableType*, size_t, size_t)) {
  RadixTransform<ValueType, RadixableType> transformer;

  RadixableType* radixable_array = reinterpret_cast<RadixableType*>(array);

  for (size_t i = 0; i < N; ++i) {
    radixable_array[i] = transformer.TransformFrom(array[i]);
  }

  RadixFunc(radixable_array, N, (8 * sizeof(RadixableType)) - 1);

  for (size_t i = 0; i < N; ++i) {
    array[i] = transformer.TransformTo(radixable_array[i]);
  }
}

}  // namespace details

/* Actual implementation of MSD Radix Sort. */
template <typename ValueType>
void RadixSortMSD(ValueType* array, size_t N) {
  using RadixableType = typename details::Radixable<sizeof(ValueType)>::type;
  details::TransformAndInvoke<RadixableType>(array, N, details::RadixSortMSD);
}

/* Actual implementation of LSD Radix Sort. */
template <size_t Base, typename ValueType>
void RadixSortLSD(ValueType* array, size_t N) {
  using RadixableType = typename details::Radixable<sizeof(ValueType)>::type;
  details::TransformAndInvoke<RadixableType>(array, N,
                                             details::RadixSortLSD<Base>);
}

#endif  // RADIX_SORT_H
