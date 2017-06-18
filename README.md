# Radix Sort
A simple implementation of MSD (Most Significant Digit) and LSD (Least Significant Digit) for the base types: signed and unsigned integers (that are 8, 16, 32, and 64 bits long), floats, and doubles.

This is not meant to be an efficient implementation of either algorithm, rather a proof-of-concept that a **_O(wn)_** sorting algorithm exists for floating point numbers (using the common IEEE format).

## Handling Unsigned Ints
This first group is by far the easiest and most straight forward. Unsigned numbers are already in radixable form.

```c++
struct RadixTransform<uint32_t, Radixable32> {
  Radixable32 TransformFrom(uint32_t value) { return value; }
  uint32_t TransformTo(Radixable32 rdxble) { return rdxble; }
};
```

## Handling Signed Ints
Now we get a little more complicated with signed numbers. Lets begin by looking at what happens if we just reinterpret the numbers as being unsigned:

    Before: 33 36 27 -35 43 -15  36  42  -1 -29
    After:  27 33 36  36 42  43 -35 -29 -15  -1

The result is ascending positive numbers, followed be ascending negative numbers. Positive numbers ascending makes sense. In the two's complement form, negative numbers have the most significant bit set, which, in the unsigned interpretation makes them valued more highly.


By flipping the high bit during the radix sort, positive numbers attain the high bit and both positives and negatives sort in ascending order.

```c++
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
```

## Handling Floating Point
Now for the most interesting group of basic types: floating point. For simplicity, I'm assuming IEEE 754 format, which is the most common, though not standard. It has several advantages that will make it easier to transform. First off the exponent bits are greater than the mantissa, and greater exponent guarentees greater value. I'm going to ignore special values such as infinity and NaN for simplicity. Lets begin by looking at what happens without any bit manipulation:

    Before: 33 36 27 -35 43 -15  36  42  -1 -29
    After:  27 33 36  36 42  43  -1 -15 -29 -35

This is similar to signed integers in that we have positive numbers followed by negative numbers. Flipping the high (sign) bit will be a key to swapping their regions. Positive numbers are still ascending (which is nice), but negative numbers are now descending (which is less nice). This is due to the fact that IEEE 754 essentially stores the absolute value and a sign bit, rather than two's complement. What would be really awesome is if we could invoke equivalent of two's complement negative operation on the negative values (so that sorting them gives positive values in ascending order), then flip the high bit on those numbers (making them negative in ascending order).

Lets see what that gives us:

    Before: 33 36 27 -35 43 -15  36  42  -1 -29
    After:  27 33 36  36 42  43 -35 -29 -15  -1

This is what signed numbers initially looked like. Flipping the high bit, we can radix floating point values.

```c++
struct RadixTransform<float, Radixable32> {
  Radixable32 TransformFrom(float value) {
    Radixable32 rdxble = *(Radixable32*)&value;
    // Two's complement negative numbers, but keep them negative.
    if ((rdxble >> 31) == 1) {
      rdxble *= -1;
      rdxble ^= (1 << 31);
    }
    // Flip the high bit on all numbers to swap the positive and negative regions.
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
```
