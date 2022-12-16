from functools import cache
from typing import Optional, Union
import math


def norm2squared(vec: list[int]) -> int:
    """
    computes the squared length of the given vector
    """
    return sum([x * x for x in vec])


class SqueezedByteArray:
    """
    FIFO - queue to store a sequence of c-bit numbers, where c varies.
    The numbers are packed to only take c bits (c may not be divisible by 8)
    """
    offset_write: int
    offset_read: int
    stored_bytes: bytearray

    def __init__(self, from_bytes=None):
        self.offset_write = 0
        self.offset_read = 0
        if from_bytes is None:
            self.stored_bytes = bytearray()
        else:
            self.stored_bytes = bytearray(from_bytes)

    def to_bytes(self) -> bytes:
        return bytes(self.stored_bytes)

    def append_number(self, *, bit_len: int, number: int):
        """
        appends a number to the end of the FIFO - queue, taking exactly bit_len bits
        number must be non-negative and use at most bit_len bits.
        Note that to read, bit_len must be known.
        """
        assert number >= 0
        assert number >> bit_len == 0
        if bit_len == 0:
            return
        if self.offset_write > 0:
            squeeze_bits = min(bit_len, 8 - self.offset_write)
            extra = number & ((1 << squeeze_bits) - 1)
            self.stored_bytes[-1] |= (extra << self.offset_write)
            bit_len -= squeeze_bits
            number >>= squeeze_bits
            self.offset_write += squeeze_bits
            if self.offset_write == 8:
                self.offset_write = 0
            if bit_len == 0:
                return
        extra_bits = bit_len % 8
        byte_len = (bit_len + 7) // 8  # rounded up
        self.stored_bytes += number.to_bytes(byte_len, byteorder='little')
        self.offset_write = extra_bits

    def get_number(self, *, bit_len: int) -> int:
        """
        reads (and removes) a bit_len - bit number from the start of the FIFO queue.
        """
        if bit_len == 0:
            return 0
        squeezed_bits = 0  # number of bits < 8 that we squeezed into the first byte for reading
        ret1 = 0  # squeezed_bits-bit value that was read from the first byte
        if self.offset_read > 0:  # Check if we should read a number of <8 bits
            squeezed_bits = 8 - self.offset_read
            if bit_len < squeezed_bits:
                # We need to actually read fewer than squeezed_bits,
                # and they are all strictly contained in the first byte.
                ret1 = self.stored_bytes[0] >> self.offset_read
                ret1 &= (1 << bit_len) - 1
                self.offset_read += bit_len
                return ret1
            # read out bits of first byte and update variables. They are all used.
            ret1 = self.stored_bytes[0] >> self.offset_read
            self.offset_read = 0
            self.stored_bytes = self.stored_bytes[1:]
            bit_len -= squeezed_bits
        # Now, pretend the start of the read is aligned to byte-boundary.
        byte_len = (bit_len + 7) // 8
        ret2 = int.from_bytes(self.stored_bytes[:byte_len], byteorder='little')
        ret2 &= (1 << bit_len) - 1
        self.offset_read = bit_len % 8
        self.stored_bytes = self.stored_bytes[(bit_len // 8):]
        return (ret2 << squeezed_bits) + ret1


def encode_ball(vec: list[int], *, n: Optional[int] = None, norm_bound: Optional[int] = None,
                norm_square_bound: Optional[int] = None, per_element_length: Optional[int] = None,
                threshold: Optional[int] = None, approximation: Optional[int] = None) -> int:
    """
    Encodes an element vec contained in the (closed) n-dimensional ball of the given radius as an integer in
    [0, BALLSIZE) where BALLSIZE is the number of elements in that ball.
    Exactly one of norm_bound, norm_square_bound, per_element_length must be given.

    n is optional (it is determined by the length of vec if not).

    if threshold or approximation is set, we may take more bits, but encode faster.
    Decoding only works with the same values for those.
    """
    if n is None:
        n = len(vec)
    else:
        assert n == len(vec)
    radius_squared: int = _get_norm_square_bound(n=n, norm_bound=norm_bound, norm_square_bound=norm_square_bound,
                                                 per_element_length=per_element_length)
    assert norm2squared(vec) <= radius_squared
    return _encode_ball(n, vec, radius_squared, threshold, approximation)


def encode_sphere(vec: list[int], *, n: Optional[int] = None, norm_bound: Optional[int] = None,
                  norm_square_bound: Optional[int] = None, per_element_length: Optional[int] = None,
                  threshold: Optional[int] = None, approximation: Optional[int] = None) -> int:
    """
    Encodes an element vec from the n-dimensional sphere of the given radius as an integer in
    [0, SPHERESIZE), where SPHERESIZE is the number of elements in the sphere.
    At most one of norm_bound, norm_square_bound, per_element_length must be given. n is optional
    (it is determined by the length of vec if not).
    If any of those are given, we check whether it matches vec.

    threshold and approximation control the tradeoff output-size vs. run-time of encoding and decoding.
    """
    if n is None:
        n = len(vec)
    else:
        assert n == len(vec)
    if norm_bound is None and norm_square_bound is None and per_element_length is None:
        radius_squared = norm2squared(vec)
    else:
        radius_squared: int = _get_norm_square_bound(n=n, norm_bound=norm_bound, norm_square_bound=norm_square_bound,
                                                     per_element_length=per_element_length)
        assert norm2squared(vec) == radius_squared
    return _encode_sphere(n, vec, radius_squared, threshold, approximation)


def decode_ball(encoding: int, *, n: int, norm_bound: Optional[int] = None, norm_square_bound: Optional[int] = None,
                per_element_length: Optional[int] = None, threshold: Optional[int] = None,
                approximation: Optional[int] = None) -> list[int]:
    """
    inverse to encode_ball
    """
    radius_squared: int = _get_norm_square_bound(n=n, norm_bound=norm_bound, norm_square_bound=norm_square_bound,
                                                 per_element_length=per_element_length)
    return _decode_ball(n=n, encoding=encoding, radius_squared=radius_squared, threshold=threshold,
                        approximation=approximation)


def decode_sphere(encoding: int, *, n: int, norm_bound: Optional[int] = None, norm_square_bound: Optional[int] = None,
                  per_element_length: Optional[int] = None, threshold: Optional[int] = None,
                  approximation: Optional[int] = None) -> list[int]:
    """
    inverse to encode_sphere
    """
    radius_squared: int = _get_norm_square_bound(n=n, norm_bound=norm_bound, norm_square_bound=norm_square_bound,
                                                 per_element_length=per_element_length)
    return _decode_sphere(n=n, encoding=encoding, radius_squared=radius_squared, threshold=threshold,
                          approximation=approximation)


def ball_volume(*, n: int, norm_bound: Optional[int] = None, norm_square_bound: Optional[int] = None,
                per_element_length: Optional[int] = None, approximation: Optional[int] = None) -> int:
    """
    Computes the number of integer points x in Z^n within a given ball around 0.
    The radius of the ball can be either given as norm_bound, norm_square_bound or per_element_length
    For norm_bound, this means ||x||^2 <= norm_bound^2
    For norm_square_bound, this means ||x||^2 <= norm_square_bound
    For per_element_length, this means ||x||^2 <= n * per_element_length^2 (This corresponds to infinity-norm-bound)

    Note that the bounds are with <=

    If approximation is set, the returned value is an upper bound.
    """
    radius_squared: int = _get_norm_square_bound(n=n, norm_bound=norm_bound, norm_square_bound=norm_square_bound,
                                                 per_element_length=per_element_length)
    return _ball_volume(n, radius_squared, approximation)


def sphere_volume(*, n: int, norm_bound: Optional[int] = None, norm_square_bound: Optional[int] = None,
                  per_element_length: Optional[int] = None, approximation: Optional[int] = None) -> int:
    """
    Computes the number of integer points x in Z^n within a given sphere around 0.
    The radius of the ball can be either given as norm_bound, norm_square_bound or per_element_length
    For norm_bound, this means ||x||^2 == norm_bound^2
    For norm_square_bound, this means ||x||^2 == norm_square_bound
    For per_element_length, this means ||x||^2 == n * per_element_length^2 (This corresponds to infinity-norm-bound)

    If approximation is set, the returned value is an upper bound
    """
    radius_squared: int = _get_norm_square_bound(n=n, norm_bound=norm_bound, norm_square_bound=norm_square_bound,
                                                 per_element_length=per_element_length)
    return _sphere_volume(n, radius_squared, approximation)


# NOTE: This function is lru - cached / memoized with unbounded memory usage.
@cache
def _sphere_volume(n: int, radius_squared: int, approximation: Optional[int]) -> int:
    """
    Actual function to compute the volume of a sphere of given squared radius in dimension n
    """

    if n == 0:
        return 1 if radius_squared == 0 else 1
    if n == 1:
        if radius_squared == 0:
            return 1  # only the 0 point
        sqrt = math.isqrt(radius_squared)
        # Check whether radius_squared is the square of an integer
        if sqrt * sqrt == radius_squared:
            return 2  # We have two points: +/-sqrt
        return 0  # If radius_squared is not a square, there are no points with that squared length.
    assert n >= 2
    left_dim = n // 2
    right_dim = n - left_dim
    return sum(
        [_sphere_volume(left_dim, r, approximation) * _sphere_volume(right_dim, radius_squared - r, approximation) for r
         in range(radius_squared + 1)])


# NOTE: This function is lru - cached / memoized with unbounded memory usage.
@cache
def _ball_volume(n: int, radius_squared: int, approximation: Optional[int]) -> int:
    """
    Actual function to compute the volume of a ball of given radius
    """
    if n == 0:
        return 1
    if n == 1:
        return 1 + 2 * math.isqrt(radius_squared)
    return sum([_sphere_volume(n, r, approximation) for r in range(radius_squared + 1)])


def _get_norm_square_bound(*, n: int, norm_bound: Optional[int] = None, norm_square_bound: Optional[int] = None,
                           per_element_length: Optional[int] = None) -> int:
    """
    The functions above can take either a norm_bound, norm_square_bound or per_element_length.
    This function converts it all to norm_square_bound
    """
    assert n >= 0
    if norm_square_bound is not None:
        assert norm_bound is None
        assert per_element_length is None
        return norm_square_bound
    elif norm_bound is not None:
        assert per_element_length is None
        return norm_bound * norm_bound
    else:
        assert per_element_length is not None
        return n * per_element_length * per_element_length


def _encode_ball(n: int, vec: list[int], radius_squared: int, threshold: Optional[int],
                 approximation: Optional[int]) -> int:
    actual_squared_norm = norm2squared(vec)
    if actual_squared_norm == 0:
        return 0
    offset = _ball_volume(n=n, radius_squared=actual_squared_norm - 1)
    return offset + _encode_sphere(n=n, vec=vec, radius_squared=actual_squared_norm, threshold=threshold,
                                   approximation=approximation)


def _decode_ball(n: int, encoding: int, *, radius_squared: int, threshold: Optional[int],
                 approximation: Optional[int]) -> list[int]:
    real_radius = 0
    for i in range(radius_squared):
        t = encoding - _sphere_volume(n=n, radius_squared=i)
        if t < 0:
            break
        encoding = t
        real_radius += 1
    return _decode_sphere(n, encoding, radius_squared=real_radius, approximation=approximation, threshold=threshold)


def _encode_sphere(n: int, vec: list[int], radius_squared: int, threshold: Optional[int],
                   approximation: Optional[int]) -> int:
    if n == 0:
        return 0
    if n == 1:
        if vec[0] >= 0:
            return 0
        else:
            return 1

    left_dim = n // 2
    right_dim = n - left_dim

    left_vec = vec[:left_dim]
    right_vec = vec[left_dim:]

    left_norm2 = norm2squared(left_vec)
    right_norm2 = norm2squared(right_vec)
    assert left_norm2 + right_norm2 == norm2squared(vec)

    left_encode = _encode_sphere(left_dim, left_vec, left_norm2, threshold, approximation)
    right_encode = _encode_sphere(right_dim, right_vec, right_norm2, threshold, approximation)

    offset = 0
    for i in range(left_norm2):
        offset += _sphere_volume(left_dim, i) * _sphere_volume(right_dim, radius_squared - i)

    result = offset + left_encode * _sphere_volume(right_dim, right_norm2) + right_encode

    return result


def _decode_sphere(n: int, encoding: int, *, radius_squared: int, threshold: Optional[int],
                   approximation: Optional[int]) -> list[int]:
    if n == 0:
        return []
    if n == 1:
        radius = math.isqrt(radius_squared)
        if encoding == 0:
            return [radius]
        else:
            assert encoding == 1 and radius > 0
            return [-radius]
    left_dim = n // 2
    right_dim = n - left_dim
    left_norm_2 = 0
    for i in range(radius_squared):
        t = encoding - _sphere_volume(left_dim, i) * _sphere_volume(right_dim, radius_squared - i)
        if t < 0:
            break
        encoding = t
        left_norm_2 += 1

    right_norm2 = radius_squared - left_norm_2
    left_encode = encoding // _sphere_volume(right_dim, right_norm2)
    right_encode = encoding - left_encode * _sphere_volume(right_dim, right_norm2)

    left_vec = _decode_sphere(n=left_dim, encoding=left_encode, radius_squared=left_norm_2, threshold=threshold,
                              approximation=approximation)
    right_vec = _decode_sphere(n=right_dim, encoding=right_encode, radius_squared=right_norm2, threshold=threshold,
                               approximation=approximation)
    return left_vec + right_vec
