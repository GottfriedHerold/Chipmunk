from functools import cache
from typing import Optional, Union
import math


def norm2squared(vec: list[int]) -> int:
    """
    computes the squared length of the given vector
    """
    return sum([x * x for x in vec])


class SizedNumber:
    """
    stores a non-negative integer together with a maximum number of bits used to represent it.
    (We need to keep track of the number of bits)
    """
    number: int
    bit_len: int

    def __init__(self, number, bit_len):
        assert number >= 0
        assert number >> bit_len == 0
        self.number = number
        self.bit_len = bit_len


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

    def append_sized_number(self, x: SizedNumber):
        self.append_number(bit_len=x.bit_len, number=x.number)

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
        self.stored_bytes += number.to_bytes(byte_len, byteorder='little', signed=False)
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
        ret2 = int.from_bytes(self.stored_bytes[:byte_len], byteorder='little', signed=False)
        ret2 &= (1 << bit_len) - 1
        self.offset_read = bit_len % 8
        self.stored_bytes = self.stored_bytes[(bit_len // 8):]
        return (ret2 << squeezed_bits) + ret1

    def read_all(self) -> int:
        ret = int.from_bytes(self.stored_bytes, byteorder='little', signed=False)
        ret >>= self.offset_read
        self.stored_bytes = bytearray()
        self.offset_read = 0
        self.offset_write = 0
        return ret

    @classmethod
    def from_single_int(cls, x: int):
        return SqueezedByteArray()


def encode_ball(vec: list[int], *, n: Optional[int] = None, norm_bound: Optional[int] = None,
                norm_square_bound: Optional[int] = None, per_element_length: Optional[int] = None,
                threshold: Optional[int] = None, approximation: Optional[int] = None) -> bytes:
    """
    Encodes an element vec contained in the (closed) n-dimensional ball of the given radius.
    If threshold and approximation are both None, the result is the bytes-representation of an integer in
    [0, BALLSIZE) where BALLSIZE is the number of elements in that ball.
    Exactly one of norm_bound, norm_square_bound, per_element_length must be given.

    n is optional (it is determined by the length of vec if not).

    if threshold or approximation is set, we may take more bits, but encode faster.
    Decoding only works with the same values for those.
    """

    # ensure input validity
    if n is None:
        n = len(vec)
    elif n != len(vec):
        raise ValueError("encode_ball called with explicit value for the dimension n, but this does not match"
                         " the length of the given vector")
    radius_squared: int = _get_norm_square_bound(n=n, norm_bound=norm_bound, norm_square_bound=norm_square_bound,
                                                 per_element_length=per_element_length)
    if norm2squared(vec) > radius_squared:
        raise ValueError("encode_ball called with vector that exceeds the maximum length given")

    # _encode_ball returns a list of numbers (together with the bit-lengths)
    ret_list = _encode_ball(n, vec, radius_squared, threshold, approximation)
    # squeeze the numbers into a single byte-array:
    ret = SqueezedByteArray()
    for x in ret_list:
        ret.append_sized_number(x)
    return ret.to_bytes()


def encode_sphere(vec: list[int], *, n: Optional[int] = None, norm_bound: Optional[int] = None,
                  norm_square_bound: Optional[int] = None, per_element_length: Optional[int] = None,
                  threshold: Optional[int] = None, approximation: Optional[int] = None) -> bytes:
    """
    Encodes an element vec from the n-dimensional sphere of the given radius.
    If both threshold and approximation are None, the result is the byte-representation of an integer in
    [0, SPHERE_SIZE), where SPHERE_SIZE is the number of elements in the sphere.
    At most one of norm_bound, norm_square_bound, per_element_length must be given. n is optional.
    (These values are determined or compared to vec, depending on whether they are given).

    threshold and approximation control a tradeoff output-size vs. run-time of encoding and decoding.
    decoding must use the same values for approximation, threshold.
    """

    # ensure input validity
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

    # encode_sphere returns a list of integers together with their length; we squeeze them into a single bytes object.
    ret_list = _encode_sphere(n, vec, radius_squared, threshold, approximation)
    ret = SqueezedByteArray()
    for x in ret_list:
        ret.append_sized_number(x)
    return ret.to_bytes()


def decode_ball(encoding: bytes, *, n: int, norm_bound: Optional[int] = None, norm_square_bound: Optional[int] = None,
                per_element_length: Optional[int] = None, threshold: Optional[int] = None,
                approximation: Optional[int] = None) -> list[int]:
    """
    inverse to encode_ball
    """
    radius_squared: int = _get_norm_square_bound(n=n, norm_bound=norm_bound, norm_square_bound=norm_square_bound,
                                                 per_element_length=per_element_length)
    return _decode_ball(n=n, encoding=SqueezedByteArray(from_bytes=encoding), radius_squared=radius_squared,
                        threshold=threshold, approximation=approximation)


def decode_sphere(encoding: bytes, *, n: int, norm_bound: Optional[int] = None, norm_square_bound: Optional[int] = None,
                  per_element_length: Optional[int] = None, threshold: Optional[int] = None,
                  approximation: Optional[int] = None) -> list[int]:
    """
    inverse to encode_sphere
    """
    radius_squared: int = _get_norm_square_bound(n=n, norm_bound=norm_bound, norm_square_bound=norm_square_bound,
                                                 per_element_length=per_element_length)
    enc = SqueezedByteArray(from_bytes=encoding)
    return _decode_sphere(n=n, encoding=enc, radius_squared=radius_squared, threshold=threshold,
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
    return _sphere_volume(n, radius_squared, approximation=approximation)


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
    if n < 0:
        raise ValueError("Negative dimension given")
    if norm_square_bound is not None:
        if norm_bound is not None or per_element_length is not None:
            raise ValueError("Exactly one of norm_bound, norm_square_bound and per_element_length must be given")
        return norm_square_bound
    elif norm_bound is not None:
        if per_element_length is not None:
            raise ValueError("Exactly one of norm_bound, norm_square_bound and per_element_length must be given")
        return norm_bound * norm_bound
    else:
        if per_element_length is None:
            raise ValueError("Exactly one of norm_bound, norm_square_bound and per_element_length must be given")
        return n * per_element_length * per_element_length


def _encode_ball(n: int, vec: list[int], radius_squared: int, threshold: Optional[int],
                 approximation: Optional[int]) -> list[SizedNumber]:
    if n == 0:
        return [SizedNumber(0, 0)]  # or empty list, doesn't matter, actually
    actual_squared_norm = norm2squared(vec)
    if threshold is None or n < threshold:
        if actual_squared_norm == 0:
            return [SizedNumber(0, 0)]
        offset = _ball_volume(n=n, radius_squared=actual_squared_norm - 1, approximation=approximation)
        ret: list[SizedNumber] = _encode_sphere(n=n, vec=vec, radius_squared=actual_squared_norm, threshold=threshold,
                                                approximation=approximation)
        assert len(ret) == 1
        ret[0].number += offset
        ret[0].bit_len = ret[0].number.bit_length()
        return ret
    else:
        return [SizedNumber(actual_squared_norm, radius_squared.bit_length())] + \
            _encode_sphere(n=n, vec=vec, radius_squared=actual_squared_norm,
                           threshold=threshold, approximation=approximation)



def _decode_ball(n: int, encoding: SqueezedByteArray, *, radius_squared: int, threshold: Optional[int],
                 approximation: Optional[int]) -> list[int]:
    if n == 0:
        return []
    if threshold is None or n < threshold:
        enc = encoding.read_all()
        if enc == 0:
            return [0] * n
        real_radius = radius_squared
        # find largest i, s.t. _ball_volume(n, radius_squared=i) is > encoding
        for i in range(radius_squared):
            if _ball_volume(n=n, radius_squared=i, approximation=approximation) > enc:
                real_radius = i
                break
        assert real_radius >= 1
        # subtract offset added during encoding
        enc -= _ball_volume(n=n, radius_squared=real_radius - 1, approximation=approximation)
        return _decode_sphere(n, enc, radius_squared=real_radius, approximation=approximation, threshold=threshold)
    else:
        real_radius = encoding.get_number(bit_len=radius_squared.bit_length())
        return _decode_sphere(n, encoding, radius_squared=real_radius, approximation=approximation, threshold=threshold)


def _encode_sphere(n: int, vec: list[int], radius_squared: int, threshold: Optional[int],
                   approximation: Optional[int]) -> list[SizedNumber]:
    if n == 0 or radius_squared == 0:  # requires 0 bits
        return [SizedNumber(0, 0)]
    if n == 1:  # can encode with a single bit
        if vec[0] >= 0:
            return [SizedNumber(number=0, bit_len=1)]
        else:
            return [SizedNumber(number=1, bit_len=1)]

    left_dim = n // 2
    right_dim = n - left_dim

    left_vec = vec[:left_dim]
    right_vec = vec[left_dim:]

    left_norm2 = norm2squared(left_vec)
    right_norm2 = norm2squared(right_vec)
    assert left_norm2 + right_norm2 == norm2squared(vec)

    left_encode = _encode_sphere(left_dim, left_vec, left_norm2, threshold, approximation)
    right_encode = _encode_sphere(right_dim, right_vec, right_norm2, threshold, approximation)

    if threshold is None or n < threshold:
        offset = 0
        for i in range(left_norm2):
            offset += _sphere_volume(left_dim, i, approximation=approximation) *\
                      _sphere_volume(right_dim, radius_squared - i, approximation=approximation)

        assert len(left_encode) == 1
        assert len(right_encode) == 1

        result = offset + left_encode[0].number * _sphere_volume(right_dim, right_norm2, approximation) + right_encode[0].number

        return [SizedNumber(result, _sphere_volume(n, radius_squared, approximation=approximation).bit_length())]
    else:
        return [SizedNumber(left_norm2, radius_squared.bit_length())] + left_encode + right_encode


def _decode_sphere(n: int, encoding: Union[SqueezedByteArray, int], *, radius_squared: int, threshold: Optional[int],
                   approximation: Optional[int]) -> list[int]:
    if n == 0:
        return []
    if radius_squared == 0:
        return [0] * n
    if n == 1:
        radius = math.isqrt(radius_squared)
        assert radius * radius == radius_squared
        if type(encoding) is SqueezedByteArray:
            enc = encoding.get_number(bit_len=1)
        else:
            enc = encoding
        if enc == 0:
            return [radius]
        else:
            assert enc == 1
            return [-radius]

    left_dim = n // 2
    right_dim = n - left_dim

    if threshold is None or n < threshold:
        enc: int
        if type(encoding) is SqueezedByteArray:
            enc = encoding.get_number(bit_len=_sphere_volume(n, radius_squared, approximation=approximation).bit_length())
        else:
            enc = encoding
        left_norm_2 = 0
        for i in range(radius_squared):
            t = enc - _sphere_volume(left_dim, i, approximation=approximation) *\
                      _sphere_volume(right_dim, radius_squared - i, approximation=approximation)
            if t < 0:
                break
            enc = t
            left_norm_2 += 1

        right_norm2 = radius_squared - left_norm_2
        left_enc = enc // _sphere_volume(right_dim, right_norm2, approximation)
        right_enc = enc - left_enc * _sphere_volume(right_dim, right_norm2, approximation)

        left_vec = _decode_sphere(n=left_dim, encoding=left_enc, radius_squared=left_norm_2, threshold=threshold,
                                  approximation=approximation)
        right_vec = _decode_sphere(n=right_dim, encoding=right_enc, radius_squared=right_norm2, threshold=threshold,
                                   approximation=approximation)
        return left_vec + right_vec
    else:
        assert type(encoding) == SqueezedByteArray
        left_norm2 = encoding.get_number(bit_len=radius_squared.bit_length())
        right_norm2 = radius_squared - left_norm2
        left_vec = _decode_sphere(n=left_dim, encoding=encoding, radius_squared=left_norm2, threshold=threshold,
                                  approximation=approximation)
        right_vec = _decode_sphere(n=right_dim, encoding=encoding, radius_squared=right_norm2, threshold=threshold,
                                   approximation=approximation)
        return left_vec + right_vec
