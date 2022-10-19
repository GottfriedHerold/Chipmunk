from functools import cache
from typing import Optional
import math

# NOTE: All bounds are assumed to be with <=, corresponding to "closed" balls (the notion closed is not appropriate/may not mean that in discrete setting, which is why we won't use it)
# The computations of volume here are exact, not (terribly) optimized for speed.

@cache
def infinity_ball(*, bound:int, n: int) -> tuple[int, float]:
    """
    computes the volume of a discrete ball wrt. the infinity norm in dimension n.
    The first returned value is the number of elements, the second the Log_2 of that number.
    """
    assert type(bound) is int
    assert type(n) is int
    assert bound >= 0
    volume = (1+(2*bound))**n
    bits = math.log2(volume)
    return volume, bits


def norm2_ball(*, n: int, norm_square_bound: Optional[int]=None, per_element_length: Optional[int]=None) -> tuple[int, float]:
    """
    computes the volume of a discrete ball wrt. the 2-norm in dimension n.
    exactly one of norm_square_bound or per_element_length must be not None. These determine the radius of the ball considered.
    For per_element_length, the radius is sqrt{n}*per_element_length and corresponds to the correct comparison with infinity_ball.

    The first returned value is the number of integral points in the ball. The second is the Log_2 of that number.
    """
    assert (norm_square_bound is None) or (per_element_length is None)
    assert not ((norm_square_bound is None) and (per_element_length is None))
    bound: int
    if norm_square_bound is not None:
        bound = norm_square_bound
    else:
        bound = per_element_length * per_element_length * n
    volume = _norm2_ball(bound=bound, n=n)
    bits = math.log2(volume)
    return volume, bits

@cache
def _norm2_ball(*, bound: int, n:int) -> int:
    assert type(bound) is int
    assert type(n) is int
    assert bound >= 0
    assert n >= 0
    if n == 0:
        return 1
    if bound == 0:
        return 1
    maxcoo = math.isqrt(bound)
    volume = 0
    volume+=_norm2_ball(bound=bound, n=n-1)
    for i in range(1, maxcoo+1):
        volume+=2 * _norm2_ball(bound= bound - i*i, n=n-1)
    return volume
