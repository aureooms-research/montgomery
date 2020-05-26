from collections import deque

def to_digits(base, x) :
    digits = []
    while x != 0:
        digits.append(x%base)
        x //= base
    digits.reverse()
    return digits

def from_digits (base, digits) :
    x = 0
    for d in digits:
        x *= base
        x += d
    return x

def _extended_euclidean_algorithm ( a , b , sa = 1 , ta = 0 , sb = 0 , tb = 1 ) :

    assert(sa * sb <= 0)
    assert(ta * tb <= 0)

    yield (a, sa, ta)
    if b == 0:
        yield (b, sb, tb)

    else:
        q, _a = a // b, a % b
        _sa = sa - q * sb
        _ta = ta - q * tb
        yield from _extended_euclidean_algorithm( b, _a, sb, tb, _sa, _ta)

def extended_euclidean_algorithm ( a , b ) :
    for step, (d,x,y) in enumerate(_extended_euclidean_algorithm(a,b)):
        print(step, d, x, y)
        assert(b == 0 or abs(x) <= b)
        assert(a == 0 or abs(y) <= a)
        assert(d == x * a + y * b)
        if step % 2 == 0:
            assert(x > 0)
            assert(y <= 0)
        else:
            assert(x <= 0)
            assert(y > 0)
        yield d, x, y

def _montgomery ( R , N ) :

    [gcd, _, M] = deque(extended_euclidean_algorithm(R,N),2)[0]
    assert(gcd == 1)
    # M = -M % R
    # assert((N*M) % R == -1 % R)

    RmodN = R % N
    R2modN = (RmodN**2) % N
    R3modN = (RmodN * R2modN) % N

    print('M', M)
    print('RmodN', RmodN)
    print('R2modN', R2modN)
    print('R3modN', R3modN)
    return [M, RmodN, R2modN, R3modN]

def _redc ( R , N , M , T ) :

    # input: Integers R = b^k > N with gcd(R, N) = 1,
    #        Integer M in [0, R − 1] such that NM ≡ −1 mod R,
    #        Integer T in the range [0, RN − 1]

    assert(R > N)
    assert(0 <= M <= R-1)
    # assert((N*M) % R == -1 % R)
    assert(0 <= T <= R*N-1)

    # output: Integer S in the range [0, N − 1] such that S ≡ TR−1 mod N

    # m ← ((T mod R)M) mod R // Can be implemented by discarding limbs
    m = ((T % R) * M) % R
    # t ← (T + mN) / R // Can be implemented with a shift
    t = (T + m*N) // R
    # // /!\ T + mN is potentially RN - 1 + (R-1) N = 2RN - N - 1 so need one
    # // extra limb for carry ?
    # assert((t*R - T) % N == 0)
    if t >= N:
        return t - N
    else:
        return t


class Montgomery ( object ) :

    def __init__(self, R, N):
        [M, RmodN, R2modN, R3modN] = _montgomery(R,N)
        self.R = R
        self.N = N
        self.M = M
        self.RmodN = RmodN
        self.R2modN = R2modN
        self.R3modN = R3modN

    def frm(self, x):
        return _redc(self.R, self.N, self.M, (x % self.N)*self.R2modN)

    def out(self, aRmodN):
        return _redc(self.R, self.N, self.M, aRmodN)

base = 67108864

R = from_digits( base, [
    0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 1, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0
] )

N = from_digits( base, [
     2097151, 67108863,
    67108863, 67108863,
    67108863, 67108863,
    67108863, 67108863,
    67108863, 67108845
])

mont = Montgomery(R,N)

x = from_digits(base,[
     2097151, 67108863,
    67108863, 67108863,
    67108863, 67108863,
    67108861, 41921216,
    63517707, 56761380
  ])

print(x)
print(mont.out(mont.frm(x)))

M = from_digits(base,[
     7615973,  3532045,
    31788409, 17660227,
    24724318, 21192272,
    56512727, 38852500,
    14128181, 60044773
  ])

print(M)
