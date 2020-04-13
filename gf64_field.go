package gf

import (
	"crypto/rand"
	"encoding/binary"
)

// We define binary field under the irreducible polynomial
// R = x^64 + x^4 + x^3 + x + 1
// R = x^64 + r
// gf64MOD is the lower part of the irreducible polynomial
// r = x^4 + x^3 + x + 1
var gf64MOD GF64 = 0x1b

// GF64 is field element of GF(64).
type GF64 uint64

// NewGF64 returns a new GF(64) element which is zero.
func NewGF64() GF64 {
	return GF64(0)
}

// randGF64 generates a new random GF(64) element.
func randGF64() GF64 {
	buf := make([]byte, 64)
	var e GF64
	for e.IsZero() {
		_, err := rand.Read(buf)
		if err != nil {
			panic(err)
		}
		e = GF64(binary.BigEndian.Uint64(buf))
	}
	return e
}

// Zero returns new field element equals to Zero
func (e GF64) Zero() GF64 {
	return GF64(0)
}

// One returns new field element equals to One
func (e GF64) One() GF64 {
	return GF64(1)
}

// IsZero returns true if GF(64) element is equal to zero
func (e GF64) IsZero() bool {
	return e.Equal(GF64(0))
}

// IsZero returns true if GF(64) element is equal to zero
func (e GF64) IsOne() bool {
	return e.Equal(GF64(1))
}

// Equal tests if given two element is equal.
func (e0 GF64) Equal(e1 GF64) bool {
	return e0 == e1
}

// Add adds two GF(64) field element.
func (e0 GF64) Add(e1 GF64) GF64 {
	return e0 ^ e1
}

// Sub subtracts two GF(64) field element.
func (e0 GF64) Sub(e1 GF64) GF64 {
	return e0.Add(e1)
}

// Neg negates the GF(64) field element.
func (e0 GF64) Neg() GF64 {
	return e0
}

// Mul multiplies two GF(64) field elements.
func (e0 GF64) Mul(e1 GF64) GF64 {
	return mul64(e0, e1)
}

// Square squares GF(64) field element.
func (e0 GF64) Square() GF64 {
	return mul64(e0, e0)
}

// Inverse inverses a GF(64) element.
func (e0 GF64) Inverse() GF64 {
	if e0.IsZero() {
		return e0
	}
	// Chain is generated with tool below.
	// https://github.com/kwantam/addchain
	// Run:
	// $ ./addchain 2^64-2
	t0 := e0          // 0
	t1 := t0.Square() // 1
	t2 := t0.Mul(t1)  // 2
	t0 = t2.Square()  // 3
	t0 = t0.Square()  // 4
	t1 = t0.Mul(t1)   // 5
	t2 = t0.Mul(t2)   // 6
	t0 = t2.Square()  // 7
	t0 = t0.Square()  // 8
	t0 = t0.Square()  // 9
	t0 = t0.Square()  // 10
	t1 = t0.Mul(t1)   // 11
	t2 = t0.Mul(t2)   // 12
	t0 = t2.Square()  // 13
	for i := 14; i < 21; i++ {
		t0 = t0.Square() // 14 - 20
	}
	t1 = t0.Mul(t1)  // 21
	t2 = t0.Mul(t2)  // 22
	t0 = t2.Square() // 23
	for i := 24; i < 39; i++ {
		t0 = t0.Square() // 24 - 38
	}
	t1 = t0.Mul(t1) // 39
	t0 = t0.Mul(t2) // 40
	for i := 41; i < 73; i++ {
		t0 = t0.Square() // 41 - 72
	}
	t0 = t0.Mul(t1) // 73
	return t0
}

// mulNaive multiplies two GF(64) field elements
// with naive multiplications method. We keep this
// function to cross test against low level
// multiplication routines.
func (e0 GF64) mulNaive(e1 GF64) GF64 {
	var result GF64 = 0
	var shifted = e0
	var i uint64
	for i = 0; i < 64; i++ {
		if (e1)&(1<<i) != 0 {
			result = result ^ shifted
		}
		if shifted&(0x8000000000000000) != 0 {
			shifted <<= 1
			shifted ^= gf64MOD
		} else {
			shifted = shifted << 1
		}
	}
	return result
}

// mulLamire multiplies two GF(64) field elements.
func (e0 GF64) mulLamire(e1 GF64) GF64 {
	e0 = lamire(e0, e1)
	return e0
}
