package gf

import (
	"crypto/rand"
	"encoding/binary"
)

// GF64 is field element of GF(64).
type GF64 uint64

// GF(64) under the irreducible polynomial R
// R = x^64 + x^4 + x^3 + x + 1
// R = x^64 + r
// gf64MOD is the lower part of the irreducible polynomial
// r = x^4 + x^3 + x + 1
var gf64MOD GF64 = 0x1b

// NewGF64 returns a new GF(64) element which is zero.
func NewGF64() GF64 {
	return GF64(0)
}

// randGF64 generates a new random GF(64) element.
func randGF64() GF64 {
	buf := make([]byte, 8)
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

// Inverse inverses a GF(64) element.
func (e0 GF64) Inverse() GF64 {
	if e0.IsZero() {
		return e0
	}
	// Chain is generated with tool below.
	// https://github.com/kwantam/addchain
	// Run:
	// $ ./addchain 2^64-2
	t0 := e0
	t1 := square64(t0)
	t2 := mul64(t0, t1)
	t0 = square64(t2)
	squareassign64(&t0)
	mulassign64(&t1, t0)
	mulassign64(&t2, t0)
	t0 = square64(t2)
	squareassign64(&t0)
	squareassign64(&t0)
	squareassign64(&t0)
	mulassign64(&t1, t0)
	mulassign64(&t2, t0)
	t0 = square64(t2)
	for i := 14; i < 21; i++ {
		squareassign64(&t0)
	}
	mulassign64(&t1, t0)
	mulassign64(&t2, t0)
	t0 = square64(t2)
	for i := 24; i < 39; i++ {
		squareassign64(&t0)
	}
	mulassign64(&t1, t0)
	mulassign64(&t0, t2)
	for i := 41; i < 73; i++ {
		squareassign64(&t0)
	}
	mulassign64(&t0, t1)
	return t0
}

// mulNaive multiplies two GF(16) element with shift and add method
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
