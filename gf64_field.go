package gf

import (
	"crypto/rand"
	"encoding/binary"
	"fmt"
)

// GF(64) under the irreducible polynomial R
// R = x^64 + x^4 + x^3 + x + 1
// R = x^64 + r
// gf64MOD is the lower part of the irreducible polynomial
// r = x^4 + x^3 + x + 1
var gf64MOD uint64 = 0x1b

// randGF64 generates a new random GF(64) element.
func randGF64() uint64 {
	buf := make([]byte, 8)
	var e uint64
	for e == 0 {
		_, err := rand.Read(buf)
		if err != nil {
			panic(err)
		}
		e = binary.BigEndian.Uint64(buf)
	}
	return e
}

func toHex(e uint64) string {
	return fmt.Sprintf("%#16.16x", e)
}

func zero() uint64 {
	return 0
}

func one() uint64 {
	return 1
}

func inverse(e0 uint64) uint64 {
	if e0 == 0 {
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

func exp(a uint64, e uint64) uint64 {
	if e == 0 {
		return 1
	}
	l := log2Ceil(int(e))
	var acc = a
	var r uint64 = 1
	for i := 0; i < l+1; i++ {
		if (e>>i)&1 == 1 {
			mulassign64(&r, acc)
		}
		squareassign64(&acc)
	}
	return r
}

// mulNaive multiplies two GF(16) element with shift and add method
func mulNaive(e0, e1 uint64) uint64 {
	var result uint64 = 0
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
