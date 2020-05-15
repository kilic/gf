// References
// A Tutorial on Reed-Solomon Coding for Fault-Tolerance in RAID-like Systems
// Plank
// http://web.eecs.utk.edu/~jplank/plank/papers/CS-96-332.pdf

package gf

import (
	"crypto/rand"
	"encoding/binary"
)

// GF16 is field element of GF(16).
type GF16 uint16

// GF(16) under the irreducible polynomial R
// R = x^16 + x^12 + x^3 + x^1 + 1
// R = x^16 + r
// gf64MOD is the lower part of the irreducible polynomial
// r = x^12 + x^3 + x^1 + 1
var gf16MOD GF16 = 0x100b

// gf16Generator generates GF(16) we use this element to generate log/antilog tables
var gf16Generator = GF16(3)

// log and antilog table of GF16
var logTable16, antiLogTable16 = func() ([65536]GF16, [65536 * 2]GF16) {
	logTable2 := [65536]GF16{}
	antiLogTable2 := [65536 * 2]GF16{}
	acc := GF16(1)
	for i := 0; i < 0xffff; i++ {
		g := acc
		logTable2[int(g)] = GF16(i)
		antiLogTable2[i] = g
		antiLogTable2[i+0xffff] = g
		acc = acc.mulNaive(gf16Generator)
	}
	return logTable2, antiLogTable2
}()

func mul16(a, b GF16) GF16 {
	if a == 0 || b == 0 {
		return 0
	}
	return antiLogTable16[int32(logTable16[a])+int32(logTable16[b])]
}

func double16(a GF16) GF16 {
	if a&(0x8000) != 0 { // must be reduced
		a <<= 1
		a ^= gf16MOD
	} else {
		a <<= 1
	}
	return a
}

func square16(a GF16) GF16 {
	if a == 0 {
		return 0
	}
	return antiLogTable16[int(logTable16[a])<<1]
}

// randGF16 generates a new random GF(16) element.
func randGF16() GF16 {
	buf := make([]byte, 2)
	var e GF16
	for e.IsZero() {
		_, err := rand.Read(buf)
		if err != nil {
			panic(err)
		}
		e = GF16(binary.BigEndian.Uint16(buf))
	}
	return e
}

// Zero returns new field element equals to Zero
func (e GF16) Zero() GF16 {
	return GF16(0)
}

// One returns new field element equals to One
func (e GF16) One() GF16 {
	return GF16(1)
}

// IsZero returns true if GF(16) element is equal to zero
func (e GF16) IsZero() bool {
	return e.Equal(GF16(0))
}

// IsZero returns true if GF(16) element is equal to zero
func (e GF16) IsOne() bool {
	return e.Equal(GF16(1))
}

// Equal tests if given two element is equal.
func (e0 GF16) Equal(e1 GF16) bool {
	return e0 == e1
}

// mulNaive multiplies two GF(16) element with shift and add method
func (e0 GF16) mulNaive(e1 GF16) GF16 {
	var result GF16 = 0
	var shifted = e0
	var i uint64
	for i = 0; i < 16; i++ {
		if (e1)&(1<<i) != 0 {
			result ^= shifted
		}
		if shifted&(0x8000) != 0 { // will overflow
			shifted <<= 1
			shifted ^= gf16MOD
		} else {
			shifted <<= 1
		}
	}
	return result
}
