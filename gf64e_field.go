// References
// Efficient Implementations of Large Finite Fields for Secure Storage Applications
// L. Luo, K. D. Bowers, A. Oprea, and L. Xu
// http://www.ccs.neu.edu/home/alina/papers/tos-gf.pdf

package gf

import (
// "fmt"
)

// GF64e is field element of GF(64) over base field GF(16).
// A field element is represented with four base field elements.
type GF64e [4]GF16

// NewGF64e returns a new element which is zero.
func NewGF64e() *GF64e {
	return &GF64e{0, 0, 0, 0}
}

// randGF64 generates a new random GF(64) element.
func randGF64e() *GF64e {
	return &GF64e{randGF16(), randGF16(), randGF16(), randGF16()}
}

// New returns new field element equals to zero
func (e *GF64e) New() *GF64e {
	return e.Zero()
}

// Set copies from input element
func (e *GF64e) Set(e2 *GF64e) {
	e[0] = e2[0]
	e[1] = e2[1]
	e[2] = e2[2]
	e[3] = e2[3]
}

// Zero returns new field element equals to zero
func (e *GF64e) Zero() *GF64e {
	return &GF64e{0, 0, 0, 0}

}

// One returns new field element equals to One
func (e *GF64e) One() *GF64e {
	return &GF64e{1, 0, 0, 0}
}

// Equal tests if given two element is equal.
func (e *GF64e) IsZero() bool {
	return e[0]|e[1]|e[2]|e[3] == 0
}

// Equal tests if given two element is equal.
func (e *GF64e) Equal(e2 *GF64e) bool {
	return e[0] == e2[0] && e[1] == e2[1] && e[2] == e2[2] && e[3] == e2[3]
}

func (e *GF64e) Mul(e2 *GF64e) {
	var c = [7]GF16{}
	for i := 0; i < 4; i++ {
		for j := 0; j < 4; j++ {
			c[i+j] ^= mul16(e[i], e2[j])
		}
	}

	// reduce with P = x^4 + x^2 + 2x + 1
	c[4] ^= c[6]
	e[3] = c[3] ^ double16(c[6]) ^ c[5]
	e[2] = c[2] ^ c[6] ^ double16(c[5]) ^ c[4]
	e[0] = c[0] ^ c[4]
	e[1] = c[1] ^ c[5] ^ double16(c[4])
}

func (e *GF64e) Add(e1 *GF64e, e2 *GF64e) {
	e[0] = e2[0] ^ e1[0]
	e[1] = e2[1] ^ e1[1]
	e[2] = e2[2] ^ e1[2]
	e[3] = e2[3] ^ e1[3]
}

func (e *GF64e) AddAssign(e2 *GF64e) {
	e[0] ^= e2[0]
	e[1] ^= e2[1]
	e[2] ^= e2[2]
	e[3] ^= e2[3]
}

func (e *GF64e) Double() {
	e[0] = double16(e[0])
	e[1] = double16(e[1])
	e[2] = double16(e[2])
	e[3] = double16(e[3])
}
