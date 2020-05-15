package gf

import (
	"errors"
)

// GF64ePoly is to represent a polynomial.
// It supports many forms of representation.
type GF64ePoly []GF64e

func NewGF64ePoly(coeffs []GF64e) GF64ePoly {
	return GF64ePoly(coeffs)
}

// NewGF64ePolyZero returns a new polynomial with coefficients equal to zero
func NewGF64ePolyZero(len int) GF64ePoly {
	return make(GF64ePoly, len)
}

// randGF64ePoly returns polynomial with random coefficients at given size
func randGF64ePoly(n uint) GF64ePoly {
	coeffs := make([]GF64e, n)
	for i := n - n; i < n; i++ {
		coeffs[i] = *randGF64e()
	}
	return NewGF64ePoly(coeffs)
}

// Degree returs degree of this polynomial.
func (p GF64ePoly) Degree() int {
	return p.Length() - 1
}

// Len returns lenght of terms of this polynomial.
func (p GF64ePoly) Length() int {
	return len(p)
}

// Clone duplicates the polynonial.
func (p GF64ePoly) Clone() GF64ePoly {
	q := make([]GF64e, p.Length())
	copy(q, p)
	return NewGF64ePoly(q)
}

// EqualInCoeff checks if two polynomial in coeffient form is equal.
func (p GF64ePoly) EqualInCoeff(q GF64ePoly) bool {
	lp := p.Length()
	lq := q.Length()
	if lp > lq {
		for i := 0; i < lq; i++ {
			if !p[i].Equal(&q[i]) {
				return false
			}
		}
		for i := lp - 1; i < lp-lq; i++ {
			if !p[i].IsZero() {
				return false
			}
		}
	} else {
		for i := 0; i < lp; i++ {
			for i := 0; i < lp; i++ {
				if !p[i].Equal(&q[i]) {
					return false
				}
			}
			for i := lq - 1; i < lp-lq; i++ {
				if !q[i].IsZero() {
					return false
				}
			}
		}
	}
	return true
}

// Eval evaluates polynomial for a given domain
// and returns a polynomial in sample form.
func (p GF64ePoly) Eval(D []GF64e) GF64ePoly {
	evals := make(GF64ePoly, len(p))
	for i := 0; i < len(D); i++ {
		evals[i] = *p.EvalSingle(&D[i])
	}
	return evals
}

// EvalSingle evaluates polynomial for single point.
func (p GF64ePoly) EvalSingle(x *GF64e) *GF64e {
	acc := NewGF64e()
	// Apply Horner formula
	for i := p.Length() - 1; i >= 0; i-- {
		// acc = acc.Mul(x).Add(p[i])
		acc.MulAssign(x)
		acc.AddAssign(&p[i])
		// acc = mul64(acc, x) ^ p[i]
	}
	return acc
}

// Add adds two polynomial and assigns the result to the first operand
// This method assumes that degree of second operand is less or equal to
// first operand
func (p GF64ePoly) Add(q GF64ePoly) {
	l := len(p)
	if len(q) > l {
		return
	}
	for i := 0; i < l; i++ {
		p[i].AddAssign(&q[i])
	}
}

// mulNaive is straight forward polynomial multiplication.
// We keep this function to cross test faster FFT multiplication function.
func (p GF64ePoly) mulNaive(q GF64ePoly) GF64ePoly {
	R := make(GF64ePoly, p.Degree()+q.Degree()+1)
	for i := 0; i < p.Length(); i++ {
		for j := 0; j < q.Length(); j++ {
			R[i+j].Mul(&p[i], &q[j])
		}
	}
	return R
}

var defaultBasis64e *Basis64e

// Default Cantor basis generator generates
// Cantor basis with last bases equals to 1
var defaultBasisGenerator64e = &GF64e{0xfa59, 0xc44a, 0xa134, 0x350e}

// Basis64 is the suite for FFT and IFFT functions
type Basis64e struct {
	n               uint
	m               uint
	bases           []GF64e
	combinations    []GF64e
	subCombinations []GF64e
}

// configureDefaultBasis64 use this to set globals
func configureDefaultBasis64e(n uint) *Basis64e {
	b, err := NewBasis64e(defaultBasisGenerator64e, n)
	if err != nil {
		panic(err)
	}
	defaultBasis64e = b
	return b
}

// NewBasis64 given generator generates Cantor bases and
// calculates all linear combinations of the basis at desired
// span size. Last element of basis expected to equal to one otherwise,
// an error is returned.
func NewBasis64e(b0 *GF64e, m uint) (*Basis64e, error) {
	N := 64
	bases := make([]GF64e, N)
	// Generates cantor bases for GF(64)
	bases[0].Set(b0)
	for i := 0; i < N-1; i++ {
		bases[i+1].Mul(&bases[i], &bases[i])
		bases[i+1].AddAssign(&bases[i])
	}
	// We expect last base to equal to 1
	if !bases[N-1].IsOne() {
		return nil, errors.New("last base expected to equal to 1")
	}
	// Slice full range bases down to what is need
	bases = bases[64-m:]
	l := uint(len(bases))
	combinations := make([]GF64e, 1<<l)
	subCombinations := make([]GF64e, 1<<(l-1))
	// Calcucate combinations
	for i := uint(0); i < l; i++ {
		a := (1 << i)
		for j := 0; j < a; j++ {
			combinations[a+j].Add(&combinations[j], &bases[l-1-i])
		}
	}
	// Collect even indexed sums as subcombinations
	for i := 0; i < len(combinations)/2; i++ {
		subCombinations[i] = combinations[2*i]
	}
	return &Basis64e{
			n:               1 << m,
			m:               m,
			bases:           bases,
			combinations:    combinations,
			subCombinations: subCombinations,
		},
		nil
}

func (basis Basis64e) fftNaive(p0 GF64ePoly) {
	p1 := p0.Eval(basis.combinations)
	copy(p0[:], p1[:])
}

// FFT calculates evaluates polynomial at combinations of this bases.
func (basis Basis64e) FFT(p GF64ePoly) error {
	if uint(p.Length()) != basis.n {
		return errors.New("fft is not applicable")
	}
	_ = basis.radixConversion(p)

	G := basis.subCombinations

	halfL := p.Length() / 2
	for i := 0; i < halfL; i++ {
		p[i+halfL].AddAssign(&p[i])
	}

	// Merging linear evaluations with base combinations
	// Step 6 in GM10 Algorithm 2.
	for i := uint(1); i < basis.m; i++ {
		d := uint(1) << (basis.m - 1 - i)
		var b uint
		for j := uint(0); j < 1<<i; j++ {
			for k := b; k < b+d; k++ {
				k2 := k + d
				p[k].Mul(&p[k2], &G[j])
				p[k2].AddAssign(&p[k])
			}
			b += (d << 1)
		}
	}
	return nil
}

func (basis Basis64e) radixConversion(p GF64ePoly) error {
	if uint(p.Length()) != basis.n {
		return errors.New("radix convertions is not applicable")
	}
	for r := uint(0); r < basis.m-1; r++ {
		for i := uint(0); i < basis.m-r-1; i++ {
			d := basis.n >> i
			d2, d4 := d>>1, d>>2
			var off uint
			for j := uint(0); j < 1<<i; j++ {
				for k := off; k < off+d4; k++ {
					p[d2+k].AddAssign(&p[d2+d4+k])
					p[d4+k].AddAssign(&p[d2+k])
				}
				off += d
			}
		}
	}
	return nil
}
