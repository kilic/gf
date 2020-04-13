package gf

import (
	"fmt"
)

// GF64Poly is to represent a polynomial.
// It supports many forms of representation.
type GF64Poly []GF64

func NewGF64Poly(coeffs []GF64) GF64Poly {
	return GF64Poly(coeffs)
}

// NewGF64PolyZero returns a new polynomial with coefficients equal to zero
func NewGF64PolyZero(len int) GF64Poly {
	return make(GF64Poly, len)
}

// randGF64Poly returns polynomial with random coefficients at given size
func randGF64Poly(n uint) GF64Poly {
	coeffs := make([]GF64, n)
	for i := n - n; i < n; i++ {
		coeffs[i] = randGF64()
	}
	return NewGF64Poly(coeffs)
}

// Degree returs degree of this polynomial.
func (p GF64Poly) Degree() int {
	return p.Length() - 1
}

// Len returns lenght of terms of this polynomial.
func (p GF64Poly) Length() int {
	return len(p)
}

// Clone duplicates the polynonial.
func (p GF64Poly) Clone() GF64Poly {
	q := make([]GF64, p.Length())
	copy(q, p)
	return NewGF64Poly(q)
}

// EqualInC checks if two polynomial in coeffient form is equal.
func (p GF64Poly) EqualInCoeff(q GF64Poly) bool {
	lp := p.Length()
	lq := q.Length()
	if lp > lq {
		for i := 0; i < lq; i++ {
			if !p[i].Equal(q[i]) {
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
				if !p[i].Equal(q[i]) {
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

// Expand pads rigth zeros if given lenght is larger than this.
func (p *GF64Poly) Expand(l int) int {
	l2 := p.Length()
	if l > l2 {
		*p = append(*p, make([]GF64, l-l2)...)
		return l
	}
	return l2
}

// Eval evaluates polynomial for a given domain
// and returns a polynomial in sample form.
func (p GF64Poly) Eval(D []GF64) GF64Poly {
	evals := make(GF64Poly, len(p))
	for i := 0; i < len(D); i++ {
		evals[i] = p.EvalSingle(D[i])
	}
	return evals
}

// EvalSingle evaluates polynomial for single point.
func (p GF64Poly) EvalSingle(x GF64) GF64 {
	acc := NewGF64()
	// Apply Horner formula
	for i := p.Length() - 1; i >= 0; i-- {
		acc = acc.Mul(x).Add(p[i])
	}
	return acc
}

// Add adds two polynomial and assigns the result to the first operand
// This method assumes that degree of second operand is less or equal to
// first operand
func (p GF64Poly) Add(q GF64Poly) {
	for i := 0; i < len(q); i++ {

	}
}

// mulNaive is straight forward polynomial multiplication.
// We keep this function to cross test faster FFT multiplication function.
func (p GF64Poly) mulNaive(q GF64Poly) GF64Poly {
	R := make(GF64Poly, p.Degree()+q.Degree()+1)
	for i := 0; i < p.Length(); i++ {
		for j := 0; j < q.Length(); j++ {
			r := p[i].Mul(q[j])
			R[i+j] = R[i+j].Add(r)
		}
	}
	return R
}

// mulSample multiplicates two polynomial sample-wise.
// It expects that given polynomials are in evaluated form
// and evaluations are in same order. Result is assigned to
// fisrt operand
func (p GF64Poly) mulSample(q GF64Poly) error {
	if p.Length() != q.Length() {
		return fmt.Errorf("sample-wise mul is not applicable")
	}
	for i := 0; i < p.Length(); i++ {
		p[i] = p[i].Mul(q[i])
	}
	return nil
}

// One may use a global basis by inializing it with configureDefaultBasis64
var defaultBasis64 *Basis64

// Default Cantor basis generator generates
// Cantor basis with last bases equals to 1
const defaultBasisGenerator = GF64(0xce41e2bee6cbe964)

// Basis64 is a suite for FFT and IFFT functions
type Basis64 struct {
	n               uint
	m               uint
	bases           []GF64
	combinations    []GF64
	subCombinations []GF64
}

// configureDefaultBasis64 use this to set globals
func configureDefaultBasis64(u GF64, n uint) *Basis64 {
	b, err := NewBasis64(u, n)
	if err != nil {
		panic(err)
	}
	defaultBasis64 = b
	return b
}

// NewBasis64 given generator generates Cantor bases and
// calculates all linear combinations of the basis at desired
// span size. Last element of basis expected to equal to one otherwise,
// an error is returned.
func NewBasis64(b0 GF64, m uint) (*Basis64, error) {
	N := 64
	bases := make([]GF64, N)
	// Generates cantor bases for GF(64)
	bases[0] = b0
	for i := 0; i < N-1; i++ {
		bi := bases[i]
		bases[i+1] = bi.Square().Add(bi)
	}
	// We expect last base to equal to 1
	if !bases[N-1].IsOne() {
		return nil, fmt.Errorf("last base expected to equal to 1")
	}
	// Slice full range bases down to what is need
	bases = bases[64-m:]
	l := uint(len(bases))
	combinations := make([]GF64, 1<<l)
	subCombinations := make([]GF64, 1<<(l-1))
	// Calcucate combinations
	for i := uint(0); i < l; i++ {
		a := (1 << i)
		for j := 0; j < a; j++ {
			combinations[a+j] = combinations[j] ^ bases[l-1-i]
		}
	}
	// Collect even indexed sums as subcombinations
	for i := 0; i < len(combinations)/2; i++ {
		subCombinations[i] = combinations[2*i]
	}
	return &Basis64{
			n:               1 << m,
			m:               m,
			bases:           bases,
			combinations:    combinations,
			subCombinations: subCombinations,
		},
		nil
}

// fftNaive calculates evaluates polynomial at combinations of this bases.
// We keep this method here to cross test faster FFT function.
func (basis Basis64) fftNaive(p0 GF64Poly) {
	p1 := p0.Eval(basis.combinations)
	copy(p0[:], p1[:])
}

// FFT calculates evaluates polynomial at combinations of this bases.
func (basis Basis64) FFT(p GF64Poly) error {
	if uint(p.Length()) != basis.n {
		return fmt.Errorf("fft is not applicable, resize polynomial to %d coefficients", basis.n)
	}
	_ = basis.radixConversion(p)

	// We will use subsets of basis combinations
	// throughout the iterations.
	// Since basis with the last element 1 is only accepted,
	// we don't need to calculate twisting operations which are
	// Step 2 and step 4 in GM10 Algorithm 2.
	G := basis.subCombinations

	// Step 1 in GM10 Algorithm 2.
	// Linear evaluation at the leafs of the recursions.
	// f(0), f(β_1) where f is linear polynomial and β_1 = β_m = 1.
	// f(0) = a_0 + 0 * a_1 = a_0
	// f(β_1) = a_0 + a_1 * β_1 = a_0 + a_1
	halfL := p.Length() / 2
	for i := 0; i < halfL; i++ {
		p[i+halfL] = p[i].Add(p[i+halfL])
	}

	// Merging linear evaluations with base combinations
	// Step 6 in GM10 Algorithm 2.
	var i, j, k, d uint
	for i = 1; i < basis.m; i++ {
		d = 1 << (basis.m - 1 - i)
		var b uint
		for j = 0; j < 1<<i; j++ {
			for k = b; k < b+d; k++ {
				k2 := k + d
				p[k] = p[k].Add(G[j].Mul(p[k2]))
				p[k2] = p[k2].Add(p[k])
			}
			b += (d << 1)
		}
	}
	return nil
}

// IFFT interpolates evaluations
// computing the FFT in reverse order
func (basis Basis64) IFFT(p GF64Poly) error {
	if uint(p.Length()) != basis.n {
		return fmt.Errorf("ifft is not applicable, resize polynomial to %d coefficients", basis.n)
	}
	G := basis.subCombinations
	var i, j, k, d uint
	for i = basis.m - 1; i > 0; i-- {
		d = 1 << (basis.m - 1 - i)
		var b uint
		for j = 0; j < 1<<i; j++ {
			for k = b; k < b+d; k++ {
				k2 := k + d
				p[k2] = p[k2].Add(p[k])
				p[k] = p[k].Add(G[j].Mul(p[k2]))
			}
			b += (d << 1)
		}
	}
	halfL := p.Length() / 2
	for i := 0; i < halfL; i++ {
		p[i+halfL] = p[i].Add(p[i+halfL])
	}
	_ = basis.iRadixConversion(p)
	return nil
}

// radixConversion is adapted from GM10 Algoritm 1,
// calculates Taylor expansion of the polinomial at (x^2 + x)
// This method assumes that input polynomial has 2^n coefficients
// Result is assigned to coefficients of the input polynomial
func (basis Basis64) radixConversion(p GF64Poly) error {
	if uint(p.Length()) != basis.n {
		return fmt.Errorf("radix convertions is not applicable, resize polynomial to %d coefficients", basis.n)
	}
	for r := uint(0); r < basis.m-1; r++ {
		for i := uint(0); i < basis.m-r-1; i++ {
			d := basis.n >> i
			d2, d4 := d>>1, d>>2
			var off uint
			for j := uint(0); j < 1<<i; j++ {
				for k := off; k < off+d4; k++ {
					p[d2+k] ^= p[d2+d4+k]
					p[d4+k] ^= p[d2+k]
				}
				off += d
			}
		}
	}
	return nil
}

// iRadixConversion given polynomial in Taylor expanded form
// calculates coefficients. It is basically calculating Taylor expantion at (x^2 + x)
// in reverse order.
func (basis Basis64) iRadixConversion(p GF64Poly) error {
	if uint(p.Length()) != basis.n {
		return fmt.Errorf("radix convertions is not applicable, resize polynomial to %d coefficients", basis.n)
	}
	for r := basis.m - 2; r != 0xffffffffffffffff; r-- {
		for i := basis.m - r - 2; i != 0xffffffffffffffff; i-- {
			d := basis.n >> i
			d2, d4 := d>>1, d>>2
			var off uint
			for j := uint(0); j < 1<<i; j++ {
				for k := off; k < off+d4; k++ {
					p[d4+k] ^= p[d2+k]
					p[d2+k] ^= p[d2+d4+k]
				}
				off += d
			}
		}
	}
	return nil
}

// Mul multiplicatates two polynomial and assigns the result to first operand.
func (basis Basis64) Mul(p0, p1 GF64Poly) error {
	if (uint(p0.Length()) != basis.n) || (uint(p0.Length()) != basis.n) {
		return fmt.Errorf("multiplication is not applicable, resize polynomials to %d coefficients", basis.n)
	}
	_ = basis.FFT(p0)
	_ = basis.FFT(p1)
	_ = p0.mulSample(p1)
	_ = basis.IFFT(p0)
	return nil
}

func (basis Basis64) Expand(p *GF64Poly) {
	p.Expand(int(basis.n))
}
