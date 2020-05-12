package gf

import (
	"errors"
	"sync"
)

// GF16Poly is to represent a polynomial.
// It supports many forms of representation.
type GF16Poly []GF16

func NewGF16Poly(coeffs []GF16) GF16Poly {
	return GF16Poly(coeffs)
}

// NewGF16PolyZero returns a new polynomial with coefficients equal to zero
func NewGF16PolyZero(len int) GF16Poly {
	return make(GF16Poly, len)
}

// randGF16Poly returns polynomial with random coefficients at given size
func randGF16Poly(n uint) GF16Poly {
	coeffs := make([]GF16, n)
	for i := n - n; i < n; i++ {
		coeffs[i] = randGF16()
	}
	return NewGF16Poly(coeffs)
}

// Degree returs degree of this polynomial.
func (p GF16Poly) Degree() int {
	return p.Length() - 1
}

// Len returns lenght of terms of this polynomial.
func (p GF16Poly) Length() int {
	return len(p)
}

// Clone duplicates the polynonial.
func (p GF16Poly) Clone() GF16Poly {
	q := make([]GF16, p.Length())
	copy(q, p)
	return NewGF16Poly(q)
}

// EqualInCoeff checks if two polynomial in coeffient form is equal.
func (p GF16Poly) EqualInCoeff(q GF16Poly) bool {
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
func (p *GF16Poly) Expand(l int) int {
	l2 := p.Length()
	if l > l2 {
		*p = append(*p, make([]GF16, l-l2)...)
		return l
	}
	return l2
}

// Eval evaluates polynomial for a given domain
// and returns a polynomial in sample form.
func (p GF16Poly) Eval(D []GF16) GF16Poly {
	evals := make(GF16Poly, len(p))
	for i := 0; i < len(D); i++ {
		evals[i] = p.EvalSingle(D[i])
	}
	return evals
}

// EvalSingle evaluates polynomial for single point.
func (p GF16Poly) EvalSingle(x GF16) GF16 {
	acc := GF16(0)
	// Apply Horner formula
	for i := p.Length() - 1; i >= 0; i-- {
		acc = mul16(acc, x) ^ p[i]
	}
	return acc
}

// Add adds two polynomial and assigns the result to the first operand
// This method assumes that degree of second operand is less or equal to
// first operand
func (p GF16Poly) Add(q GF16Poly) {
	l := len(p)
	if len(q) > l {
		return
	}
	for i := 0; i < l; i++ {
		p[i] ^= q[i]
	}
}

// mulNaive is straight forward polynomial multiplication.
// We keep this function to cross test faster FFT multiplication function.
func (p GF16Poly) mulNaive(q GF16Poly) GF16Poly {
	R := make(GF16Poly, p.Degree()+q.Degree()+1)
	for i := 0; i < p.Length(); i++ {
		for j := 0; j < q.Length(); j++ {
			R[i+j] ^= mul16(p[i], q[j])
		}
	}
	return R
}

// mulSample multiplicates two polynomial sample-wise.
// It expects that given polynomials are in evaluated form
// and evaluations are in same order. Result is assigned to
// fisrt operand
func (p GF16Poly) mulSample(q GF16Poly) error {
	if p.Length() != q.Length() {
		return errors.New("sample-wise mul is not applicable")
	}
	for i := 0; i < p.Length(); i++ {
		p[i] = mul16(p[i], q[i])
	}
	return nil
}

// One may use a global basis by inializing it with configureDefaultBasis16
var defaultBasis16 *Basis16

// Default Cantor basis generator generates
// Cantor basis with last bases equals to 1
const defaultBasisGenerator16 = GF16(0x2000)

// Basis16 is a suite for FFT and IFFT functions
type Basis16 struct {
	n               uint
	m               uint
	bases           []GF16
	combinations    []GF16
	subCombinations []GF16
}

// configureDefaultBasis16 use this to set globals
func configureDefaultBasis16(n uint) *Basis16 {
	b, err := NewBasis16(defaultBasisGenerator16, n)
	if err != nil {
		panic(err)
	}
	defaultBasis16 = b
	return b
}

// NewBasis16 given generator generates Cantor bases and
// calculates all linear combinations of the basis at desired
// span size. Last element of basis expected to equal to one otherwise,
// an error is returned.
func NewBasis16(b0 GF16, m uint) (*Basis16, error) {
	N := 16
	bases := make([]GF16, N)
	// Generates cantor bases for GF(16)
	bases[0] = b0
	for i := 0; i < N-1; i++ {
		bi := bases[i]
		bases[i+1] = square16(bi) ^ bi
	}
	// We expect last base to equal to 1
	if !bases[N-1].IsOne() {
		return nil, errors.New("last base expected to equal to 1")
	}
	// Slice full range bases down to what is need
	bases = bases[16-m:]
	l := uint(len(bases))
	combinations := make([]GF16, 1<<l)
	subCombinations := make([]GF16, 1<<(l-1))
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
	return &Basis16{
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
func (basis Basis16) fftNaive(p0 GF16Poly) {
	p1 := p0.Eval(basis.combinations)
	copy(p0[:], p1[:])
}

// FFT calculates evaluates polynomial at combinations of this bases.
func (basis Basis16) FFT(p GF16Poly) error {
	if uint(p.Length()) != basis.n {
		return errors.New("fft is not applicable")
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
		p[i+halfL] ^= p[i]
	}

	// Merging linear evaluations with base combinations
	// Step 6 in GM10 Algorithm 2.
	// var i, j, k uint
	for i := uint(1); i < basis.m; i++ {
		d := uint(1) << (basis.m - 1 - i)
		var b uint
		for j := uint(0); j < 1<<i; j++ {
			for k := b; k < b+d; k++ {
				k2 := k + d
				p[k] ^= mul16(G[j], p[k2])
				p[k2] ^= p[k]
			}
			b += (d << 1)
		}
	}
	return nil
}

// lFFT stands for lazy FFT, lFFt skips radix conversion phase
func (basis *Basis16) lFFT(p GF16Poly) error {
	if uint(p.Length()) != basis.n {
		return errors.New("fft is not applicable")
	}

	G := basis.subCombinations

	halfL := p.Length() / 2
	for i := 0; i < halfL; i++ {
		p[i+halfL] ^= p[i]
	}

	var i, j, k, d uint
	for i = 1; i < basis.m; i++ {
		d = 1 << (basis.m - 1 - i)
		var b uint
		for j = 0; j < 1<<i; j++ {
			for k = b; k < b+d; k++ {
				k2 := k + d
				p[k] ^= mul16(G[j], p[k2])
				p[k2] ^= p[k]
			}
			b += (d << 1)
		}
	}
	return nil
}

// clFFT stands for concurrent lazy FFT, clFFt skips radix conversion phase
// and implements simple concurrent processing
func (basis Basis16) clFFT(p GF16Poly) error {
	if uint(p.Length()) != basis.n {
		return errors.New("fft is not applicable")
	}
	G := basis.subCombinations
	runHalfJ := func(p GF16Poly, j0, j1, d, b uint, wg *sync.WaitGroup) {
		for j := j0; j < j1; j++ {
			for k := b; k < b+d; k++ {
				k2 := k + d
				p[k] ^= mul16(G[j], p[k2])
				p[k2] ^= p[k]
			}
			b += (d << 1)
		}
		wg.Done()
	}
	halfL := p.Length() / 2
	for i := 0; i < halfL; i++ {
		p[i+halfL] ^= p[i]
	}
	var wg sync.WaitGroup
	for i := uint(1); i < basis.m; i++ {
		d := uint(1) << (basis.m - 1 - i)
		j := uint(1) << i
		j2 := j / 2
		wg.Add(2)
		go runHalfJ(p, 0, j2, d, 0, &wg)
		go runHalfJ(p, j2, j, d, 1<<(basis.m-1), &wg)
		wg.Wait()
	}
	return nil
}

// IFFT interpolates evaluations
// computing the FFT in reverse order
func (basis Basis16) IFFT(p GF16Poly) error {
	if uint(p.Length()) != basis.n {
		return errors.New("ifft is not applicable")
	}
	G := basis.subCombinations
	var i, j, k, d uint
	for i = basis.m - 1; i > 0; i-- {
		d = 1 << (basis.m - 1 - i)
		var b uint
		for j = 0; j < 1<<i; j++ {
			for k = b; k < b+d; k++ {
				k2 := k + d
				p[k2] ^= p[k]
				p[k] ^= mul16(G[j], p[k2])
			}
			b += (d << 1)
		}
	}
	halfL := p.Length() / 2
	for i := 0; i < halfL; i++ {
		p[i+halfL] ^= p[i]
	}
	_ = basis.iRadixConversion(p)
	return nil
}

func (basis Basis16) lIFFT(p GF16Poly) error {
	if uint(p.Length()) != basis.n {
		return errors.New("ifft is not applicable")
	}
	G := basis.subCombinations
	var i, j, k, d uint
	for i = basis.m - 1; i > 0; i-- {
		d = 1 << (basis.m - 1 - i)
		var b uint
		for j = 0; j < 1<<i; j++ {
			for k = b; k < b+d; k++ {
				k2 := k + d
				p[k2] ^= p[k]
				p[k] ^= mul16(G[j], p[k2])
			}
			b += (d << 1)
		}
	}
	halfL := p.Length() / 2
	for i := 0; i < halfL; i++ {
		p[i+halfL] ^= p[i]
	}
	return nil
}

func (basis Basis16) clIFFT(p GF16Poly) error {
	if uint(p.Length()) != basis.n {
		return errors.New("ifft is not applicable")
	}
	G := basis.subCombinations
	runHalfJ := func(p GF16Poly, j0, j1, d, b uint, wg *sync.WaitGroup) {
		for j := j0; j < j1; j++ {
			for k := b; k < b+d; k++ {
				k2 := k + d
				p[k2] ^= p[k]
				p[k] ^= mul16(G[j], p[k2])
			}
			b += (d << 1)
		}
		wg.Done()
	}
	var wg sync.WaitGroup
	for i := basis.m - 1; i > 0; i-- {
		d := uint(1) << (basis.m - 1 - i)
		j := uint(1) << i
		j2 := j / 2
		wg.Add(2)
		go runHalfJ(p, 0, j2, d, 0, &wg)
		go runHalfJ(p, j2, j, d, 1<<(basis.m-1), &wg)
		wg.Wait()
	}
	halfL := p.Length() / 2
	for i := 0; i < halfL; i++ {
		p[i+halfL] ^= p[i]
	}
	return nil
}

// radixConversion is adapted from GM10 Algoritm 1,
// calculates Taylor expansion of the polinomial at (x^2 + x)
// This method assumes that input polynomial has 2^n coefficients
// Result is assigned to coefficients of the input polynomial
func (basis Basis16) radixConversion(p GF16Poly) error {
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
func (basis Basis16) iRadixConversion(p GF16Poly) error {
	if uint(p.Length()) != basis.n {
		return errors.New("radix convertions is not applicable")
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
func (basis Basis16) Mul(p0, p1 GF16Poly) error {
	if (uint(p0.Length()) != basis.n) || (uint(p0.Length()) != basis.n) {
		return errors.New("multiplication is not applicable")
	}
	_ = basis.FFT(p0)
	_ = basis.FFT(p1)
	_ = p0.mulSample(p1)
	_ = basis.IFFT(p0)
	return nil
}

func (basis Basis16) Expand(p *GF16Poly) {
	p.Expand(int(basis.n))
}
