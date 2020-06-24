package gf

import (
	"errors"
	"fmt"
	"sync"
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
func randGF64Poly(n int) GF64Poly {
	coeffs := make([]GF64, n)
	for i := 0; i < n; i++ {
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

// EqualInCoeff checks if two polynomial in coeffient form is equal.
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

// expand pads rigth zeros if given lenght is larger than this.
func (p *GF64Poly) expand(l int) int {
	l2 := p.Length()
	if l > l2 {
		*p = append(*p, make([]GF64, l-l2)...)
		return l
	}
	return l2
}

// Eval evaluates polynomial for a given domain
// and returns a polynomial in sample form.
func (p GF64Poly) Eval(n int, D []GF64) GF64Poly {
	evals := make(GF64Poly, len(p))
	for i := 0; i < len(D); i++ {
		evals[i] = p.EvalSingle(n, D[i])
	}
	return evals
}

// EvalSingle evaluates polynomial for single point.
// func (p GF64Poly) EvalSingle(x GF64) GF64 {
func (p GF64Poly) EvalSingle(n int, x GF64) GF64 {
	if n < p.Length() {
		n = p.Length()
	}
	acc := NewGF64()
	// Apply Horner formula
	for i := n - 1; i >= 0; i-- {
		// acc = acc.Mul(x).Add(p[i])
		acc = mul64(acc, x) ^ p[i]
	}
	return acc
}

// Add adds two polynomial and assigns the result to the first operand
// This method assumes that degree of second operand is less or equal to
// first operand
func (p GF64Poly) Add(q GF64Poly) {
	l := len(p)
	if len(q) > l {
		return
	}
	for i := 0; i < l; i++ {
		p[i] ^= q[i]
	}
}

// mul implements straight forward polynomial multiplication.
// We keep this function to cross test faster FFT multiplication function.
func (p GF64Poly) mul(q GF64Poly) GF64Poly {
	R := make(GF64Poly, p.Degree()+q.Degree()+1)
	for i := 0; i < p.Length(); i++ {
		for j := 0; j < q.Length(); j++ {
			R[i+j] ^= mul64(p[i], q[j])
		}
	}
	return R
}

func (p GF64Poly) mulSample(n int, q GF64Poly) error {
	if n < p.Length() || n < q.Length() {
		return errors.New("sample-wise mul is not applicable")
	}
	for i := 0; i < n; i++ {
		mulassign64(&p[i], q[i])
	}
	return nil
}

func (p GF64Poly) debug(desc string) {
	fmt.Println(desc, len(p))
	for i := 0; i < len(p); i++ {
		fmt.Printf("%#16.16x\n", p[i])
	}
}

// One may use a global basis by inializing it with configureDefaultBasis64
var defaultBasis64 *Basis64

// Default Cantor basis generator generates
// Cantor basis with last bases equals to 1
const defaultBasisGenerator64 = GF64(0xce41e2bee6cbe964)

// Basis64 is the suite for FFT and IFFT functions
type Basis64 struct {
	n               int
	m               int
	combinations    []GF64
	subCombinations []GF64
}

// configureDefaultBasis64 use this to set globals
func configureDefaultBasis64(n int) *Basis64 {
	b, err := NewBasis64(defaultBasisGenerator64, n)
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
func NewBasis64(b0 GF64, m int) (*Basis64, error) {
	N := 64
	bases := make([]GF64, N)
	// Generates cantor bases for GF(64)
	bases[0] = b0
	for i := 0; i < N-1; i++ {
		bases[i+1] = mul64(bases[i], bases[i]) ^ bases[i]
	}
	// We expect last base to equal to 1
	if !bases[N-1].IsOne() {
		return nil, errors.New("last base expected to equal to 1")
	}
	// Slice full range bases down to what is need
	bases = bases[64-m:]
	l := len(bases)
	combinations := make([]GF64, 1<<l)
	subCombinations := make([]GF64, 1<<(l-1))
	// Calcucate combinations
	for i := 0; i < l; i++ {
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
			combinations:    combinations,
			subCombinations: subCombinations,
		},
		nil
}

// fftNaive calculates evaluates polynomial at combinations of this bases.
// We keep this method here to cross test faster FFT function.
func (basis Basis64) fftNaive(m int, p0 GF64Poly) {
	p1 := p0.Eval(1<<m, basis.combinations)
	copy(p0[:], p1[:])
}

func (basis *Basis64) fft(m int, p GF64Poly) error {
	if uint(p.Length()) > 1<<m {
		return errors.New("fft is not applicable")
	}
	_ = basis.radixConversion(m, p)

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
	halfL := 1 << (m - 1)
	for i := 0; i < halfL; i++ {
		p[i+halfL] ^= p[i]
	}

	// Merging linear evaluations with base combinations
	// Step 6 in GM10 Algorithm 2.
	for i := 1; i < m; i++ {
		d := 1 << (m - 1 - i)
		var b int
		for j := 0; j < 1<<i; j++ {
			for k := b; k < b+d; k++ {
				// butterfly(&p[k], &p[k+d], G[j])
				k2 := k + d
				p[k] ^= mul64(G[j], p[k2])
				p[k2] ^= p[k]
			}
			b += (d << 1)
		}
	}
	return nil
}

// lFFT stands for lazy FFT, lFFt skips radix conversion phase
func (basis *Basis64) lFFT(p GF64Poly) error {
	if p.Length() != basis.n {
		return errors.New("fft is not applicable")
	}
	G := basis.subCombinations
	halfL := p.Length() / 2
	for i := 0; i < halfL; i++ {
		p[i+halfL] ^= p[i]
	}
	for i := 1; i < basis.m; i++ {
		d := 1 << (basis.m - 1 - i)
		var b int
		for j := 0; j < 1<<i; j++ {
			for k := b; k < b+d; k++ {
				butterfly(&p[k], &p[k+d], G[j])
			}
			b += (d << 1)
		}
	}
	return nil
}

// clFFT stands for concurrent lazy FFT, clFFt skips radix conversion phase
// and implements simple concurrent processing
func (basis *Basis64) clFFT(p GF64Poly) error {
	if p.Length() != basis.n {
		return errors.New("fft is not applicable")
	}
	G := basis.subCombinations
	runHalfJ := func(p GF64Poly, j0, j1, d, b int, wg *sync.WaitGroup) {
		for j := j0; j < j1; j++ {
			for k := b; k < b+d; k++ {
				butterfly(&p[k], &p[k+d], G[j])
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
	for i := 1; i < basis.m; i++ {
		d := 1 << (basis.m - 1 - i)
		j := 1 << i
		j2 := j / 2
		wg.Add(2)
		go runHalfJ(p, 0, j2, d, 0, &wg)
		go runHalfJ(p, j2, j, d, 1<<(basis.m-1), &wg)
		wg.Wait()
	}
	return nil
}

func (basis *Basis64) ifft(m int, p GF64Poly) error {
	if p.Length() > 1<<m {
		return errors.New("fft is not applicable")
	}
	G := basis.subCombinations
	var i, j, k, d int
	for i = m - 1; i > 0; i-- {
		d = 1 << (m - 1 - i)
		var b int
		for j = 0; j < 1<<i; j++ {
			for k = b; k < b+d; k++ {
				ibutterfly(&p[k], &p[k+d], G[j])
			}
			b += (d << 1)
		}
	}
	halfL := 1 << (m - 1)
	for i := 0; i < halfL; i++ {
		p[i+halfL] ^= p[i]
	}
	_ = basis.iRadixConversion(m, p)
	return nil
}

func (basis *Basis64) lIFFT(p GF64Poly) error {
	if p.Length() != basis.n {
		return errors.New("ifft is not applicable")
	}
	G := basis.subCombinations
	var i, j, k, d int
	for i = basis.m - 1; i > 0; i-- {
		d = 1 << (basis.m - 1 - i)
		var b int
		for j = 0; j < 1<<i; j++ {
			for k = b; k < b+d; k++ {
				ibutterfly(&p[k], &p[k+d], G[j])
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

func (basis *Basis64) clIFFT(p GF64Poly) error {
	if p.Length() != basis.n {
		return errors.New("ifft is not applicable")
	}
	G := basis.subCombinations
	runHalfJ := func(p GF64Poly, j0, j1, d, b int, wg *sync.WaitGroup) {
		for j := j0; j < j1; j++ {
			for k := b; k < b+d; k++ {
				ibutterfly(&p[k], &p[k+d], G[j])
			}
			b += (d << 1)
		}
		wg.Done()
	}
	var wg sync.WaitGroup
	for i := basis.m - 1; i > 0; i-- {
		d := 1 << (basis.m - 1 - i)
		j := 1 << i
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
func (basis *Basis64) radixConversion(m int, p GF64Poly) error {
	n := 1 << m
	if p.Length() < n {
		return errors.New("radix convertions is not applicable")
	}
	for r := 0; r < m-1; r++ {
		for i := 0; i < m-r-1; i++ {
			d := n >> i
			d2, d4 := d>>1, d>>2
			var off int
			for j := 0; j < 1<<i; j++ {
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
func (basis Basis64) iRadixConversion(m int, p GF64Poly) error {
	n := 1 << m
	if p.Length() < n {
		return errors.New("radix convertions is not applicable")
	}
	for r := m - 2; r >= 0; r-- {
		for i := m - r - 2; i >= 0; i-- {
			d := n >> i
			d2, d4 := d>>1, d>>2
			var off int
			for j := 0; j < 1<<i; j++ {
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

func (basis *Basis64) mul(m int, p0, p1 *GF64Poly) error {
	n := 1 << m
	if p0.Length() < n || p1.Length() < n {
		return errors.New("multiplication is not applicable")
	}
	// if n < 4 {
	// 	p0.mul(p1)
	// 	return nil
	// }
	// apply expand here

	fmt.Printf("%p\n", &p0)
	p0.expand(1 << (m + 1))
	p1.expand(1 << (m + 1))
	fmt.Printf("%p\n", &p0)

	// fmt.Println(len(p0))
	// fmt.Println(len(p1))
	// if err := basis.fft(m+1, p0); err != nil {
	// 	return err
	// }
	// if err := basis.fft(m+1, p1); err != nil {
	// 	return err
	// }
	// if err := p0.mulSample(1<<(m+1), p1); err != nil {
	// 	return err
	// }
	// if err := basis.ifft(m+1, p0); err != nil {
	// 	return err
	// }

	return nil
}

func (basis *Basis64) z(roots []GF64) (GF64Poly, error) {
	l := len(roots)
	p := make([]GF64Poly, l)
	if len(roots)&1 != 0 {
		return nil, errors.New("expect even number of roots for simplicity")
	}
	for i := 0; i < l; i++ {
		p[i] = NewGF64Poly([]GF64{roots[i], GF64(1)})
	}
	for i := 0; i < l>>2; i++ {
		fmt.Println("xxx", i)
		p[i].mul(p[i+1])

	}

	return nil, nil
}
