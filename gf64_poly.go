package gf

import (
	"errors"
	"fmt"
)

type basis struct {
	n               int
	m               int
	bases           []uint64
	combinations    []uint64
	subCombinations []uint64
}

// defaultBasisGenerator generates Cantor basis with last bases equals to 1
const defaultBasisGenerator uint64 = 0xce41e2bee6cbe964

var defaultBasis *basis

func initDefaultBasis(m int) {
	var err error
	defaultBasis, err = newBasis(defaultBasisGenerator, m)
	if err != nil {
		panic(err)
	}
}

// newBasis given generator generates Cantor bases and
// calculates all linear combinations of the basis at desired
// span size. Last element of basis expected to equal to one otherwise,
// an error is returned.
func newBasis(b0 uint64, m int) (*basis, error) {
	N := 64
	bases := make([]uint64, N)
	// Generates cantor bases for GF(64)
	bases[0] = b0
	for i := 0; i < N-1; i++ {
		bases[i+1] = mul64(bases[i], bases[i]) ^ bases[i]
	}
	// We expect last base to equal to 1
	if bases[N-1] != 1 {
		return nil, errors.New("last base expected to equal to 1")
	}
	// Slice full range bases down to what is need
	bases = bases[64-m:]
	l := len(bases)
	combinations := make([]uint64, 1<<l)
	subCombinations := make([]uint64, 1<<(l-1))
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
	return &basis{
			n:               1 << m,
			m:               m,
			bases:           bases,
			combinations:    combinations,
			subCombinations: subCombinations,
		},
		nil
}

// Poly is to represent a polynomial.
// It supports coefficient and evaluation forms of representation.
type poly struct {
	a []uint64
}

func newPoly(coeffs []uint64) *poly {
	return &poly{coeffs}
}

func newEmptyPoly(len int) *poly {
	return &poly{make([]uint64, len)}
}

func randPoly(n int) *poly {
	coeffs := make([]uint64, n)
	for i := 0; i < n; i++ {
		coeffs[i] = randGF64()
	}
	return newPoly(coeffs)
}

func (p *poly) degree() int {
	return p.length() - 1
}

func (p *poly) length() int {
	return len(p.a)
}

func (p *poly) m() int {
	return log2Ceil(p.length())
}

func (p *poly) clone() *poly {
	q := make([]uint64, p.length())
	copy(q, p.a)
	return newPoly(q)
}

func (p *poly) equalInCoeff(q *poly) bool {
	lp := p.length()
	lq := q.length()
	if lp > lq {
		for i := 0; i < lq; i++ {
			if p.a[i] != q.a[i] {
				return false
			}
		}
		for i := lp - 1; i < lp-lq; i++ {
			if p.a[i] != 0 {
				return false
			}
		}
	} else {
		for i := 0; i < lp; i++ {
			for i := 0; i < lp; i++ {
				if p.a[i] != q.a[i] {
					return false
				}
			}
			for i := lq - 1; i < lp-lq; i++ {
				if q.a[i] != 0 {
					return false
				}
			}
		}
	}
	return true
}

func (p *poly) eval(D []uint64) *poly {
	evals := make([]uint64, len(p.a))
	for i := 0; i < len(D); i++ {
		evals[i] = p.evalSingle(D[i])
	}
	return &poly{evals}
}

func (p *poly) evalSingle(x uint64) uint64 {
	var acc uint64
	for i := p.length() - 1; i >= 0; i-- {
		acc = mul64(acc, x) ^ p.a[i]
	}
	return acc
}

func (p *poly) expand(l int) int {
	l2 := p.length()
	if l > l2 {
		p.a = append(p.a, make([]uint64, l-l2)...)
		return l
	}
	return l2
}

func (p *poly) trimZeros() {
	l := 0
	for i := p.length() - 1; i >= 0; i-- {
		if p.a[i] == 0 {
			l++
		} else {
			break
		}
	}
	p.a = p.a[:p.length()-l]
}

func (p *poly) substitute(k uint64) *poly {
	n := p.length()
	acc := k
	for i := 1; i < n; i++ {
		mulassign64(&p.a[i], acc)
		mulassign64(&acc, k)
	}
	return p
}

func (p *poly) add(q *poly) {
	l := len(p.a)
	if l > len(q.a) {
		l = len(q.a)
	}
	for i := 0; i < l; i++ {
		p.a[i] ^= q.a[i]
	}
}

func (p *poly) radixConversion() error {
	m := p.m()
	n := p.length()
	if n != 1<<m {
		return fmt.Errorf("radix conversion expects polynomial length as power of two: %d, %d", m, n)
	}
	for r := 0; r < m-1; r++ {
		for i := 0; i < m-r-1; i++ {
			d := n >> i
			d2, d4 := d>>1, d>>2
			var off int
			for j := 0; j < 1<<i; j++ {
				for k := off; k < off+d4; k++ {
					p.a[d2+k] ^= p.a[d2+d4+k]
					p.a[d4+k] ^= p.a[d2+k]
				}
				off += d
			}
		}
	}
	return nil
}

func (p *poly) iRadixConversion() error {
	m := p.m()
	n := p.length()
	if n != 1<<m {
		return fmt.Errorf("inverse radix conversion expects polynomial length as power of two: %d, %d", m, n)
	}
	for r := m - 2; r >= 0; r-- {
		for i := m - r - 2; i >= 0; i-- {
			d := n >> i
			d2, d4 := d>>1, d>>2
			var off int
			for j := 0; j < 1<<i; j++ {
				for k := off; k < off+d4; k++ {
					p.a[d4+k] ^= p.a[d2+k]
					p.a[d2+k] ^= p.a[d2+d4+k]
				}
				off += d
			}
		}
	}
	return nil
}

// fftNaive evaluates polynomial at combinations of default bases.
func (p *poly) fftNaive() {
	copy(p.a[:], p.eval(defaultBasis.combinations).a[:])
}

func (p *poly) fft() (*poly, error) {
	m := p.m()
	n := p.length()
	if n != 1<<m {
		return nil, fmt.Errorf("fft operation expects polynomial length as power of two, %d, %d", m, n)
	}
	// Applies taylor expantion
	_ = p.radixConversion()

	// We will use subsets of basis combinations
	// throughout the iterations.
	// Since basis with the last element 1 is only accepted,
	// we don't need to calculate twisting operations which are
	// step 2 and step 4 in GM10 Algorithm 2.
	G := defaultBasis.subCombinations

	// Step 1 in GM10 Algorithm 2.
	// Linear evaluation at the leafs of the recursions.
	// f(0), f(β_1) where f is linear polynomial and β_1 = β_m = 1.
	// f(0) = a_0 + 0 * a_1 = a_0
	// f(β_1) = a_0 + a_1 * β_1 = a_0 + a_1
	halfL := 1 << (m - 1)
	for i := 0; i < halfL; i++ {
		p.a[i+halfL] ^= p.a[i]
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
				p.a[k] ^= mul64(G[j], p.a[k2])
				p.a[k2] ^= p.a[k]
			}
			b += (d << 1)
		}
	}
	return p, nil
}

// lfft stands for lazy fft, lfft skips radix conversion phase
func (p *poly) lfft() (*poly, error) {
	m := p.m()
	n := p.length()
	if n != 1<<m {
		return nil, fmt.Errorf("fft operation expects polynomial length as power of two, %d, %d", m, n)
	}
	G := defaultBasis.subCombinations
	halfL := 1 << (m - 1)
	for i := 0; i < halfL; i++ {
		p.a[i+halfL] ^= p.a[i]
	}
	for i := 1; i < m; i++ {
		d := 1 << (m - 1 - i)
		var b int
		for j := 0; j < 1<<i; j++ {
			for k := b; k < b+d; k++ {
				butterfly(&p.a[k], &p.a[k+d], G[j])
			}
			b += (d << 1)
		}
	}
	return p, nil
}

func (p *poly) ifft() (*poly, error) {
	m := p.m()
	n := p.length()
	if n != 1<<m {
		return nil, fmt.Errorf("ifft operation expects polynomial length as power of two, %d, %d", m, n)
	}
	G := defaultBasis.subCombinations
	var i, j, k, d int
	for i = m - 1; i > 0; i-- {
		d = 1 << (m - 1 - i)
		var b int
		for j = 0; j < 1<<i; j++ {
			for k = b; k < b+d; k++ {
				ibutterfly(&p.a[k], &p.a[k+d], G[j])
			}
			b += (d << 1)
		}
	}
	halfL := 1 << (m - 1)
	for i := 0; i < halfL; i++ {
		p.a[i+halfL] ^= p.a[i]
	}
	_ = p.iRadixConversion()
	return p, nil
}

func (p *poly) lifft() (*poly, error) {
	m := p.m()
	n := p.length()
	if n != 1<<m {
		return nil, fmt.Errorf("ifft operation expects polynomial length as power of two, %d, %d", m, n)
	}
	G := defaultBasis.subCombinations
	var i, j, k, d int
	for i = m - 1; i > 0; i-- {
		d = 1 << (m - 1 - i)
		var b int
		for j = 0; j < 1<<i; j++ {
			for k = b; k < b+d; k++ {
				ibutterfly(&p.a[k], &p.a[k+d], G[j])
			}
			b += (d << 1)
		}
	}
	halfL := 1 << (m - 1)
	for i := 0; i < halfL; i++ {
		p.a[i+halfL] ^= p.a[i]
	}
	return p, nil
}

func (p *poly) mulN(q *poly) {
	r := p.degree() + q.degree() + 1
	R := make([]uint64, r)
	for i := 0; i < p.length(); i++ {
		for j := 0; j < q.length(); j++ {
			R[i+j] ^= mul64(p.a[i], q.a[j])
		}
	}
	p.expand(r)
	for i := 0; i < r; i++ {
		p.a[i] = R[i]
	}
}

func (p *poly) muls2(q *poly) {
	p.expand(3)
	p.a[2] = mul64(p.a[1], q.a[1])
	p.a[1] = mul64(p.a[1], q.a[0]) ^ mul64(p.a[0], q.a[1])
	p.a[0] = mul64(p.a[0], q.a[0])
}

func (p *poly) muls3(q *poly) {
	p.expand(5)
	p.a[4] = mul64(p.a[2], q.a[2])
	p.a[3] = mul64(p.a[1], q.a[2]) ^ mul64(p.a[2], q.a[1])
	p.a[2] = mul64(p.a[0], q.a[2]) ^ mul64(p.a[1], q.a[1]) ^ mul64(p.a[2], q.a[0])
	p.a[1] = mul64(p.a[1], q.a[0]) ^ mul64(p.a[0], q.a[1])
	p.a[0] = mul64(p.a[0], q.a[0])
}

func (p *poly) muls4(q *poly) {
	p.expand(7)
	p.a[6] = mul64(p.a[3], q.a[3])
	p.a[5] = mul64(p.a[3], q.a[2]) ^ mul64(p.a[2], q.a[3])
	p.a[4] = mul64(p.a[3], q.a[1]) ^ mul64(p.a[1], q.a[3]) ^ mul64(p.a[2], q.a[2])
	p.a[3] = mul64(p.a[3], q.a[0]) ^ mul64(p.a[0], q.a[3]) ^ mul64(p.a[1], q.a[2]) ^ mul64(p.a[2], q.a[1])
	p.a[2] = mul64(p.a[1], q.a[1]) ^ mul64(p.a[0], q.a[2]) ^ mul64(p.a[2], q.a[0])
	p.a[1] = mul64(p.a[1], q.a[0]) ^ mul64(p.a[0], q.a[1])
	p.a[0] = mul64(p.a[0], q.a[0])
}

func (p *poly) mulSample(q *poly) (*poly, error) {
	// TODO expand?
	n := p.length()
	if n != q.length() {
		return nil, errors.New("expect equal sized polynomials")
	}
	for i := 0; i < n; i++ {
		mulassign64(&p.a[i], q.a[i])
	}
	return p, nil
}

func (p *poly) mul(q *poly) (*poly, error) {
	q1 := q.clone()
	n := p.length()
	if n == 1 {
		p.a[0] = mul64(p.a[0], q.a[0])
	} else if n == 2 {
		p.muls2(q1)
	} else if n == 3 {
		p.muls3(q1)
	} else if n == 4 {
		p.muls4(q1)
	} else {
		m := p.m()
		n := 1 << m
		p.expand(2 * n)
		q1.expand(2 * n)
		if _, err := p.fft(); err != nil {
			return nil, err
		}
		if _, err := q1.fft(); err != nil {
			return nil, err
		}
		if _, err := p.mulSample(q1); err != nil {
			return nil, err
		}
		if _, err := p.ifft(); err != nil {
			return nil, err
		}
	}
	return p, nil
}

func (p *poly) invSample() (*poly, error) {
	n := p.length()
	tA := []uint64{}
	j, setFirst := 0, false
	for i := 0; i < n; i++ {
		if p.a[i] != 0 {
			if !setFirst {
				tA = append(tA, p.a[i])
				setFirst = true
			} else {
				tA = append(tA, mul64(p.a[i], tA[j-1]))
			}
			j++
		}
	}

	tB := make([]uint64, j)
	tB[j-1] = inverse(tA[j-1])
	j--
	for i := n - 1; i >= 0; i-- {
		if p.a[i] != 0 {
			tB[j-1] = mul64(p.a[i], tB[j])
			j--
			if j == 0 {
				break
			}
		}
	}

	if j != 0 {
		panic("bad multi inv implementation")
	}

	setFirst = false
	for i := 0; i < n; i++ {
		if p.a[i] != 0 {
			if !setFirst {
				p.a[i] = tB[j]
				setFirst = true
			} else {
				p.a[i] = mul64(tA[j-1], tB[j])
			}
			j++
		}
	}
	return p, nil
}

// only expect perfect division
func (p *poly) div(q *poly) (*poly, error) {
	m := p.m()
	n := 1 << m
	q1 := q.clone()
	p.expand(2 * n)
	q1.expand(2 * n)
	if _, err := p.fft(); err != nil {
		return nil, err
	}
	if _, err := q1.fft(); err != nil {
		return nil, err
	}
	if _, err := q1.invSample(); err != nil {
		return nil, err
	}
	if _, err := p.mulSample(q1); err != nil {
		return nil, err
	}
	if _, err := p.ifft(); err != nil {
		return nil, err
	}

	return p, nil
}

func (p *poly) debug(desc string) {
	fmt.Println(desc, len(p.a))
	for i := 0; i < len(p.a); i++ {
		fmt.Println(toHex(p.a[i]))
	}
}

func z(roots []uint64) (*poly, error) {
	l := len(roots)
	m := log2Ceil(l)
	n := 1 << m
	p := make([]poly, n)
	for i := 0; i < l; i++ {
		p[i] = poly{[]uint64{roots[i], uint64(1)}}
	}
	// fill with 1 * x ^ 0 + 0 * x ^ 1
	// this is obviously suboptimal
	// but let's leave it for a while in sake of simplicity
	for i := l; i < n; i++ {
		p[i] = poly{[]uint64{1, 0}}
	}
	for i := 0; i < m; i++ {
		z := 1 << (m - i - 1)
		for j := 0; j < z; j++ {
			if _, err := p[j].mul(&p[j+z]); err != nil {
				return nil, err
			}
		}
	}
	return &p[0], nil
}
