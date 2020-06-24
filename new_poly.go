package gf

import (
	"errors"
	"fmt"
)

var basis *Basis64

func configureDefaultBasis(m int) {
	var err error
	basis, err = NewBasis64(defaultBasisGenerator64, m)
	if err != nil {
		panic(err)
	}
}

type Poly struct {
	a []GF64
}

func newPoly(coeffs []GF64) *Poly {
	return &Poly{coeffs}
}

func newEmptyPoly(len int) *Poly {
	return &Poly{make([]GF64, len)}
}

func randPoly(n int) *Poly {
	coeffs := make([]GF64, n)
	for i := 0; i < n; i++ {
		coeffs[i] = randGF64()
	}
	return newPoly(coeffs)
}

func (p *Poly) degree() int {
	return p.length() - 1
}

func (p *Poly) length() int {
	return len(p.a)
}

func (p *Poly) clone() *Poly {
	q := make([]GF64, p.length())
	copy(q, p.a)
	return newPoly(q)
}

func (p *Poly) equalInCoeff(q *Poly) bool {
	lp := p.length()
	lq := q.length()
	if lp > lq {
		for i := 0; i < lq; i++ {
			if !p.a[i].Equal(q.a[i]) {
				return false
			}
		}
		for i := lp - 1; i < lp-lq; i++ {
			if !p.a[i].IsZero() {
				return false
			}
		}
	} else {
		for i := 0; i < lp; i++ {
			for i := 0; i < lp; i++ {
				if !p.a[i].Equal(q.a[i]) {
					return false
				}
			}
			for i := lq - 1; i < lp-lq; i++ {
				if !q.a[i].IsZero() {
					return false
				}
			}
		}
	}
	return true
}

func (p *Poly) eval(n int, D []GF64) Poly {
	evals := make([]GF64, len(p.a))
	for i := 0; i < len(D); i++ {
		evals[i] = p.evalSingle(n, D[i])
	}
	return Poly{evals}
}

func (p *Poly) evalSingle(n int, x GF64) GF64 {
	if n > p.length() {
		n = p.length()
	}
	acc := NewGF64()
	for i := n - 1; i >= 0; i-- {
		acc = mul64(acc, x) ^ p.a[i]
	}
	return acc
}

func (p *Poly) expand(l int) int {
	l2 := p.length()
	if l > l2 {
		p.a = append(p.a, make([]GF64, l-l2)...)
		return l
	}
	return l2
}

func (p *Poly) kx(n int, k GF64) {
	acc := k
	for i := 1; i < n; i++ {
		mulassign64(&p.a[i], acc)
		mulassign64(&acc, k)
	}
}

func (p *Poly) add(q Poly) {
	l := len(p.a)
	if len(q.a) > l {
		return
	}
	for i := 0; i < l; i++ {
		p.a[i] ^= q.a[i]
	}
}

func (p *Poly) radixConversion(m int) error {
	n := 1 << m
	if p.length() < n {
		return errors.New("radix convertions is not applicable")
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

func (p *Poly) iRadixConversion(m int) error {
	n := 1 << m
	if p.length() < n {
		return errors.New("radix convertions is not applicable")
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

func (p *Poly) fftNaive(m int) {
	copy(p.a[:], p.eval(1<<m, basis.combinations).a[:])
}

func (p *Poly) fft(m int) error {
	if p.length() > 1<<m {
		return errors.New("fft is not applicable")
	}
	_ = p.radixConversion(m)
	G := basis.subCombinations
	halfL := 1 << (m - 1)
	for i := 0; i < halfL; i++ {
		p.a[i+halfL] ^= p.a[i]
	}
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
	return nil
}

func (p *Poly) ifft(m int) error {
	if p.length() > 1<<m {
		return errors.New("fft is not applicable")
	}
	G := basis.subCombinations
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
	_ = p.iRadixConversion(m)
	return nil
}

func (a *Poly) mulN(b *Poly) {
	r := a.degree() + b.degree() + 1
	R := make([]GF64, r)
	for i := 0; i < a.length(); i++ {
		for j := 0; j < b.length(); j++ {
			R[i+j] ^= mul64(a.a[i], b.a[j])
		}
	}
	a.expand(r)
	for i := 0; i < r; i++ {
		a.a[i] = R[i]
	}
}

func (a *Poly) muls2(b *Poly) {
	a.expand(3)
	a.a[2] = mul64(a.a[1], b.a[1])
	a.a[1] = mul64(a.a[1], b.a[0]) ^ mul64(a.a[0], b.a[1])
	a.a[0] = mul64(a.a[0], b.a[0])
}

func (a *Poly) muls3(b *Poly) {
	a.expand(5)
	a.a[4] = mul64(a.a[2], b.a[2])
	a.a[3] = mul64(a.a[1], b.a[2]) ^ mul64(a.a[2], b.a[1])
	a.a[2] = mul64(a.a[0], b.a[2]) ^ mul64(a.a[1], b.a[1]) ^ mul64(a.a[2], b.a[0])
	a.a[1] = mul64(a.a[1], b.a[0]) ^ mul64(a.a[0], b.a[1])
	a.a[0] = mul64(a.a[0], b.a[0])
}

func (a *Poly) muls4(b *Poly) {
	a.expand(7)
	a.a[6] = mul64(a.a[3], b.a[3])
	a.a[5] = mul64(a.a[3], b.a[2]) ^ mul64(a.a[2], b.a[3])
	a.a[4] = mul64(a.a[3], b.a[1]) ^ mul64(a.a[1], b.a[3]) ^ mul64(a.a[2], b.a[2])
	a.a[3] = mul64(a.a[3], b.a[0]) ^ mul64(a.a[0], b.a[3]) ^ mul64(a.a[1], b.a[2]) ^ mul64(a.a[2], b.a[1])
	a.a[2] = mul64(a.a[1], b.a[1]) ^ mul64(a.a[0], b.a[2]) ^ mul64(a.a[2], b.a[0])
	a.a[1] = mul64(a.a[1], b.a[0]) ^ mul64(a.a[0], b.a[1])
	a.a[0] = mul64(a.a[0], b.a[0])
}

func (a *Poly) mulSample(n int, b *Poly) error {
	// TODO expand?
	if n < a.length() || n < b.length() {
		return errors.New("sample-wise mul is not applicable")
	}
	for i := 0; i < n; i++ {
		mulassign64(&a.a[i], b.a[i])
	}
	return nil
}

func (a *Poly) mul(m int, b *Poly) error {
	// TODO: should be non descructive for polyn b
	n := a.length()
	if n == 1 {
		a.a[0] = mul64(a.a[0], b.a[0])
	} else if n == 2 {
		a.muls2(b)
	} else if n == 3 {
		a.muls3(b)
	} else if n == 4 {
		a.muls4(b)
	} else {
		// TODO: use log?
		n := 1 << m
		a.expand(2 * n)
		b.expand(2 * n)
		if err := a.fft(m + 1); err != nil {
			return err
		}
		if err := b.fft(m + 1); err != nil {
			return err
		}
		if err := a.mulSample(1<<(m+1), b); err != nil {
			return err
		}
		if err := a.ifft(m + 1); err != nil {
			return err
		}
	}
	return nil
}

func (f *Poly) invSample(n int) error {
	if f.length() < n {
		return errors.New("poly is too short")
	}

	tA := []GF64{}
	j, setFirst := 0, false
	for i := 0; i < n; i++ {
		if !f.a[i].IsZero() {
			if !setFirst {
				tA = append(tA, f.a[i])
				setFirst = true
			} else {
				tA = append(tA, mul64(f.a[i], tA[j-1]))
			}
			j++
		}
	}

	tB := make([]GF64, j)
	tB[j-1] = tA[j-1].Inverse()
	j--
	for i := n - 1; i >= 0; i-- {
		if !f.a[i].IsZero() {
			tB[j-1] = mul64(f.a[i], tB[j])
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
		if !f.a[i].IsZero() {
			if !setFirst {
				f.a[i] = tB[j]
				setFirst = true
			} else {
				f.a[i] = mul64(tA[j-1], tB[j])
			}
			j++
		}
	}
	return nil
}

func (f *Poly) invSampleNaive(n int) {
	for i := 0; i < n; i++ {
		f.a[i] = f.a[i].Inverse()
	}
}

func (a *Poly) div(m int, b *Poly) error {
	n := 1 << m
	a.expand(2 * n)
	b.expand(2 * n)
	if err := a.fft(m + 1); err != nil {
		return err
	}
	if err := b.fft(m + 1); err != nil {
		return err
	}
	if err := b.invSample(1 << (m + 1)); err != nil {
		return err
	}
	if err := a.mulSample(1<<(m+1), b); err != nil {
		return err
	}
	if err := a.ifft(m + 1); err != nil {
		return err
	}

	return nil
}

func (p *Poly) debug(desc string) {
	fmt.Println(desc, len(p.a))
	for i := 0; i < len(p.a); i++ {
		fmt.Println(p.a[i].hex())
	}
}

func z(m int, roots []GF64) (*Poly, error) {
	l := len(roots)
	if l&1 != 0 {
		return nil, errors.New("expect even number of roots for simplicity")
	}
	p := make([]Poly, l)
	for i := 0; i < l; i++ {
		p[i] = *newPoly([]GF64{roots[i], GF64(1)})
	}
	for i := 0; i < m; i++ {
		z := 1 << (m - i - 1)
		for j := 0; j < z; j++ {
			err := p[j].mul(i+1, &p[j+z])
			if err != nil {
				return nil, err
			}
		}
	}
	return &p[0], nil
}

func log2(m int) int {
	r := 0
	for {
		m = m >> 1
		if m == 0 {
			return r
		}
		r += 1
	}
}

// func (f *Poly) inv2(n int) error {
// 	if f.length() < n {
// 		return errors.New("poly is too short")
// 	}
// 	tA := make([]GF64, n)
// 	tA[0] = f.a[0]
// 	for i := 1; i < n; i++ {
// 		tA[i] = mul64(f.a[i], tA[i-1])
// 	}

// 	tB := make([]GF64, n)
// 	tB[n-1] = tA[n-1].Inverse()
// 	for i := n - 1; i > 0; i-- {
// 		fmt.Println(i, f.a[i].hex(), tB[i].hex())
// 		tB[i-1] = mul64(f.a[i], tB[i])
// 	}

// 	fmt.Println("TB")
// 	for i := 0; i < len(tA); i++ {
// 		fmt.Println(tB[i].hex())
// 	}

// 	f.a[0] = tB[0]
// 	for i := 1; i < n; i++ {
// 		f.a[i] = mul64(tA[i-1], tB[i])
// 	}
// 	return nil
// }

// func (basis *Basis64) z(roots []GF64) (GF64Poly, error) {
// 	l := len(roots)
// 	p := make([]GF64Poly, l)
// 	if len(roots)&1 != 0 {
// 		return nil, errors.New("expect even number of roots for simplicity")
// 	}
// 	for i := 0; i < l; i++ {
// 		p[i] = NewGF64Poly([]GF64{roots[i], GF64(1)})
// 	}
// 	for i := 0; i < l>>2; i++ {
// 		fmt.Println("xxx", i)
// 		p[i].mul(p[i+1])

// 	}

// 	return nil, nil
// }
