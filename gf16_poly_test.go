package gf

import (
	"fmt"
	"testing"
)

func (p GF16Poly) debug(desc string) {
	fmt.Println(desc, len(p))
	for i := 0; i < len(p); i++ {
		fmt.Printf("%#4.4x\n", p[i])
	}
}

func TestBasisConversion16(t *testing.T) {
	// Extension degree is 16
	n := 16
	// Find cantor bases where βm = 1.
	a := GF16(0x2000)
	bases := make([]GF16, n)
	bases[0] = a
	for j := 0; j < n-1; j++ {
		bases[j+1] = mul16(bases[j], bases[j]) ^ bases[j]
	}
	if !bases[n-1].IsOne() {
		// never expected though
		t.Fatalf("given initial β0 last basis expected to be 1")
	}
	// Twist initial set.
	// γi = βi * βm ^ -1
	// βm = 1
	// γi = βi
	// δi = γi ^ 2 − γi
	m := n - 1
	twisted := make([]GF16, m)
	for j := 0; j < m; j++ {
		twisted[j] = mul16(bases[j], bases[j]) ^ bases[j]
	}
	if !twisted[m-1].IsOne() {
		// never expected though
		t.Fatalf("last element of twisted basis expected to be 1")
	}
	for i := 0; i < m; i++ {
		if !twisted[i].Equal(bases[i+1]) {
			t.Fatal("bad twisting", i)
		}
	}
}

func TestRadixConversation16(t *testing.T) {
	m := uint(8)
	b := configureDefaultBasis16(m)
	coeffs0 := make([]GF16, b.n)
	for i := 0; i < int(b.n); i++ {
		coeffs0[i] = GF16(1 << uint(i))
	}
	f0 := NewGF16Poly(coeffs0)
	f1 := f0.Clone()
	_ = b.radixConversion(f0)
	_ = b.iRadixConversion(f0)
	if !f0.EqualInCoeff(f1) {
		t.Fatalf("radix conversion fails")
	}
}

func TestFFT16(t *testing.T) {
	m := uint(8)
	b := configureDefaultBasis16(m)
	f0 := randGF16Poly(b.n)
	f1 := f0.Clone()
	b.fftNaive(f0)
	_ = b.FFT(f1)
	if !f0.EqualInCoeff(f1) {
		t.Fatalf("fft16 fails")
	}
}

func TestLazyFFT16(t *testing.T) {
	m := uint(8)
	b := configureDefaultBasis16(m)
	f0 := randGF16Poly(b.n)
	f1 := f0.Clone()
	_ = b.lFFT(f1)
	_ = b.lIFFT(f1)
	if !f0.EqualInCoeff(f1) {
		t.Fatalf("lazy fft16 fails")
	}
}

func TestPolyMultiplication16(t *testing.T) {
	m := uint(4)
	b := configureDefaultBasis16(m)
	f0 := randGF16Poly(b.n/2 - 1)
	f1 := randGF16Poly(b.n/2 - 1)
	r0 := f0.mulNaive(f1)
	b.Expand(&f0)
	b.Expand(&f1)
	_ = b.Mul(f0, f1)
	if !r0.EqualInCoeff(f0) {
		t.Fatalf("polynomial multiplication fails")
	}
}
