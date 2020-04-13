package gf

import (
	"fmt"
	"testing"
)

func (p GF64Poly) debug(desc string) {
	fmt.Println(desc, len(p))
	for i := 0; i < len(p); i++ {
		fmt.Printf("%#16.16x\n", p[i])
	}
}

func TestBasisConversion(t *testing.T) {
	// Extension degree is 64
	n := 64
	// Find cantor bases where βm = 1.
	a := GF64(14862409554843134308)
	bases := make([]GF64, n)
	bases[0] = a
	for j := 0; j < n-1; j++ {
		a = a.Square().Add(a)
		bases[j+1] = a
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
	twisted := make([]GF64, m)
	for j := 0; j < m; j++ {
		twisted[j] = bases[j].Square().Add(bases[j])
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

func TestRadixConversation(t *testing.T) {
	m := uint(16)
	b := configureDefaultBasis64(defaultBasisGenerator, m)
	coeffs0 := make([]GF64, b.n)
	for i := 0; i < int(b.n); i++ {
		coeffs0[i] = GF64(1 << uint(i))
	}
	f0 := NewGF64Poly(coeffs0)
	f1 := f0.Clone()
	_ = b.radixConversion(f0)
	_ = b.iRadixConversion(f0)
	for i := 0; i < int(b.n); i++ {
		if !f0[i].Equal(f1[i]) {
			t.Errorf("")
		}
	}
}

func TestFFT(t *testing.T) {
	m := uint(8)
	b := configureDefaultBasis64(defaultBasisGenerator, m)
	f0 := randGF64Poly(b.n)
	f1 := f0.Clone()
	b.fftNaive(f0)
	_ = b.FFT(f1)
	for i := 0; i < int(b.n); i++ {
		if !f0[i].Equal(f1[i]) {
			t.Fatalf("%d %#16.16x, %#16.16x\n", f0[i], f1[i], i)
		}
	}
}

func TestPolyMultiplication(t *testing.T) {
	m := uint(4)
	b := configureDefaultBasis64(defaultBasisGenerator, m)
	f0 := randGF64Poly(b.n/2 - 1)
	f1 := randGF64Poly(b.n/2 - 1)
	r0 := f0.mulNaive(f1)
	b.Expand(&f0)
	b.Expand(&f1)
	_ = b.Mul(f0, f1)
	if !r0.EqualInCoeff(f0) {
		r0.debug("r0")
		f0.debug("f0")
		t.Fatalf("polynomial multiplication fails")
	}
}

func BenchmarkRadix(t *testing.B) {
	m := uint(16)
	b := configureDefaultBasis64(defaultBasisGenerator, m)
	f0 := randGF64Poly(b.n)
	t.ResetTimer()
	for i := 0; i < t.N; i++ {
		_ = b.radixConversion(f0)
	}
}

func BenchmarkFFT(t *testing.B) {
	m := uint(16)
	b := configureDefaultBasis64(defaultBasisGenerator, m)
	f0 := randGF64Poly(b.n)
	t.ResetTimer()
	for i := 0; i < t.N; i++ {
		_ = b.FFT(f0)
	}
}
