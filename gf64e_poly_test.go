package gf

import (
	"testing"
)

func TestFFT64e(t *testing.T) {
	m := uint(8)
	b := configureDefaultBasis64e(m)
	f0 := randGF64ePoly(b.n)
	f1 := f0.Clone()
	b.fftNaive(f0)
	_ = b.FFT(f1)
	if !f0.EqualInCoeff(f1) {
		t.Fatalf("fft fails / expected for now")
	}
}

func BenchmarkFFT64e(t *testing.B) {
	m := polyLen
	b := configureDefaultBasis64e(m)
	f0 := randGF64ePoly(b.n)
	t.ResetTimer()
	for i := 0; i < t.N; i++ {
		_ = b.FFT(f0)
	}
}
