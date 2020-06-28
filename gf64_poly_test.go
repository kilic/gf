package gf

import (
	"testing"
)

func TestRadixConversation(t *testing.T) {
	m := 8
	n := 1 << m
	initDefaultBasis(m)
	coeffs0 := make([]uint64, n)
	for i := 0; i < int(n); i++ {
		coeffs0[i] = uint64(1 << i)
	}
	f0 := newPoly(coeffs0)
	f1 := f0.clone()
	if err := f0.radixConversion(); err != nil {
		t.Fatal(err)
	}
	if err := f0.iRadixConversion(); err != nil {
		t.Fatal(err)
	}
	if !f0.equalInCoeff(f1) {
		t.Fatal("radix conversion failed")
	}
}

func TestFFT(t *testing.T) {
	m := 8
	n := 1 << m
	initDefaultBasis(m)
	f0 := randPoly(n)
	f1 := f0.clone()
	f0.fftNaive()
	if _, err := f1.fft(); err != nil {
		t.Fatal(err)
	}
	if !f0.equalInCoeff(f1) {
		t.Fatalf("fft failed")
	}
	coeffs := make([]uint64, n)
	a0 := uint64(0xff)
	coeffs[0] = a0
	f0 = newPoly(coeffs)
	if _, err := f0.fft(); err != nil {
		t.Fatal(err)
	}
	for i := 0; i < 1<<8; i++ {
		if f0.a[i] != a0 {
			t.Fatalf("fft failed")
		}
	}
}

func TestMulN(t *testing.T) {
	A := randPoly(4)
	B := randPoly(4)
	A.mulN(B)
	if A.length() != 7 {
		t.Fatal("bad naive mul result degree")
	}
	A0 := randPoly(2)
	B0 := randPoly(2)
	A1 := A0.clone()
	B1 := B0.clone()
	A0.mulN(B0)
	A1.muls2(B1)
	if A0.length() != A1.length() {
		t.Fatal("bad degree, s2")
	}
	if !A0.equalInCoeff(A1) {
		t.Fatal("bad mul result, s2")
	}
	A0 = randPoly(3)
	B0 = randPoly(3)
	A1 = A0.clone()
	B1 = B0.clone()
	A0.mulN(B0)
	A1.muls3(B1)
	if A0.length() != A1.length() {
		t.Fatal("bad degree, s3")
	}
	if !A0.equalInCoeff(A1) {
		t.Fatal("bad mul result, s3")
	}
	A0 = randPoly(4)
	B0 = randPoly(4)
	A1 = A0.clone()
	B1 = B0.clone()
	A0.mulN(B0)
	A1.muls4(B1)
	if A0.length() != A1.length() {
		t.Fatal("bad degree, s4")
	}
	if !A0.equalInCoeff(A1) {
		t.Fatal("bad mul result, s4")
	}
}

func TestPolyMultiplication(t *testing.T) {
	initDefaultBasis(16)
	A0 := randPoly(1 << 8)
	B0 := randPoly(1 << 8)
	A1 := A0.clone()
	B1 := B0.clone()
	if _, err := A0.mul(B0); err != nil {
		t.Fatal(err)
	}
	A1.mulN(B1)
	if !A1.equalInCoeff(A0) {
		t.Error("polynomial multiplicaiton failed")
	}
	A0 = randPoly(1 << 15)
	B0 = randPoly(1 << 15)
	C0 := A0.clone()
	if _, err := C0.mul(B0.clone()); err != nil {
		t.Fatal(err)
	}
	for i := 0; i < 100; i++ {
		k := randGF64()
		// a0 := A0.evalSingle(1<<15, k)
		a0 := A0.evalSingle(k)
		// b0 := B0.evalSingle(1<<15, k)
		b0 := B0.evalSingle(k)
		// c0 := C0.evalSingle(1<<16, k)
		c0 := C0.evalSingle(k)
		c1 := mul64(a0, b0)
		if c0 != c1 {
			t.Error("polynomial multiplication failed")
		}
	}
}

func TestPolySubstitution(t *testing.T) {
	initDefaultBasis(16)
	m := 8
	A0 := randPoly(1 << m)
	A1 := A0.clone()
	k := randGF64()
	A1.substitute(k)
	i := randGF64()
	ik := mul64(i, k)
	// e0 := A0.evalSingle(1<<m, ik)
	e0 := A0.evalSingle(ik)
	// e1 := A1.evalSingle(1<<m, i)
	e1 := A1.evalSingle(i)
	if e0 != e1 {
		t.Fatal("kx transformation failed")
	}
}

func TestPolyDiv(t *testing.T) {
	initDefaultBasis(16)
	for i := 0; i < 100; i++ {
		m := 3
		n := 1 << m
		roots := make([]uint64, n)
		for i := 0; i < n; i++ {
			roots[i] = randGF64()
		}
		a0, _ := z(roots)
		a1 := a0.clone()
		b, _ := z(roots[:n/2])
		q, _ := z(roots[n/2:])
		if _, err := a1.div(b); err != nil {
			t.Fatal(err)
		}
		if !a1.equalInCoeff(q) {
			t.Fatal(i)
		}
	}
	for i := 0; i < 100; i++ {
		a0 := randPoly(1 << 8)
		b0 := randPoly(1 << 8)
		a1 := a0.clone()
		b1 := b0.clone()
		a1.mul(b0)
		a1.div(b1)
		a1.trimZeros()
		if a1.length() != a0.length() {
			t.Fatal("division for quotient failed")
		}
		if !a1.equalInCoeff(a0) {
			t.Fatal("division for quotient failed")
		}
	}
}

func TestPolySampleInv(t *testing.T) {
	initDefaultBasis(16)
	A0 := newPoly([]uint64{10, 11, 12, 13})
	A1 := A0.clone()
	if _, err := A1.invSample(); err != nil {
		t.Fatal(err)
	}
	for i := 0; i < A0.length(); i++ {
		c := mul64(A0.a[i], A1.a[i])
		if c != 1 {
			t.Fatal("multi inversion failed", i)
		}
	}
	A0 = newPoly([]uint64{0, 0, 10, 11, 12, 13, 0, 0})
	A1 = A0.clone()
	if _, err := A1.invSample(); err != nil {
		t.Fatal(err)
	}
	for i := 0; i < A0.length(); i++ {
		if A0.a[i] == 0 {
			if A1.a[i] != 0 {
				t.Fatal("multi inversion failed", i)
			}
			continue
		}
		c := mul64(A0.a[i], A1.a[i])
		if c != 1 {
			t.Fatal("multi inversion failed", i)
		}
	}
	A0 = newPoly([]uint64{21, 20, 10, 0, 0, 13, 30, 31})
	A1 = A0.clone()
	if _, err := A1.invSample(); err != nil {
		t.Fatal(err)
	}
	for i := 0; i < A0.length(); i++ {
		if A0.a[i] == 0 {
			if A1.a[i] != 0 {
				t.Fatal("multi inversion failed", i)
			}
			continue
		}
		c := mul64(A0.a[i], A1.a[i])
		if c != 1 {
			t.Fatal("multi inversion failed", i)
		}
	}
	A0 = randPoly(1 << 8)
	A1 = A0.clone()
	if _, err := A0.invSample(); err != nil {
		t.Fatal(err)
	}
	for i := 0; i < A0.length(); i++ {
		if A0.a[i] == 0 {
			if A1.a[i] != 0 {
				t.Fatal("multi inversion failed", i)
			}
			continue
		}
		c := mul64(A0.a[i], A1.a[i])
		if c != 1 {
			t.Fatal("multi inversion failed", i)
		}
	}
	if _, err := A0.invSample(); err != nil {
		t.Fatal(err)
	}
	if !A0.equalInCoeff(A1) {
		t.Fatal("multi inversion failed")
	}
}

func TestZPoly(t *testing.T) {
	m := 16
	initDefaultBasis(m)
	roots := randPoly(1 << 4).a
	Z, err := z(roots)
	if err != nil {
		t.Fatal(err)
	}
	for i := 0; i < len(roots); i++ {
		// e := Z.evalSingle(len(roots)+1, roots[i])
		e := Z.evalSingle(roots[i])
		if e != 0 {
			t.Fatalf("evaluation at root must be zero")
		}
	}
}

func BenchmarkRadixConversion(t *testing.B) {
	m := polyLen
	n := 1 << m
	initDefaultBasis(m)
	f0 := randPoly(n)
	t.ResetTimer()
	for i := 0; i < t.N; i++ {
		if err := f0.radixConversion(); err != nil {
			t.Fatal(err)
		}
	}
}

func BenchmarkFFT(t *testing.B) {
	m := polyLen
	n := 1 << m
	initDefaultBasis(m)
	f0 := randPoly(n)
	t.ResetTimer()
	for i := 0; i < t.N; i++ {
		if _, err := f0.fft(); err != nil {
			t.Fatal(err)
		}
	}
}

func BenchmarkZPoly(t *testing.B) {
	m := polyLen
	initDefaultBasis(m)
	roots := randPoly(1 << (m - 1))
	t.ResetTimer()
	for i := 0; i < t.N; i++ {
		if _, err := z(roots.a); err != nil {
			t.Fatal(err)
		}
	}
}
