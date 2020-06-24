package gf

import (
	"fmt"
	"testing"
)

// func (f *Poly) invSampleNaive(n int) {
// 	for i := 0; i < n; i++ {
// 		f.a[i] = f.a[i].Inverse()
// 	}
// }

func TestRadixConversation(t *testing.T) {
	m := 4
	configureDefaultBasis(m)
	coeffs0 := make([]GF64, basis.n)
	for i := 0; i < int(basis.n); i++ {
		coeffs0[i] = GF64(1 << i)
	}
	f0 := newPoly(coeffs0)
	f1 := f0.clone()
	if err := f0.radixConversion(m); err != nil {
		t.Fatal(err)
	}
	if err := f0.iRadixConversion(m); err != nil {
		t.Fatal(err)
	}
	if !f0.equalInCoeff(f1) {
		t.Fatal("bad radix failed")
	}
}

func TestFFT(t *testing.T) {
	m := 8
	configureDefaultBasis(m)
	f0 := randPoly(basis.n)
	f1 := f0.clone()
	f0.fftNaive(m)
	if err := f1.fft(m); err != nil {
		t.Fatal(err)
	}
	if !f0.equalInCoeff(f1) {
		t.Fatalf("fft failed")
	}
	coeffs := make([]GF64, basis.n)
	a0 := GF64(0xff)
	coeffs[0] = a0
	f0 = newPoly(coeffs)
	if err := f0.fft(m); err != nil {
		t.Fatal(err)
	}
	for i := 0; i < 1<<8; i++ {
		if !f0.a[i].Equal(a0) {
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
	configureDefaultBasis(16)
	A0 := randPoly(1 << 8)
	B0 := randPoly(1 << 8)
	A1 := A0.clone()
	B1 := B0.clone()
	if err := A0.mul(8, B0); err != nil {
		t.Fatal(err)
	}
	A1.mulN(B1)
	if !A1.equalInCoeff(A0) {
		t.Error("polynomial multiplicaiton failed")
	}
	A0 = randPoly(1 << 15)
	B0 = randPoly(1 << 15)
	C0 := A0.clone()
	if err := C0.mul(15, B0.clone()); err != nil {
		t.Fatal(err)
	}
	for i := 0; i < 100; i++ {
		k := randGF64()
		a0 := A0.evalSingle(1<<15, k)
		b0 := B0.evalSingle(1<<15, k)
		c0 := C0.evalSingle(1<<16, k)
		c1 := mul64(a0, b0)
		if !c0.Equal(c1) {
			t.Error("polynomial multiplication failed")
		}
	}
}

func TestPolyKX(t *testing.T) {
	configureDefaultBasis(16)
	m := 8
	n := 1 << m
	A0 := randPoly(1 << m)
	A1 := A0.clone()
	k := randGF64()
	A1.kx(n, k)
	i := randGF64()
	ik := mul64(i, k)
	e0 := A0.evalSingle(1<<m, ik)
	e1 := A1.evalSingle(1<<m, i)
	fmt.Println(e0.Equal(e1))
	if !e0.Equal(e1) {
		t.Fatal("kx transformation failed")
	}
}

func TestPolyDivision(t *testing.T) {
	configureDefaultBasis(16)

	// m := 3
	// A0 := randPoly(1 << m)
	// A1 := A0.clone()
	// if err := A1.fft(m); err != nil {
	// 	t.Fatal(err)
	// }
	// A1.debug("fft")
	// if err := A1.invSample(1 << m); err != nil {
	// 	t.Fatal(err)
	// }
	// A1.debug("inv sample")
	// if err := A1.ifft(m); err != nil {
	// 	t.Fatal(err)
	// }
	// A1.debug("ifft")

	// // k := randGF64()
	// // k := basis.combinations[2]
	// k := randGF64()
	// e0 := A0.evalSingle(1<<5, k)
	// e1 := A1.evalSingle(1<<5, k)
	// z := mul64(e0, e1)
	// fmt.Println("")
	// fmt.Println(z.hex())
	// fmt.Println(e0.hex())
	// fmt.Println(e1.hex())

	for i := 0; i < 1000000; i++ {
		m := 3
		n := 1 << m
		roots := make([]GF64, n)
		for i := 0; i < n; i++ {
			roots[i] = randGF64()
		}
		a0, _ := z(m, roots)
		a1 := a0.clone()
		b, _ := z(m-1, roots[:n/2])
		q, _ := z(m-1, roots[n/2:])
		if err := a1.div(m+1, b.clone()); err != nil {
			t.Fatal(err)
		}

		// fmt.Println()
		// a0.debug("a0")
		// a1.debug("a1")
		// q.debug("q")
		// fmt.Println(a1.equalInCoeff(q))
		if !a1.equalInCoeff(q) {
			t.Fatal(i)
		}
	}

	// a0, _ := z(2, []GF64{5, 6, 7, 8})
	// a1 := a0.clone()
	// b, _ := z(1, []GF64{6, 7})
	// q, _ := z(1, []GF64{5, 8})
	// if err := a1.div(3, b.clone()); err != nil {
	// 	t.Fatal(err)
	// }
	// a0.debug("a0")
	// a1.debug("a1")
	// q.debug("q")

	// a0, _ := z(3, []GF64{95, 96, 97, 98, 99, 100, 101, 102})
	// a1 := a0.clone()
	// b, _ := z(2, []GF64{95, 96, 97, 98})
	// q, _ := z(2, []GF64{99, 100, 101, 102})
	// if err := a1.div(4, b.clone()); err != nil {
	// 	t.Fatal(err)
	// }
	// a0.debug("a0")
	// a1.debug("a1")
	// q.debug("q")
}

func TestPolySampleInv(t *testing.T) {
	configureDefaultBasis(16)
	A0 := newPoly([]GF64{10, 11, 12, 13})
	A1 := A0.clone()
	_ = A1.invSample(A0.length())
	for i := 0; i < A0.length(); i++ {
		c := mul64(A0.a[i], A1.a[i])
		if !c.IsOne() {
			t.Fatal("multi inversion failed", i)
		}
	}
	A0 = newPoly([]GF64{0, 0, 10, 11, 12, 13, 0, 0})
	A1 = A0.clone()
	_ = A1.invSample(A0.length())
	for i := 0; i < A0.length(); i++ {
		if A0.a[i].IsZero() {
			if !A1.a[i].IsZero() {
				t.Fatal("multi inversion failed", i)
			}
			continue
		}
		c := mul64(A0.a[i], A1.a[i])
		if !c.IsOne() {
			t.Fatal("multi inversion failed", i)
		}
	}
	A0 = newPoly([]GF64{21, 20, 10, 0, 0, 13, 30, 31})
	A1 = A0.clone()
	_ = A1.invSample(A0.length())
	for i := 0; i < A0.length(); i++ {
		if A0.a[i].IsZero() {
			if !A1.a[i].IsZero() {
				t.Fatal("multi inversion failed", i)
			}
			continue
		}
		c := mul64(A0.a[i], A1.a[i])
		if !c.IsOne() {
			t.Fatal("multi inversion failed", i)
		}
	}
	A0 = randPoly(1 << 8)
	A1 = A0.clone()
	if err := A0.invSample(1 << 8); err != nil {
		t.Fatal(err)
	}
	for i := 0; i < A0.length(); i++ {
		if A0.a[i].IsZero() {
			if !A1.a[i].IsZero() {
				t.Fatal("multi inversion failed", i)
			}
			continue
		}
		c := mul64(A0.a[i], A1.a[i])
		if !c.IsOne() {
			t.Fatal("multi inversion failed", i)
		}
	}
	_ = A0.invSample(1 << 8)
	if !A0.equalInCoeff(A1) {
		t.Fatal("multi inversion failed")
	}
	A0 = randPoly(1 << 8)
	A1 = A0.clone()
	if err := A0.invSample(1 << 8); err != nil {
		t.Fatal(err)
	}
	A1.invSampleNaive(1 << 8)
	if !A0.equalInCoeff(A1) {
		t.Fatal("multi inversion failed")
	}

}

func TestZPoly(t *testing.T) {
	m := 16
	configureDefaultBasis(m)
	mz := 10
	roots := randPoly(1 << mz).a
	Z, err := z(mz, roots)
	if err != nil {
		t.Fatal(err)
	}
	for i := 0; i < len(roots); i++ {
		e := Z.evalSingle(len(roots)+1, roots[i])
		if !e.IsZero() {
			t.Fatalf("evaluation at root must be zero")
		}
	}
}

func BenchmarkRadixConversion(t *testing.B) {
	m := polyLen
	configureDefaultBasis(m)
	f0 := randPoly(basis.n)
	t.ResetTimer()
	for i := 0; i < t.N; i++ {
		if err := f0.radixConversion(m); err != nil {
			t.Fatal(err)
		}
	}
}

func BenchmarkFFT(t *testing.B) {
	m := polyLen
	configureDefaultBasis(m)
	f0 := randPoly(basis.n)
	t.ResetTimer()
	for i := 0; i < t.N; i++ {
		if err := f0.fft(m); err != nil {
			t.Fatal(err)
		}
	}
}

func BenchmarkZPoly(t *testing.B) {
	m := polyLen
	configureDefaultBasis(m)
	roots := randPoly(1 << (m - 1))
	t.ResetTimer()
	for i := 0; i < t.N; i++ {
		if _, err := z(m-1, roots.a); err != nil {
			t.Fatal(err)
		}
	}
}

func BenchmarkPolyInversion(t *testing.B) {
	m := polyLen
	configureDefaultBasis(m)
	f0 := randPoly(basis.n)
	t.ResetTimer()
	for i := 0; i < t.N; i++ {
		if err := f0.invSample(basis.n); err != nil {
			t.Fatal(err)
		}
	}
}

// func TestBasisX(t *testing.T) {
// 	m := 4
// 	configureDefaultBasis(m)
// 	combs := map[GF64]int{}
// 	for i := 0; i < len(basis.combinations); i++ {
// 		// fmt.Printf("%16.16x %d\n", basis.combinations[i], i)
// 		combs[basis.combinations[i]] = i
// 	}
// 	l := len(basis.combinations)
// 	x := 0
// 	for i := 0; i < l; i++ {
// 		// y := mul64(basis.combinations[3], basis.combinations[((i+8)+l)%l])
// 		// x := mul64(basis.combinations[3], basis.combinations[i])
// 		// fmt.Printf("%16.16x, %d\t%d\n", x, combs[x], combs[x]^combs[y])
// 		r := mul64(basis.combinations[x], basis.combinations[i])
// 		fmt.Printf("%d\t%d\t%d\t%d\n", x, i, combs[r], x^i)
// 	}
// }

// func TestLog2(t *testing.T) {

// 	fmt.Println(0, log2(0))
// 	fmt.Println(1, log2(1))
// 	fmt.Println(2, log2(2))
// 	fmt.Println(3, log2(3))
// 	fmt.Println(4, log2(4))
// 	fmt.Println(5, log2(5))
// 	fmt.Println(6, log2(6))
// 	fmt.Println(7, log2(7))
// 	fmt.Println(8, log2(8))
// 	fmt.Println(9, log2(9))
// }

// func TestMulBasis(t *testing.T) {
// 	m := 4
// 	configureDefaultBasis(m)

// 	combs := map[int]GF64{}
// 	for i := 0; i < len(basis.combinations); i++ {
// 		// fmt.Printf("%16.16x %d\n", basis.combinations[i], i)
// 		combs[i] = basis.combinations[i]
// 	}

// 	x := 2
// 	y := 7
// 	R := mul64(basis.combinations[x], basis.combinations[y])
// 	Z := combs[(y+x)%16]

// 	fmt.Printf("R: %16.16x\n", R)
// 	fmt.Printf("R: %16.16x\n", Z)

// }

// func TestBasisX(t *testing.T) {
// 	m := 4
// 	configureDefaultBasis(m)
// 	combs := map[GF64]int{}
// 	for i := 0; i < len(basis.combinations); i++ {
// 		fmt.Printf("%16.16x %d\n", basis.combinations[i], i)
// 		combs[basis.combinations[i]] = i
// 	}
// 	fmt.Println("xxx")
// 	l := len(basis.combinations)
// 	// for i := 0; i < l; i++ {
// 	// 	y := mul64(basis.combinations[2], basis.combinations[((i+8)+l)%l])
// 	// 	x := mul64(basis.combinations[2], basis.combinations[i])
// 	// 	fmt.Printf("%16.16x, %d\t%d\n", x, combs[x], combs[x]^combs[y])
// 	// 	// fmt.Printf("%16.16x, %d\t%d\n", x, combs[x], combs[y])
// 	// }
// 	R0 := mul64(basis.combinations[3], basis.combinations[(9+8)%l])
// 	R1 := mul64(basis.combinations[3], basis.combinations[9])
// 	fmt.Printf("%16.16x, %16.16x, %d\t%d\t%d\n", R0, R1, combs[R0], combs[R1], combs[R0]^combs[R1])
// }
