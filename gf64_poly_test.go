package gf

// import (
// 	"fmt"
// 	"testing"
// )

// func TestBasisConversion64(t *testing.T) {
// 	// Extension degree is 64
// 	n := 64
// 	// Find cantor bases where βm = 1.
// 	bases := make([]GF64, n)
// 	bases[0] = defaultBasisGenerator64
// 	for j := 0; j < n-1; j++ {
// 		bases[j+1] = mul64(bases[j], bases[j]) ^ bases[j]
// 	}
// 	if !bases[n-1].IsOne() {
// 		// never expected though
// 		t.Fatalf("given initial β0 last basis expected to be 1")
// 	}
// 	// Twist initial set.
// 	// γi = βi * βm ^ -1
// 	// βm = 1
// 	// γi = βi
// 	// δi = γi ^ 2 − γi
// 	m := n - 1
// 	twisted := make([]GF64, m)
// 	for j := 0; j < m; j++ {
// 		twisted[j] = mul64(bases[j], bases[j]) ^ bases[j]
// 	}
// 	if !twisted[m-1].IsOne() {
// 		// never expected though
// 		t.Fatalf("last element of twisted basis expected to be 1")
// 	}
// 	for i := 0; i < m; i++ {
// 		if !twisted[i].Equal(bases[i+1]) {
// 			t.Fatal("bad twisting", i)
// 		}
// 	}
// }

// func TestRadixConversation64(t *testing.T) {
// 	m := 16
// 	b := configureDefaultBasis64(m)
// 	coeffs0 := make([]GF64, b.n)
// 	for i := 0; i < int(b.n); i++ {
// 		coeffs0[i] = GF64(1 << i)
// 	}
// 	f0 := NewGF64Poly(coeffs0)
// 	f1 := f0.Clone()
// 	_ = b.radixConversion(m, f0)
// 	_ = b.iRadixConversion(m, f0)
// 	for i := 0; i < int(b.n); i++ {
// 		if !f0[i].Equal(f1[i]) {
// 			t.Errorf("")
// 		}
// 	}
// }

// func TestFFT64(t *testing.T) {
// 	m := 8
// 	b := configureDefaultBasis64(m)
// 	f0 := randGF64Poly(b.n)
// 	f1 := f0.Clone()
// 	b.fftNaive(8, f0)
// 	_ = b.fft(8, f1)
// 	if !f0.EqualInCoeff(f1) {
// 		t.Fatalf("fft fails")
// 	}
// }

// // func TestPolyMultiplication64(t *testing.T) {
// // 	m := 8
// // 	b := configureDefaultBasis64(m)
// // 	f0 := randGF64Poly(b.n/2 - 1)
// // 	f1 := randGF64Poly(b.n/2 - 1)
// // 	r0 := f0.mul(f1)
// // 	b.expand(&f0)
// // 	b.expand(&f1)
// // 	err := b.mul(m, f0, f1)
// // 	if err != nil {
// // 		t.Fatal(err)
// // 	}
// // 	if !r0.EqualInCoeff(f0) {
// // 		t.Fatalf("polynomial multiplication fails")
// // 	}
// // }

// func TestPolyMulX(t *testing.T) {
// 	m := 8
// 	b := configureDefaultBasis64(m)
// 	f0 := randGF64Poly(8)
// 	f1 := randGF64Poly(8)
// 	// r0 := f0.mul(f1)

// 	f0.debug("f0")
// 	fmt.Printf("%p\n", &f0)
// 	err := b.mul(3, &f0, &f1)
// 	if err != nil {
// 		t.Fatal(err)
// 	}
// 	f0.debug("f0")
// 	fmt.Printf("%p\n", &f0)
// 	// r0.debug("r0")

// 	// if !r0.EqualInCoeff(f0) {
// 	// 	t.Fatalf("polynomial multiplication fails")
// 	// }
// }

// func BenchmarkFFT64(t *testing.B) {
// 	m := polyLen
// 	b := configureDefaultBasis64(m)
// 	f0 := randGF64Poly(b.n)
// 	t.ResetTimer()
// 	for i := 0; i < t.N; i++ {
// 		_ = b.fft(m, f0)
// 	}
// }

// // func TestSubBasis(t *testing.T) {

// // }

// // func TestZPoly(t *testing.T) {

// // 	var err error
// // 	// b1, _ := NewBasis64(defaultBasisGenerator64, 2)
// // 	b3, _ := NewBasis64(defaultBasisGenerator64, 3)
// // 	b4, _ := NewBasis64(defaultBasisGenerator64, 4)

// // 	// p0 := NewGF64Poly([]GF64{GF64(1), GF64(1)})
// // 	// p1 := NewGF64Poly([]GF64{GF64(2), GF64(1)})
// // 	// p2 := NewGF64Poly([]GF64{GF64(3), GF64(1)})
// // 	// p3 := NewGF64Poly([]GF64{GF64(4), GF64(1)})
// // 	// p4 := NewGF64Poly([]GF64{GF64(5), GF64(1)})
// // 	// p5 := NewGF64Poly([]GF64{GF64(6), GF64(1)})
// // 	// p6 := NewGF64Poly([]GF64{GF64(7), GF64(1)})
// // 	// p7 := NewGF64Poly([]GF64{GF64(8), GF64(1)})

// // 	p0 := NewGF64Poly([]GF64{b4.combinations[1], GF64(1)})
// // 	p1 := NewGF64Poly([]GF64{b4.combinations[2], GF64(1)})
// // 	p2 := NewGF64Poly([]GF64{b4.combinations[3], GF64(1)})
// // 	p3 := NewGF64Poly([]GF64{b4.combinations[4], GF64(1)})
// // 	p4 := NewGF64Poly([]GF64{b4.combinations[5], GF64(1)})
// // 	p5 := NewGF64Poly([]GF64{b4.combinations[6], GF64(1)})
// // 	p6 := NewGF64Poly([]GF64{b4.combinations[7], GF64(1)})
// // 	p7 := NewGF64Poly([]GF64{b4.combinations[8], GF64(1)})

// // 	p0 = p0.mulNaive(p1)
// // 	p1 = p2.mulNaive(p3)
// // 	p2 = p4.mulNaive(p5)
// // 	p3 = p6.mulNaive(p7)

// // 	{
// // 		r := p0.mulNaive(p1).mulNaive(p2).mulNaive(p3)
// // 		r.debug("r")
// // 	}

// // 	p0.Expand(8)
// // 	p1.Expand(8)
// // 	p2.Expand(8)
// // 	p3.Expand(8)

// // 	err = b3.Mul(p0, p1)
// // 	if err != nil {
// // 		t.Fatal(err)
// // 	}
// // 	err = b3.Mul(p2, p3)
// // 	if err != nil {
// // 		t.Fatal(err)
// // 	}
// // 	p0.Expand(16)
// // 	p2.Expand(16)
// // 	err = b4.Mul(p0, p2)
// // 	if err != nil {
// // 		t.Fatal(err)
// // 	}
// // 	p0.debug("p0")
// // }

// // func TestSubspace(t *testing.T) {

// // 	n := 64
// // 	m := 4
// // 	// Find cantor bases where βm = 1.
// // 	cantor := make([]GF64, n)
// // 	cantor[0] = defaultBasisGenerator64
// // 	for j := 0; j < n-1; j++ {
// // 		cantor[j+1] = mul64(cantor[j], cantor[j]) ^ cantor[j]
// // 	}
// // 	if !cantor[n-1].IsOne() {
// // 		// never expected though
// // 		t.Fatalf("given initial β0 last basis expected to be 1")
// // 	}
// // 	cantor = cantor[64-m:]
// // 	l := uint(len(cantor))
// // 	combinations := make([]GF64, 1<<l)
// // 	// subCombinations := make([]GF64, 1<<(l-1))
// // 	// Calcucate combinations
// // 	for i := uint(0); i < l; i++ {
// // 		a := (1 << i)
// // 		for j := 0; j < a; j++ {
// // 			combinations[a+j] = combinations[j] ^ cantor[l-1-i]
// // 		}
// // 	}

// // 	// fmt.Println(len(combinations))
// // 	for i := 0; i < len(combinations); i++ {
// // 		fmt.Println(combinations[i])
// // 	}
// // 	fmt.Println(combinations[2] ^ combinations[7])
// // 	fmt.Println(combinations[9])
// // }

// // func TestLazyFFT64(t *testing.T) {
// // 	m := uint(16)
// // 	b := configureDefaultBasis64(m)
// // 	f0 := randGF64Poly(b.n)
// // 	f1 := f0.Clone()
// // 	_ = b.lFFT(f1)
// // 	_ = b.lIFFT(f1)
// // 	if !f0.EqualInCoeff(f1) {
// // 		t.Fatalf("lazy fft fails")
// // 	}
// // }

// // func TestConcurrentFFT64(t *testing.T) {
// // 	m := uint(16)
// // 	b := configureDefaultBasis64(m)
// // 	f0 := randGF64Poly(b.n)
// // 	f1 := f0.Clone()
// // 	_ = b.clFFT(f1)
// // 	_ = b.clIFFT(f1)
// // 	if !f0.EqualInCoeff(f1) {
// // 		t.Fatalf("lazy fft fails")
// // 	}
// // }

// // func BenchmarkRadix64(t *testing.B) {
// // 	m := polyLen
// // 	b := configureDefaultBasis64(m)
// // 	f0 := randGF64Poly(b.n)
// // 	t.ResetTimer()
// // 	for i := 0; i < t.N; i++ {
// // 		_ = b.radixConversion(f0)
// // 	}
// // }

// // func BenchmarkLazyParallelFFT64(t *testing.B) {
// // 	m := polyLen
// // 	b := configureDefaultBasis64(m)
// // 	f0 := randGF64Poly(b.n)
// // 	t.ResetTimer()
// // 	for i := 0; i < t.N; i++ {
// // 		_ = b.clFFT(f0)
// // 	}
// // }

// // func BenchmarkLazyFFT64(t *testing.B) {
// // 	m := polyLen
// // 	b := configureDefaultBasis64(m)
// // 	f0 := randGF64Poly(b.n)
// // 	t.ResetTimer()
// // 	for i := 0; i < t.N; i++ {
// // 		_ = b.lFFT(f0)
// // 	}
// // }
