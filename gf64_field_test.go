package gf

import (
	"testing"
)

func TestGF64MultiplicationCrossAgainstNaive(t *testing.T) {
	for i := 0; i < 1000; i++ {
		a, b := randGF64(), randGF64()
		c0 := mul64(a, b)
		c1 := a.mulNaive(b)
		if !c0.Equal(c1) {
			t.Fatalf("a, b, c0, c1 %x, %x, %x, %x", a, b, c0, c1)
		}
	}
}

func TestGF64MultiplicationProperties(t *testing.T) {
	zero, one := GF64(0), GF64(1)
	a := randGF64()
	c0 := mul64(a, zero)
	if !c0.IsZero() {
		t.Fatalf("a * 0 == 0")
	}
	c0 = square64(zero)
	if !c0.IsZero() {
		t.Fatalf("0 ^ 2 == 0")
	}
	c1 := mul64(a, one)
	if !c1.Equal(a) {
		t.Fatalf("a * 1 == a")
	}
	c1 = square64(one)
	if !c1.IsOne() {
		t.Fatalf("1 ^ 2 == 1")
	}
	for i := 0; i < 1000; i++ {
		a := randGF64()
		b := randGF64()
		c0 = mul64(a, b)
		c1 = mul64(b, a)
		if !c0.Equal(c1) {
			t.Fatalf("a * b == b * a")
		}
		c := randGF64()
		c0 = mul64(mul64(a, b), c)
		c1 = mul64(mul64(a, c), b)
		if !c0.Equal(c1) {
			t.Fatalf("(a * b) * c == (a * c) * b")
		}
		c0 = mul64(a, a)
		c1 = square64(a)
		if !c0.Equal(c1) {
			t.Fatal("a^2 == a * a", i)
		}
	}
}

func TestGF64Inverse(t *testing.T) {
	one := GF64(1)
	for i := 0; i < 1000; i++ {
		a := randGF64()
		ai := a.Inverse()
		b := mul64(a, ai)
		if !b.Equal(one) {
			t.Fatalf("a * a ^ -1 == 1")
		}
		ione := one.Inverse()
		if !ione.Equal(one) {
			t.Fatalf("1 ^ -1 == 1")
		}
	}
}

func BenchmarkGF64Mul(t *testing.B) {
	r0, r1 := randGF64(), randGF64()
	t.ResetTimer()
	for i := 0; i < t.N; i++ {
		r0 = mul64(r0, r1)
	}
}

func BenchmarkGF64MulLamire(t *testing.B) {
	r0, r1 := randGF64(), randGF64()
	t.ResetTimer()
	for i := 0; i < t.N; i++ {
		r0 = lamire(r0, r1)
	}
}

func BenchmarkGF64MulAssign(t *testing.B) {
	r0, r1 := randGF64(), randGF64()
	t.ResetTimer()
	for i := 0; i < t.N; i++ {
		mulassign64(&r0, r1)
	}
}

func BenchmarkGF64Butterfly(t *testing.B) {
	r0, r1, g := randGF64(), randGF64(), randGF64()
	t.ResetTimer()
	for i := 0; i < t.N; i++ {
		butterfly(&r0, &r1, g)
	}
}
