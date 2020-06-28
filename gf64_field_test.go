package gf

import (
	"testing"
)

func TestGF64MultiplicationCrossAgainstNaive(t *testing.T) {
	for i := 0; i < 1000; i++ {
		a, b := randGF64(), randGF64()
		c0 := mul64(a, b)
		c1 := mulNaive(a, b)
		if c0 != c1 {
			t.Fatalf("a, b, c0, c1 %x, %x, %x, %x", a, b, c0, c1)
		}
	}
}

func TestGF64MultiplicationProperties(t *testing.T) {
	var zero uint64 = 0
	var one uint64 = 1
	a := randGF64()
	c0 := mul64(a, zero)
	if c0 != 0 {
		t.Fatalf("a * 0 == 0")
	}
	c0 = square64(zero)
	if c0 != 0 {
		t.Fatalf("0 ^ 2 == 0")
	}
	c1 := mul64(a, one)
	if c1 != a {
		t.Fatalf("a * 1 == a")
	}
	c1 = square64(one)
	if c1 != 1 {
		t.Fatalf("1 ^ 2 == 1")
	}
	for i := 0; i < 1000; i++ {
		a := randGF64()
		b := randGF64()
		c0 = mul64(a, b)
		c1 = mul64(b, a)
		if c0 != c1 {
			t.Fatalf("a * b == b * a")
		}
		c := randGF64()
		c0 = mul64(mul64(a, b), c)
		c1 = mul64(mul64(a, c), b)
		if c0 != c1 {
			t.Fatalf("(a * b) * c == (a * c) * b")
		}
		c0 = mul64(a, a)
		c1 = square64(a)
		if c0 != c1 {
			t.Fatal("a^2 == a * a", i)
		}
	}
}

func TestGF64Inverse(t *testing.T) {
	var one uint64 = 1
	if one != inverse(one) {
		t.Fatalf("1 ^ -1 == 1")
	}
	for i := 0; i < 1000; i++ {
		a := randGF64()
		ai := inverse(a)
		b := mul64(a, ai)
		if b != one {
			t.Fatalf("a * a ^ -1 == 1")
		}
	}
}

func TestGF64Exp(t *testing.T) {
	a := randGF64()
	r := exp(a, 1)
	if r != a {
		t.Fatalf("a ^ 1 == a")
	}
	r = exp(a, 0)
	if r != 1 {
		t.Fatalf("a ^ 0 == 1")
	}
	b := mul64(a, a)
	b = mul64(b, a)
	r = exp(a, 3)
	if b != r {
		t.Fatalf("a * a * a == a ^ 3")
	}
	b = mul64(b, a)
	r = exp(a, 4)
	if b != r {
		t.Fatalf("a * a * a * a == a ^ 4")
	}
	var e0 uint64 = 1000
	var e1 uint64 = 2001
	r0 := exp(a, e0)
	r1 := exp(a, e1)
	r2 := exp(a, e0+e1)
	if r2 != mul64(r0, r1) {
		t.Fatalf("a ^ e0 * a ^ e1 == a ^ (e0 + e1)")
	}
}

func BenchmarkGF64Inverse(t *testing.B) {
	r0 := randGF64()
	t.ResetTimer()
	for i := 0; i < t.N; i++ {
		r0 = inverse(r0)
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
