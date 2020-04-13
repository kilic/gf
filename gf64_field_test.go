package gf

import (
	// "fmt"
	"testing"
)

func TestGF64MultiplicationCrossAgainstNaive(t *testing.T) {
	for i := 0; i < 1000; i++ {
		a, b := randGF64(), randGF64()
		c0 := a.Mul(b)
		c1 := a.mulNaive(b)
		if !c0.Equal(c1) {
			t.Fatalf("a, b, c0, c1 %x, %x, %x, %x", a, b, c0, c1)
		}
	}
}

func TestGF64MultiplicationProperties(t *testing.T) {
	zero, one := GF64(0), GF64(1)
	for i := 0; i < 1000; i++ {
		a := randGF64()
		c0 := a.Mul(zero)
		if !c0.IsZero() {
			t.Fatalf("a * 0 == 0")
		}
		c1 := a.Mul(one)
		if !c1.Equal(a) {
			t.Fatalf("a * 1 == a")
		}
		b := randGF64()
		c0 = a.Mul(b)
		c1 = b.Mul(a)
		if !c0.Equal(c1) {
			t.Fatalf("a * b == b * a")
		}
		c := randGF64()
		c0 = a.Mul(b).Mul(c)
		c1 = a.Mul(c).Mul(b)
		if !c0.Equal(c1) {
			t.Fatalf("(a * b) * c == (a * c) * b")
		}
	}
}

func TestGF64Inverse(t *testing.T) {
	one := GF64(1)
	for i := 0; i < 1000; i++ {
		a := randGF64()
		ai := a.Inverse()
		b := a.Mul(ai)
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
		r0.Mul(r1)
	}
}
