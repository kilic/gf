package gf

import (
	"testing"
)

func TestGF16MultiplicationCrossAgainstNaive(t *testing.T) {
	for i := 0; i < 1000; i++ {
		a, b := randGF16(), randGF16()
		c0 := mul16(a, b)
		c1 := a.mulNaive(b)
		if !c0.Equal(c1) {
			t.Fatalf("a, b, c0, c1 %x, %x, %x, %x", a, b, c0, c1)
		}
	}
}

func TestGF16MultiplicationProperties(t *testing.T) {
	zero, one := GF16(0), GF16(1)
	c0 := mul16(randGF16(), zero)
	if !c0.IsZero() {
		t.Fatalf("a * 0 == 0")
	}
	a := randGF16()
	c1 := mul16(a, one)
	if !c1.Equal(a) {
		t.Fatalf("a * 1 == a")
	}
	for i := 0; i < 1000; i++ {
		a := randGF16()
		b := randGF16()
		c0 = mul16(a, b)
		c1 = mul16(b, a)
		if !c0.Equal(c1) {
			t.Fatalf("a * b == b * a")
		}
		c := randGF16()
		c0 = mul16(a, mul16(b, c))
		c1 = mul16(b, mul16(c, a))
		if !c0.Equal(c1) {
			t.Fatalf("(a * b) * c == (a * c) * b")
		}
	}
}
