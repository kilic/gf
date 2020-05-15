package gf

import (
	"fmt"
	"testing"
)

func (e *GF64e) debug(desc string) {
	fmt.Println(desc)
	for i := 0; i < len(e); i++ {
		fmt.Printf("%#4.4x\n", e[i])
	}
}

func TestGF64eMultiplicationProperties(t *testing.T) {
	zero, one := new(GF64e).Zero(), new(GF64e).One()
	a := randGF64e()
	a.Mul(zero)
	if !a.IsZero() {
		t.Fatal("a * 0 == 0")
	}

	c := new(GF64e)
	c.Set(a)
	a.Mul(one)
	if !c.Equal(a) {
		t.Fatal("a * 1 == a")
	}
	for i := 0; i < 1; i++ {
		a, b := randGF64e(), randGF64e()
		c0, c1 := new(GF64e), new(GF64e)
		c0.Set(a)
		c1.Set(b)
		c1.Mul(a)
		c0.Mul(b)
		if !c0.Equal(c1) {
			c0.debug("c0")
			c1.debug("c1")
			t.Fatal("a * b == b * a", i)
		}
		a, b = randGF64e(), randGF64e()
		c := randGF64e()
		c0.Set(a)
		c1.Set(a)
		c0.Mul(b)
		c1.Mul(c)
		c0.Mul(c)
		c1.Mul(b)
		if !c0.Equal(c1) {
			t.Fatal("(a * b) * c == (a * c) * b")
		}
	}
}
