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
	a.MulAssign(zero)
	if !a.IsZero() {
		t.Fatal("a * 0 == 0")
	}

	c := new(GF64e)
	c.Set(a)
	a.MulAssign(one)
	if !c.Equal(a) {
		t.Fatal("a * 1 == a")
	}
	for i := 0; i < 1; i++ {
		a, b := randGF64e(), randGF64e()
		c0, c1 := new(GF64e), new(GF64e)
		c0.Set(a)
		c1.Set(b)
		c1.MulAssign(a)
		c0.MulAssign(b)
		if !c0.Equal(c1) {
			c0.debug("c0")
			c1.debug("c1")
			t.Fatal("a * b == b * a", i)
		}
		a, b = randGF64e(), randGF64e()
		c := randGF64e()
		c0.Set(a)
		c1.Set(a)
		c0.MulAssign(b)
		c1.MulAssign(c)
		c0.MulAssign(c)
		c1.MulAssign(b)
		if !c0.Equal(c1) {
			t.Fatal("(a * b) * c == (a * c) * b")
		}
	}
}
func BenchmarkGF64eMul(t *testing.B) {
	r0, r1 := randGF64e(), randGF64e()
	t.ResetTimer()
	for i := 0; i < t.N; i++ {
		r0.MulAssign(r1)
	}
}
