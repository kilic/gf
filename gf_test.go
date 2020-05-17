package gf

import (
	"flag"
	"testing"
)

var polyLen uint

func TestMain(m *testing.M) {
	pl := flag.Int("pl", 16, "lenght of polynomial")
	flag.Parse()
	polyLen = uint(*pl)
	m.Run()
}
