package gf

import (
	"flag"
	"testing"
)

var polyLen uint

func TestMain(m *testing.M) {
	pl := flag.Int("pl", 8, "lenght of polynomial")
	flag.Parse()
	polyLen = uint(*pl)
	m.Run()
}
