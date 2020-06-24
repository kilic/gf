package gf

import (
	"flag"
	"os"
	"testing"
)

var polyLen int

func TestMain(m *testing.M) {
	pl := flag.Int("pl", 16, "lenght of polynomial")
	flag.Parse()
	polyLen = *pl
	v := m.Run()
	os.Exit(v)
}
