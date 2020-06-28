package gf

import (
	"testing"
)

func TestRS(t *testing.T) {

	// Set the basis
	m := 16
	initDefaultBasis(m)

	// Generate some data
	data := randPoly(1 << 14)

	// Encode the data
	encodedData, err := encode(data, 2)

	if err != nil {
		t.Fatal(err)
	}

	// Make some missing points
	erasureData := encodedData.clone()
	missing := []uint64{2, 3}
	for _, i := range missing {
		erasureData.a[i] = 0
	}

	recoveredData, err := recover(erasureData, missing)
	if err != nil {
		t.Fatal(err)
	}

	// recoveredData.debug("recovered")
	// data.debug("data")
	if !recoveredData.equalInCoeff(data) {
		t.Fatal("rs recovery failed")
	}
}

func BenchmarkRSEncoding(t *testing.B) {
	m := polyLen
	initDefaultBasis(m)
	data := randPoly(1 << (m - 2))
	for i := 0; i < t.N; i++ {
		_, err := encode(data, 2)
		if err != nil {
			t.Fatal(err)
		}
	}
}

func BenchmarkRSDecoding(t *testing.B) {
	m := polyLen
	initDefaultBasis(m)
	data := randPoly(1 << (m - 2))
	encodedData, err := encode(data, 2)
	if err != nil {
		t.Fatal(err)
	}
	erasureData := encodedData.clone()
	missing := []uint64{2, 3}
	for _, i := range missing {
		erasureData.a[i] = 0
	}
	for i := 0; i < t.N; i++ {
		_, err := recover(erasureData, missing)
		if err != nil {
			t.Fatal(err)
		}
	}
}
