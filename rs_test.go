package gf

import (
	"testing"
)

func TestRS(t *testing.T) {

	// Set the basis
	initDefaultBasis(16)

	m := 14

	// Generate some data
	data := randPoly(1 << m)

	// Encode the data
	encodedData, err := encode(data, 2)

	if err != nil {
		t.Fatal(err)
	}

	// Make some missing points
	erasureData := encodedData.clone()

	missingDataSize := 1 << m
	missing := make([]uint64, missingDataSize)
	for i := 0; i < missingDataSize; i++ {
		missing[i] = 2 + uint64(i)
		erasureData.a[missing[i]] = 0
	}
	for _, i := range missing {
		erasureData.a[i] = 0
	}

	recoveredData, err := recover(erasureData, missing)
	if err != nil {
		t.Fatal(err)
	}

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
	// bench with half of the data is missing
	missingDataSize := 1 << (m - 2)
	missing := make([]uint64, missingDataSize)
	for i := 0; i < missingDataSize; i++ {
		missing[i] = 2 + uint64(i)
		erasureData.a[missing[i]] = 0
	}
	for _, i := range missing {
		erasureData.a[i] = 0
	}

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
