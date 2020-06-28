package gf

import (
	"errors"
	"fmt"
)

func encode(data *poly, factor int) (*poly, error) {
	encodedData := data.clone()
	encodedData.expand(data.length() * factor)
	if _, err := encodedData.fft(); err != nil {
		return nil, err
	}
	return encodedData, nil
}

func recover(erasureData *poly, missing []uint64) (*poly, error) {

	m := erasureData.m()
	n := 1 << m

	I := make([]uint64, len(missing))
	for k := 0; k < len(missing); k++ {
		i := missing[k]
		if i >= uint64(erasureData.length()) {
			return nil, fmt.Errorf("missing data index %d exceeds erasure data length %d", i, erasureData.length())
		}
		if erasureData.a[i] != 0 {
			return nil, fmt.Errorf("erasure data at %d expected to be zero", missing[k])

		}
		I[k] = defaultBasis.combinations[i]
	}

	// Z(x)
	Zx, err := z(I)
	if err != nil {
		return nil, err
	}

	// Z'(x) = eval(Z(x))
	Zx.expand(n)
	Zxval, err := Zx.clone().fft()
	if err != nil {
		return nil, err
	}

	// DZ'(x) = Z'(x) * E'(x)
	DZxval, err := Zxval.mulSample(erasureData)
	if err != nil {
		return nil, err
	}

	// DZ(x) = interpolate(DZ'(x))
	DZx, err := DZxval.ifft()
	if err != nil {
		return nil, err
	}

	// pick random k
	// TODO: check if k is picked wisely
	k := randGF64()
	// Z(k * x)
	Zkx := Zx.substitute(k)
	// DZ(k * x)
	DZxk := DZx.substitute(k)

	// D(k * x) = DZ(k * x) / Z(k * x)
	Dxk, err := DZxk.div(Zkx)
	if err != nil {
		return nil, err
	}
	Dxk.trimZeros()
	if Dxk.length() > 1<<(m-1) {
		return nil, errors.New("perfect division is expected")
	}
	// D(x)
	Dx := Dxk.substitute(inverse(k))
	return Dx, nil

}
