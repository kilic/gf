package gf

//go:noescape
func mul64(a, b GF64) GF64

//go:noescape
func mulassign64(a *GF64, b GF64)

//go:noescape
func square64(a GF64) GF64

//go:noescape
func squareassign64(a *GF64)

//go:noescape
func butterfly(k, k2 *GF64, G GF64)

//go:noescape
func ibutterfly(k, k2 *GF64, G GF64)

//go:noescape
func lamire(a, b GF64) GF64
