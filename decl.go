package gf

//go:noescape
func mul64(a, b uint64) uint64

//go:noescape
func mulassign64(a *uint64, b uint64)

//go:noescape
func square64(a uint64) uint64

//go:noescape
func squareassign64(a *uint64)

//go:noescape
func butterfly(k, k2 *uint64, G uint64)

//go:noescape
func ibutterfly(k, k2 *uint64, G uint64)

//go:noescape
func lamire(a, b uint64) uint64
