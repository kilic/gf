package gf

func log2Floor(n int) int {
	if n == 0 {
		panic("log 0")
	}
	r := 0
	for {
		n = n >> 1
		if n == 0 {
			return r
		}
		r += 1
	}
}

func log2Ceil(n int) int {
	if n == 0 {
		panic("log 0")
	}
	r := 0
	n = n - 1
	for {
		if n == 0 {
			return r
		}
		n = n >> 1
		r += 1
	}
}
