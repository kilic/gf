## gf

__Work in progress.__

`gf` implements 64 bit binary field and polynomial operations.

## Run Tests and Benchmarks

``` bash
# run tests
go test -v

# run benchmarks
go test -run ^$ -bench=.
```

## References

* [Additive Fast Fourier Transforms over Finite Fields](http://www.math.clemson.edu/~sgao/papers/GM10.pdf)
* [McBits Revisited](https://tungchou.github.io/papers/mcbits_revisited.pdf)
* [Reversing a Finite Field Multiplication Optimization](https://blog.quarkslab.com/reversing-a-finite-field-multiplication-optimization.html#gk2010)
* [Faster 64-bit universal hashing using carry-less multiplications](https://arxiv.org/pdf/1503.03465.pdf)
