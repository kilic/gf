#include "textflag.h"

DATA P<>+0(SB)/8, $0x1b
GLOBL P<>(SB), RODATA|NOPTR, $8


TEXT ·mul64(SB), NOSPLIT, $0-24
  // move inputs to xmm registers
  MOVQ a+0(FP), X0
  MOVQ b+8(FP), X1
  // cH, cL = a * b
  PCLMULQDQ $0x00, X1, X0
  // get irreducible polynomial
  MOVQ P<>+0(SB), X2
  // dH, dL = cH * P
  VPCLMULQDQ $0x10, X0, X2, X3
  //  0, eL = dH * P
  PCLMULQDQ $0x10, X3, X2
  // cL + dL
  PXOR X0, X3
  // r = (cL + dL) + eL
  PXOR X3, X2
  // return result
  MOVQ X2, ret+16(FP)
  RET


TEXT ·mulassign64(SB), NOSPLIT, $0-16
  MOVQ a+0(FP), SI
  MOVQ (SI), X0
  MOVQ b+8(FP), X1
  PCLMULQDQ $0x00, X1, X0
  MOVQ P<>+0(SB), X2
  VPCLMULQDQ $0x10, X0, X2, X3
  PCLMULQDQ $0x10, X3, X2
  PXOR X0, X3
  PXOR X3, X2
  MOVQ X2, (SI)
  RET

TEXT ·square64(SB), NOSPLIT, $0-16
  MOVQ a+0(FP), X0
  PCLMULQDQ $0x00, X0, X0
  MOVQ P<>+0(SB), X1
  VPCLMULQDQ $0x10, X0, X1, X2
  PCLMULQDQ $0x10, X2, X1
  PXOR X0, X2
  PXOR X2, X1
  MOVQ X1, ret+8(FP)
  RET

TEXT ·squareassign64(SB), NOSPLIT, $0-0
  MOVQ a+0(FP), SI
  MOVQ (SI), X0
  PCLMULQDQ $0x00, X0, X0
  MOVQ P<>+0(SB), X1
  VPCLMULQDQ $0x10, X0, X1, X2
  PCLMULQDQ $0x10, X2, X1
  PXOR X0, X2
  PXOR X2, X1
  MOVQ X1, (SI)
  RET


TEXT ·butterfly(SB), NOSPLIT, $0-24
  MOVQ k2+8(FP), SI
  MOVQ (SI), X0
  MOVQ G+16(FP), X1
  PCLMULQDQ $0x00, X1, X0
  MOVQ P<>+0(SB), X2
  VPCLMULQDQ $0x10, X0, X2, X3
  PCLMULQDQ $0x10, X3, X2
  PXOR X0, X3
  PXOR X3, X2

  MOVQ X2, R8
  MOVQ k+0(FP), DX
  XORQ (DX), R8
  MOVQ R8, (DX)
  XORQ (SI), R8
  MOVQ R8, (SI)
  RET


TEXT ·ibutterfly(SB), NOSPLIT, $0-24
  MOVQ k2+8(FP), BX
  MOVQ k+0(FP), CX
  MOVQ (BX), R8 
  XORQ (CX), R8
  MOVQ R8, (BX)

  MOVQ R8, X0
  MOVQ G+16(FP), X1
  PCLMULQDQ $0x00, X1, X0
  MOVQ P<>+0(SB), X2
  VPCLMULQDQ $0x10, X0, X2, X3
  PCLMULQDQ $0x10, X3, X2
  PXOR X0, X3
  PXOR X3, X2

  MOVQ X2, R8
  XORQ (CX), R8
  MOVQ R8, (CX)
  RET



DATA LAMIRE<>+0(SB)/1, $0
DATA LAMIRE<>+1(SB)/1, $27
DATA LAMIRE<>+2(SB)/1, $54
DATA LAMIRE<>+3(SB)/1, $45
DATA LAMIRE<>+4(SB)/1, $108
DATA LAMIRE<>+5(SB)/1, $119
DATA LAMIRE<>+6(SB)/1, $90
DATA LAMIRE<>+7(SB)/1, $65
DATA LAMIRE<>+8(SB)/1, $216
DATA LAMIRE<>+9(SB)/1, $195
DATA LAMIRE<>+10(SB)/1, $238
DATA LAMIRE<>+11(SB)/1, $245
DATA LAMIRE<>+12(SB)/1, $180
DATA LAMIRE<>+13(SB)/1, $175
DATA LAMIRE<>+14(SB)/1, $130
DATA LAMIRE<>+15(SB)/1, $153
GLOBL LAMIRE<>(SB), RODATA|NOPTR, $128


// "Faster 64-bit universal hashing using carry-less multiplications"
// Daniel Lamire, Owen Kaser
// Algorithm 3
// https://arxiv.org/pdf/1503.03465.pdf
TEXT ·lamire(SB), NOSPLIT, $0-24
  MOVQ a+0(FP), X0
  MOVQ b+8(FP), X1
  PCLMULQDQ $0x00, X1, X0
  MOVQ P<>+0(SB), X1
  PCLMULQDQ $0x10, X0, X1
  VPSRLDQ $8, X1, X2
  VMOVDQU LAMIRE<>+0(SB), X3
  VPSHUFB X2, X3, X2
  PXOR X1, X0
  PXOR X2, X0
  MOVQ X0, ret+16(FP)
  RET
