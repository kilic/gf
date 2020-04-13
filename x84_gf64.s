#include "textflag.h"

DATA P<>+0(SB)/8, $0x1b
GLOBL P<>(SB), RODATA|NOPTR, $8


TEXT ·mul64(SB), NOSPLIT, $0-24
  // move inputs to xmm registerss
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
  // move result to first operand
  MOVQ X2, ret+16(FP)
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
