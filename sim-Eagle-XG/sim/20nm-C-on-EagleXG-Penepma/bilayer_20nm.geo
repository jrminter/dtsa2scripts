XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
C
C  Thin film on substrate
C
0000000000000000000000000000000000000000000000000000000000000000
SURFACE (   1)   Plane Z=0
INDICES=( 0, 0, 0, 1, 0)
0000000000000000000000000000000000000000000000000000000000000000
SURFACE (   2)   Plane Z=-20 nm
INDICES=( 0, 0, 0, 1, 0)
Z-SHIFT=(-2.000000000000000E-06,   0)
0000000000000000000000000000000000000000000000000000000000000000
SURFACE (   3)   Plane Z=-0.1 cm
INDICES=( 0, 0, 0, 1, 0)
Z-SHIFT=(-1.000000000000000E-01,   0)
0000000000000000000000000000000000000000000000000000000000000000
SURFACE (   4)   Cylinder, 1 cm radius
INDICES=( 1, 1, 0, 0,-1)
X-SCALE=(+1.000000000000000E+00,   0)               (DEFAULT=1.0)
Y-SCALE=(+1.000000000000000E+00,   0)               (DEFAULT=1.0)
0000000000000000000000000000000000000000000000000000000000000000
MODULE  (   1)  Thin film
MATERIAL(   1)
SURFACE (   1), SIDE POINTER=(-1)
SURFACE (   2), SIDE POINTER=( 1)
SURFACE (   4), SIDE POINTER=(-1)
0000000000000000000000000000000000000000000000000000000000000000
MODULE  (   2)  Substrate
MATERIAL(   2)
SURFACE (   2), SIDE POINTER=(-1)
SURFACE (   3), SIDE POINTER=( 1)
SURFACE (   4), SIDE POINTER=(-1)
0000000000000000000000000000000000000000000000000000000000000000
END      0000000000000000000000000000000000000000000000000000000
