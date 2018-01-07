//=============================================================================
//              Hilbert-curve (a space-filling Peano curve) library
//=============================================================================
#ifndef HILBERT_H
#define HILBERT_H

typedef unsigned int coord_t;	// char, short, int for upto 8, 16, 32 bits per word

#ifdef __cplusplus
extern "C" {
#endif

// Functions: LinetoAxes
//            AxestoLine
// Purpose:   Serial Hilbert length  <---->   multidimensional Axes position.
//   Space  = n-dimensional hypercube of side R = 2^b
//            Number of cells = N = R^n = 2^(n*b)
//
//   Line   = serial number of cell along Hilbert curve through hypercube
//          = extended integer of n*b bits ranging from 0 to N-1, stored as vector of n unsigned b-bit integers with [0] high.
//
//   Axes   = Geometrical position of cell
//          = n b-bit integers representing coordinates.
//
// Example:   side R = 16, dimension n = 2, number of cells = N = 256.
//            Line = 9, stored in base-16 words as
//                   Line[0] = 0 (high),   Line[1] = 9 (low),
//            corresponds to position (2,3) as in diagram, stored as
//                   Axes[0] = 2,   Axes[1] = 3.
// 
//        |
//     15 |    @---@   @---@   @---@   @---@   @---@   @---@   @---@   @---@
//        |    |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |
//        |    @   @---@   @   @   @---@   @   @   @---@   @   @   @---@   @
//        |    |           |   |           |   |           |   |           |
//        |    @---@   @---@   @---@   @---@   @---@   @---@   @---@   @---@
//        |        |   |           |   |           |   |           |   |    
//        |    @---@   @---@---@---@   @---@   @---@   @---@---@---@   @---@
//        |    |                           |   |                           |
//        |    @   @---@---@   @---@---@   @   @   @---@---@   @---@---@   @
//        |    |   |       |   |       |   |   |   |       |   |       |   |
// Axes[1]|    @---@   @---@   @---@   @---@   @---@   @---@   @---@   @---@
//        |            |           |                   |           |        
//        |    @---@   @---@   @---@   @---@   @---@   @---@   @---@   @---@
//        |    |   |       |   |       |   |   |   |       |   |       |   |
//        |    @   @---@---@   @---@---@   @---@   @---@---@   @---@---@   @
//        |    |                                                           |
//        |    @---@   @---@---@   @---@---@   @---@---@   @---@---@   @---@
//        |        |   |       |   |       |   |       |   |       |   |    
//        |    @---@   @---@   @---@   @---@   @---@   @---@   @---@   @---@
//        |    |           |           |           |           |           |
//        |    @   @---@   @   @---@   @---@   @---@   @---@   @   @---@   @
//        |    |   |   |   |   |   |       |   |       |   |   |   |   |   |
//        |    @---@   @---@   @   @---@---@   @---@---@   @   @---@   @---@
//        |                    |                           |                
//      3 |    5---6   9---@   @   @---@---@   @---@---@   @   @---@   @---@
//        |    |   |   |   |   |   |       |   |       |   |   |   |   |   |
//      2 |    4   7---8   @   @---@   @---@   @---@   @---@   @   @---@   @
//        |    |           |           |           |           |           |
//      1 |    3---2   @---@   @---@   @---@   @---@   @---@   @---@   @---@
//        |        |   |       |   |       |   |       |   |       |   |    
//      0 |    0---1   @---@---@   @---@---@   @---@---@   @---@---@   @--255
//        |
//         -------------------------------------------------------------------
//             0   1   2   3          ---> Axes[0]                         15
//
// Notes: (1) Unit change in Line yields single unit change in Axes position:
//            the Hilbert curve is maximally local.
//        (2) CPU proportional to total number of bits, = b * n.
//-----------------------------------------------------------------------------
	extern void LinetoAxes(coord_t*, coord_t*, int, int);	// multi-dimensional geometrical axes, linear serial #, # bits, dimension
	extern void AxestoLine(coord_t*, coord_t*, int, int);	// linear serial #, multi-dimensional geometrical axes, # bits, dimension


// Functions: LinetoTranspose
//            TransposetoLine
// Purpose:   Recover Hilbert integer by bit-transposition
// Example:   b=5 bits for each of n=3 coordinates
//               15-bit Hilbert integer = A B C D E a b c d e 1 2 3 4 5
//                                        X[0]..... X[1]..... X[2].....
//            transposed to
//               X[0](high) = A D b e 3
//               X[1]       = B E c 1 4
//               X[2](low)  = C a d 2 5
//                            high  low
//-----------------------------------------------------------------------------
	static void LinetoTranspose(coord_t*, coord_t*, int, int);	// transpose, Hilbert integer, # bits, dimension
	static void TransposetoLine(coord_t*, coord_t*, int, int);	// Hilbert integer, transpose, # bits, dimension


// Functions: TransposetoAxes
//            AxestoTranspose
// Purpose:   Transform in-place between Hilbert transpose and geometrical axes
// Example:   b=5 bits for each of n=3 coordinates.
//            15-bit Hilbert integer = A B C D E F G H I J K L M N O  is stored as its Transpose
//                   X[0] = A D G J M                X[2]|
//                   X[1] = B E H K N    <------->       | /X[1]
//                   X[2] = C F I L O               axes |/
//                          high  low                    0------ X[0]
//            Axes are stored conventially as b-bit integers.
//-----------------------------------------------------------------------------
	static void TransposetoAxes(coord_t*, int, int);	// position, # bits, dimension
	static void AxestoTranspose(coord_t*, int, int);	// position, # bits, dimension
#ifdef __cplusplus
};
#endif

#endif