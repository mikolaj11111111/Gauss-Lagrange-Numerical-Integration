# Gauss-Legrange Numerical Integration

## Description
Implements Gauss-Legrange quadrature for numerical integration with composite subdivision. Compares accuracy of different node counts (2, 3, 4 nodes) and subdivision levels for complex functions.

## Features
- **Gauss-Legendre Quadrature** - High-precision numerical integration
- **Composite Integration** - Subdivisions for improved accuracy  
- **Multiple Node Orders** - 2, 3, and 4-point Gauss quadrature
- **Error Analysis** - Relative error calculation vs analytical values
- **Convergence Study** - Analysis of accuracy improvement with subdivisions
- **Data Export** - Results saved for plotting and analysis

## Mathematical Background
- **Gauss-Legrange**: Optimal quadrature using Legrange polynomial roots
- **Composite Method**: ∫f(x)dx = Σ∫f(x)dx over subintervals
- **Coordinate Transform**: Maps [a,b] to standard interval [-1,1]
- **Convergence**: Higher node count and more subdivisions → better accuracy

## Test Functions
1. **f₁(x) = x² · sin³(x)** on [1.0, 4.8], analytical = -10.9001
2. **f₂(x) = e^(x²) · (1-x)** on [-1.5, 3.2], analytical = -9364.62

## Usage
Program automatically tests subdivisions {1, 2, 4, 8, 16, 32} with 2, 3, and 4-node quadrature, showing convergence behavior and exporting results to "convergence_composite.txt".
