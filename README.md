# Multi-Step Methods for ODEs

Implementation of various multi-step methods for solving Ordinary Differential Equations (ODEs).

## Project Structure
```
├── src/
│   ├── methods/          # Implementation of numerical methods
│   │   ├── AB_AM_1.m    # First-order Adams-Bashforth-Adams-Moulton
│   │   ├── AB_AM_2.m    # Second-order Adams-Bashforth-Adams-Moulton
│   │   ├── BDF.m        # Backward Differentiation Formula
│   │   ├── FE.m         # Forward Euler
│   │   └── Prob4.m      # Problem 4 specific implementation
│   └── main.m           # Main driver script
```

## Methods Implemented
1. Adams-Bashforth-Adams-Moulton (ABAM) Methods
   - First-order ABAM predictor-corrector
   - Second-order ABAM predictor-corrector
2. Backward Differentiation Formula (BDF)
3. Forward Euler Method
4. Special handling for advection problems

## Problem Types
1. Linear ODE System
2. Predator-Prey System (Lotka-Volterra)
3. Van der Pol Oscillator
4. Advection Equation

## Author
Faranak Rajabi  
