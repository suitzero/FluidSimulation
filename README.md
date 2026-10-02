# Fluid Simulation

This project is a 2D fluid simulation built in JavaScript that leverages Geometric Algebra for its mathematical foundation. It features a custom `Multivector` class for performing Euclidean geometric operations, such as geometric, outer, and inner products. The simulation models fluid dynamics through Semi-Lagrangian advection, explicit diffusion, and pressure projection to maintain incompressibility on a grid.

## Running Tests

You can run the test suite for the Geometric Algebra multivector operations using the following command:

```bash
node tests/test_ga.js
```

## Open Tasks

- Add a `package.json` file to manage project metadata and execution scripts.
- Add unit tests for `fluid_simulation/simulation.js` (currently only `ga.js` is covered by tests).
- Refactor the hardcoded grid dimensions (`NX`, `NY`) and simulation parameters in `fluid_simulation/simulation.js` into configurable options.
