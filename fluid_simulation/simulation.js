// Basic 2D Fluid Simulation using Geometric Algebra
const Multivector = require('../ga/ga.js');

// Grid Parameters
const NX = 20; // Number of cells in X (rows)
const NY = 20; // Number of cells in Y (columns)
const DX = 1.0 / NX; // Cell width
const DY = 1.0 / NY; // Cell height

// Fields Initialization
// Velocity field (vel[r][c] is a Multivector for velocity at cell (r,c))
let vel = Array(NX).fill(null).map(() => Array(NY).fill(null).map(() => new Multivector()));
// Pressure field (pressure[r][c] is a Multivector for pressure (scalar) at cell (r,c))
let pressure = Array(NX).fill(null).map(() => Array(NY).fill(null).map(() => new Multivector({ "1": 0 }))); // Initialize with zero pressure

// Example initial condition: set a small initial velocity at the center
if (NX > 5 && NY > 5) {
    const centerR = Math.floor(NX / 2);
    const centerC = Math.floor(NY / 2);
    vel[centerR][centerC] = new Multivector({ "e1": 0.5, "e2": 0.3 });
    // Add some pressure variation for testing gradient
    if (centerR + 1 < NX) pressure[centerR + 1][centerC] = new Multivector({ "1": 1.0 });
    if (centerR - 1 >= 0) pressure[centerR - 1][centerC] = new Multivector({ "1": -1.0 });
}

console.log("\nFluid simulation placeholder loaded.");
console.log("\nMultivector class available:", Multivector ? "Yes" : "No");
console.log(`Grid: ${NX}x${NY}, Cell size: ${DX.toFixed(3)}x${DY.toFixed(3)}`);

// Helper to safely get a component from a multivector in the field, handling boundaries (clamping)
function getComponent(field, r, c, bladeName) {
    const R = Math.max(0, Math.min(r, NX - 1)); // Clamp row index
    const C = Math.max(0, Math.min(c, NY - 1)); // Clamp col index
    return field[R][C].coefficients.get(bladeName) || 0;
}

// Finite difference for d(component)/dx at (r, c) using central differences
function Dx(field, r, c, bladeName, dx) {
    // For r=0, use forward difference. For r=NX-1, use backward difference. Otherwise central.
    if (r === 0 && NX > 1) { // Forward difference at left boundary
        return (getComponent(field, r + 1, c, bladeName) - getComponent(field, r, c, bladeName)) / dx;
    } else if (r === NX - 1 && NX > 1) { // Backward difference at right boundary
        return (getComponent(field, r, c, bladeName) - getComponent(field, r - 1, c, bladeName)) / dx;
    } else if (NX <=1) { // Not enough points for difference
        return 0;
    }
    // Central difference
    return (getComponent(field, r + 1, c, bladeName) - getComponent(field, r - 1, c, bladeName)) / (2 * dx);
}

// Finite difference for d(component)/dy at (r, c) using central differences
function Dy(field, r, c, bladeName, dy) {
    // For c=0, use forward difference. For c=NY-1, use backward difference. Otherwise central.
    if (c === 0 && NY > 1) { // Forward difference at top boundary
        return (getComponent(field, r, c + 1, bladeName) - getComponent(field, r, c, bladeName)) / dy;
    } else if (c === NY - 1 && NY > 1) { // Backward difference at bottom boundary
        return (getComponent(field, r, c, bladeName) - getComponent(field, r, c - 1, bladeName)) / dy;
    } else if (NY <= 1) { // Not enough points for difference
        return 0;
    }
    // Central difference
    return (getComponent(field, r, c + 1, bladeName) - getComponent(field, r, c - 1, bladeName)) / (2 * dy);
}

/**
 * Computes ∇v = (∂_x v_x + ∂_y v_y) + (∂_x v_y - ∂_y v_x)e_12 at cell (r,c).
 * This is the geometric product of the del operator ∇ = e1 ∂_x + e2 ∂_y with v = v_x e1 + v_y e2.
 * ∇v = (e1 ∂_x + e2 ∂_y)(v_x e1 + v_y e2)
 *    = e1 ∂_x (v_x e1) + e1 ∂_x (v_y e2) + e2 ∂_y (v_x e1) + e2 ∂_y (v_y e2)
 *    = (∂_x v_x) e1e1 + (∂_x v_y) e1e2 + (∂_y v_x) e2e1 + (∂_y v_y) e2e2
 *    = (∂_x v_x) + (∂_x v_y) e12 - (∂_y v_x) e12 + (∂_y v_y)
 *    = (∂_x v_x + ∂_y v_y) [scalar part, divergence]
 *    + (∂_x v_y - ∂_y v_x) e12 [bivector part, 2D curl component]
 * @param {Multivector[][]} velocityField - The 2D array of velocity multivectors.
 * @param {number} r - Row index.
 * @param {number} c - Column index.
 * @param {number} dx - Cell width.
 * @param {number} dy - Cell height.
 * @returns {Multivector} ∇v as a multivector (scalar part is divergence, e12 part is curl).
 */
function computeGradV(velocityField, r, c, dx, dy) {
    const dvx_dx = Dx(velocityField, r, c, "e1", dx);
    const dvy_dy = Dy(velocityField, r, c, "e2", dy);
    const dvx_dy = Dy(velocityField, r, c, "e1", dy); // Note: dvx/dy
    const dvy_dx = Dx(velocityField, r, c, "e2", dx); // Note: dvy/dx

    const divergence = dvx_dx + dvy_dy;       // Scalar part of ∇v
    const curl_component_e12 = dvx_dy - dvy_dx; // e12 part of ∇v

    return new Multivector({ "1": divergence, "e12": curl_component_e12 });
}

/**
 * Computes ∇p = (∂_x p)e1 + (∂_y p)e2 at cell (r,c) for a scalar pressure field p.
 * @param {Multivector[][]} pressureField - The 2D array of pressure multivectors (each holding a scalar in "1").
 * @param {number} r - Row index.
 * @param {number} c - Column index.
 * @param {number} dx - Cell width.
 * @param {number} dy - Cell height.
 * @returns {Multivector} ∇p as a vector multivector.
 */
function computePressureGradient(pressureField, r, c, dx, dy) {
    const dp_dx = Dx(pressureField, r, c, "1", dx);
    const dp_dy = Dy(pressureField, r, c, "1", dy);
    return new Multivector({ "e1": dp_dx, "e2": dp_dy });
}

/**
 * Computes ∇²v = (∇²v_x)e1 + (∇²v_y)e2 for a vector field v.
 * Where ∇²f = (d²f/dx² + d²f/dy²) is the Laplacian of a scalar function f.
 * d²f/dx² ≈ (f(i+1) - 2f(i) + f(i-1)) / dx² (Central difference for 2nd derivative)
 * @param {Multivector[][]} vectorField - The 2D array of vector multivectors.
 * @param {number} r - Row index.
 * @param {number} c - Column index.
 * @param {number} dx - Cell width.
 * @param {number} dy - Cell height.
 * @returns {Multivector} ∇²v as a vector multivector.
 */
function computeVectorLaplacian(vectorField, r, c, dx, dy) {
    // Laplacian for v_x component
    const v_x_center   = getComponent(vectorField, r, c, "e1");
    const v_x_plus_dx  = getComponent(vectorField, r + 1, c, "e1");
    const v_x_minus_dx = getComponent(vectorField, r - 1, c, "e1");
    const v_x_plus_dy  = getComponent(vectorField, r, c + 1, "e1");
    const v_x_minus_dy = getComponent(vectorField, r, c - 1, "e1");

    const d2vx_dx2 = (v_x_plus_dx - 2 * v_x_center + v_x_minus_dx) / (dx * dx);
    const d2vx_dy2 = (v_x_plus_dy - 2 * v_x_center + v_x_minus_dy) / (dy * dy);
    const laplacian_vx = d2vx_dx2 + d2vx_dy2;

    // Laplacian for v_y component
    const v_y_center   = getComponent(vectorField, r, c, "e2");
    const v_y_plus_dx  = getComponent(vectorField, r + 1, c, "e2");
    const v_y_minus_dx = getComponent(vectorField, r - 1, c, "e2");
    const v_y_plus_dy  = getComponent(vectorField, r, c + 1, "e2");
    const v_y_minus_dy = getComponent(vectorField, r, c - 1, "e2");

    const d2vy_dx2 = (v_y_plus_dx - 2 * v_y_center + v_y_minus_dx) / (dx * dx);
    const d2vy_dy2 = (v_y_plus_dy - 2 * v_y_center + v_y_minus_dy) / (dy * dy);
    const laplacian_vy = d2vy_dx2 + d2vy_dy2;

    return new Multivector({ "e1": laplacian_vx, "e2": laplacian_vy });
}


// Helper to perform bilinear interpolation of a multivector field
function interpolateBilinear(field, r_float, c_float) {
    r_float = Math.max(0, Math.min(r_float, NX - 1));
    c_float = Math.max(0, Math.min(c_float, NY - 1));

    const r0 = Math.floor(r_float);
    const r1 = Math.min(r0 + 1, NX - 1);
    const c0 = Math.floor(c_float);
    const c1 = Math.min(c0 + 1, NY - 1);

    const s1 = r_float - r0;
    const s0 = 1.0 - s1;
    const t1 = c_float - c0;
    const t0 = 1.0 - t1;

    const v00 = field[r0][c0];
    const v10 = field[r1][c0];
    const v01 = field[r0][c1];
    const v11 = field[r1][c1];

    const blades = new Set([
        ...v00.coefficients.keys(),
        ...v10.coefficients.keys(),
        ...v01.coefficients.keys(),
        ...v11.coefficients.keys()
    ]);

    const coeffs = {};
    for (const blade of blades) {
        const val00 = v00.coefficients.get(blade) || 0;
        const val10 = v10.coefficients.get(blade) || 0;
        const val01 = v01.coefficients.get(blade) || 0;
        const val11 = v11.coefficients.get(blade) || 0;

        const interpolated = s0 * (t0 * val00 + t1 * val01) + s1 * (t0 * val10 + t1 * val11);
        if (Math.abs(interpolated) > 1e-10) {
            coeffs[blade] = interpolated;
        }
    }

    return new Multivector(coeffs);
}

// Semi-Lagrangian Advection
function advect(field, velocityField, dt) {
    const newField = Array(NX).fill(null).map(() => Array(NY).fill(null));

    for (let r = 0; r < NX; r++) {
        for (let c = 0; c < NY; c++) {
            const v_x = getComponent(velocityField, r, c, "e1");
            const v_y = getComponent(velocityField, r, c, "e2");

            // Backtrace
            const prev_r = r - v_x * dt / DX;
            const prev_c = c - v_y * dt / DY;

            newField[r][c] = interpolateBilinear(field, prev_r, prev_c);
        }
    }
    return newField;
}


// Explicit Diffusion (Viscosity)
function diffuse(field, nu, dt) {
    const newField = Array(NX).fill(null).map(() => Array(NY).fill(null));

    for (let r = 0; r < NX; r++) {
        for (let c = 0; c < NY; c++) {
            const laplacian = computeVectorLaplacian(field, r, c, DX, DY);

            const currentCoeffs = new Map(field[r][c].coefficients);
            for (const [blade, coeff] of laplacian.coefficients) {
                currentCoeffs.set(blade, (currentCoeffs.get(blade) || 0) + nu * dt * coeff);
            }

            // Clean up small values
            const finalCoeffs = {};
            for (const [blade, coeff] of currentCoeffs) {
                if (Math.abs(coeff) > 1e-10) {
                    finalCoeffs[blade] = coeff;
                }
            }
            newField[r][c] = new Multivector(finalCoeffs);
        }
    }
    return newField;
}


// Computes divergence of vector field (e1 and e2 parts)
function computeDivergence(vectorField, r, c, dx, dy) {
    const dvx_dx = Dx(vectorField, r, c, "e1", dx);
    const dvy_dy = Dy(vectorField, r, c, "e2", dy);
    return dvx_dx + dvy_dy;
}

// Pressure Projection using Jacobi iteration for Poisson equation: ∇²p = ∇ · v
function project(velocityField, iterations = 20) {
    // Calculate divergence of velocity
    const div = Array(NX).fill(null).map(() => Array(NY).fill(0));
    for (let r = 0; r < NX; r++) {
        for (let c = 0; c < NY; c++) {
            div[r][c] = computeDivergence(velocityField, r, c, DX, DY);
        }
    }

    // Initialize pressure to 0
    let p = Array(NX).fill(null).map(() => Array(NY).fill(0));
    let p_new = Array(NX).fill(null).map(() => Array(NY).fill(0));

    // Jacobi iteration
    for (let iter = 0; iter < iterations; iter++) {
        for (let r = 0; r < NX; r++) {
            for (let c = 0; c < NY; c++) {
                const p_left = r > 0 ? p[r - 1][c] : 0;
                const p_right = r < NX - 1 ? p[r + 1][c] : 0;
                const p_up = c > 0 ? p[r][c - 1] : 0;
                const p_down = c < NY - 1 ? p[r][c + 1] : 0;

                // p[r][c] = (div[r][c] * dx * dy - (p_left + p_right + p_up + p_down)) / -4
                // Let's assume dx == dy for simplicity, or we can write the exact form:
                // d2p/dx2 + d2p/dy2 = div
                // (p(r+1,c) - 2p(r,c) + p(r-1,c))/dx^2 + (p(r,c+1) - 2p(r,c) + p(r,c-1))/dy^2 = div(r,c)
                const denom = 2 / (DX * DX) + 2 / (DY * DY);
                const sum = (p_right + p_left) / (DX * DX) + (p_down + p_up) / (DY * DY) - div[r][c];
                p_new[r][c] = sum / denom;
            }
        }
        // Swap p and p_new
        const temp = p;
        p = p_new;
        p_new = temp;
    }

    // Subtract pressure gradient from velocity
    const newVelocity = Array(NX).fill(null).map(() => Array(NY).fill(null));

    // Create a MultiVector field for pressure to reuse computePressureGradient
    const pressureField = Array(NX).fill(null).map((_, r) =>
        Array(NY).fill(null).map((_, c) => new Multivector({"1": p[r][c]}))
    );

    for (let r = 0; r < NX; r++) {
        for (let c = 0; c < NY; c++) {
            const gradP = computePressureGradient(pressureField, r, c, DX, DY);
            newVelocity[r][c] = velocityField[r][c].subtract(gradP);
        }
    }

    return { projectedVelocity: newVelocity, pressureField };
}



function main_simulation_loop() {
    console.log("\n--- Simulation loop tick ---");
    const dt = 0.01;
    const nu = 0.001; // Kinematic viscosity

    // Step 1: Advection
    vel = advect(vel, vel, dt);

    // Step 2: Diffusion
    vel = diffuse(vel, nu, dt);

    // Step 3: Projection (Enforce incompressibility)
    const result = project(vel, 20);
    vel = result.projectedVelocity;
    pressure = result.pressureField;

    if (NX > 5 && NY > 5) {
      const r_test = Math.floor(NX/2);
      const c_test = Math.floor(NY/2);
      console.log(`Velocity at (${r_test},${c_test}): ${vel[r_test][c_test].toString()}`);
      console.log(`Pressure at (${r_test},${c_test}): ${pressure[r_test][c_test].toString()}`);
    }
}

main_simulation_loop();

// Example usage:
// const vec = new Multivector({"e1": 1, "e2": 2});
// console.log("\nVector:", vec.toString());
// const scalarField = new Multivector({"1": 10});
// console.log("\nScalar field value:", scalarField.toString());

if (typeof module !== 'undefined') {
  module.exports = {
    getComponent,
    Dx,
    Dy,
    computeGradV,
    computePressureGradient,
    computeVectorLaplacian,
    interpolateBilinear,
    advect,
    diffuse,
    computeDivergence,
    project,
    NX,
    NY,
    DX,
    DY
  };
}
