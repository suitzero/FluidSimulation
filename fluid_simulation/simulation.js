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

console.log("Fluid simulation placeholder loaded.");
console.log("Multivector class available:", Multivector ? "Yes" : "No");
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

function main_simulation_loop() {
    console.log("\n--- Simulation loop tick (testing derivatives) ---");
    // TODO: Implement advection, diffusion, pressure projection, velocity update

    if (NX > 5 && NY > 5) {
      const r_test = Math.floor(NX/2);
      const c_test = Math.floor(NY/2);

      console.log(`Initial velocity at (${r_test},${c_test}): ${vel[r_test][c_test].toString()}`);

      const gradV_test = computeGradV(vel, r_test, c_test, DX, DY);
      console.log(`GradV at (${r_test},${c_test}): ${gradV_test.toString()}`);
      // Expected for vel[r_test][c_test] = 0.5e1 + 0.3e2 and zeros elsewhere:
      // dvx_dx = (0 - 0) / (2*DX) = 0
      // dvy_dy = (0 - 0) / (2*DY) = 0
      // dvx_dy = (0 - 0) / (2*DY) = 0  (for v_x component, centered at 0.5)
      // dvy_dx = (0 - 0) / (2*DX) = 0  (for v_y component, centered at 0.3)
      // This will be zero if neighbors are zero. Let's try non-centered point or more complex field.
      // For the Dx/Dy, I've added boundary conditions, so at center it will be (0-0)/(2*DX) if neighbors are default Multivector()
      // The current initial condition is a single point of non-zero velocity.
      // So, Dx(vel, r_test, c_test, "e1", DX) will be (vel[r_test+1][c_test].e1 - vel[r_test-1][c_test].e1) / (2*DX) = (0-0)/(2*DX) = 0
      // And Dx(vel, r_test, c_test, "e2", DX) will be (vel[r_test+1][c_test].e2 - vel[r_test-1][c_test].e2) / (2*DX) = (0-0)/(2*DX) = 0
      // Similarly for Dy. So gradV will be 0.

      const laplacianV_test = computeVectorLaplacian(vel, r_test, c_test, DX, DY);
      console.log(`LaplacianV at (${r_test},${c_test}): ${laplacianV_test.toString()}`);
      // d2vx_dx2 = (0 - 2*0.5 + 0) / DX^2 = -1 / DX^2
      // d2vx_dy2 = (0 - 2*0.5 + 0) / DY^2 = -1 / DY^2
      // laplacian_vx = -1/DX^2 - 1/DY^2
      // d2vy_dx2 = (0 - 2*0.3 + 0) / DX^2 = -0.6 / DX^2
      // d2vy_dy2 = (0 - 2*0.3 + 0) / DY^2 = -0.6 / DY^2
      // laplacian_vy = -0.6/DX^2 - 0.6/DY^2
      // lap_vx = (-1/(0.05*0.05)) * 2 = -400 * 2 = -800 (if DX=DY=0.05)
      // lap_vy = (-0.6/(0.05*0.05)) * 2 = -240 * 2 = -480 (if DX=DY=0.05)

      const gradP_test = computePressureGradient(pressure, r_test, c_test, DX, DY);
      console.log(`GradP at (${r_test},${c_test}): ${gradP_test.toString()}`);
      // pressure[centerR+1][centerC] = 1, pressure[centerR-1][centerC] = -1. pressure[centerR][centerC] = 0.
      // dp_dx = ( P(r+1,c) - P(r-1,c) ) / (2*DX) = (1 - (-1)) / (2*DX) = 2 / (2*DX) = 1/DX = 1/0.05 = 20
      // dp_dy = ( P(r,c+1) - P(r,c-1) ) / (2*DY) = (0 - 0) / (2*DY) = 0
      // Expected: 20e1
    }
}
main_simulation_loop();

// Example usage:
// const vec = new Multivector({"e1": 1, "e2": 2});
// console.log("Vector:", vec.toString());
// const scalarField = new Multivector({"1": 10});
// console.log("Scalar field value:", scalarField.toString());
