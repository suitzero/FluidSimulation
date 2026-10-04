const assert = require('assert');
const Multivector = require('../ga/ga.js');
const sim = require('../fluid_simulation/simulation.js');

function runTests() {
  console.log("Running simulation tests...");
  let testsPassed = 0;
  let testsTotal = 0;

  function runTest(name, testFn) {
    testsTotal++;
    try {
      testFn();
      console.log(`  [PASS] ${name}`);
      testsPassed++;
    } catch (e) {
      console.log(`  [FAIL] ${name}`);
      console.log(e);
    }
  }

  function createField(valFn) {
    let field = Array(sim.NX).fill(null).map(() => Array(sim.NY).fill(null));
    for(let r=0; r<sim.NX; r++){
      for(let c=0; c<sim.NY; c++){
        field[r][c] = valFn(r,c);
      }
    }
    return field;
  }

  runTest("Dx should compute forward difference at r=0", () => {
    const field = createField((r, c) => new Multivector({"1": r}));
    const res = sim.Dx(field, 0, 0, "1", sim.DX);
    assert.ok(Math.abs(res - 1/sim.DX) < 1e-6);
  });

  runTest("Dx should compute central difference in middle", () => {
    const field = createField((r, c) => new Multivector({"1": r*r}));
    const r = 2;
    const res = sim.Dx(field, r, 0, "1", sim.DX);
    assert.ok(Math.abs(res - 2*r/sim.DX) < 1e-6);
  });

  runTest("Dy should compute central difference in middle", () => {
    const field = createField((r, c) => new Multivector({"1": c*c}));
    const c = 2;
    const res = sim.Dy(field, 0, c, "1", sim.DY);
    assert.ok(Math.abs(res - 2*c/sim.DY) < 1e-6);
  });

  runTest("computeDivergence should compute correct value", () => {
    const field = createField((r, c) => new Multivector({"e1": r, "e2": c}));
    const div = sim.computeDivergence(field, 2, 2, sim.DX, sim.DY);
    assert.ok(Math.abs(div - (1/sim.DX + 1/sim.DY)) < 1e-6);
  });

  runTest("computePressureGradient should compute grad P", () => {
    const field = createField((r, c) => new Multivector({"1": r + 2*c}));
    const grad = sim.computePressureGradient(field, 2, 2, sim.DX, sim.DY);
    assert.ok(Math.abs(grad.coefficients.get("e1") - 1/sim.DX) < 1e-6);
    assert.ok(Math.abs(grad.coefficients.get("e2") - 2/sim.DY) < 1e-6);
  });
  
  runTest("computeVectorLaplacian should compute laplacian", () => {
    const field = createField((r, c) => new Multivector({"e1": r*r, "e2": c*c}));
    const lap = sim.computeVectorLaplacian(field, 2, 2, sim.DX, sim.DY);
    assert.ok(Math.abs(lap.coefficients.get("e1") - 2/(sim.DX*sim.DX)) < 1e-6);
    assert.ok(Math.abs(lap.coefficients.get("e2") - 2/(sim.DY*sim.DY)) < 1e-6);
  });

  runTest("advect should translate field based on velocity", () => {
    let field = createField((r, c) => new Multivector({"1": (r===5 && c===5) ? 1 : 0}));
    let vel = createField((r,c) => new Multivector({"e1": 1, "e2": 0}));
    let advected = sim.advect(field, vel, sim.DX);
    assert.ok(Math.abs(advected[6][5].coefficients.get("1") - 1) < 1e-6);
  });

  runTest("diffuse should apply viscosity", () => {
    let field = createField((r, c) => new Multivector({"e1": (r===2 && c===2) ? 1 : 0, "e2": 0}));
    let diffused = sim.diffuse(field, 1.0, 1.0);
    let val_center = diffused[2][2].coefficients.get("e1");
    assert.ok(val_center < 1);
    let val_right = diffused[3][2].coefficients.get("e1");
    assert.ok(val_right > 0);
  });

  runTest("project should calculate pressure to counteract divergence", () => {
    // If vel field has uniform divergence (like v = r e1), div = 1/dx
    // Project will compute a pressure field that causes a gradient that reduces vel.
    // We'll just verify that project returns an object with a pressure field and a projected velocity
    let vel = createField((r, c) => new Multivector({"e1": (r===10 && c===10) ? 1 : 0}));
    let result = sim.project(vel, 2); 
    
    assert.ok(result.projectedVelocity);
    assert.ok(result.pressureField);
    
    // Original velocity was strictly e1. Projected will have different values.
    // We just check it didn't explode or turn to NaN.
    let pv = result.projectedVelocity[10][10].coefficients.get("e1") || 0;
    assert.ok(!isNaN(pv));
  });

  console.log(`\nSimulation tests passed: ${testsPassed} / ${testsTotal}`);
  if (testsPassed < testsTotal) {
    process.exit(1);
  }
}

runTests();
