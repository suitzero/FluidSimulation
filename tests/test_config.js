const assert = require('assert');
const sim = require('../fluid_simulation/simulation.js');

function runConfigTests() {
  console.log("\nRunning config tests...");
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
      console.error(e);
    }
  }

  runTest("Default configuration should initialize to 20x20", () => {
    assert.strictEqual(sim.NX, 20);
    assert.strictEqual(sim.NY, 20);
    assert.strictEqual(sim.config.NX, 20);
    assert.strictEqual(sim.config.NY, 20);
    assert.strictEqual(sim.DX, 1.0 / 20);
    assert.strictEqual(sim.DY, 1.0 / 20);
  });

  runTest("setConfig should correctly update dimensions and cell sizes", () => {
    sim.setConfig({ NX: 10, NY: 8 });
    assert.strictEqual(sim.NX, 10);
    assert.strictEqual(sim.NY, 8);
    assert.strictEqual(sim.DX, 1.0 / 10);
    assert.strictEqual(sim.DY, 1.0 / 8);
    assert.strictEqual(sim.config.NX, 10);
    assert.strictEqual(sim.config.NY, 8);
  });

  runTest("setConfig should reject negative values", () => {
    assert.throws(() => {
      sim.setConfig({ NX: -5 });
    }, /NX must be a positive integer/);
  });
  
  runTest("setConfig should reject zero", () => {
    assert.throws(() => {
      sim.setConfig({ NY: 0 });
    }, /NY must be a positive integer/);
  });

  runTest("setConfig should reject non-integers", () => {
    assert.throws(() => {
      sim.setConfig({ NX: 10.5 });
    }, /NX must be a positive integer/);
  });

  console.log(`\nConfig tests passed: ${testsPassed} / ${testsTotal}`);
  if (testsPassed < testsTotal) {
    process.exit(1);
  }
}

runConfigTests();
