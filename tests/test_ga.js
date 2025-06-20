const assert = require('assert');
const Multivector = require('../ga/ga.js');

// Test Suite for Multivector operations in 2D Geometric Algebra
console.log("Running Multivector tests...");

// Setup: Define basic multivectors
const s = (val) => new Multivector({ "1": val });
const e1 = new Multivector({ "e1": 1 });
const e2 = new Multivector({ "e2": 1 });
const e12 = new Multivector({ "e12": 1 });

// toString() Tests
assert.strictEqual(new Multivector({ "e1": 2, "1": 1, "e12": -1 }).toString(), "1 + 2e1 - e12", "Test: toString() order and formatting");
assert.strictEqual(new Multivector().toString(), "0", "Test: toString() for zero multivector");
assert.strictEqual(s(5).toString(), "5", "Test: toString() for scalar");
assert.strictEqual(new Multivector({ "e1": 1 }).toString(), "e1", "Test: toString() for e1");
assert.strictEqual(new Multivector({ "e1": -1 }).toString(), "-e1", "Test: toString() for -e1");
assert.strictEqual(new Multivector({ "e12": 1, "e1": -2, "1": 3}).toString(), "3 - 2e1 + e12", "Test: toString() comprehensive");


// Addition Tests
const m1_add = new Multivector({ "1": 1, "e1": 2 });
const m2_add = new Multivector({ "e1": 1, "e2": 3 });
assert.strictEqual(m1_add.add(m2_add).toString(), "1 + 3e1 + 3e2", "Test: Addition m1+m2");
assert.strictEqual(e1.add(e1).toString(), "2e1", "Test: Addition e1+e1");
const m_add_c1 = new Multivector({"1": 1, "e1": 1});
const m_add_c2 = new Multivector({"1": -1, "e1": -1});
assert.strictEqual(m_add_c1.add(m_add_c2).toString(), "0", "Test: Addition resulting in zero");

// Subtraction Tests
const m1_sub = new Multivector({ "1": 1, "e1": 2 }); // Same as m1_add
const m2_sub = new Multivector({ "e1": 1, "e2": 3 }); // Same as m2_add
assert.strictEqual(m1_sub.subtract(m2_sub).toString(), "1 + e1 - 3e2", "Test: Subtraction m1-m2");
assert.strictEqual(e1.subtract(e1).toString(), "0", "Test: Subtraction e1-e1 (zero result)");
const m_sub_c1 = new Multivector({"1": 1, "e1": 1, "e2": 1});
const m_sub_c2 = new Multivector({"e1": 1, "e2": -1, "e12": 5});
assert.strictEqual(m_sub_c1.subtract(m_sub_c2).toString(), "1 + 2e2 - 5e12", "Test: Subtraction comprehensive");

// Geometric Product Tests (2D: e1*e1=1, e2*e2=1, e1*e2=e12, e2*e1=-e12)
assert.strictEqual(e1.geometricProduct(e1).toString(), "1", "Test: GP e1*e1");
assert.strictEqual(e2.geometricProduct(e2).toString(), "1", "Test: GP e2*e2");
assert.strictEqual(e1.geometricProduct(e2).toString(), "e12", "Test: GP e1*e2");
assert.strictEqual(e2.geometricProduct(e1).toString(), "-e12", "Test: GP e2*e1");
assert.strictEqual(e12.geometricProduct(e1).toString(), "-e2", "Test: GP e12*e1"); // e1e2e1 = -e1e1e2 = -e2
assert.strictEqual(e1.geometricProduct(e12).toString(), "e2", "Test: GP e1*e12"); // e1e1e2 = e2
assert.strictEqual(e12.geometricProduct(e12).toString(), "-1", "Test: GP e12*e12"); // e1e2e1e2 = -1
const m3_gp = new Multivector({ "1": 1, "e1": 1 }); // 1 + e1
const m4_gp = new Multivector({ "e2": 1, "e12": 1 }); // e2 + e12
// (1+e1)(e2+e12) = 1*e2 + 1*e12 + e1*e2 + e1*e12 = e2 + e12 + e12 + e2 = 2e2 + 2e12
assert.strictEqual(m3_gp.geometricProduct(m4_gp).toString(), "2e2 + 2e12", "Test: GP (1+e1)(e2+e12)");
const m5_gp = new Multivector({ "e1": 1, "e2": 1 }); // e1 + e2
const m6_gp = new Multivector({ "e1": 1, "e2": -1 }); // e1 - e2
// (e1+e2)(e1-e2) = e1*e1 + e2*e1 - e1*e2 - e2*e2 = 1 - e12 - e12 - 1 = -2e12
assert.strictEqual(m5_gp.geometricProduct(m6_gp).toString(), "-2e12", "Test: GP (e1+e2)(e1-e2)");


// Outer Product Tests
assert.strictEqual(e1.outerProduct(e1).toString(), "0", "Test: OP e1^e1");
assert.strictEqual(e1.outerProduct(e2).toString(), "e12", "Test: OP e1^e2");
assert.strictEqual(e2.outerProduct(e1).toString(), "-e12", "Test: OP e2^e1");
assert.strictEqual(e1.outerProduct(s(5)).toString(), "5e1", "Test: OP e1 ^ 5");
assert.strictEqual(s(5).outerProduct(e1).toString(), "5e1", "Test: OP 5 ^ e1");
assert.strictEqual(e1.outerProduct(e12).toString(), "0", "Test: OP e1^e12");
assert.strictEqual(e12.outerProduct(e1).toString(), "0", "Test: OP e12^e1");
const m3_op = new Multivector({ "1": 1, "e1": 1 }); // 1 + e1
const m4_op = new Multivector({ "e2": 1, "e12": 1 }); // e2 + e12
// (1+e1)^(e2+e12) = 1^e2 + 1^e12 + e1^e2 + e1^e12
// = e2 + e12 + e12 + 0 = e2 + 2e12
assert.strictEqual(m3_op.outerProduct(m4_op).toString(), "e2 + 2e12", "Test: OP (1+e1)^(e2+e12)");
const m5_op = new Multivector({ "1": 2, "e12": 3 }); // 2 + 3e12
const m6_op = new Multivector({ "e1": 4, "e2": 5 }); // 4e1 + 5e2
// (2+3e12)^(4e1+5e2) = 2^(4e1) + 2^(5e2) + 3e12^(4e1) + 3e12^(5e2)
// = 8e1 + 10e2 + 0 + 0 = 8e1 + 10e2
assert.strictEqual(m5_op.outerProduct(m6_op).toString(), "8e1 + 10e2", "Test: OP (2+3e12)^(4e1+5e2)");


// Inner Product Tests (Symmetric: (AB)_|gA-gB|)
assert.strictEqual(e1.innerProduct(e1).toString(), "1", "Test: IP e1.e1");
assert.strictEqual(e1.innerProduct(e2).toString(), "0", "Test: IP e1.e2");
assert.strictEqual(e12.innerProduct(e1).toString(), "-e2", "Test: IP e12.e1");
assert.strictEqual(e1.innerProduct(e12).toString(), "e2", "Test: IP e1.e12");
assert.strictEqual(e12.innerProduct(e12).toString(), "-1", "Test: IP e12.e12");
assert.strictEqual(s(5).innerProduct(e1).toString(), "5e1", "Test: IP 5.e1"); // GP is 5e1 (grade 1), |0-1|=1. So 5e1.
assert.strictEqual(e1.innerProduct(s(5)).toString(), "5e1", "Test: IP e1.5"); // GP is 5e1 (grade 1), |1-0|=1. So 5e1.
const m3_ip = new Multivector({ "1": 1, "e1": 1 }); // 1 + e1
const m4_ip = new Multivector({ "e2": 1, "e12": 1 }); // e2 + e12
// (1+e1).(e2+e12) = 1.e2 + 1.e12 + e1.e2 + e1.e12
// 1.e2 = <1*e2>_|0-1|=1 = <e2>_1 = e2
// 1.e12 = <1*e12>_|0-2|=2 = <e12>_2 = e12
// e1.e2 = <e1*e2>_|1-1|=0 = <e12>_0 = 0
// e1.e12 = <e1*e12>_|1-2|=1 = <e2>_1 = e2
// Sum = e2 + e12 + 0 + e2 = 2e2 + e12 (Sorted by toString: 2e2 + e12)
assert.strictEqual(m3_ip.innerProduct(m4_ip).toString(), "2e2 + e12", "Test: IP (1+e1).(e2+e12)");

// Test cases from prompt for IP e1.s(5) and s(5).e1
// e1 . 5: blade1="e1", coeff1=1; blade2="1", coeff2=5
// gpResult = _multiplyBasisBlades("e1", "1") = { coefficient: 1, blade: "e1" }
// grade1 = _getBladeGrade("e1") = 1
// grade2 = _getBladeGrade("1") = 0
// gradeGpBlade = _getBladeGrade("e1") = 1
// gradeGpBlade === Math.abs(grade1 - grade2) -> 1 === Math.abs(1 - 0) -> 1 === 1. True.
// productCoeff = coeff1 * coeff2 * gpResult.coefficient = 1 * 5 * 1 = 5
// newCoefficients.set("e1", 5)
// Result: 5e1. This matches the test output.

// s(5) . e1: blade1="1", coeff1=5; blade2="e1", coeff2=1
// gpResult = _multiplyBasisBlades("1", "e1") = { coefficient: 1, blade: "e1" }
// grade1 = _getBladeGrade("1") = 0
// grade2 = _getBladeGrade("e1") = 1
// gradeGpBlade = _getBladeGrade("e1") = 1
// gradeGpBlade === Math.abs(grade1 - grade2) -> 1 === Math.abs(0 - 1) -> 1 === 1. True.
// productCoeff = coeff1 * coeff2 * gpResult.coefficient = 5 * 1 * 1 = 5
// newCoefficients.set("e1", 5)
// Result: 5e1. The test output for s(5).innerProduct(e1) was "0", which was incorrect.
// My manual trace shows it should be 5e1. Correcting the assert.
assert.strictEqual(s(5).innerProduct(e1).toString(), "5e1", "Test: IP 5.e1 (Corrected)");


console.log("All GA tests passed!");
