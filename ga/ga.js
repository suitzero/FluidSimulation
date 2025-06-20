/**
 * @file Geometric Algebra Library Core
 * Provides the Multivector class for performing geometric algebra operations.
 * Assumes an orthonormal basis with e_i*e_i = 1 for all basis vectors e_i.
 * The current implementation is primarily suited for Euclidean geometric algebras.
 * Extending to Projective (PGA) or Conformal (CGA) geometric algebras would require
 * adjustments to the metric (e.g., e_i*e_i can be 0 or -1) in `_multiplyBasisBlades`.
 */

/**
 * Represents a multivector in a geometric algebra.
 * Coefficients are stored in a Map where keys are basis blade strings
 * (e.g., "1" for scalar, "e1", "e2", "e12" for vectors and bivectors)
 * and values are their numerical coefficients.
 * Assumes an orthonormal basis where e_i * e_i = 1 for underlying vector basis.
 * Blade names are strings like "1", "e1", "e2", "e12", where "e" is a prefix
 * and indices are sorted (e.g., "e12" not "e21").
 * @class Multivector
 */
class Multivector {
    /**
     * Creates an instance of a Multivector.
     * @param {Object<string, number>} [initialCoefficients={}] - An object mapping basis blade strings
     *                                                            to their coefficients.
     *                                                            Example: { "1": 1.0, "e1": 2.0, "e12": 3.0 }.
     *                                                            If empty or not provided, a zero multivector is created.
     */
    constructor(initialCoefficients = {}) {
        this.coefficients = new Map(Object.entries(initialCoefficients));
    }

    /**
     * Adds another multivector to this multivector.
     * @param {Multivector} other - The multivector to add.
     * @returns {Multivector} A new multivector representing the sum.
     */
    add(other) {
        const newCoefficients = new Map(this.coefficients);
        for (const [blade, coeff] of other.coefficients) {
            newCoefficients.set(blade, (newCoefficients.get(blade) || 0) + coeff);
        }
        // Remove zero coefficients for a cleaner representation and to simplify other operations.
        for (const [blade, coeff] of newCoefficients) {
            if (coeff === 0) {
                newCoefficients.delete(blade);
            }
        }
        return new Multivector(Object.fromEntries(newCoefficients));
    }

    /**
     * Subtracts another multivector from this multivector.
     * @param {Multivector} other - The multivector to subtract.
     * @returns {Multivector} A new multivector representing the difference.
     */
    subtract(other) {
        const newCoefficients = new Map(this.coefficients);
        for (const [blade, coeff] of other.coefficients) {
            newCoefficients.set(blade, (newCoefficients.get(blade) || 0) - coeff);
        }
        // Remove zero coefficients.
        for (const [blade, coeff] of newCoefficients) {
            if (coeff === 0) {
                newCoefficients.delete(blade);
            }
        }
        return new Multivector(Object.fromEntries(newCoefficients));
    }

    /**
     * Generates a string representation of the multivector.
     * Terms are sorted: scalar first, then by blade length (grade), then alphabetically.
     * Example: "1 + 2e1 - e12". A zero multivector is represented as "0".
     * @returns {string} A string representation of the multivector.
     */
    toString() {
        if (this.coefficients.size === 0) {
            return "0";
        }
        const terms = [];
        // Sort blades: scalar first, then by length (grade), then alphabetically for consistent output.
        const sortedBlades = Array.from(this.coefficients.keys()).sort((a, b) => {
            if (a === "1") return -1; // Scalar "1" always comes first.
            if (b === "1") return 1;
            if (a.length !== b.length) {
                return a.length - b.length; // Shorter blades (lower grade) first.
            }
            return a.localeCompare(b); // Alphabetical for same grade.
        });

        for (const blade of sortedBlades) {
            const coeff = this.coefficients.get(blade);
            // Coefficients that become zero (e.g. after addition/subtraction) should already be removed,
            // but this check is a safeguard.
            if (coeff === 0) continue;

            if (blade === "1") {
                terms.push(`${coeff}`);
            } else {
                // Handle coefficient display: 1 -> "", -1 -> "-", N -> N
                terms.push(`${coeff === 1 ? "" : coeff === -1 ? "-" : coeff}${blade}`);
            }
        }
        // Join terms with " + " and replace " + -" with " - " for cleaner output.
        return terms.join(" + ").replace(/ \+ -/g, " - ");
    }

    /**
     * Helper to compute the geometric product of two basis blades (e.g., "e1", "e2").
     * Assumes e_i*e_i = 1 (orthonormal Euclidean basis).
     * Handles sorting of indices and sign changes due to anti-commutation.
     * For example, e2 * e1 = -e12. e1 * e1 = 1.
     * @param {string} blade1 - String for the first basis blade (e.g., "1", "e1", "e12").
     * @param {string} blade2 - String for the second basis blade.
     * @returns {{coefficient: number, blade: string}} The resulting blade string (canonical form) and its coefficient multiplier (sign).
     * @private
     */
    _multiplyBasisBlades(blade1, blade2) {
        if (blade1 === "1") return { coefficient: 1, blade: blade2 };
        if (blade2 === "1") return { coefficient: 1, blade: blade1 };

        // Remove 'e' prefix to get indices, e.g., "e12" -> "12"
        let chars1 = blade1.substring(1);
        let chars2 = blade2.substring(1);

        let combinedChars = (chars1 + chars2).split(''); // E.g., "12" + "2" -> ['1','2','2']

        let sign = 1;

        // Sort characters to form canonical blade name and determine sign from swaps.
        // This is a bubble-sort style pass; for longer blade names, a more efficient sort might be considered.
        let i = 0;
        while (i < combinedChars.length) {
            let j = i + 1;
            while (j < combinedChars.length) {
                if (combinedChars[i] > combinedChars[j]) {
                    sign *= -1; // Each swap flips the sign.
                    [combinedChars[i], combinedChars[j]] = [combinedChars[j], combinedChars[i]];
                }
                j++;
            }
            i++;
        }

        // Annihilation: remove pairs of identical indices (e.g., e_i * e_i = 1).
        // This loop restarts the scan (i=0) after each annihilation because removal can create new adjacent pairs.
        i = 0;
        while (i < combinedChars.length - 1) {
            if (combinedChars[i] === combinedChars[i+1]) {
                combinedChars.splice(i, 2); // Remove the pair.
                // The sign change for bringing identical indices together (e.g. e2e1e2 -> e1e2e2)
                // is already handled by the sorting pass. Annihilation e_i*e_i=1 itself doesn't add a sign.
                i = 0; // Restart scan for further pairs (e.g. e1e2e1e2 -> e1e1e2e2 -> e2e2 -> 1)
            } else {
                i++;
            }
        }

        const finalBladeName = combinedChars.length > 0 ? 'e' + combinedChars.join('') : '1'; // If all indices annihilated, result is scalar "1".
        return { coefficient: sign, blade: finalBladeName };
    }

    /**
     * Computes the geometric product of this multivector with another.
     * @param {Multivector} other - The multivector to multiply with.
     * @returns {Multivector} A new multivector representing the geometric product.
     */
    geometricProduct(other) {
        const newCoefficients = new Map();

        for (const [blade1, coeff1] of this.coefficients) {
            for (const [blade2, coeff2] of other.coefficients) {
                const { coefficient: productSign, blade: newBlade } = this._multiplyBasisBlades(blade1, blade2);
                const totalProductCoeff = coeff1 * coeff2 * productSign;

                newCoefficients.set(newBlade, (newCoefficients.get(newBlade) || 0) + totalProductCoeff);
            }
        }

        for (const [blade, coeff] of newCoefficients) {
            if (coeff === 0) {
                newCoefficients.delete(blade);
            }
        }

        return new Multivector(Object.fromEntries(newCoefficients));
    }

    /**
     * Helper to get the grade of a basis blade string.
     * Grade is 0 for scalar "1", 1 for "ex", 2 for "exy", etc.
     * @param {string} blade - The basis blade string (e.g., "1", "e1", "e12").
     * @returns {number} The grade of the blade.
     * @private
     */
    _getBladeGrade(blade) {
        if (blade === '1') return 0;
        // Assumes blade is 'e' followed by unique, sorted digits, e.g., "e1", "e12".
        // The grade is the number of basis vectors in the blade, which is length of digits string.
        return blade.length - 1; // e.g. "e1" is length 2, grade 1. "e12" is length 3, grade 2.
    }

    /**
     * Helper to compute the outer (wedge) product of two basis blades.
     * Returns { coefficient: 0, blade: '1' } if blades share any common basis vector index.
     * Otherwise, concatenates indices, sorts them (adjusting sign for swaps), and returns the new blade.
     * @param {string} blade1 - String for the first basis blade.
     * @param {string} blade2 - String for the second basis blade.
     * @returns {{coefficient: number, blade: string}} The resulting blade and its coefficient (0, 1, or -1).
     * @private
     */
    _outerProductBasisBlades(blade1, blade2) {
        if (blade1 === '1') return { coefficient: 1, blade: blade2 }; // Scalar multiplication behavior for wedge product.
        if (blade2 === '1') return { coefficient: 1, blade: blade1 };

        const chars1 = blade1.substring(1).split('');
        const chars2 = blade2.substring(1).split('');

        // Check for common elements: if any index is shared, the outer product is zero.
        const set1 = new Set(chars1);
        for (const char of chars2) {
            if (set1.has(char)) {
                return { coefficient: 0, blade: '1' }; // Represent zero term, blade name is arbitrary but '1' is simple.
            }
        }

        // No common elements, combine and determine sign by sorting.
        const combinedChars = chars1.concat(chars2);
        let sign = 1;

        // Sort combinedChars to form canonical blade and count swaps for sign.
        for (let i = 0; i < combinedChars.length; i++) {
            for (let j = i + 1; j < combinedChars.length; j++) {
                if (combinedChars[i] > combinedChars[j]) {
                    sign *= -1;
                    [combinedChars[i], combinedChars[j]] = [combinedChars[j], combinedChars[i]]; // Swap for canonical name.
                }
            }
        }

        // Should not be empty if inputs are not "1" and no common elements.
        const finalBladeName = combinedChars.length > 0 ? 'e' + combinedChars.join('') : '1';
        return { coefficient: sign, blade: finalBladeName };
    }

    /**
     * Computes the outer (wedge) product of this multivector with another.
     * @param {Multivector} other - The multivector for the outer product.
     * @returns {Multivector} A new multivector representing the outer product.
     */
    outerProduct(other) {
        const newCoefficients = new Map();
        for (const [blade1, coeff1] of this.coefficients) {
            for (const [blade2, coeff2] of other.coefficients) {
                const { coefficient: productSign, blade: newBlade } = this._outerProductBasisBlades(blade1, blade2);
                // Only add if the product is not zero (outer product of blades with common factors is zero).
                if (productSign !== 0) {
                    const totalProductCoeff = coeff1 * coeff2 * productSign;
                    newCoefficients.set(newBlade, (newCoefficients.get(newBlade) || 0) + totalProductCoeff);
                }
            }
        }

        for (const [blade, coeff] of newCoefficients) {
            if (coeff === 0) {
                newCoefficients.delete(blade);
            }
        }
        return new Multivector(Object.fromEntries(newCoefficients));
    }

    /**
     * Computes the symmetric inner product of this multivector with another.
     * The symmetric inner product A . B is defined as the part of the geometric product AB
     * that has a grade equal to |grade(A) - grade(B)|.
     * @param {Multivector} other - The multivector for the inner product.
     * @returns {Multivector} A new multivector representing the symmetric inner product.
     */
    innerProduct(other) {
        const newCoefficients = new Map();
        for (const [blade1, coeff1] of this.coefficients) {
            for (const [blade2, coeff2] of other.coefficients) {
                // Calculate the geometric product of the individual basis blades.
                const gpResult = this._multiplyBasisBlades(blade1, blade2);

                const grade1 = this._getBladeGrade(blade1);
                const grade2 = this._getBladeGrade(blade2);
                const gradeGpBlade = this._getBladeGrade(gpResult.blade);

                // Filter condition for symmetric inner product.
                if (gradeGpBlade === Math.abs(grade1 - grade2)) {
                    const totalProductCoeff = coeff1 * coeff2 * gpResult.coefficient;
                    newCoefficients.set(gpResult.blade, (newCoefficients.get(gpResult.blade) || 0) + totalProductCoeff);
                }
            }
        }

        for (const [blade, coeff] of newCoefficients) {
            if (coeff === 0) {
                newCoefficients.delete(blade);
            }
        }
        return new Multivector(Object.fromEntries(newCoefficients));
    }

    /**
     * Extracts terms of a specific grade from this multivector.
     * For example, `mv.getGrade(0)` returns the scalar part, `mv.getGrade(1)` returns the vector part.
     * @param {number} grade - The grade to extract (e.g., 0 for scalar, 1 for vector, 2 for bivector).
     * @returns {Multivector} A new multivector containing only the terms of the specified grade.
     *                        Returns a zero multivector if no terms of that grade exist.
     */
    getGrade(grade) {
        const newCoefficients = new Map();
        for (const [blade, coeff] of this.coefficients) {
            if (this._getBladeGrade(blade) === grade) {
                newCoefficients.set(blade, coeff);
            }
        }
        return new Multivector(Object.fromEntries(newCoefficients));
    }
}

module.exports = Multivector;
