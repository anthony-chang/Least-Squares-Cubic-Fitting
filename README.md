# Least-Squares-Cubic-Fitting
 Solve for the general polynomial (cubic) equation using the method of least squares. Created for an AP Physics C project in grade 12.

 This program fits a cubic function to the input data coordinate points. Let y be the vector with the y-coordinates. The outputs at each step:
 * the Vandermonde matrix from the x-coordinates, V
 * the transpose of the Vandermonde matrix, V<sup>T</sup>
 * V<sup>T</sup>V
 * (V<sup>T</sup>V)<sup>-1</sup>
 * V<sup>T</sup>y
 * (V<sup>T</sup>V)<sup>-1</sup>V<sup>T</sup>y, or the coefficient vector
