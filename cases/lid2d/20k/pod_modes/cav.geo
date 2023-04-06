// Number of points and progressions
Nx1 = 25; Rx1 = 1.05;  // Cavity x-direction
Ny1 = 25; Ry1 = 1.05;   // Cavity y-direction

// Local densities??
lc1 = 1.0;
lc2 = 1.0;


// Boundary points
Point(1) = {-1.0,  1.0, 0, lc2};
Point(2) = {-1.0,  0.0, 0, lc2};
Point(3) = {-1.0, -1.0, 0, lc2};
Point(4) = { 0.0,  1.0, 0, lc2};
Point(5) = { 0.0,  0.0, 0, lc2};
Point(6) = { 0.0, -1.0, 0, lc2};
Point(7) = { 1.0,  1.0, 0, lc2};
Point(8) = { 1.0,  0.0, 0, lc2};
Point(9) = { 1.0, -1.0, 0, lc2};


Line(1) = {1, 2};  Transfinite Curve {1} = Ny1 Using Progression Ry1;  // Upper left cavity
Line(2) = {2, 3};  Transfinite Curve {2} = Ny1 Using Progression 1/Ry1;  // Lower left cavity
Line(3) = {4, 5};  Transfinite Curve {3} = Ny1 Using Progression Ry1;  // Upper mid cavity
Line(4) = {5, 6};  Transfinite Curve {4} = Ny1 Using Progression 1/Ry1;  // Lower mid cavity
Line(5) = {7, 8};  Transfinite Curve {5} = Ny1 Using Progression Ry1;  // Upper right cavity
Line(6) = {8, 9};  Transfinite Curve {6} = Ny1 Using Progression 1/Ry1;  // Lower right cavity
Line(7) = {1, 4};  Transfinite Curve {7} = Nx1 Using Progression Rx1;  // Top left cavity
Line(8) = {4, 7};  Transfinite Curve {8} = Nx1 Using Progression 1/Rx1;  // Top right cavity 
Line(9) = {2, 5};  Transfinite Curve {9} = Nx1 Using Progression Rx1;  // Mid left cavity
Line(10) = {5, 8};  Transfinite Curve {10} = Nx1 Using Progression 1/Rx1;  //  Mid right cavity
Line(11) = {3, 6};  Transfinite Curve {11} = Nx1 Using Progression Rx1;  // Bottom left cavity
Line(12) = {6, 9};  Transfinite Curve {12} = Nx1 Using Progression 1/Rx1;  // Bottom right cavity


Curve Loop(1) = {1, 9, -3, -7};
Plane Surface(1) = {1};
Curve Loop(2) = {2, 11, -4, -9};
Plane Surface(2) = {2};
Curve Loop(3) = {3, 10, -5, -8};
Plane Surface(3) = {3};
Curve Loop(4) = {4, 12, -6, -10};
Plane Surface(4) = {4};

Transfinite Surface {1};
Transfinite Surface {2};
Transfinite Surface {3};
Transfinite Surface {4};

Recombine Surface {1};
Recombine Surface {2};
Recombine Surface {3};
Recombine Surface {4};

Physical Line("Wall") = {1, 2, 11, 12, 6, 5};
Physical Line("Roof") = {7, 8};


Physical Surface("Fluid") = {1, 2, 3, 4};
