//+
SetFactory("OpenCASCADE");
R = 0.01;
Rsmall = 0.0095;
Z = 0.001;
Cylinder(1) = {0, 0, 0, Z, 0, 0, R, 2*Pi};

Cylinder(2) = {0, 0, 0.0, Z, 0, 0, Rsmall, 2*Pi};
Box(3) = {0.0, -0.02, -0.01, Z, 2*R, 2*R};

//+
BooleanDifference(4) = { Volume{2}; Delete; }{ Volume{3}; Delete; };
BooleanDifference(5) = { Volume{1}; Delete; }{ Volume{4}; Delete; };
//+
Physical Surface("endo", 9) = {4};
Physical Surface("endo", 9) += {5};
//+
Physical Surface("epi", 10) = {1};
//+
