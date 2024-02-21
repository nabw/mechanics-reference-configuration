//+
SetFactory("OpenCASCADE");
Cylinder(1) = {0, 0, 0, 0.001, 0, 0, 0.01, 2*Pi};

Cylinder(2) = {0, 0, 0.0028, 0.001, 0, 0, 0.007, 2*Pi};

//+
BooleanDifference{ Volume{1}; Delete; }{ Volume{2}; Delete; }
//+
Physical Surface("endo", 9) = {4};
//+
Physical Surface("epi", 10) = {5};
