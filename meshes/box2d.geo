lc = 0.005;
Point(1) = {0.0,0.0,0.0,lc};
Point(2) = {0.16,0.0,0.0,lc};
Point(3) = {0.16,0.04,0.0,lc};
Point(4) = {0,0.04,0.0,lc};
Point(5) = {0.02,0.02,0.0,lc};
//Point(6) = {0.0,0.02,0.0,lc};
//Point(7) = {0.08,0.02,0.0,lc};
Line(1) = {4,3};
Line(2) = {3,2};
Line(3) = {2,1};
Line(4) = {1,4};
//Line(7) = {6,7};
Line Loop(9) = {2,3,4,1};
Plane Surface(6) = {9};
//tmp[] = Extrude {0,0.0,0.01} {
//  Surface{6};
//};
//Physical Volume(1) = tmp[1];

Field[1] = Distance;
Field[1].NodesList = {5};

Field[2] = Threshold;
Field[2].IField = 1;
Field[2].LcMin = lc / 10;
Field[2].LcMax = lc;
Field[2].DistMin = 0.005;
Field[2].DistMax = 0.025;

Background Field = 2;

Mesh.CharacteristicLengthExtendFromBoundary = 0;
Mesh.CharacteristicLengthFromPoints = 0;
Mesh.CharacteristicLengthFromCurvature = 0;


