%after loading in the plan into CERR, do:

global planC;

aorta = 20;
infvenacav = 21;


%STRUCTURE MASKS: 0/1: 1 if voxel is in structure.
maskAorta = getUniformStr(aorta);
maskInfVenaCava = getUniformStr(infvenacav);



%DOSES

%dA = getUniformDose(doseNum, scanNum, structNum, planC)
d1a = getUniformDose(1,1,aorta,planC);
d1v = getUniformDose(1,1,infvenacav,planC);
%entire dose in region of interest (structure number 16)
dr = getUniformDose(1,1,16,planC);
%entire dose in patient
dp = getUniformDose(1,1,1,planC);




%now make up a velocity field in aorta and inf vena cava

%first aorta:
%define a u,v, and w for the vector field on the same grid as all
%the other data: 512x512x168

%u = 0*dr;v = u;w=u;

nx=512;
ny=512;
nz=168;

%for now say w is the main direction we want:
%u = .1*rand(nx,ny,nz);
%v = .1*rand(nx,ny,nz);
%w = 1+.2*rand(nx,ny,nz);

%rough better estimate is to go in this ration of z:y (or, w:v as
%the case may be:  1.45:0.65. this is the slope that the aorta
%descends back into the body in our region of interest - in
%direction from naval back towards back of neck)

u = .1*rand(nx,ny,nz);
v = 0.45 + .1*rand(nx,ny,nz);
w = -1.65 + .1*rand(nx,ny,nz);

%for inf vena cava, negate those .45 and -1.65.

%now zero out anything that is not in the aorta:
u = u.*maskAorta;
v = v.*maskAorta;
w = w.*maskAorta;

%x y z positions for plotting 
%x = 1:1:nx;
%y = 1:1:ny;
%z = 1:1:nz;
%[X,Y,Z] = meshgrid(x,y,z);


%trying to quiver3 plot on these data killed the system. instead
%let's load in each velocity component as a "dose cube" into cerr
%so I can at least view the components one at a time.

fractionGroupID = 'velocityAortaW';
assocScanNum    = 1;
showIMDose(w,fractionGroupID,assocScanNum)



%for inf vena cava, negate those
u = .1*rand(nx,ny,nz);
v = -0.45 + .1*rand(nx,ny,nz);
w = 1.65 + .1*rand(nx,ny,nz);
%now zero out anything that is not in the inf vena cava:
u = u.*maskInfVenaCava;
v = v.*maskInfVenaCava;
w = w.*maskInfVenaCava;

fractionGroupID = 'velocityInfVenaCavaW';
assocScanNum    = 1;
showIMDose(w,fractionGroupID,assocScanNum)