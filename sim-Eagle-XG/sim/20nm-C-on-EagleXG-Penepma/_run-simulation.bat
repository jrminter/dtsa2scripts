echo Running Penepma Simulation 20 nm C on Eagle XG glass 15 kV
echo Using composition from midpoints in US7935649 B2

date /T
time /T

penepma.exe < "20-nm-C-coated-EagleXG.in"

echo Simulation done...
date /T
time /T
pause
