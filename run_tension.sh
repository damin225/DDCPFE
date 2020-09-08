# get the current working directory
cwd=$(pwd)
echo pwd: $cwd

## Run the tension test
ABAQUS="/home/damin/abaqus/Commands/abaqus"
rm *.lck
ulimit -s unlimited
$ABAQUS job=Direct_9grain user=umat.f -inter
