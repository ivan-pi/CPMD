

find src -name "*.inc" -exec basename {} \; > srclist.txt
find src -name "*.F90" -exec basename {} \; >> srclist.txt
find src -name "*.f90" -exec basename {} \; >> srclist.txt

# Create CMakeLists.txt and add the required content
{
  echo "add_library(CPMD_core"
  sed 's/^/    /' srclist.txt  # This adds 4 spaces of indentation
  echo ")"
} > src/CMakeLists.txt


find src -name "*.cu" -exec basename {} \; > culist.txt
{
  echo "add_library(CPMD_cuda"
  sed 's/^/    /' culist.txt  # This adds 4 spaces of indentation
  echo ")"
} >> src/CMakeLists.txt