#include "fem3d0/fem3d0_driver_entry.hpp"

int main(int argc, char **argv)
{
  return run_fem3d0_case_driver(argc, argv, fem3d::CaseId::half, "fem3d0_half");
}
