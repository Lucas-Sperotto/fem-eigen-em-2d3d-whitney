#include "fem3d1/fem3d1_driver_entry.hpp"

int main(int argc, char **argv)
{
  return run_fem3d1_case_driver(argc, argv, fem3d::CaseId::cyl, "fem3d1_cyl");
}
