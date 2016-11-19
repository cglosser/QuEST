#include "bloch_system.h"

BlochSystem::BlochSystem(std::vector<QuantumDot> ds)
{
  dots.swap(ds);
}
