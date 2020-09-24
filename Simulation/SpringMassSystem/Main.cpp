#include "SpringMassSystem.h"

int main()
{
	SpringMassSystem<2> spring_mass_system;
	spring_mass_system.Set_Particles({ 1.0 }, { {1.0,2.0 },{2.0,3.0 } });
	return 0;
}
