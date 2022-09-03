/*
Current Simplifications:
No Obstacles
Incompressible Flow
External Boundary Conditions
Simple 3D Visualisation (haven't figured out yet)
*/
#include <Mesh.hpp>

int main() {
	// Set Constants and Parameters (Or use separate header file to define constants)
	float X = 10;
	float Y = X;
	float Z = X;
	float dx = X / 10;
	float dy = Y / 10;
	float dz = Z / 10;
	float dt = 0.5;
	float sim_T = 20;
	// Initialise Mesh
	Mesh First_Mesh = Mesh(X, Y, Z, dx, dy, dz);
	First_Mesh.Project(dt);
	// Set boundary Conditions
	First_Mesh.AddVelocityUniform( glm::vec3(0.5, 0.5, 0.5) );
	for (int t = 0; t < int(sim_T / dt); t++) {
		First_Mesh.Advect(dt);
		First_Mesh.Project(dt);
	}
	// Set Simulation Parameters
	//Start Loop

	//Diffuse

	//Project (Conserve Density)

	//Advect

	//Project

	//Data Visualisation
};