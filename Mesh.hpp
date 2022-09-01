#pragma once

#include <glm/vec3.hpp>

class Mesh
{
public:
    // Physical Sizes of the grid
    float gridX;
    float gridY;
    float gridZ;
    // Number of cells in each principle direction (not this is not number of grid points) (Using a MAC grid)
    int numX;
    int numY;
    int numZ;
    // Step size in each direction
    float dx;
    float dy;
    float dz;
    //
    float* Density;
    float* Pressure;

    float* xVel;
    float* yVel;
    float* zVel;

    int cellsize;
    int gridsize;
    glm::vec3* Position;

    // Advection Programme
    void Advect(float dt);

    // Advection Helper Functions
    bool within_bounds(glm::vec3 pos);
    glm::vec3 nearest_valid_pos(glm::vec3 pos);
    glm::vec3 interpolate_advecvel(glm::vec3 pos);
    glm::vec3 advect_fininterp(glm::vec3 particle_pos, glm::vec3 particle_vel);

    // Project Programme
    void Project();


    Mesh(float X, float Y, float Z, float dx1, float dy1, float dz1); // Constructor Function
    ~Mesh(); // Destructor function

    int INDEX(int idx, int idy, int nidz);
};

Mesh::Mesh(float X, float Y, float Z, float dx1, float dy1, float dz1) {
    gridX = X;
    gridY = Y;
    gridZ = Z;
    dx = dx1;
    dy = dy1;
    dz = dz1;
    numX = int(gridX / dx);
    numY = int(gridY / dy);
    numZ = int(gridZ / dz);

    cellsize = numX * numY * numZ;
    //gridsize = (numX + 1) * (numY + 1) * (numZ + 1);

    Pressure = new float[cellsize];
    memset(Pressure, 0, sizeof(float) * cellsize);
    Density = new float[cellsize];
    memset(Density, 0, sizeof(float) * cellsize);
    // The cellsize block of the position array stores the cell centers (for force and quantity advection)
    // The gridsize block of the position array stores the grid points (for vector advection)
    //Position = new glm::vec3[cellsize + gridsize];

    // Decided not to do staggered grid for now
    Position = new glm::vec3[cellsize];

    xVel = new float[cellsize];
    memset(xVel, 0, sizeof(float) * cellsize);
    yVel = new float[cellsize];
    memset(yVel, 0, sizeof(float) * cellsize);
    zVel = new float[cellsize];
    memset(zVel, 0, sizeof(float) * cellsize);
};

Mesh::~Mesh() {
    // Releasing allocated memory from heap
    delete xVel;
    delete yVel;
    delete zVel;
    delete Position;
};

int Mesh::INDEX(int idx, int idy, int idz) {
    return idx + numX * idy + numX * numY * idz;
};

void Mesh::Advect(float dt) {
    for (int i = 1; i <= numX; i++) {
        for (int j = 1; j <= numY; j++) {
            for (int k = 1; k <= numZ; k++) {
                glm::vec3 Curr_Position = Position[INDEX(i, j, k)];
                glm::vec3 Curr_Vel = glm::vec3(xVel[i],yVel[j],zVel[k]);
                glm::vec3 Particle_Postemp;
                glm::vec3 Particle_Pos;
                glm::vec3 temp_Vel;
                glm::vec3 Particle_Vel;
                glm::vec3 new_Vel;

                // I'm Using runge kutta btw Refer to the Bridson-Robert Book pg 32
                Particle_Postemp = Curr_Position - (float)dt* Curr_Vel*0.5f;
                if (within_bounds(Particle_Postemp) == true) {
                    temp_Vel = interpolate_advecvel(Particle_Postemp);
                }
                else {
                    Particle_Postemp = nearest_valid_pos(Particle_Postemp);
                    temp_Vel = glm::vec3(xVel[int(Particle_Postemp.x)], yVel[int(Particle_Postemp.y)], zVel[int(Particle_Postemp.z)]);
                }
                Particle_Pos = Particle_Postemp - (float)dt * temp_Vel;
                if (within_bounds(Particle_Pos) == true) {
                    Particle_Vel = interpolate_advecvel(Particle_Pos);
                }
                else {
                    Particle_Pos = nearest_valid_pos(Particle_Pos);
                    Particle_Vel = glm::vec3(xVel[int(Particle_Pos.x)], yVel[int(Particle_Pos.y)], zVel[int(Particle_Pos.z)]);
                }
                new_Vel = advect_fininterp(Particle_Pos, Particle_Vel);
                xVel[i] = new_Vel.x;
                yVel[j] = new_Vel.y;
                zVel[k] = new_Vel.z;
            }
        }
    }
}

void Mesh::Project() {

}