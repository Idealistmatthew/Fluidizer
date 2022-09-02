#pragma once

#include <glm/vec3.hpp>
#include <glm/glm.hpp>

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
    glm::vec3* Velocity;

    // Conversion physical units and cell units
    glm::vec3 cell_to_dist_vector;
    glm::vec3 dist_to_cell_vector;

    int cellsize;
    int gridsize;
    glm::vec3* Position;

    // Helper functions to build initial conditions
    void AddVelocityUniform(glm::vec3 AddVel);

    // Advection Programme
    void Advect(float dt);

    // Advection Helper Functions
    bool within_bounds(glm::vec3 pos);
    glm::vec3 nearest_valid_pos(glm::vec3 pos);
    glm::vec3 interpolate_advecvel(glm::vec3 pos);

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
    cell_to_dist_vector = glm::vec3(dx, dy, dz);
    dist_to_cell_vector = glm::vec3(1 / dx, 1 / dy, 1 / dz);

    Pressure = new float[cellsize];
    memset(Pressure, 0, sizeof(float) * cellsize);
    Density = new float[cellsize];
    memset(Density, 0, sizeof(float) * cellsize);
    // The cellsize block of the position array stores the cell centers (for force and quantity advection)
    // The gridsize block of the position array stores the grid points (for vector advection)
    //Position = new glm::vec3[cellsize + gridsize];

    // Decided not to do staggered grid for now
    Position = new glm::vec3[cellsize];
    for (int i = 0; i < numX; i++) {
        for (int j = 0; j < numY; j++) {
            for (int k = 0; k < numZ; k++) {
                Position[INDEX(i, j, k)] = glm::vec3(i * dx, j * dy, k * dz);
            }
        }
    }

    Velocity = new glm::vec3[cellsize];
    for (int i = 0; i < cellsize; i++) {
        Velocity[i] = glm::vec3(0, 0, 0);
    }
    
    // Keeping for Staggered Grid implementation maybe
    /*xvel = new float[cellsize];
    memset(xvel, 0, sizeof(float) * cellsize);
    yvel = new float[cellsize];
    memset(yvel, 0, sizeof(float) * cellsize);
    zvel = new float[cellsize];
    memset(zvel, 0, sizeof(float) * cellsize);*/
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

void Mesh::AddVelocityUniform(glm::vec3 AddVel) {
    for (int i = 0; i < cellsize; i++) {
        Position[i] += AddVel;
    }
}

void Mesh::Advect(float dt) {
    for (int i = 0; i < numX; i++) {
        for (int j = 0; j < numY; j++) {
            for (int k = 0; k < numZ; k++) {
                glm::vec3 Curr_Position = Position[INDEX(i,j,k)];
                glm::vec3 Curr_Vel = Velocity[INDEX(i,j,k)];
                glm::vec3 Particle_Postemp;
                glm::vec3 Particle_Pos;
                glm::vec3 temp_Vel;
                glm::vec3 Particle_Vel;

                // I'm Using runge kutta btw Refer to the Bridson-Robert Book pg 32
                Particle_Postemp = Curr_Position - (float)dt* (Curr_Vel*dist_to_cell_vector)*0.5f;
                if (within_bounds(Particle_Postemp) == false) {
                    Particle_Postemp = nearest_valid_pos(Particle_Postemp);
                }
                temp_Vel = interpolate_advecvel(Particle_Postemp);
                Particle_Pos = Particle_Postemp - (float)dt * temp_Vel* dist_to_cell_vector ;
                if (within_bounds(Particle_Pos) == false) {
                    Particle_Pos = nearest_valid_pos(Particle_Pos);
                }
                Particle_Vel = interpolate_advecvel(Particle_Pos);
                Velocity[INDEX(i, j, k)] = Particle_Vel;
            }
        }
    }
}

bool Mesh::within_bounds(glm::vec3 Particle_Pos) {
    return Particle_Pos.x > numX || Particle_Pos.x < 0 || Particle_Pos.y > numY || Particle_Pos.y < 0 || Particle_Pos.z > numZ || Particle_Pos.z < 0;
}

glm::vec3 Mesh::nearest_valid_pos(glm::vec3 Pos) {
    glm::vec3 valid_pos;
    float x_pos;
    float y_pos;
    float z_pos;
    if (Pos.x > numX ) {
        x_pos = numX;
    }
    else if (Pos.x < 0) {
        x_pos = 0;
    }
    else {
        x_pos = Pos.x;
    }
    if (Pos.y > numY) {
        y_pos = numY;
    }
    else if (Pos.y < 0) {
        y_pos = 0;
    }
    else {
        y_pos = Pos.y;
    }
    if (Pos.z > numZ) {
        z_pos = numZ;
    }
    else if (Pos.z < 0) {
        z_pos = 0;
    }
    else {
        z_pos = Pos.z;
    }
    valid_pos = glm::vec3(x_pos, y_pos, z_pos);
    return valid_pos;
}

glm::vec3 Mesh::interpolate_advecvel(glm::vec3 pos) {
    glm::vec3 bottom_pos = glm::vec3(floor(pos.x), floor(pos.y), floor(pos.z));
    glm::vec3 top_pos = glm::vec3(ceil(pos.x), ceil(pos.y), ceil(pos.z));
    float inter_dist = glm::length(pos - bottom_pos);
    float cell_diag_dist = glm::length(top_pos - bottom_pos);
    glm::vec3 bottom_vel = Velocity[INDEX(int(bottom_pos.x), int(bottom_pos.y), int(bottom_pos.z))];
    glm::vec3 top_vel = Velocity[INDEX(int(top_pos.x), int(top_pos.y), int(top_pos.z))];
    glm::vec3 interp_vel = bottom_vel + (inter_dist / cell_diag_dist) * top_vel;
    return interp_vel;
}

void Mesh::Project() {

}