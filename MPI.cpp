#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <cmath>
#include <unistd.h>
#include <cfloat>

#define OCEAN_SIZE 100
#define MREQUIN 10
#define MPOISSON 5
#define PREP 0.1
#define EMPTY -1
#define MAX_ANIMALS 1000
#define HUNGER_LIMIT 10

#define REPPOISSON 1.0         
#define ATTRSHARK 1.0       
#define ATTRSHARK_CLOSEST 2.0 
#define VISIBILITY_RANGE 20 

typedef struct {
    int type; // 0 pour poisson, 1 pour requin
    float x, y; 
    float vx, vy; 
    float ax, ay; 
    int hunger; 
} Animal;

Animal ocean[MAX_ANIMALS];

void initializeOcean() {
    srand(time(NULL));
    for (int i = 0; i < MAX_ANIMALS; i++) {
        ocean[i].type = rand() % 2;
        ocean[i].x = rand() % OCEAN_SIZE;
        ocean[i].y = rand() % OCEAN_SIZE;
        ocean[i].vx = ocean[i].vy = 0;
        ocean[i].ax = ocean[i].ay = 0;
        ocean[i].hunger = 0;
    }
}

void updateForces(Animal* a) {
    float force_x = 0.0, force_y = 0.0;

    force_x += 0.1; 

    for (int i = 0; i < MAX_ANIMALS; i++) {
        float distance = sqrt(pow(ocean[i].x - a->x, 2) + pow(ocean[i].y - a->y, 2));

        if (a->type == 1 && ocean[i].type == 0) {
            if (distance < VISIBILITY_RANGE) {
                force_x += ATTRSHARK_CLOSEST / pow(distance, 2);
                force_y += ATTRSHARK_CLOSEST / pow(distance, 2);
            } else {
                force_x += ATTRSHARK / distance;
                force_y += ATTRSHARK / distance;
            }
        }

        if (a->type == 0 && ocean[i].type == 1) {
            force_x -= REPPOISSON / distance;
            force_y -= REPPOISSON / distance;
        }
    }

    a->ax = force_x / (a->type == 0 ? MPOISSON : MREQUIN);
    a->ay = force_y / (a->type == 0 ? MPOISSON : MREQUIN);
}


void handleCollisionsAndReproduction() {
    for (int i = 0; i < MAX_ANIMALS; i++) {
        for (int j = i + 1; j < MAX_ANIMALS; j++) {
            float distance = sqrt(pow(ocean[i].x - ocean[j].x, 2) + pow(ocean[i].y - ocean[j].y, 2));
            if (distance < 1.0) { 
                if (ocean[i].type == ocean[j].type) {
                    float temp_vx = ocean[i].vx;
                    float temp_vy = ocean[i].vy;
                    ocean[i].vx = ocean[j].vx;
                    ocean[i].vy = ocean[j].vy;
                    ocean[j].vx = temp_vx;
                    ocean[j].vy = temp_vy;
                }
                else {
                    // Shark eats fish
                    if (ocean[i].type == 1) {
                        ocean[j].type = EMPTY; 
                        ocean[i].hunger = 0; 
                    }
                    else {
                        ocean[i].type = EMPTY; 
                        ocean[j].hunger = 0; 
                    }
                }
            }
        }
    }

    for (int i = 0; i < MAX_ANIMALS; i++) {
        if (ocean[i].type != EMPTY) {
            if (rand() < PREP * RAND_MAX) {
                for (int k = 0; k < MAX_ANIMALS; k++) {
                    if (ocean[k].type == EMPTY) {
                        ocean[k].type = ocean[i].type;
                        ocean[k].x = ocean[i].x + ((rand() % 3) - 1);
                        ocean[k].y = ocean[i].y + ((rand() % 3) - 1);
                        ocean[k].vx = ocean[k].vy = 0;
                        ocean[k].ax = ocean[k].ay = 0;
                        ocean[k].hunger = 0;
                        break;
                    }
                }
            }

            if (ocean[i].type == 1) {
                ocean[i].hunger++;
                if (ocean[i].hunger > HUNGER_LIMIT) {
                    ocean[i].type = EMPTY; 
                }
            }
        }
    }
}

void updatePosition(Animal* a, float timeStep) {
    a->vx += a->ax * timeStep;
    a->vy += a->ay * timeStep;
    a->x += a->vx * timeStep;
    a->y += a->vy * timeStep;
}

void exchangeAnimals(int world_rank, int world_size, int count, Animal* buffer, int* num_received) {
    MPI_Status status;

    if (world_size > 1) {
        if (world_rank == 0) {
            MPI_Send(buffer, count * sizeof(Animal), MPI_BYTE, 1, 0, MPI_COMM_WORLD);
            MPI_Probe(1, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            MPI_Get_count(&status, MPI_BYTE, num_received);
            MPI_Recv(buffer, *num_received, MPI_BYTE, 1, 0, MPI_COMM_WORLD, &status);
        } else if (world_rank == 1) {
            MPI_Recv(buffer, MAX_ANIMALS * sizeof(Animal), MPI_BYTE, 0, 0, MPI_COMM_WORLD, &status);
            MPI_Send(buffer, count * sizeof(Animal), MPI_BYTE, 0, 0, MPI_COMM_WORLD);
        }
    }
}

void processReceivedAnimals(Animal* buffer, int numAnimals) {
    for (int i = 0; i < numAnimals; i++) {
        for (int j = 0; j < MAX_ANIMALS; j++) {
            if (ocean[j].type == EMPTY) {
                ocean[j] = buffer[i];
                break;
            }
        }
    }
}

void printOcean(Animal ocean[], int oceanSize) {
    printf("  ");
    for (int j = 0; j < oceanSize; ++j) {
        printf("%2d", j);
    }
    printf("\n");

    for (int i = 0; i < oceanSize; ++i) {
        printf("%2d", i);

        for (int j = 0; j < oceanSize; ++j) {
            char displayChar = '.'; 

            for (int k = 0; k < MAX_ANIMALS; ++k) {
                if (ocean[k].type != EMPTY && (int)ocean[k].x == j && (int)ocean[k].y == i) {
                    displayChar = (ocean[k].type == 0) ? 'P' : 'R';  // 'P' for fish, 'R' for shark
                    break;
                }
            }
            printf(" %c", displayChar);
        }
        printf("\n");
    }

     printf("\nPress Enter to continue to the next round...\n");
     getchar(); 
}

int main(int argc, char** argv) {

    MPI_Init(&argc, &argv);

    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    int domain_size = OCEAN_SIZE / sqrt(world_size);
    int start_x = (world_rank % (int)sqrt(world_size)) * domain_size;
    int start_y = (world_rank / (int)sqrt(world_size)) * domain_size;

    int num_received = 0;

    initializeOcean();

    float timeStep = 1.0;
    int count = 0;
    Animal buffer[MAX_ANIMALS];

    for (int step = 0; step < 1000; step++) {
        count = 0;
        for (int i = 0; i < MAX_ANIMALS; i++) {
            if (ocean[i].x >= start_x && ocean[i].x < start_x + domain_size &&
                ocean[i].y >= start_y && ocean[i].y < start_y + domain_size) {
                updateForces(&ocean[i]);
                updatePosition(&ocean[i], timeStep);

            }
            if (ocean[i].x < start_x || ocean[i].x >= start_x + domain_size ||
                ocean[i].y < start_y || ocean[i].y >= start_y + domain_size) {
                buffer[count++] = ocean[i];
                ocean[i].type = EMPTY;
            }
        }

        handleCollisionsAndReproduction();

        MPI_Barrier(MPI_COMM_WORLD);
        exchangeAnimals(world_rank, world_size, count, buffer, &num_received);
        MPI_Barrier(MPI_COMM_WORLD);

        int numAnimalsReceived = (world_rank == 0) ? num_received / sizeof(Animal) : count;
        processReceivedAnimals(buffer, numAnimalsReceived);
        printOcean(ocean, OCEAN_SIZE);
        printf("\n");
        printf("\n");
    }

    MPI_Finalize();
    return 0;
}
