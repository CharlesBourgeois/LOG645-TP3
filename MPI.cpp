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
    float x_center;
    float y_center;
    float total_force_x;
    float total_force_y;
} SubdomainForce;

typedef struct {
    int type; // 0 pour poisson, 1 pour requin
    float x, y; 
    float vx, vy; 
    float ax, ay; 
    int hunger; 
} Animal;

Animal ocean[MAX_ANIMALS];

void initializeLocalOcean(Animal* local_ocean, int* local_count, int start_x, int start_y, int subdomain_size, int world_rank, int world_size) {
    srand((unsigned int)(time(NULL)) ^ (world_rank + 1));
    *local_count = 0;

    int animals_per_process = MAX_ANIMALS / world_size;
    if (world_rank == world_size - 1) {
        animals_per_process += MAX_ANIMALS % world_size;
    }

    for (int i = 0; i < animals_per_process; i++) {
        local_ocean[i].type = rand() % 2;
        local_ocean[i].x = start_x + (rand() % subdomain_size);
        local_ocean[i].y = start_y + (rand() % subdomain_size);
        local_ocean[i].vx = local_ocean[i].vy = 0;
        local_ocean[i].ax = local_ocean[i].ay = 0;
        local_ocean[i].hunger = 0;
        (*local_count)++;
    }
}

void updateForces(Animal* a, Animal* local_ocean, int local_count) {
    float force_x = 0.0, force_y = 0.0;
    force_x += 0.1; 

    for (int i = 0; i < local_count; i++) {
        float distance = sqrt(pow(local_ocean[i].x - a->x, 2) + pow(local_ocean[i].y - a->y, 2));

        if (a->type == 1 && local_ocean[i].type == 0) {
            if (distance < VISIBILITY_RANGE) {
                force_x += ATTRSHARK_CLOSEST / pow(distance, 2);
                force_y += ATTRSHARK_CLOSEST / pow(distance, 2);
            } else {
                force_x += ATTRSHARK / distance;
                force_y += ATTRSHARK / distance;
            }
        }

        if (a->type == 0 && local_ocean[i].type == 1) {
            force_x -= REPPOISSON / distance;
            force_y -= REPPOISSON / distance;
        }
    }

    a->ax = force_x / (a->type == 0 ? MPOISSON : MREQUIN);
    a->ay = force_y / (a->type == 0 ? MPOISSON : MREQUIN);
}

void updateLocalForces(Animal* local_ocean, int local_count) {
    for (int i = 0; i < local_count; i++) {
        updateForces(&local_ocean[i], local_ocean, local_count);
    }
}

void handleLocalCollisionsAndReproduction(Animal* local_ocean, int* local_count) {
    for (int i = 0; i < *local_count; i++) {
        for (int j = i + 1; j < *local_count; j++) {
            float distance = sqrt(pow(local_ocean[i].x - local_ocean[j].x, 2) + 
                                  pow(local_ocean[i].y - local_ocean[j].y, 2));            
            if (distance < 1.0) { 
                if (local_ocean[i].type == local_ocean[j].type) {
                    float temp_vx = local_ocean[i].vx;
                    float temp_vy = local_ocean[i].vy;
                    local_ocean[i].vx = local_ocean[j].vx;
                    local_ocean[i].vy = local_ocean[j].vy;
                    local_ocean[j].vx = temp_vx;
                    local_ocean[j].vy = temp_vy;
                }
                else {
                    // Shark eats fish
                    if (local_ocean[i].type == 1) {
                        local_ocean[j].type = EMPTY; 
                        local_ocean[i].hunger = 0; 
                    }
                    else {
                        local_ocean[i].type = EMPTY; 
                        local_ocean[j].hunger = 0; 
                    }
                }
            }
        }
    }

    for (int i = 0; i < *local_count; i++) {
        if (local_ocean[i].type != EMPTY) {
            if ((double)rand() / RAND_MAX < PREP) {
                for (int k = 0; k < *local_count; k++) {
                    if (local_ocean[k].type == EMPTY) {
                        local_ocean[k].type = local_ocean[i].type;
                        local_ocean[k].x = local_ocean[i].x + ((rand() % 3) - 1);
                        local_ocean[k].y = local_ocean[i].y + ((rand() % 3) - 1);
                        local_ocean[k].vx = local_ocean[k].vy = 0;
                        local_ocean[k].ax = local_ocean[k].ay = 0;
                        local_ocean[k].hunger = 0;
                        break;
                    }
                }
            }

            if (local_ocean[i].type == 1) {
                local_ocean[i].hunger++;
                if (local_ocean[i].hunger > HUNGER_LIMIT) {
                    local_ocean[i].type = EMPTY; 
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

void processReceivedAnimals(Animal* local_ocean, Animal* buffer, int numAnimals, int* local_count, int world_size) {
    int max_animals_per_process = MAX_ANIMALS / world_size;
    for (int i = 0; i < numAnimals; i++) {
        if (*local_count < max_animals_per_process) {
            local_ocean[*local_count] = buffer[i];
            (*local_count)++;
        } else {
            int replace_index = rand() % max_animals_per_process;
            local_ocean[replace_index] = buffer[i];
        }
    }
}

void printOcean(Animal* local_ocean, int local_count, int oceanSize, int world_rank, int world_size) {
    Animal* all_ocean = NULL;

    if (world_rank == 0) {
        all_ocean = (Animal*)malloc(MAX_ANIMALS * sizeof(Animal));
        if (all_ocean == NULL) {
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
    }

    int gather_result = MPI_Gather(local_ocean, local_count * sizeof(Animal), MPI_BYTE,
                                   all_ocean, local_count * sizeof(Animal), MPI_BYTE,
                                   0, MPI_COMM_WORLD);

    if (gather_result != MPI_SUCCESS) {
        if (world_rank == 0) {
            free(all_ocean);
        }
        MPI_Abort(MPI_COMM_WORLD, gather_result);
    }

    if (world_rank == 0) {
        char display[oceanSize][oceanSize];
        memset(display, '.', sizeof(display)); 


        for (int i = 0; i < MAX_ANIMALS; i++) {
            if (all_ocean[i].type != EMPTY) {
                int x = (int)all_ocean[i].x;
                int y = (int)all_ocean[i].y;
                if (x >= 0 && x < oceanSize && y >= 0 && y < oceanSize) {
                    display[y][x] = (all_ocean[i].type == 0) ? 'P' : 'R';
                }
            }
        }

        for (int i = 0; i < oceanSize; i++) {
            for (int j = 0; j < oceanSize; j++) {
                printf(" %c", display[i][j]);
            }
            printf("\n");
        }

        free(all_ocean);
    }

    printf("\nPress Enter to continue to the next round...\n");
    getchar(); 
}

int exchangeAnimals(int world_rank, int world_size, Animal* buffer, int count, Animal* local_ocean, int* local_count) {
    MPI_Status status;
    int numAnimalsReceived = 0;

    int next_rank = (world_rank + 1) % world_size;
    int prev_rank = (world_rank == 0) ? world_size - 1 : world_rank - 1;

    MPI_Sendrecv(&count, 1, MPI_INT, next_rank, 0, &numAnimalsReceived, 1, MPI_INT, prev_rank, 0, MPI_COMM_WORLD, &status);

    Animal* recv_buffer = new Animal[numAnimalsReceived];

    MPI_Sendrecv(buffer, count * sizeof(Animal), MPI_BYTE, next_rank, 0, 
                 recv_buffer, numAnimalsReceived * sizeof(Animal), MPI_BYTE, prev_rank, 0, 
                 MPI_COMM_WORLD, &status);

    for (int i = 0; i < numAnimalsReceived; i++) {
        if (*local_count < MAX_ANIMALS / world_size) {
            local_ocean[*local_count] = recv_buffer[i];
            (*local_count)++;
        }
    }

    delete[] recv_buffer;
    return numAnimalsReceived;
}

int prepareExchange(Animal* local_ocean, Animal* buffer, int* local_count, int start_x, int start_y, int subdomain_size) {
    int count = 0;
    
    for (int i = 0; i < *local_count; i++) {
        int end_x = start_x + subdomain_size;
        int end_y = start_y + subdomain_size;
        
        if (local_ocean[i].x < start_x || local_ocean[i].x >= end_x ||
            local_ocean[i].y < start_y || local_ocean[i].y >= end_y) {
            buffer[count] = local_ocean[i];
            count++;
            
            local_ocean[i].type = EMPTY;
        }
    }

    *local_count -= count;
    
    return count;
}


int main(int argc, char** argv) {

    MPI_Init(&argc, &argv);

    int world_size, world_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    int subdomain_size = OCEAN_SIZE / sqrt(world_size);
    int start_x = (world_rank % (int)sqrt(world_size)) * subdomain_size;
    int start_y = (world_rank / (int)sqrt(world_size)) * subdomain_size;

    Animal local_ocean[MAX_ANIMALS / world_size]; 
    int local_count = 0;

    initializeLocalOcean(local_ocean, &local_count, start_x, start_y, subdomain_size, world_rank, world_size);
            
    float timeStep = 1.0;

    for (int step = 0; step < 1000; step++) {
        updateLocalForces(local_ocean, local_count);
        for (int i = 0; i < local_count; i++) {
            updatePosition(&local_ocean[i], timeStep);
        } 

        handleLocalCollisionsAndReproduction(local_ocean, &local_count);

        Animal buffer[MAX_ANIMALS]; 
        int count = prepareExchange(local_ocean, buffer, &local_count, start_x, start_y, subdomain_size);
        
        int numAnimalsReceived = exchangeAnimals(world_rank, world_size, buffer, count, local_ocean, &local_count);
        processReceivedAnimals(local_ocean, buffer, numAnimalsReceived, &local_count, world_size); 
        
        MPI_Barrier(MPI_COMM_WORLD);

        printOcean(local_ocean, local_count, OCEAN_SIZE, world_rank, world_size);
        
        MPI_Barrier(MPI_COMM_WORLD);
    }

    MPI_Finalize();
    return 0;
}
