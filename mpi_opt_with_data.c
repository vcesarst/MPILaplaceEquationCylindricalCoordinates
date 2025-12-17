#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <mpi.h>

#define PI 3.14159265358979323846
#define TOL 1e-6
#define ITMAX 1000000
#define NR 400
#define NTHETA 400

//Estrutura de dados padrão que cada processo receberá, assumindo vetorização dos dados
typedef struct {
    int start_i;      
    int end_i;        
    int local_nr;     
    int ntheta;       
    double *r_vals;   
    double *theta;    
    double *V;        
    double *V_new;   
    double dr;
    double dtheta;
    double R;
} LocalGrid;

//Função para acesso aos índices dos vetores de dados
int idx(int i, int j, int ntheta) {
    return i * ntheta + j;
}

//Função para liberação de memória
void free_local_grid(LocalGrid *grid) {
    if (grid->r_vals) free(grid->r_vals);
    if (grid->theta) free(grid->theta);
    if (grid->V) free(grid->V);
    if (grid->V_new) free(grid->V_new);
}

//Discretização realizadana componente radial da malha
void discretize_r(double *r, int n, double R) {
    double dr = R / (n - 1);
    for(int i = 0; i < n; i++) {
        r[i] = i * dr;
    }
    r[0] = 0.0;
}

//discretização na componente angular
void discretize_theta(double *theta, int m) {
    double dtheta = 2.0 * PI / m;
    for(int j = 0; j < m; j++) {
        theta[j] = j * dtheta;
    }
}

// Inicialição da malha local, com o processo master distribuindo para todos os processos dentro do comunicador os elementos inicializados na estrutura de dados. Assume divisão com resto e seu subsequente tratamento
void initialize_local_grid(LocalGrid *grid, int rank, int size, double R, int nr, int ntheta, int use_good_guess) {
    int base_points = nr / size;
    int remainder = nr % size;
    
    if (rank < remainder) {
        grid->local_nr = base_points + 1;
        grid->start_i = rank * grid->local_nr;
    } else {
        grid->local_nr = base_points;
        grid->start_i = remainder * (base_points + 1) + (rank - remainder) * base_points;
    }
    grid->end_i = grid->start_i + grid->local_nr - 1;
    
    grid->ntheta = ntheta;
    grid->R = R;
    
    grid->r_vals = (double*)malloc(grid->local_nr * sizeof(double));
    grid->theta = (double*)malloc(ntheta * sizeof(double));
    grid->V = (double*)calloc(grid->local_nr * ntheta, sizeof(double));
    grid->V_new = (double*)calloc(grid->local_nr * ntheta, sizeof(double));
    
    double *global_r = (double*)malloc(nr * sizeof(double));
    if (rank == 0) {
        discretize_r(global_r, nr, R);
    }
    MPI_Bcast(global_r, nr, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    for(int i = 0; i < grid->local_nr; i++) {
        int global_i = grid->start_i + i;
        grid->r_vals[i] = global_r[global_i];
    }
    
    discretize_theta(grid->theta, ntheta);
    grid->dr = global_r[1] - global_r[0];
    grid->dtheta = grid->theta[1] - grid->theta[0];
    
    for(int i = 0; i < grid->local_nr; i++) {
        for(int j = 0; j < ntheta; j++) {
            int index = idx(i, j, ntheta);
            if (use_good_guess) {
                grid->V[index] = 100.0 * (grid->r_vals[i] / R) * cos(grid->theta[j]);
            }
            grid->V_new[index] = grid->V[index];
        }
    }
    
    free(global_r);
}
//Aplicação da condição de contorno na casca cilíndrica: 100*cos(theta) = V(R,theta)
void apply_boundary_conditions(LocalGrid *grid, int rank, int size) {
    if (rank == size - 1) {
        int last_i = grid->local_nr - 1;
        for(int j = 0; j < grid->ntheta; j++) {
            int index = idx(last_i, j, grid->ntheta);
            grid->V[index] = 100.0 * cos(grid->theta[j]);
            grid->V_new[index] = grid->V[index];
        }
    }
}
//Comunicação de dados entre os processos para que um determinado processos receba os dados de seu antecessos e sucessor, e envie seus dados para os mesmos. É realizada uma regra de envio: primeiro os ranks pares enviam os dados e os ímpares recebem. Depois os pares recebem enquanto os ímpares enviam. É um tipo de comunicação que evita deadlock.
void exchange_boundary_data(LocalGrid *grid, int rank, int size, 
                           double *left_boundary, double *right_boundary) {
    MPI_Status status;
    
    memset(left_boundary, 0, grid->ntheta * sizeof(double));
    memset(right_boundary, 0, grid->ntheta * sizeof(double));
    
    if (size > 1) {
        if (rank % 2 == 0) {
            if (rank < size - 1) {
                MPI_Send(&grid->V[idx(grid->local_nr - 1, 0, grid->ntheta)], 
                        grid->ntheta, MPI_DOUBLE, rank + 1, 100, MPI_COMM_WORLD);
            }
            if (rank > 0) {
                MPI_Send(&grid->V[idx(0, 0, grid->ntheta)], 
                        grid->ntheta, MPI_DOUBLE, rank - 1, 200, MPI_COMM_WORLD);
            }
            
            if (rank > 0) {
                MPI_Recv(left_boundary, grid->ntheta, MPI_DOUBLE, rank - 1, 
                        300, MPI_COMM_WORLD, &status);
            }
            if (rank < size - 1) {
                MPI_Recv(right_boundary, grid->ntheta, MPI_DOUBLE, rank + 1, 
                        400, MPI_COMM_WORLD, &status);
            }
        } else {
            if (rank > 0) {
                MPI_Recv(left_boundary, grid->ntheta, MPI_DOUBLE, rank - 1, 
                        100, MPI_COMM_WORLD, &status);
            }
            if (rank < size - 1) {
                MPI_Recv(right_boundary, grid->ntheta, MPI_DOUBLE, rank + 1, 
                        200, MPI_COMM_WORLD, &status);
            }
            
            if (rank < size - 1) {
                MPI_Send(&grid->V[idx(grid->local_nr - 1, 0, grid->ntheta)], 
                        grid->ntheta, MPI_DOUBLE, rank + 1, 300, MPI_COMM_WORLD);
            }
            if (rank > 0) {
                MPI_Send(&grid->V[idx(0, 0, grid->ntheta)], 
                        grid->ntheta, MPI_DOUBLE, rank - 1, 400, MPI_COMM_WORLD);
            }
        }
    }
}


// Função para realizar o passo de Jacobi para resolução do sistema de equações para os potenciais.
double jacobi_iteration_mpi(LocalGrid *grid, int rank, int size) {
    double local_error_sum = 0.0;
    int local_interior_points = 0;
    
    double dr2 = grid->dr * grid->dr;
    double dtheta2 = grid->dtheta * grid->dtheta;
    
    double *left_boundary = (double*)malloc(grid->ntheta * sizeof(double));
    double *right_boundary = (double*)malloc(grid->ntheta * sizeof(double));
    
    if (!left_boundary || !right_boundary) {
        return 0.0;
    }
    
    exchange_boundary_data(grid, rank, size, left_boundary, right_boundary);
    
    for(int i = 0; i < grid->local_nr; i++) {
        double r_i = grid->r_vals[i];
        
        for(int j = 0; j < grid->ntheta; j++) {
            int idx_center = idx(i, j, grid->ntheta);
            double V_old = grid->V[idx_center];
            double V_new = V_old;
            
            if (rank == 0 && i == 0 && r_i < 1e-12) {
                double sum = 0.0;
                for(int k = 0; k < grid->ntheta; k++) {
                    sum += grid->V[idx(1, k, grid->ntheta)];
                }
                V_new = sum / grid->ntheta;
            }
            else if (rank == size - 1 && i == grid->local_nr - 1) {
                V_new = V_old;
            }
            else {
                int j_plus = (j + 1) % grid->ntheta;
                int j_minus = (j == 0) ? grid->ntheta - 1 : j - 1;
                
                double V_r_plus, V_r_minus;
                
                if (i == grid->local_nr - 1) {
                    V_r_plus = (rank < size - 1) ? right_boundary[j] : V_old;
                } else {
                    V_r_plus = grid->V[idx(i + 1, j, grid->ntheta)];
                }
                
                if (i == 0) {
                    V_r_minus = (rank > 0) ? left_boundary[j] : V_old;
                } else {
                    V_r_minus = grid->V[idx(i - 1, j, grid->ntheta)];
                }
                
                double V_theta_plus = grid->V[idx(i, j_plus, grid->ntheta)];
                double V_theta_minus = grid->V[idx(i, j_minus, grid->ntheta)];
                
                double r_safe = (r_i < 1e-12) ? grid->dr : r_i;
                double r2 = r_safe * r_safe;
                
                double numerator = (V_r_plus + V_r_minus) / dr2 + 
                                 (V_r_plus - V_r_minus) / (2.0 * r_safe * grid->dr) +
                                 (V_theta_plus + V_theta_minus) / (r2 * dtheta2);
                
                double denominator = 2.0 / dr2 + 2.0 / (r2 * dtheta2);
                
                if (denominator > 1e-15) {
                    V_new = numerator / denominator;
                }
            }
            
            grid->V_new[idx_center] = V_new;
            
            if (!(rank == size - 1 && i == grid->local_nr - 1)) {
                double diff = fabs(V_new - V_old);
                local_error_sum += diff;
                local_interior_points++;
            }
        }
    }
    
    double *temp = grid->V;
    grid->V = grid->V_new;
    grid->V_new = temp;
    
    apply_boundary_conditions(grid, rank, size);
    
    free(left_boundary);
    free(right_boundary);
    
    return (local_interior_points > 0) ? local_error_sum / local_interior_points : 0.0;
}

//Dados salvos por rank para então serem concatenados via bash. Essa técnica evita que os processos enviem de maneira errada os dados finais para o nó master.
void save_complete_data(LocalGrid *grid, int rank, int size, int np) {
    // Cada processo salva seus dados locais
    char filename[100];
    sprintf(filename, "potential_data_%dprocs_rank%02d.dat", np, rank);
    
    FILE *fp = fopen(filename, "w");
    if (fp) {
        fprintf(fp, "# r theta V rank\n");
        for(int i = 0; i < grid->local_nr; i++) {
            for(int j = 0; j < grid->ntheta; j++) {
                int index = idx(i, j, grid->ntheta);
                fprintf(fp, "%.6f %.6f %.6f %d\n", 
                        grid->r_vals[i], grid->theta[j], grid->V[index], rank);
            }
        }
        fclose(fp);
    }
}

int main(int argc, char *argv[]) {
    int rank, size;
    double total_start, total_end;
    
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    total_start = MPI_Wtime();
    
    if (rank == 0) {
        printf("=== INICIANDO EXECUÇÃO MPI ===\n");
        printf("Processos: %d, Grid: %dx%d\n", size, NR, NTHETA);
        printf("Tolerancia: 1e-8, Max iter: %d\n", ITMAX);
        printf("Horario de inicio: %s", ctime(&(time_t){total_start}));
    }
    
    LocalGrid grid;
    int use_good_guess = 0;
    
    initialize_local_grid(&grid, rank, size, 5.0, NR, NTHETA, use_good_guess);
    apply_boundary_conditions(&grid, rank, size);
    
    // Loop principal com output reduzido
    int it = 0;
    double avg_error = 1.0;
    double solve_start = MPI_Wtime();
    
    while(avg_error > TOL && it < ITMAX) {
        double local_error = jacobi_iteration_mpi(&grid, rank, size);
        
        double global_error_sum;
        MPI_Allreduce(&local_error, &global_error_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        avg_error = global_error_sum / size;
        it++;
        
        // Output a cada 10000 iterações apenas
        if (it % 10000 == 0 && rank == 0) {
            printf("Iter %d, Erro: %.2e\n", it, avg_error);
        }
    }
    
    double solve_end = MPI_Wtime();
    total_end = MPI_Wtime();
    
    // Salvar dados completos
    save_complete_data(&grid, rank, size, size);
    
    // Output final apenas pelo rank 0
    if (rank == 0) {
        printf("\n=== RESULTADOS FINAIS ===\n");
        printf("Processos: %d\n", size);
        printf("Iteracoes: %d\n", it);
        printf("Erro final: %.2e\n", avg_error);
        printf("Tempo total: %.2f segundos\n", total_end - total_start);
        printf("Tempo de solucao: %.2f segundos\n", solve_end - solve_start);
        printf("Status: %s\n", (avg_error <= TOL) ? "CONVERGIU" : "MAX ITERACOES");
        printf("Horario de termino: %s", ctime(&(time_t){total_end}));
    }
    
    free_local_grid(&grid);
    MPI_Finalize();
    return 0;
}
