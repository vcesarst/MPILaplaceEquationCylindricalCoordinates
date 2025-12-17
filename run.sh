#!/bin/bash

# Script para executar experimento MPI com 1-8 cores
# Salva tempos em arquivo .dat

# Nome do executável
EXEC="laplace_mpi.x"
SOURCE_CODE="mpi_opt_with_data.c"
OUTPUT_FILE="mpi_timings.dat"
LOG_DIR="nohup_logs"

# Compilar o código
echo "Compilando código MPI..."
mpicc -O2 -mtune=native -o $EXEC $SOURCE_CODE -lm

if [ $? -ne 0 ]; then
    echo "Erro na compilação!"
    exit 1
fi

# Criar diretório para logs
mkdir -p $LOG_DIR

# Cabeçalho do arquivo de tempos
echo "# Cores Tempo_Total(s) Tempo_Solucao(s) Iteracoes" > $OUTPUT_FILE

# Executar para 1 a 8 cores
for cores in {1..8}; do
    echo "Executando com $cores core(s)..."
    
    # Nome do arquivo de log
    LOG_FILE="$LOG_DIR/nohup_${cores}cores.log"
    
    # Comando MPI - ajuste conforme seu ambiente
    MPI_CMD="mpirun -np $cores ./$EXEC"
    
    # Executar com nohup e capturar output
    echo "=== Execução com $cores core(s) ===" >> $LOG_FILE
    echo "Início: $(date)" >> $LOG_FILE
    
    # Executar e capturar tempo
    start_time=$(date +%s.%N)
    
    nohup $MPI_CMD >> $LOG_FILE 2>&1
    
    end_time=$(date +%s.%N)
    execution_time=$(echo "$end_time - $start_time" | bc)
    
    echo "Término: $(date)" >> $LOG_FILE
    echo "Tempo de execução do shell: ${execution_time}s" >> $LOG_FILE
    
    # Extrair informações dos logs
    total_time=$(grep "Tempo total" $LOG_FILE | awk '{print $4}')
    solve_time=$(grep "Tempo de solucao" $LOG_FILE | awk '{print $4}')
    iterations=$(grep "Iteracoes" $LOG_FILE | tail -1 | awk '{print $3}')
    
    # Se não encontrou nos logs, usar tempo do shell
    if [ -z "$total_time" ]; then
        total_time=$execution_time
        solve_time=$execution_time
    fi
    
    if [ -z "$iterations" ]; then
        iterations="N/A"
    fi
    
    # Salvar no arquivo .dat
    echo "$cores $total_time $solve_time $iterations" >> $OUTPUT_FILE
    
    echo "Finalizado $cores core(s) - Tempo: ${total_time}s"
    echo "----------------------------------------"
    
    # Pequena pausa entre execuções
    sleep 2
done

echo "Experimento concluído!"
echo "Resultados salvos em: $OUTPUT_FILE"
echo "Logs salvos em: $LOG_DIR/"
