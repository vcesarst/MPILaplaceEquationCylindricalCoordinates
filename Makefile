# Makefile para solução paralela da equação de Laplace com MPI, compilado com GNU

# Compilador e flags
CC = mpicc
CFLAGS = -O2 -mtune=native
TARGET = laplace_mpi.x
SOURCE = mpi_opt_with_data.c

# Scripts
RUN_SCRIPT = run.sh
CONCAT_SCRIPT = concat_robust.sh

# Arquivos gerados
DAT_FILES = potential_data_*.dat
RESULTS = resultado_final.dat
BINARY = $(TARGET)

# Targets principais
.PHONY: all compile run clean full

# Target padrão: apenas compilar
all: compile

# Compilar o programa
compile: $(SOURCE)
	$(CC) $(CFLAGS) -o $(TARGET) $(SOURCE) -lm

# Executar o programa (chama run.sh e depois concat_robust.sh)
run: compile
	@echo "Executando o programa MPI..."
	@if [ -f "$(RUN_SCRIPT)" ]; then \
		chmod +x $(RUN_SCRIPT) && ./$(RUN_SCRIPT); \
	else \
		echo "Erro: $(RUN_SCRIPT) não encontrado!"; \
		exit 1; \
	fi
	@echo "Concatenando arquivos .dat..."
	@if [ -f "$(CONCAT_SCRIPT)" ]; then \
		chmod +x $(CONCAT_SCRIPT) && ./$(CONCAT_SCRIPT); \
	else \
		echo "Erro: $(CONCAT_SCRIPT) não encontrado!"; \
		exit 1; \
	fi
	@echo "Limpando arquivos .dat intermediários..."
	@if ls $(DAT_FILES) 1> /dev/null 2>&1; then \
		rm -f $(DAT_FILES); \
		echo "Arquivos .dat removidos."; \
	else \
		echo "Nenhum arquivo .dat encontrado para limpar."; \
	fi

# Executar tudo em sequência (compilar + executar + concatenar + limpar)
full: run

# Limpar tudo exceto o código fonte e scripts
clean:
	@echo "Limpando arquivos gerados..."
	@if [ -f "$(TARGET)" ]; then \
		rm -f $(TARGET); \
		echo "Arquivo executável $(TARGET) removido."; \
	fi
	@if ls $(DAT_FILES) 1> /dev/null 2>&1; then \
		rm -f $(DAT_FILES); \
		echo "Arquivos .dat removidos."; \
	fi
	@if [ -f "$(RESULTS)" ]; then \
		rm -f $(RESULTS); \
		echo "Arquivo $(RESULTS) removido."; \
	fi
	@echo "Limpeza concluída. Arquivos .c e .sh mantidos."

# Help
help:
	@echo "Opções disponíveis:"
	@echo "  make compile      - Compila o programa MPI (gera $(TARGET))"
	@echo "  make run          - Compila, executa run.sh, concatena e limpa .dat"
	@echo "  make full         - Mesmo que 'make run'"
	@echo "  make clean        - Remove todos os arquivos gerados"
	@echo "  make help         - Mostra esta ajuda"
	@echo ""
	@echo "Arquivos mantidos após clean:"
	@echo "  - $(SOURCE)"
	@echo "  - $(RUN_SCRIPT)"
	@echo "  - $(CONCAT_SCRIPT)"
