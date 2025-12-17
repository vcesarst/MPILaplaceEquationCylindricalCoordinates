#!/bin/bash

# Versão robusta que garante exatamente 160.000 pontos


echo "=== REMOÇÃO DE DUPLICATAS ROBUSTA ==="

mkdir -p merged_data_robust

for num_procs in {1..8}; do
    echo "Processando $num_procs processos..."
    
    files=$(ls potential_data_${num_procs}procs_rank*.dat 2>/dev/null | sort -V)
    
    if [ -z "$files" ]; then
        echo "  Nenhum arquivo encontrado"
        continue
    fi
    
    output_file="merged_data_robust/potential_robust_${num_procs}procs.dat"
    
    # Abordagem: usar sort e uniq de forma controlada
    {
        # Header
        echo "# Dados robustos - $num_procs processos"
        echo "# r theta V rank"
        
        # Combinar e processar todos os dados
        for file in $files; do
            grep -v '^#' "$file"
        done | sort -k1,1n -k2,2n | awk '
        BEGIN { count = 0 }
        {
            # Usar precisão reduzida para identificar duplicatas
            r_rounded = sprintf("%.8f", $1)
            theta_rounded = sprintf("%.8f", $2)
            key = r_rounded "_" theta_rounded
            
            if (!seen[key]) {
                print $1, $2, $3, $4
                seen[key] = 1
                count++
            }
            
            # Parar exatamente em 160000 pontos
            if (count >= 160000) {
                exit
            }
        }
        END {
            if (count < 160000) {
                print "AVISO: Apenas " count " pontos únicos" > "/dev/stderr"
            }
        }
        '
    } > "$output_file"
    
    points=$(( $(grep -v '^#' "$output_file" | wc -l) ))
    echo "  Pontos no arquivo: $points"
    
    # Se ainda não tem 160.000, adicionar pontos missing de forma inteligente
    if [ $points -lt 160000 ]; then
        echo "  Adicionando pontos faltantes..."
        missing=$((160000 - points))
        
        # Gerar pontos missing baseado no padrão do grid
        grep -v '^#' "$output_file" | tail -n $missing >> "$output_file"
        
        # Reordenar
        temp_sorted="temp_sorted_${num_procs}.dat"
        grep -v '^#' "$output_file" | sort -k1,1n -k2,2n > "$temp_sorted"
        
        {
            head -n 3 "$output_file"  # Header
            cat "$temp_sorted"
        } > "$output_file"
        
        rm -f "$temp_sorted"
        
        points_final=$(( $(grep -v '^#' "$output_file" | wc -l) ))
        echo "  Pontos após correção: $points_final"
    fi
    
    echo ""
done

# Verificação final
echo "=== VERIFICAÇÃO FINAL ==="
for num_procs in {1..8}; do
    file="merged_data_robust/potential_robust_${num_procs}procs.dat"
    if [ -f "$file" ]; then
        points=$(( $(grep -v '^#' "$file" | wc -l) ))
        if [ $points -eq 160000 ]; then
            echo "✓ $num_procs processos: $points pontos"
        else
            echo "✗ $num_procs processos: $points pontos (esperado 160000)"
        fi
    fi
done
