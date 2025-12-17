PARA EXECUTAR O PROGRAMA:

Conforme explicitado no trabalho escrito, os dados são escritos por rank e são concatenados no fim de modo a evitar deadlock e sobrecarga no processo master. Com este objetivo, somado à eliminação das duplicata, o script de concatenação robusta foi criado. Esse script roda após o script de execução variando o número decores de 1 a 8, assumindo que o programa será rodado para avaliação em um computador local, e não de um ambiente de cluster. 

Via Makefile:

make compile : compila o programa simplesmente

make run : compila e roda o programa, concatenando os arquivos .dat de cada rank

make clean :  limpa os dados

make full : realiza todos os processos acima sequencialmente.
