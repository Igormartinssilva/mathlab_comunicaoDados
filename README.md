# mathlab_comunicaoDados
trabalho final, coisar para se ver:
 - se o ldpc soft e hard funcionam e qual a ordem do BER do, sem codificação, soft e hard.
    - sem codificação: não tem correção de bits.
    - codificação hard: trunca o bit, ou ta certo ou errado, nao analisa a probabilidade dele estar errado.
    - codificação soft: analiza a probabilidade do bit estar certo ou errado e corrige ele.
 - ver se a modulação qpsk e 64 qam estão corretas.
    - qpsk: modula somente a fase.
    - qam: modula a fase e amplitude, permite transmitir mais bits por simbolo.
 - ver a matriz de paridade e outros tipos de matriz.
    - Matriz de paridade: Deve-se usar matriz esparsa quase ciclica
    - tamanho do bloco ciclico deve ser z = 54, pois 54 é divisor de n=1296 e k = 43.
    - mariz de deslocamento:  é usada em códigos LDPC esparsos quase cíclicos (QC-LDPC) para definir
       como os blocos circulantes (ou submatrizes) devem ser organizados na matriz de paridade H. Cada 
       elemento da matriz de deslocamento determina o grau de deslocamento (em colunas) para o bloco circulante associado.
 - tentar usar awgn para gerar o ruido do sinal
 - fazer tudo sem ruido
