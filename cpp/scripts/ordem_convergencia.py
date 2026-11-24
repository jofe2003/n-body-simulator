"""
Estima a ordem de convergência de um simulador refinando o passo de tempo.

Roda o executável com dt, dt/2, dt/4, ... (mesmo tempo final) e compara as
posições finais de passos consecutivos:

    p ≈ log2( |y(h) - y(h/2)| / |y(h/2) - y(h/4)| )

Uso (de dentro de cpp/):
    python scripts/ordem_convergencia.py ./rk4_nbody   entradas/figura8.txt
    python scripts/ordem_convergencia.py ./euler_nbody entradas/figura8.txt

Argumentos opcionais: [refinamentos=6] [fator=1]. O fator multiplica o dt
inicial (e divide o número de passos): com passos muito pequenos o RK4 atinge
o limite de arredondamento do double antes de mostrar a ordem 4.
"""
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np


def ler_entrada(caminho):
    """Retorna (cabeçalho, linhas dos corpos), ignorando comentários."""
    linhas = []
    for linha in Path(caminho).read_text().splitlines():
        linha = linha.split('#')[0].strip()
        if linha:
            linhas.append(linha)
    n, dt, passos, _ = linhas[0].split()
    return int(n), float(dt), int(passos), linhas[1:]


def posicao_final(executavel, n, dt, passos, corpos, pasta):
    entrada = Path(pasta) / 'entrada.txt'
    saida = Path(pasta) / 'saida.txt'
    # salvar_a_cada = passos -> grava apenas o estado inicial e o final
    entrada.write_text(f'{n} {dt!r} {passos} {passos}\n' + '\n'.join(corpos) + '\n')
    subprocess.run([executavel, str(entrada), str(saida)], check=True,
                   stdout=subprocess.DEVNULL)
    return np.loadtxt(saida)[-1]


def main():
    if len(sys.argv) < 3:
        print(__doc__)
        sys.exit(1)

    executavel = str(Path(sys.argv[1]).resolve())
    refinamentos = int(sys.argv[3]) if len(sys.argv) > 3 else 6
    fator = int(sys.argv[4]) if len(sys.argv) > 4 else 1
    n, dt, passos, corpos = ler_entrada(sys.argv[2])
    dt, passos = dt * fator, passos // fator

    with tempfile.TemporaryDirectory() as pasta:
        finais = [posicao_final(executavel, n, dt / 2**k, passos * 2**k, corpos, pasta)
                  for k in range(refinamentos)]

    print(f'{"dt":>14} {"|y(h)-y(h/2)|":>16} {"ordem p":>9}')
    diferencas = [np.linalg.norm(finais[k] - finais[k + 1]) for k in range(refinamentos - 1)]
    for k, d in enumerate(diferencas):
        ordem = f'{np.log2(diferencas[k - 1] / d):9.3f}' if k > 0 else ' ' * 9
        print(f'{dt / 2**k:14.4e} {d:16.6e} {ordem}')


if __name__ == '__main__':
    main()
