"""
Compara a convergência de Euler e RK4 em um gráfico log-log (erro x passo).

Uso (de dentro de cpp/, depois de `make`):
    python scripts/comparar_metodos.py [entrada] [saida.png]
    (padrão: entradas/figura8.txt e ../assets/ordem-convergencia.png)
"""
import sys
import tempfile
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from ordem_convergencia import ler_entrada, posicao_final

CPP = Path(__file__).resolve().parent.parent
FUNDO = '#0b0f1a'
METODOS = [('Euler semi-implícito', CPP / 'euler_nbody', '#f06292', 1),
           ('Runge-Kutta 4', CPP / 'rk4_nbody', '#4dd0e1', 4)]
REFINAMENTOS = 6
FATOR = 8   # começa com dt 8x maior que o do arquivo (antes do limite do double)


def main():
    entrada = Path(sys.argv[1]) if len(sys.argv) > 1 else CPP / 'entradas' / 'figura8.txt'
    saida = Path(sys.argv[2]) if len(sys.argv) > 2 else CPP.parent / 'assets' / 'ordem-convergencia.png'
    n, dt, passos, corpos = ler_entrada(entrada)
    dt, passos = dt * FATOR, passos // FATOR
    dts = np.array([dt / 2**k for k in range(REFINAMENTOS - 1)])

    fig, ax = plt.subplots(figsize=(7, 4.5), dpi=110, facecolor=FUNDO)
    ax.set_facecolor(FUNDO)
    with tempfile.TemporaryDirectory() as pasta:
        for nome, exe, cor, ordem in METODOS:
            finais = [posicao_final(str(exe), n, dt / 2**k, passos * 2**k, corpos, pasta)
                      for k in range(REFINAMENTOS)]
            erros = np.array([np.linalg.norm(finais[k] - finais[k + 1])
                              for k in range(REFINAMENTOS - 1)])
            p = np.polyfit(np.log2(dts), np.log2(erros), 1)[0]
            ax.loglog(dts, erros, 'o-', color=cor, lw=2, label=f'{nome}: ordem medida ≈ {p:.2f}')
            ref = erros[0] * (dts / dts[0]) ** ordem
            ax.loglog(dts, ref, '--', color=cor, alpha=0.4, lw=1)

    ax.set_xlabel('passo de tempo Δt (s)', color='white')
    ax.set_ylabel('|y(Δt) − y(Δt/2)|  (m)', color='white')
    ax.set_title('Convergência na órbita em figura 8 (tracejado: inclinação teórica)',
                 color='white', fontsize=11)
    ax.tick_params(colors='#9aa4b2', which='both')
    for lado in ax.spines.values():
        lado.set_color('#3a4252')
    ax.grid(True, which='major', color='#232a38', lw=0.6)
    ax.legend(frameon=False, labelcolor='white')
    fig.tight_layout()
    saida.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(saida, facecolor=FUNDO)
    print(f'Gráfico salvo em {saida}')


if __name__ == '__main__':
    main()
