"""
Animação 3D das órbitas calculadas por RK44.py.

Lê dados/dados.txt (coluna 0 = tempo; depois x, y, z, vx, vy, vz de cada corpo)
e salva um GIF. Como a estrela parte com velocidade vz = 5 km/s, os planetas a
acompanham e as órbitas desenham uma hélice.

Uso:  python plotter.py [saida.gif]      (padrão: ../assets/tres-corpos-rk44.gif)
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import animation

BASE = Path(__file__).resolve().parent
SAIDA = Path(sys.argv[1]) if len(sys.argv) > 1 else BASE.parent / 'assets' / 'tres-corpos-rk44.gif'

TEMPO_FINAL = 1.6e8   # segundos mostrados na animação (~5 órbitas)
QUADROS = 160
FPS = 25
UA = 1.495978707e11

NOMES = ['Estrela', 'Planeta A', 'Planeta B']
CORES = ['#ffcc33', '#4dd0e1', '#f06292']
FUNDO = '#0b0f1a'


def main():
    data = np.loadtxt(BASE / 'dados' / 'dados.txt')
    data = data[data[:, 0] <= TEMPO_FINAL]
    tempo = data[:, 0]
    corpos = [data[:, 1 + 6*i: 4 + 6*i] / UA for i in range(3)]
    indices = np.linspace(1, len(data) - 1, QUADROS).astype(int)

    fig = plt.figure(figsize=(6, 6), dpi=100, facecolor=FUNDO)
    ax = fig.add_subplot(projection='3d', facecolor=FUNDO)
    fig.subplots_adjust(0, 0, 1, 1)
    ax.set_axis_off()

    todos = np.vstack(corpos)
    ax.set_xlim(-1.1, 1.1)
    ax.set_ylim(-1.1, 1.1)
    ax.set_zlim(todos[:, 2].min() - 0.2, todos[:, 2].max() + 0.2)
    ax.set_box_aspect((1, 1, 1.6))

    rastros = [ax.plot([], [], [], color=c, lw=1.2, alpha=0.75)[0] for c in CORES]
    pontos = [ax.plot([], [], [], 'o', color=c, ms=ms, mec='white', mew=0.4)[0]
              for c, ms in zip(CORES, (11, 6, 6))]
    titulo = fig.text(0.5, 0.95, 'Problema de 3 corpos — Runge-Kutta 4 (Python)',
                      color='white', ha='center', fontsize=12, weight='bold')
    relogio = fig.text(0.5, 0.91, '', color='#9aa4b2', ha='center', fontsize=10)
    fig.legend(pontos, NOMES, loc='lower center', ncol=3, frameon=False,
               labelcolor='white', fontsize=9)

    def atualizar(q):
        n = indices[q]
        for corpo, rastro, ponto in zip(corpos, rastros, pontos):
            rastro.set_data_3d(corpo[:n, 0], corpo[:n, 1], corpo[:n, 2])
            ponto.set_data_3d(corpo[n-1:n, 0], corpo[n-1:n, 1], corpo[n-1:n, 2])
        ax.view_init(elev=18, azim=-60 + 90 * q / QUADROS)
        relogio.set_text(f't = {tempo[n] / 3.15576e7:4.2f} anos')
        return rastros + pontos

    anim = animation.FuncAnimation(fig, atualizar, frames=QUADROS, blit=False)
    SAIDA.parent.mkdir(parents=True, exist_ok=True)
    anim.save(SAIDA, writer=animation.PillowWriter(fps=FPS),
              savefig_kwargs={'facecolor': FUNDO})
    print(f'GIF salvo em {SAIDA}')


if __name__ == '__main__':
    main()
