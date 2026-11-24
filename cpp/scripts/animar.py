"""
Gera um GIF animado a partir da saída de rk4_nbody / euler_nbody.

A saída tem uma linha por instante, com "x y z" de cada corpo em sequência,
então o número de corpos é deduzido do número de colunas.

Exemplos (de dentro de cpp/):
    python scripts/animar.py output.txt figura8.gif --2d --nomes "A,B,C"
    python scripts/animar.py output.txt solar.gif --nomes "Sol,Mercúrio,Vênus,Terra,Marte"
"""
import argparse

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import animation

CORES = ['#ffcc33', '#4dd0e1', '#f06292', '#81c784', '#ff8a65', '#b39ddb', '#e0e0e0']
FUNDO = '#0b0f1a'
UA = 1.495978707e11


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('entrada', help='arquivo de saída do simulador')
    p.add_argument('gif', help='GIF a ser gerado')
    p.add_argument('--titulo', default='Simulação de N corpos')
    p.add_argument('--nomes', help='nomes dos corpos separados por vírgula')
    p.add_argument('--2d', dest='plano', action='store_true', help='vista de cima (plano xy)')
    p.add_argument('--quadros', type=int, default=150)
    p.add_argument('--fps', type=int, default=25)
    p.add_argument('--rastro', type=float, default=1.0,
                   help='fração da trajetória mantida como rastro (1 = tudo)')
    p.add_argument('--tamanhos', help='tamanho do marcador de cada corpo, separados por vírgula')
    args = p.parse_args()

    dados = np.loadtxt(args.entrada) / UA
    n = dados.shape[1] // 3
    corpos = [dados[:, 3*i:3*i + 3] for i in range(n)]
    nomes = args.nomes.split(',') if args.nomes else [f'Corpo {i + 1}' for i in range(n)]
    tamanhos = [float(t) for t in args.tamanhos.split(',')] if args.tamanhos else [7] * n
    cores = [CORES[i % len(CORES)] for i in range(n)]
    indices = np.linspace(1, len(dados), args.quadros).astype(int)
    rastro = max(2, int(args.rastro * len(dados)))

    fig = plt.figure(figsize=(6, 6), dpi=100, facecolor=FUNDO)
    if args.plano:
        ax = fig.add_axes([0.04, 0.08, 0.92, 0.82], facecolor=FUNDO)
        ax.set_aspect('equal')
        xs, ys = dados[:, 0::3], dados[:, 1::3]
        cx, cy = (xs.max() + xs.min()) / 2, (ys.max() + ys.min()) / 2
        meia = max(np.ptp(xs), np.ptp(ys)) / 2 * 1.1
        ax.set_xlim(cx - meia, cx + meia)
        ax.set_ylim(cy - meia, cy + meia)
        ax.set_axis_off()
        rastros = [ax.plot([], [], color=c, lw=1.4, alpha=0.8)[0] for c in cores]
        pontos = [ax.plot([], [], 'o', color=c, ms=s, mec='white', mew=0.4)[0]
                  for c, s in zip(cores, tamanhos)]
    else:
        ax = fig.add_subplot(projection='3d', facecolor=FUNDO)
        fig.subplots_adjust(0, 0, 1, 1)
        ax.set_axis_off()
        lim = np.abs(np.delete(dados, np.s_[2::3], axis=1)).max() * 1.02
        ax.set_xlim(-lim, lim)
        ax.set_ylim(-lim, lim)
        ax.set_zlim(-lim / 2, lim / 2)
        ax.set_box_aspect((1, 1, 0.5), zoom=1.45)
        rastros = [ax.plot([], [], [], color=c, lw=1.2, alpha=0.8)[0] for c in cores]
        pontos = [ax.plot([], [], [], 'o', color=c, ms=s, mec='white', mew=0.4)[0]
                  for c, s in zip(cores, tamanhos)]

    fig.text(0.5, 0.95, args.titulo, color='white', ha='center', fontsize=12, weight='bold')
    fig.legend(pontos, nomes, loc='lower center', ncol=min(n, 5), frameon=False,
               labelcolor='white', fontsize=9)

    def atualizar(q):
        fim = indices[q]
        ini = max(0, fim - rastro)
        for corpo, linha, ponto in zip(corpos, rastros, pontos):
            if args.plano:
                linha.set_data(corpo[ini:fim, 0], corpo[ini:fim, 1])
                ponto.set_data(corpo[fim-1:fim, 0], corpo[fim-1:fim, 1])
            else:
                linha.set_data_3d(corpo[ini:fim, 0], corpo[ini:fim, 1], corpo[ini:fim, 2])
                ponto.set_data_3d(corpo[fim-1:fim, 0], corpo[fim-1:fim, 1], corpo[fim-1:fim, 2])
        if not args.plano:
            ax.view_init(elev=25, azim=-70 + 60 * q / args.quadros)
        return rastros + pontos

    anim = animation.FuncAnimation(fig, atualizar, frames=args.quadros, blit=False)
    anim.save(args.gif, writer=animation.PillowWriter(fps=args.fps),
              savefig_kwargs={'facecolor': FUNDO})
    print(f'GIF salvo em {args.gif}')


if __name__ == '__main__':
    main()
