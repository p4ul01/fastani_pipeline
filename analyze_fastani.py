#!/usr/bin/env python3
"""
analyze_fastani.py
==================
Analisa saida do FastANI e gera:
  - fastani_matrix.csv         : matriz ANI (lower triangular, simetrica)
  - fastani_heatmap.png        : heatmap + dendrograma (clustermap) com
                                 barras de cor por grupo biologico
  - fastani_dendrogram.png     : dendrograma retangular isolado com anotacoes
  - fastani_summary.csv        : estatisticas por genoma (mean, min, max ANI)
  - possible_outliers.csv      : pares com ANI < 95% (possivel contaminacao)

Formato de saida do FastANI (tabular):
  genome1.fna  genome2.fna  ANI_value  count_orfrags  total_orfrags

"""
from __future__ import annotations

import argparse
import os
import re
import sys
import warnings
from pathlib import Path

# Silencia warnings
warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=DeprecationWarning)

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from scipy.spatial.distance import squareform


# ===================================================================
# Normalizacao de nomes (hifen/underscore/espaco) 
# =====================================================================
def _canon_name(name: str) -> str:
    """Normaliza nome: remove extensao, converte hifen/espaco/underscore para underscore unico."""
    name = str(name).strip()
    name = re.sub(r"\.(fna|fasta|fa)$", "", name, flags=re.IGNORECASE)
    name = re.sub(r"[\s\-_]+", "_", name)
    return name.lower()


def _canon_name_loose(name: str) -> str:
    """remove separadores."""
    name = str(name).strip()
    name = re.sub(r"\.(fna|fasta|fa)$", "", name, flags=re.IGNORECASE)
    name = re.sub(r"[\s\-_]+", "", name)
    return name.lower()


def load_groups(groups_csv: str) -> dict:
    """Carrega grupos do CSV. Retorna dict nome -> grupo."""
    if not os.path.exists(groups_csv):
        return {}
    df = pd.read_csv(groups_csv, comment="#")
    df = df.dropna(subset=["Genome"])

    gmap = {}
    for _, row in df.iterrows():
        name = str(row["Genome"]).strip()
        name = re.sub(r"\.(fna|fasta|fa)$", "", name)
        grupo = str(row.get("Grupo", "unknown")).strip().lower()

        if not grupo or grupo == "unknown":
            ptype = str(row.get("Pathotype_original", "")).strip().lower()
            if ptype:
                grupo = ptype

        if "probiotic" in grupo:
            grupo = "probiotico"
        elif ("patogen" in grupo or "pathogen" in grupo
              or any(m in grupo for m in ["ehec","stec","etec","expec","upec",
                      "epec","eiec","eaec","dec","apec","nmec","abu"])):
            grupo = "patogenico"
        elif "comensal" in grupo or "commensal" in grupo:
            grupo = "comensal"
        else:
            grupo = "unknown" 

        # Armazena com multiplas chaves para matching flexivel
        gmap[name] = grupo
        gmap[name.replace("_", " ")] = grupo
        gmap[name.replace(" ", "_")] = grupo
        gmap[name.replace("-", "_")] = grupo
        gmap[name.replace("_", "-")] = grupo
        gmap[_canon_name(name)] = grupo
        gmap[_canon_name_loose(name)] = grupo
    return gmap


def get_group(name: str, gmap: dict) -> str:
    """Resolve grupo de um genoma pelo nome, com multiplas tentativas de matching."""
    name = str(name).strip()
    if not name:
        return "unknown"

    # Tentativa 1: nome exato
    if name in gmap:
        return gmap[name]

    # Tentativa 2: sem extensao
    name2 = re.sub(r"\.(fna|fasta|fa)$", "", name, flags=re.IGNORECASE)
    if name2 in gmap:
        return gmap[name2]

    # Tentativa 3: transformacoes comuns
    for transform in [
        lambda s: s.replace("_", " "),
        lambda s: s.replace(" ", "_"),
        lambda s: s.replace("-", "_"),
        lambda s: s.replace("_", "-"),
        _canon_name,
        _canon_name_loose,
    ]:
        n2 = transform(name)
        if n2 in gmap:
            return gmap[n2]

    # Tentativa 4: matching parcial (substring)
    canon = _canon_name(name)
    for key, val in gmap.items():
        if len(key) >= 3 and (canon == key or canon in key or key in canon):
            return val

    return "unknown"


# =====================================================================
# Parse do arquivo de saida do FastANI
# ====================================================================
def parse_fastani_output(input_file: str) -> pd.DataFrame:
    """
    Le saida tabular do FastANI:
        genome1  genome2  ANI  count_orfrags  total_orfrags

    Retorna DataFrame com colunas: genome_a, genome_b, ani, orthologous_fraction
    """
    rows = []
    with open(input_file) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 5:
                continue
            g1 = os.path.basename(parts[0])
            g2 = os.path.basename(parts[1])
            g1 = re.sub(r"\.(fna|fasta|fa)$", "", g1, flags=re.IGNORECASE)
            g2 = re.sub(r"\.(fna|fasta|fa)$", "", g2, flags=re.IGNORECASE)
            try:
                ani = float(parts[2])
                n_orf = int(parts[3])
                n_tot = int(parts[4])
            except (ValueError, IndexError):
                continue
            frac = n_orf / n_tot if n_tot > 0 else 0.0
            rows.append({
                "genome_a": g1,
                "genome_b": g2,
                "ani": ani,
                "orthologous_fraction": frac,
            })

    if not rows:
        print("ERRO: nenhuma comparacao valida encontrada no arquivo de entrada",
              file=sys.stderr)
        sys.exit(1)

    return pd.DataFrame(rows)


def build_matrix(df: pd.DataFrame) -> pd.DataFrame:
    """
    Constroi matriz simetrica N x N de ANI a partir do DataFrame de pares.
    Diagonal = 100.0 (mesmo genoma).
    """
    genomes = sorted(set(df["genome_a"]).union(set(df["genome_b"])))
    n = len(genomes)
    idx = {g: i for i, g in enumerate(genomes)}

    mat = np.full((n, n), np.nan)
    # Diagonal
    for i in range(n):
        mat[i, i] = 100.0

    # Preenche simetrico (FastANI reporta ambas as direcoes, mas por seguranca
    # copiamos A->B para B->A se faltar)
    for _, row in df.iterrows():
        i, j = idx[row["genome_a"]], idx[row["genome_b"]]
        if not np.isnan(mat[i, j]):
            # Media se ja existe (FastANI as vezes da valores levemente
            # diferentes nas duas direcoes)
            mat[i, j] = (mat[i, j] + row["ani"]) / 2
        else:
            mat[i, j] = row["ani"]
        mat[j, i] = mat[i, j]

    return pd.DataFrame(mat, index=genomes, columns=genomes)


# ==================================================================
# Plotagem: heatmap + dendrograma com barras de grupo
# ==================================================================
PALETTE_GROUP = {
    "probiotico": "#2ca02c",   # verde
    "patogenico": "#d62728",   # vermelho
    "comensal":   "#1f77b4",   # azul
    "unknown":    "#7f7f7f",   # cinza
}


def plot_clustermap(matrix: pd.DataFrame, gmap: dict, output_path: str):
    """
    Heatmap + dendrograma usando seaborn.clustermap.
    Barras de cor lateral indicam o grupo biologico de cada genoma.
    """
    # Distancia = 100 - ANI (ANI alto = distancia baixa)
    dist = 100.0 - matrix.values
    np.fill_diagonal(dist, 0.0)
    # Garante simetria exata
    dist = (dist + dist.T) / 2

    # Condensada para linkage
    condensed = squareform(dist, checks=False)
    Z = linkage(condensed, method="average")

    # Cores por grupo
    groups = [get_group(g, gmap) for g in matrix.index]
    row_colors = pd.Series(
        [PALETTE_GROUP.get(g, PALETTE_GROUP["unknown"]) for g in groups],
        index=matrix.index,
        name="Grupo"
    )

    # Clustermap: ANI em escala 95-100 para destacar diferencas
    fig = sns.clustermap(
        matrix,
        row_linkage=Z, col_linkage=Z,
        row_colors=row_colors, col_colors=row_colors,
        cmap="viridis", vmin=95, vmax=100, # pode ser alterado a vontade
        figsize=(16, 16),
        xticklabels=False, yticklabels=False,
        dendrogram_ratio=(0.15, 0.15),
        cbar_pos=(0.02, 0.82, 0.03, 0.15),
        cbar_kws={"label": "ANI (%)", "shrink": 0.6},
    )

    # Legenda dos grupos
    handles = [mpatches.Patch(color=c, label=g)
               for g, c in PALETTE_GROUP.items() if g in set(groups)]
    fig.ax_heatmap.legend(
        handles=handles, title="Grupo biologico",
        loc="upper left", bbox_to_anchor=(1.02, 1.0),
        frameon=True, fontsize=9,
    )

    fig.fig.suptitle(
        f"FastANI all-vs-all ({len(matrix)} genomes)\n"
        f"Probiotico={groups.count('probiotico')} | "
        f"Patogenico={groups.count('patogenico')} | "
        f"Comensal={groups.count('comensal')} | "
        f"Unknown={groups.count('unknown')}",
        y=1.01, fontsize=12,
    )

    fig.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close("all")
    print(f"  Heatmap: {output_path}")


def plot_dendrogram(matrix: pd.DataFrame, gmap: dict, output_path: str):
    """
    Dendrograma retangular isolado, com barras de cor nas folhas.
    """
    dist = 100.0 - matrix.values
    np.fill_diagonal(dist, 0.0)
    dist = (dist + dist.T) / 2
    condensed = squareform(dist, checks=False)
    Z = linkage(condensed, method="average")

    fig, ax = plt.subplots(figsize=(20, 8))

    # Calcula coordenadas das folhas para posicionar as barras de cor
    ddata = dendrogram(
        Z,
        labels=matrix.index.tolist(),
        leaf_rotation=90,
        leaf_font_size=6,
        ax=ax,
        color_threshold=2.0,  # ANI < 98% = cor diferente
        above_threshold_color="#888888",
    )

    leaf_labels = [lbl.get_text() for lbl in ax.get_xticklabels()]
    groups = [get_group(g, gmap) for g in leaf_labels]
    color_map = [PALETTE_GROUP.get(g, PALETTE_GROUP["unknown"]) for g in groups]

    for tick, color in zip(ax.get_xticklabels(), color_map):
        tick.set_color(color)
        tick.set_fontsize(6)

    ax.set_ylabel("Distance (100 - ANI)")
    ax.set_title(f"FastANI dendrogram ({len(matrix)} genomes)")

    # Legenda
    handles = [mpatches.Patch(color=c, label=g)
               for g, c in PALETTE_GROUP.items() if g in set(groups)]
    ax.legend(handles=handles, title="Grupo", loc="upper right",
              frameon=True, fontsize=9)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    plt.tight_layout()
    fig.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close("all")
    print(f"  Dendrograma: {output_path}")


def generate_summary(matrix: pd.DataFrame, gmap: dict, output_path: str):
    """
    Estatisticas por genoma: ANI medio, minimo (excluindo self), maximo,
    e numero de genomas com ANI >= 99% (proximos).
    """
    rows = []
    for genome in matrix.index:
        vals = matrix.loc[genome].drop(genome).values  # exclui self (100.0)
        if len(vals) == 0:
            continue
        # Filtra NaN (comparacoes que falharam)
        vals = vals[~np.isnan(vals)]
        if len(vals) == 0:
            continue
        n_close = int((vals >= 99.0).sum())
        rows.append({
            "genome": genome,
            "group": get_group(genome, gmap),
            "mean_ani": float(np.mean(vals)),
            "min_ani": float(np.min(vals)),
            "max_ani": float(np.max(vals)),
            "median_ani": float(np.median(vals)),
            "n_genomes_ani_ge_99": n_close,
            "n_comparisons": len(vals),
        })

    summary = pd.DataFrame(rows).sort_values(
        ["group", "genome"], ascending=[True, True])
    summary.to_csv(output_path, index=False)
    print(f"  Summary: {output_path}")
    return summary


def find_outliers(df_pairs: pd.DataFrame, threshold: float = 95.0,
                  output_path: str | None = None) -> pd.DataFrame:
    """
    Lista pares com ANI < threshold (default 95%) - possivel contaminacao
    ou especie diferente.
    """
    outliers = df_pairs[df_pairs["ani"] < threshold].copy()
    outliers = outliers.sort_values("ani")
    if output_path and len(outliers) > 0:
        outliers.to_csv(output_path, index=False)
        print(f"  Outliers (ANI < {threshold}%): {output_path} ({len(outliers)} pares)")
    elif output_path:
        # Cria arquivo vazio com header para registro
        outliers.to_csv(output_path, index=False)
        print(f"  Outliers (ANI < {threshold}%): nenhum par encontrado")
    return outliers


# ====================================================================
# Main
# ===================================================================
def main():
    ap = argparse.ArgumentParser(
        description="Analise de saida do FastANI (matriz + heatmap + dendrograma)")
    ap.add_argument("--input", required=True,
                    help="Arquivo de saida do FastANI (.out)")
    ap.add_argument("--groups", default="",
                    help="CSV com Genome,Grupo[,Pathotype_original]")
    ap.add_argument("--output-dir", default=".",
                    help="Diretorio de saida")
    ap.add_argument("--ani-min", type=float, default=95.0,
                    help="ANI minimo para considerar mesma especie (default 95)")
    args = ap.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"  Lendo: {args.input}")
    df_pairs = parse_fastani_output(args.input)
    n_unique_pairs = len(df_pairs.drop_duplicates(subset=["genome_a", "genome_b"]))
    print(f"  Comparacoes: {len(df_pairs)} linhas ({n_unique_pairs} pares unicos)")

    # Constroi matriz
    matrix = build_matrix(df_pairs)
    n = len(matrix)
    print(f"  Matriz: {n}x{n} genomas")

    # Salva matriz como CSV
    matrix_path = output_dir / "fastani_matrix.csv"
    matrix.to_csv(matrix_path)
    print(f"  Matriz CSV: {matrix_path}")

    # Carrega grupos
    gmap = load_groups(args.groups) if args.groups else {}
    if gmap:
        n_with_group = sum(1 for g in matrix.index if get_group(g, gmap) != "unknown")
        n_without = n - n_with_group
        print(f"  Grupos: {n_with_group}/{n} genomas mapeados")
        if n_without > 0:
            print("  Genomas SEM grupo encontrado:")
            for g in matrix.index:
                if get_group(g, gmap) == "unknown":
                    print(f"      - {g}")

    # Heatmap + dendrograma (clustermap)
    print(f"  Gerando heatmap + dendrograma...")
    plot_clustermap(matrix, gmap, str(output_dir / "fastani_heatmap.png"))

    # Dendrograma retangular isolado
    print(f"  Gerando dendrograma retangular...")
    plot_dendrogram(matrix, gmap, str(output_dir / "fastani_dendrogram.png"))

    # Summary por genoma
    summary = generate_summary(matrix, gmap,
                                str(output_dir / "fastani_summary.csv"))

    # Outliers (possiveis contaminacoes)
    find_outliers(df_pairs, threshold=args.ani_min,
                  output_path=str(output_dir / "possible_outliers.csv"))

    # Estatisticas globais
    ani_vals = df_pairs["ani"].values
    print()
    print(f"  ANI global: min={ani_vals.min():.2f}%, "
          f"median={np.median(ani_vals):.2f}%, max={ani_vals.max():.2f}%")
    print(f"  Pares com ANI < {args.ani_min}%: "
          f"{(df_pairs['ani'] < args.ani_min).sum()}")
    print()
    print(f"  OK: resultados em {output_dir}/")


if __name__ == "__main__":
    main()
