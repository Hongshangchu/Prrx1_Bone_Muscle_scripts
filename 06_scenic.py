# ============================================================================
# 06_scenic.py —— Aim2c: pySCENIC regulon 分析 (断点续跑版)
# 依赖: conda 环境 scenic (python3.10 + dask2023.8.1 + numpy1.26 + setuptools<70)
# 重要: 可执行代码必须在 main() 内 (Windows spawn); 已完成的阶段自动跳过
# ============================================================================
import os
import pandas as pd
import numpy as np
import loompy
import subprocess
from arboreto.algo import grnboost2


def main():
    ## 1. 表达矩阵
    ex_matrix = pd.read_csv("results/msc_counts.csv", index_col=0)
    tf_names = [l.strip() for l in open("scenic/mm_mgi_tfs.txt")]
    print("表达矩阵:", ex_matrix.shape, "| TF数:", len(tf_names))

    ## 2. GRNBoost2 (已有结果则跳过, 支持断点续跑)
    if os.path.exists("results/scenic_adjacencies.csv"):
        print("检测到 scenic_adjacencies.csv, 跳过 GRNBoost2")
    else:
        adjacencies = grnboost2(ex_matrix.T, tf_names=tf_names, verbose=True)
        adjacencies.to_csv("results/scenic_adjacencies.csv", index=False)
        print("GRNBoost2 完成:", adjacencies.shape)

    ## 3. loom
    if not os.path.exists("results/msc.loom"):
        row_attrs = {"Gene": ex_matrix.index.values}
        col_attrs = {"CellID": ex_matrix.columns.values}
        loompy.create("results/msc.loom",
                      ex_matrix.values.astype(float), row_attrs, col_attrs)

    feather = ("scenic/mm10__refseq-r80__10kb_up_and_down_tss.mc9nr."
               "genes_vs_motifs.rankings.feather")

    ## 4. cisTarget 剪枝 (跳过若已存在)
    if os.path.exists("results/scenic_regulons.csv"):
        print("检测到 scenic_regulons.csv, 跳过剪枝")
    else:
        subprocess.run(["pyscenic", "ctx",   # pyscenic>=0.11 剪枝子命令为 ctx
            "results/scenic_adjacencies.csv", feather,
            "--annotations_fname", "scenic/motifs-v9-nr.mgi-m0.001-o0.0.tbl",
            "--expression_mtx_fname", "results/msc.loom",   # loom 行列语义明确, 避免 CSV 方向歧义
            "-o", "results/scenic_regulons.csv", "--num_workers", "8"],
            check=True)

    ## 5. Regulon AUC (跳过若已存在)
    if os.path.exists("results/scenic_aucell.loom"):
        print("检测到 scenic_aucell.loom, 跳过 AUC")
    else:
        subprocess.run(["pyscenic", "aucell",
            "results/msc.loom", "results/scenic_regulons.csv",
            "-o", "results/scenic_aucell.loom", "--num_workers", "8"],
            check=True)
    print("SCENIC 全部完成")


if __name__ == "__main__":
    main()
