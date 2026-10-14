# Prrx1+ 骨折修复多组学分析项目

**主题**: Artificial Intelligence and Multi-omics Decoding of Inflammatory Mechanisms
and Therapeutic Interventions in Musculoskeletal Soft-Tissue Injury

整合 GSE195940 / GSE164573 / GSE138689 + 自有
Prrx1CreERT;R26tdTomato 骨折 d7 切片 + 稳态骨髓 Prrx1+ Bulk RNA-seq (n=3)。

## 目录结构

```
prrx1_fracture_multiomics/
├── README.md
├── data/
│   ├── README_download.md     <- 全部下载直链 (GEO/NicheNet/SCENIC/Enrichr)
│   ├── GSM_sample_map.csv     <- 样本编号与实验条件映射
│   ├── bulk/prrx1_bulk_counts.csv   <- 自有 Bulk count (放入)
│   ├── functional/qPCR_rescue.csv   <- rescue 实验数据 (放入)
│   ├── GSE195940/ GSE164573/ GSE138689/
│   └── DSigDB/                <- DSigDB_mm.rds + L1000_down/up.rds
├── nichenetr/                 <- NicheNet v2 网络 (4 个 rds)
├── scenic/                    <- SCENIC feather + motifs + mm_tfs.txt
├── scripts/                   <- 00-11 + run_all.R
└── results/{figures,tables}/
```

## 运行流程

```bash
# 0. 依赖安装 (一次)
Rscript scripts/00_setup.R
# 1. 数据下载 (国内建议浏览器手动下, 见 data/README_download.md)
Rscript scripts/00b_download_geo.R
Rscript scripts/00c_download_resources.R
# 2. 自动解包整理 RAW.tar -> raw/<cond>/ 目录 (无需联网)
Rscript scripts/00d_organize_geo.R
#    核对 results/tables/geo_sample_organization.csv, 放入自有 Bulk 数据
# 3. 一键分析
Rscript scripts/run_all.R
# 4. SCENIC (Python)
pip install pyscenic arboreto loompy && python scripts/06_scenic.py
```

## 关键注意

1. **celltype 注释**: `02` 运行后须根据 `FigS1_marker_dotplot.pdf` 人工完成
   `seurat_clusters -> celltype` 映射, 再继续后续脚本。
2. **Prx1Cre vs Prrx1CreERT**: 公共数据为胚胎谱系标记, 自有数据为成年诱导
   标记 —— 保留 `cre` 元数据列加以区分, 是本文核心创新点。
3. **Bulk n=3** 统计效力有限, 仅用于身份鉴定与反卷积参考。
4. **网络问题**: R 下载失败但浏览器能开 -> 浏览器下载后放对应目录,
   脚本检测到文件存在会自动跳过。代理: Sys.setenv(https_proxy="http://127.0.0.1:7890")

## 更新日志
- v1.0: 初版 (Aim1-4 全套)
- v1.1: 修复 R 4.4 安装源 (nichenetr/CellChat->GitHub, BisqueRNA->CRAN,
  tradeSeq/slingshot/fgsea->Bioconductor); 移除弃用的 loomR
- v1.2: 修复 NicheNet 文件名 (21122021 日期后缀); DSigDB 改 Enrichr API;
  mm_tfs 增加 jsdelivr 镜像; 下载脚本全部增加"已存在跳过"
- v1.3 (最终版): 10_drug_repurposing 增加小鼠->人源 toupper 转换与
  离线 L1000 逆转分析 (LINCS_L1000_Chem_Pert_down/up)

- v1.4: 新增 00d_organize_geo.R (自动解包 RAW.tar + 按关键词分类到 raw/<cond>/
  标准 10X 目录); 01 的 list.files 改递归兼容

- v1.6: 02 的 IntegrateLayers 替换为 harmony::RunHarmony
  (修复 JoinLayers+SCT 下 contrasts 单水平因子错误)
- v1.7: 完成 30 个 cluster 四重验证注释 (见 02 脚本内注释块);
  重要发现: 环境 RNA 极强, Ptprc 亦不可信, 方法学已记录
- v1.8: 新增 01b_soupx.R (SoupX 环境RNA校正), 校正后重跑 02;
  02 自动优先使用 results/soupx_objects.rds
- v1.9: 环境RNA校正改用 DecontX (celda) 替代 SoupX —— SoupX 依赖空液滴
  原始矩阵, 本数据只有过滤后矩阵时其内部分布估计退化(NA); DecontX 仅需
  矩阵+聚类标签, 更适配本场景。01b_soupx.R 已移除, 新增 01b_decontx.R
- v2.0: DecontX 后完成最终 29-cluster 注释 (见 02 内注释块, 依据
  cluster_marker_pct.csv 定量表); 03/04 细胞类型列表加回 FCP
- v2.2: 样本输入定稿为 9 个; 新增 02b_annotate.R (外部注释定稿);
  02 不再自动应用历史注释。衰老臂改用 BM_SS vs BM_d7 + 自有数据。
- v2.3: 9样本终版 (54,606细胞/23群) 注释定稿并写入 02b_annotate.R。
  亮点: 双Chondro群(纤维软骨12/成熟软骨17) + 教科书级Osteo(21) 构成分化序列
- v2.4: 全链路移除问题数据集 (下载/解包/质控/DecontX/重建等全部相关脚本与描述);
  删除 01c_callus_from_authors.R 与 01d_final9.R, 02 输入优先级简化为 decontx > qc;
  全部图形统一经 fig_utils.R 导出 PDF + 600 dpi TIFF
