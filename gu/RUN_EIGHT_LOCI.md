在 WSL Bash 中运行下列命令。每一步成功完成后再运行下一步；分析命令使用前台模式，避免 `final` 在分析完成前启动。

```bash
cd /mnt/d/scripts/gu

# 1. 完整生成八个位点的区域匹配、候选树和 marker 注释。
# 相同输入可复用 raw comparison；派生证据会重新生成。
./gu.sh phyml --loci /mnt/d/files/gu.37.bed --grch 37 \
  --loci-flank 100kb --plot-phy TRUE --jobs 4 --foreground TRUE

# 2. 补齐八条染色体的 IBDmix；已完成且参数相同的结果会复用。
./gu.sh ibdmix --chr 1,3,6,8,9,12,19,X --grch 37 --target 1kg \
  --target-dir /mnt/d/data.BIG/refGen/1kg/37/pfile/chr \
  --jobs 4 --foreground TRUE

# 3. 为 TRACE 准备完整染色体 ARG 输入并构树。
./arg.sh prep_gen --dir-gen /mnt/d/data.BIG/refGen/1kg/37 \
  --method needle --format trace --chr 1,3,6,8,9,12,19,X \
  --threads 8 --jobs 1

./arg.sh build --dir-gen /mnt/d/data.BIG/refGen/1kg/37 \
  --method needle --format trace --chr 1,3,6,8,9,12,19,X \
  --threads 8 --jobs 1

# 4. 在这些染色体上运行 TRACE。
# 替换已有 TRACE 输出，使这些染色体统一采用本次 ARG backend。
./gu.sh trace --chr 1,3,6,8,9,12,19,X --grch 37 --target 1kg \
  --target-dir /mnt/d/data.BIG/refGen/1kg/37/pfile/chr \
  --arg-dir /mnt/d/data.BIG/refGen/1kg/37/arg/trace/needle \
  --replace-trace TRUE --jobs 4 --foreground TRUE

# 5. 整合所有已完成结果并生成逐谱系、逐单倍型的 Shiny 报告。
./gu.sh final --grch 37

# 6. 启动 Shiny，打开 PhyML evidence。
./gu.sh shiny
```

若先查看 PhyML 和现有 IBDmix/TRACE，可在第 1 步后直接执行第 5–6 步。补齐 TRACE 后再次运行 `final` 并重启 Shiny。无需手工调用报告脚本。

ARG 必须覆盖完整染色体，不能只对八个小区间建树。第 3 步通常最耗时；已存在且通过校验的同参数产物可复用。Needle 的样本选择、时间标定和非 PAR X 处理会影响 TRACE 的证据解释。TRACE 的 Ghost/Unknown 不负责确定 Neanderthal 或 Denisovan 来源。

Shiny 地址默认为 http://127.0.0.1:3838 。

- 谱系汇总：每个位点 × Neanderthal/Denisovan 一行，包括无候选或不可评估的结果。
- 单倍型明细：所有通过筛选的 H 类型，以及明确标记的序列比较对象；显示指定 SNP 和最佳 tag SNP。
- 方法证据：任意重叠人数、覆盖 ≥80% 人数、候选人数，TRACE 另列同一染色体拷贝的支持数。
- 候选树：随选定谱系切换，展示 bootstrap、purity 与 sensitivity。

`partial_overlap` 表示存在片段重叠但未达到80%覆盖，`not_run` 表示方法未完成，`not_evaluable_scope_or_panel` 表示候选不在可确认的运行区间或样本集合内，`exploratory` 表示该方法结果仅供探索。这些状态不能当作生物学阴性。
