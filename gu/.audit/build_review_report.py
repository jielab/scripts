from pathlib import Path
import csv, json

root = Path('D:/analysis/gu/threshold_review/20260906')
current = root / 'legacy_compatible'
original = Path('D:/analysis/gu/phyml/1kg')

def rows(p):
    with p.open(encoding='utf-8') as f:
        return list(csv.DictReader(f, delimiter='\t'))

order = [('1','rs12405323','Neanderthal'), ('3','rs35044562','Neanderthal'),
         ('6','rs9349320','Denisovan'), ('8','rs1041791','Neanderthal'),
         ('9','rs8176719','Neanderthal'), ('12','rs417915','Denisovan'),
         ('19','rs7412','Neanderthal'), ('X','rs6525167','Neanderthal')]
lines = ['# GU PhyML 八个位点复核与旧阈值兼容验证', '',
         '复核日期：2026-09-06。旧表实际文件为 `D:/files/gu.b37.xlsx`（Sheet1，B1:F16）；共五个性状位点，三个对照不在该表内。', '',
         '已修改代码默认阈值，并在独立目录完成八个位点的证据层重算与候选树分析。原始 `phyml/1kg`、IBDmix、`gu.sqlite` 与 normalize 输出保留，便于前后核对；本目录的新结果尚未替换 Shiny 所读的旧结果。区域发现阶段完全复用原结果，8×8 个区域压缩表的 SHA-256 均一致。', '',
         '**结论：可以恢复多个旧信号的候选资格，但不能声称旧 Excel 的五个渐渗结论已全部复现。** 无风险标签的完整区域分析与旧的 risk-LD SNP 子集分析不是同一个统计问题；阈值相同也不会保证树拓扑相同。旧Excel也没有列出bootstrap值，其Inferred archaic lineage不能等同于通过当前新增的候选树确认层。', '',
         '|位点|原严格配置的区域判定|兼容配置：旧表对应谱系|常见候选拷贝|该谱系候选树 bootstrap|该谱系判定|候选携带旧风险等位基因|',
         '|---|---|---|---:|---:|---|---|']
for chrom, rs, lineage in order:
    d = current / f'chr{chrom}_{rs}'
    old = rows(original / d.name / 'final/evidence_loci.tsv')[0]
    new = rows(d / 'final/evidence_loci.tsv')[0]
    assert new['candidate_tree_status'] != 'pending', d
    hs = [r for r in rows(d / 'final/evidence_haplotypes.tsv') if r['diagnostic_lineage'] == lineage and r['diagnostic_candidate_pass']=='1']
    tree = next((r for r in rows(d / 'final/evidence_trees.tsv') if r['expected_lineage']==lineage), {})
    risk = next((r for r in rows(d / 'final/gwas_risk_links.tsv') if r.get('lineage')==lineage), {})
    n = sum(int(r['n']) for r in hs)
    call = ('not_evaluable' if new['introgression_call']=='not_evaluable' else 'no_candidate') if not hs else ('supported_candidate' if tree.get('candidate_clade_pass')=='1' else 'candidate_sequence_only')
    rr = f"{risk['risk_allele']}: {risk['n_candidate_risk_copies']}/{risk['n_candidate_anchor_called_copies']}" if risk and n else '—'
    lines.append(f"|chr{chrom} {rs}|{old['introgression_call']}|{lineage}|{n}|{tree.get('candidate_clade_bootstrap') or '—'}|{call}|{rr}|")

lines += ['', '拷贝数不是人数；表中是通过兼容阈值的常见完整区域 H types。旧 Excel 第5行统计所有携带 lead 风险等位基因的拷贝，所以不能要求它等于本表的候选拷贝数。chr6 按旧表的 Denisovan 谱系展示；该区域另有 Neanderthal 候选，单行 locus 摘要可能优先选择数量更多的谱系，风险核对表逐谱系保留结果。', '',
'''**各位点的解读**

- chr1：原来45个候选 H types、692 copies，候选树 bootstrap=100，但纯度45/57=78.95%，刚好低于新增80%门槛。607名序列候选携带者全部与同谱系 IBDmix 在候选范围重叠。兼容配置使用 n>10 的常见 H types，并取消旧代码没有的纯度硬门槛，恢复树支持候选。风险 A 的总拷贝737与旧表完全一致；支持保留 Neanderthal 渐渗候选。
- chr3：旧严格配置394 copies、351人，候选树73%，351人均有 IBDmix 重叠。兼容设置改用4个常见 H types、363 copies后，同版本的新focused tree bootstrap=66%，低于保留的70%确认门槛，因此自动判定为candidate_sequence_only。未为通过阳性对照而降低门槛或挑选随机重复；原全古参考区域树bootstrap=99%且既有IBDmix和文献支持仍在，不能将本次focused筛选降级解释为没有Neanderthal渐渗。100次bootstrap本身有抽样波动，种子记录在stats文件中。
- chr6：旧风险 G 总拷贝55完全复现。修改前 Denisovan 候选48 copies，全部携带 G，48人均有 Denisovan IBDmix 重叠。古参考在 rs9349320 却是 T，故原 `named_anchor_matches_lineage=0` 不表示“与 GWAS 风险不一致”。新的逐谱系风险表解决这个语义问题。区域可能有不止一种 archaic-like 单倍型，应把风险 G 的 Denisovan 候选与其他候选分开解释。
- chr8：风险 C 的978 copies与旧表完全一致；旧表5个 diagnostic matches在当前常见风险 H types中仍可看到。原流程因5/7=71.43%<80%拒绝候选。全古参考的常见区域树也已有10个H types、747 copies与三个 Neanderthal 聚类（bootstrap=71）。兼容配置恢复候选并重建了focused tree。但当前 IBDmix 在这一核心候选区没有同谱系重叠，故这是有系统发育支持、仍缺独立 tract caller 佐证的候选，不能称为已确认渐渗。
- chr12：风险 A 的1250 copies完全复现。常见风险H types匹配27/34=79.41%，原80%门槛将其排除；匹配34/34的两个H types则各只有1 copy。兼容配置恢复12个常见H types、1016 copies，其中1013携带A。但是全古参考focused tree没有合格的独立 Denisovan–candidate 边；仅含Denisovan的常见区域树可以给出100%，而加入其他古参考后谱系归属发生歧义。因此“Denisovan-like序列”有依据，“已证实Denisovan渐渗”尚未复现。原27/34组1010名携带者中127人与Denisovan IBDmix在该候选范围重叠，不能用normalize中大量Neanderthal重叠替代Denisovan支持。
- chrX：rs6525167实际存在于GRCh37 X:66423003 G>A，保存ID为 `X:66423003:G:A`，不是rsID。现已通过build+坐标+REF/ALT精确映射修复。当前只分析1233名男性，A=837 copies；旧3414不能直接与male-only分母相比。源pvar INFO给出全人群AC=2577，且2577+837恰为3414，与旧计算中男性按两份计数的情形一致，但未取得旧运行VCF，不能把这一推测当成已核实的来源。风险A也被INFO/AA标注为祖先状态，单个位点的 archaic match不能证明渐渗。兼容配置能保留更多区域候选，但候选树未建立可靠的Neanderthal归属；旧结论仍未完整复现。
- ABO：原严格配置检出H91、6个东亚copies，7/7 diagnostic matches、bootstrap99%，6人均有IBDmix重叠。但其观察到的diagnostic span仅15,735 bp，用旧ILS模型p约0.20，未通过<0.1；且n=6不符合旧常见H阈值。兼容配置不再保留这一正式候选。它适合作为对新增流程灵敏度的检查，但不能先假定整个ABO区域绝无古老共享谱系。
- APOE：没有通过derived+outgroup条件的diagnostic位点，仍为not_evaluable。它没有提供渐渗支持，但“不可评估”并不是统计证明无渐渗。

**与真实 bald GWAS 的核对（GRCh38）**

|标记|GRCh38位置|EA|BETA|原始P|
|---|---|---|---:|---:|
|rs1041791|8:108099842|C|0.0407958|5.79696e-38|
|rs417915|12:52448290|A|0.0143388|6.24799e-08|
|rs6525167|X:67203161（文件CHR=23）|A|0.227762|文件为0，数值精度限制|

这三个风险方向与Excel一致。rs417915原始P略高于常用5e-8；COJO在附近选择12:52471808 A>T，bJ=0.0233744、pJ=1.69382e-16。chr8附近也有多个条件信号，包括8:108090690 G>A、pJ=4.394e-09。原三个rsID对应位置均不在当前jma的独立标记表；该jma文件仅含chr1–22，没有X，所以不能据此说X“条件分析后不显著”。这些结果支持保留区域先行设计，并在末端另做LD/条件信号或共定位核对。单个lead的风险一致性还不是同一因果信号的证明。

**实现与阈值依据**

旧GitHub版本默认simple模式：`pick_lineage_simple`只要求ILS p<0.1并按谱系合并的risk-allele匹配比例排序，没有80%匹配率下限；`make_phy`用n>10的常见H types。YRI上限0.05、risk-vs-nonrisk allele frequency delta>0.5是旧diagnostic定义的一部分。

新默认profile为`legacy_compatible`：n>=11；diagnostic至少1个，比例门槛0；取消额外50kb gap、80% purity和50% sensitivity门槛；对观察到的diagnostic span应用旧ILS p<0.1。零diagnostic匹配不能通过，短片段由ILS门槛排除。这个统一的0.53 cM/Mb旧模型没有使用各位点的局部重组图，尤其不能将X或Denisovan结果的p值当作经过单独校准的最终ILS检验。原`strict`参数可通过 `PHYML_THRESHOLD_PROFILE=strict`恢复，显式数值环境变量优先于profile。HKY85、100 bootstrap，以及新流程的bootstrap>=70树确认层保留。

没有把旧risk enrichment delta硬塞进前端：这需要先知道GWAS风险组，会违反区域先行约束。没有取消新流程的archaic-derived检查，没有更改区域SNP集合，也没有把male-only X改回男女混合或男性双拷贝。旧Excel的61/156/26/89/157个LD选择SNP与当前完整区域SNP不同，不能直接比较匹配百分比及diagnostic数量。

因此本次是有明确出处的阈值兼容与实测复核，不能包装成五个阳性结论的逐值重放。下一步若要解决chr12/X，应核对旧LD SNP清单、祖先等位基因及古参考的缺失模式，并检查完整区域内的重组分段；不应继续为了得到阳性而任意改阈值。

**验证和文件**

- 6项回归测试通过，包括旧ILS数值与稳定尾概率、显式参数覆盖、strict profile、gap开关、rsID alias的build/allele防错、不同谱系的独立风险核对。
- chr1和chr8使用strict profile重算后，每个haplotype×lineage的pass/fail与修改前一致。
- 8个位点共64个第一阶段region压缩表SHA-256一致，发现结果没有被风险信息重新筛选。
- 8个位点的证据准备、候选树重算、汇总与风险核对已执行。最终树均使用与原分析相同的GU环境PhyML 3.3.20260528、HKY85、100次bootstrap，单核每棵树、最多四棵同时运行。系统默认旧版PhyML的试运行已排除，其产物另标为.system-version，不用于本报告。随机种子保存在各stats文件中。
- 代码：`D:/scripts/gu/f/phyml_thresholds.py`、`phyml_evidence.py`、`phyml_anchor.py`、`phyml_risk.py`、`phyml.sh`。精确标记及风险说明：`f/phyml_markers.tsv`。
- 每个位点的新结果在本报告同级 `legacy_compatible/chr*/final/`；重点看 `evidence_loci.tsv`、`evidence_trees.tsv`、`evidence_haplotypes.tsv`、`evidence_parameters.json` 和新增 `gwas_risk_links.tsv`。
- 原始GWAS：`D:/data.BIG/gwas/main/common/bald/gwas/bald.gz`；条件分析：同目录 `bald.jma.cojo`。

方法参考：[旧分析代码](https://github.com/zhaoran7/Bald/tree/main/archaic%20introgression)；[Zeberg与Pääbo，Nature 2020](https://www.nature.com/articles/s41586-020-2818-3)。该文章结合系统树与片段长度/ILS论证；bootstrap与ILS筛选概率都不是渐渗后验概率。
''']
(root / 'gu.phyml.review.md').write_text('\n'.join(lines), encoding='utf-8')
print(root / 'gu.phyml.review.md')
