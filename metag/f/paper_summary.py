"""Render a ranked summary from the published tables; never infer truth from counts."""
from collections import Counter, defaultdict
import csv
import heapq
import html
from paper_layout import NAMES


def records(path):
    with path.open(encoding='utf-8', newline='') as f:
        yield from csv.DictReader(f, delimiter='\t')


def table(headers, values):
    def cells(tag, row):
        return '<tr>' + ''.join(f'<{tag}>{html.escape(str(v))}</{tag}>' for v in row) + '</tr>'
    return '<table><thead>'+cells('th', headers)+'</thead><tbody>'+''.join(cells('td', r) for r in values)+'</tbody></table>'


STYLE = ('<meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">'
         '<style>body{font:16px/1.55 system-ui;max-width:78rem;margin:2rem auto;padding:0 1rem;color:#163449}'
         'section{overflow:auto}table{border-collapse:collapse;margin:1em 0}td,th{border-bottom:1px solid #ced8db;'
         'padding:.3em .65em;text-align:left;white-space:nowrap}th{background:#e9f3f1}a{color:#006056}'
         'h2{margin-top:1.7em}li{margin:.35em 0}</style>')


def target_panel_qc(folder, samples):
    """Show panel opportunity for human and every plant species called in this report."""
    taxa = {'9606': ('human', 'Homo sapiens')}
    calls = Counter()
    for r in records(folder/'target_species_evidence.tsv'):
        if r['group'] == 'plant' or r['taxid'] == '9606':
            tid = r['taxid']
            taxa[tid] = (r['group'], r['name'])
            calls[r['sample'], tid] += int(r['pairs'])
    panels = {(r['sample'], r['taxid']): r for r in records(folder/'panels.tsv') if r['taxid'] in taxa}
    omitted = {(r['sample'], r['taxid']): r for r in records(folder/'panel_omitted.tsv') if r['taxid'] in taxa}
    counts = defaultdict(Counter)
    fields = ('sample', 'input_pairs', 'group', 'taxid', 'name', 'in_panel', 'eligible',
              'kraken_pairs', 'anchor_windows_capped', 'split_min_windows', 'exact_pairs', 'status')
    with (folder/'target_panel_qc.tsv').open('w', encoding='utf-8', newline='') as out:
        writer = csv.writer(out, delimiter='\t', lineterminator='\n')
        writer.writerow(fields)
        for sample, row in samples.items():
            dataset = sample.split('/', 1)[0]
            total = int(row['pairs'])
            threshold = row.get('rank_min_windows', 'NA')
            for tid, (group, name) in sorted(taxa.items(), key=lambda item: int(item[0])):
                key = sample, tid
                panel = panels.get(key)
                missing = omitted.get(key, {})
                n = calls[key]
                eligible = bool(panel and int(panel['near_relative_and_background']))
                anchors = int(panel['agreeing_nuclear_windows_capped']) if panel else None
                if not total:
                    status = 'empty_input'
                elif n:
                    supported = eligible and (threshold == 'NA' or anchors >= int(threshold))
                    status = 'exact_call' if supported else 'check_exact_call'
                elif panel is None:
                    status = missing.get('reason', 'not_in_panel')
                elif not eligible:
                    status = 'no_eligible_comparison'
                elif anchors == 0:
                    status = 'no_agreeing_windows'
                elif threshold != 'NA' and anchors < int(threshold):
                    status = 'below_window_threshold'
                else:
                    status = 'windows_without_exact_call'
                writer.writerow([sample, total, group, tid, name, int(panel is not None), int(eligible),
                                 (panel or missing).get('kraken_pairs', 'NA'),
                                 anchors if anchors is not None else 'NA', threshold, n, status])
                c = counts[dataset, group, tid, name]
                if total:
                    c['nonempty'] += 1
                    c['panel'] += panel is not None
                    c['eligible'] += eligible
                    c['positive'] += n > 0
                    c['pairs'] += n
    columns = ('nonempty', 'panel', 'eligible', 'positive', 'pairs')
    summary = [(*key, *(c[k] for k in columns)) for key, c in sorted(counts.items())]
    with (folder/'target_panel_by_dataset.tsv').open('w', encoding='utf-8', newline='') as out:
        writer = csv.writer(out, delimiter='\t', lineterminator='\n')
        writer.writerow(('dataset', 'group', 'taxid', 'name', 'nonempty_samples', 'panel_samples',
                         'eligible_samples', 'exact_call_samples', 'exact_pairs'))
        writer.writerows(summary)
    return summary


def render(folder):
    samples = {}
    for row in records(folder/'samples.tsv'):
        name = row['sample']
        if not name or name in samples:
            raise ValueError('Missing or duplicate sample: ' + name)
        samples[name] = row
    if not samples:
        raise ValueError('Empty report')
    def read(kind):
        for row in records(folder/(kind+'.tsv')):
            name = row['sample']
            if name not in samples:
                raise ValueError(f'{kind}: unknown sample {name}')
            yield name, row
    totals = defaultdict(Counter)
    group, rank, exact, initial = (defaultdict(Counter) for _ in range(4))
    plants, panels, omitted, quality = defaultdict(list), {}, Counter(), defaultdict(Counter)
    panel_size, plant_comparison = Counter(), Counter()
    target_support = Counter()
    for name, row in read('target_species_evidence'):
        target_support[name, row['group'], row['evidence']] += int(row['pairs'])
    for name, row in read('taxa'):
        n, g, level = int(row['pairs']), row['group'], row['rank']
        if n < 0:
            raise ValueError('Negative bin count')
        group[name][g] += n
        rank[name][level] += n
        if level == 'species':
            exact[name][g] += n
            if g == 'plant':
                heapq.heappush(plants[name], (n, row['name'], row['taxid']))
                if len(plants[name]) > 10:
                    heapq.heappop(plants[name])
    for kind, field, counter in [('groups','group',group), ('ranks','rank',rank)]:
        declared = defaultdict(Counter)
        for name, row in read(kind):
            if row[field] in declared[name] or int(row['total']) != int(samples[name]['pairs']):
                raise ValueError('Inconsistent ' + kind)
            declared[name][row[field]] = int(row['pairs'])
        if any(declared[s] != counter[s] for s in samples):
            raise ValueError(kind + ' disagrees with taxa')
    final_groups = defaultdict(Counter)
    for name, row in read('transitions'):
        n = int(row['pairs'])
        initial[name][row['kraken_group']] += n
        final_groups[name][row['final_group']] += n
    for kind in ('evidence', 'taxon_evidence'):
        for name, row in read(kind):
            totals[name][kind] += int(row['pairs'])
    for name, row in read('panels'):
        panel_size[name] += 1
        if row['role'] in ('plant_candidate', 'plant_relative'):
            plant_comparison[name, 'selected'] += 1
            plant_comparison[name, 'eligible'] += int(row['near_relative_and_background'])
        if row['taxid'] == '9606':
            panels[name] = row['agreeing_nuclear_windows_capped']
    for name, row in read('panel_omitted'):
        omitted[name] += 1
    qc = (folder/'fastq_qc.tsv').is_file()
    if qc:
        for name, row in read('fastq_qc'):
            acc = quality[name, row['mate']]
            bases = int(row['bases'])
            acc['bases'] += bases
            acc['pairs'] += int(row['pairs'])
            for key in ('q30_fraction','gc_fraction','n_fraction','mean_phred'):
                acc[key] += bases * float(row[key])
    overview, targets, details, human_note, qrows = [], [], [], [], []
    for s, row in samples.items():
        n = int(row['pairs'])
        sp = sum(exact[s].values())
        unc = group[s]['unclassified']
        higher = n - sp - unc
        if (sum(group[s].values()) != n or final_groups[s] != group[s] or
                any(totals[s][k] != n for k in ('evidence','taxon_evidence')) or
                (sp,higher,unc) != tuple(int(row[k]) for k in ('species_pairs','higher_rank_pairs','unclassified_pairs')) or
                any(group[s][g] != int(row[g+'_pairs']) for g in ('human','plant','virus'))):
            raise ValueError('Inconsistent totals for ' + s)
        def amount(x):
            return f'{x:,} ({100*x/n:.3f}%)' if n else '0'
        overview.append((s,f'{n:,}',amount(sp),amount(higher),amount(unc),row['integrity'],row['truth_pairs']))
        for g in ('human','plant','virus'):
            targets.append((s,g,f'{initial[s][g]:,}',f'{group[s][g]:,}',f'{exact[s][g]:,}',
                            f'{group[s][g]-exact[s][g]:,}'))
        human_note.append((s, panel_size[s], panels.get(s,'未纳入'), plant_comparison[s,'eligible'],
                           plant_comparison[s,'selected'], omitted[s]))
        details.extend((s,t,nm,f'{cnt:,}') for cnt,nm,t in sorted(plants[s],reverse=True))
        if qc:
            for mate in ('1','2'):
                acc = quality[s,mate]
                if acc['pairs'] != n:
                    raise ValueError('FASTQ QC total differs: ' + s)
                b = acc['bases']
                qrows.append((s,mate,f'{b:,}',*[(fmt % (factor*acc[k]/b) if b else 'NA')
                    for k,factor,fmt in [('q30_fraction',100,'%.2f%%'),('gc_fraction',100,'%.2f%%'),
                                         ('n_fraction',100,'%.4f%%'),('mean_phred',1,'%.2f')]]))
    def section(title, headers, data):
        return '<section><h2>'+title+'</h2>'+table(headers,data)+'</section>'
    page=['<!doctype html><html lang="zh"><title>METAG ranked 结果</title>'+STYLE+'<h1>宏基因组双端分类结果</h1>',
          '<p>每对 R1/R2 只计一次。分类层级百分比以输入读对为分母，样本检出频率以非空样本为分母。</p>',
          '<h2>论文图表</h2><p>各队列和全部样本的汇总见下表。汇总百分比按读对加权，样本中位数和四分位数按样本等权计算。</p>',
          '<p><a href="table1_classification.tsv">表 1 · 分类层级</a>　'
          '<a href="table2_target_groups.tsv">表 2 · 人源、植物、病毒</a>　'
          '<a href="table3_plant_species.tsv">表 3 · 全部植物物种</a>　'
          '<a href="sample_targets.tsv">逐样本重点类群计数</a></p>',
          *[f'<figure><a href="{stem}.pdf"><img style="max-width:100%;height:auto" '
            f'src="{stem}.png" alt="{caption}"></a><figcaption>{caption}　'
            f'<a href="{stem}.pdf">PDF</a> · <a href="{stem}.png">PNG</a></figcaption></figure>'
            for stem, caption in [
                ('fig1_classification','图 1　队列分类层级与逐样本种级分配率。'),
                ('fig2_target_groups','图 2　人源、植物、病毒的初筛、最终类群及最终种级计数，每百万输入读对。'),
                ('fig3_plant_detection','图 3　植物种级检出频率。至少一对种级 reads 即计为检出，分母为有输入的已完成样本。'
                 '最多展示十二种，按检出样本数、种级读对数、taxid 排序。')]],
          section('分类层级',['样本','总 pairs','物种','较高层级','未分类','完整性检查','真值 pairs'],overview),
          section('重点类群',['样本','类群','全库初筛','最终归入该组','其中到种','其中较高层级'],targets),
          section('参考面板',['样本','已选参考数','人源一致核窗口数（封顶）','有同属竞争和背景的候选及近缘植物',
                           '已选候选及近缘植物','未纳入 taxid 数（所有角色）'],human_note)]
    page.append(section('人和已检出植物的面板覆盖情况',
                        ['队列','类群','taxid','物种','非空样本','纳入面板','可做近缘竞争比较','种级阳性样本','种级读对'],
                        target_panel_qc(folder, samples)))
    versions = Counter((r['database'], r['code'], r.get('rank_min_windows', 'NA')) for r in samples.values())
    page.append(section('各样本 split 记录中的运行版本',
                        ['数据库标识（前 12 位）','split 方法标识（前 12 位）','核窗口门槛','样本数'],
                        [(db[:12], code[:12], minimum, n) for (db, code, minimum), n in sorted(versions.items())]))
    if len(versions) > 1:
        page.append('<p>这批样本的数据库、方法标识或门槛存在差异。合并解释前请在 samples.tsv 中核对各样本的运行记录。</p>')
    page.append('<p>图 3 的分母是全部非空样本。不同样本的植物面板由初筛结果和预算决定。'
                '符合近缘竞争条件的样本数描述检测机会，无法据此估计人群阳性率或分析灵敏度。'
                '窗口数为通过双端 Kraken 与核比对一致性检查的不同 10 kb 窗口数，达到门槛后封顶。'
                'split 门槛取自各样本的运行记录，历史汇总缺少此项时记为 NA。</p>'
                '<p><a href="target_panel_by_dataset.tsv">下载队列分母</a>　'
                '<a href="target_panel_qc.tsv">下载逐样本面板、窗口与种级读对</a></p>')
    if target_support:
        page.append(section('重点类群的种级证据',['样本','类群','判定依据','pairs'],
                            ((*key, f'{n:,}') for key, n in sorted(target_support.items()))))
        page.append('<p>nuclear 和 refined_nuclear 表示通过核序列竞争比对。'
                    'both_mates_kraken 表示两端的 Kraken 证据支持同一物种。'
                    '窗口数在达到设定门槛后停止累计。</p>')
    if details:
        page.append(section('植物种级调用（每样本最多十种）',['样本','taxid','参考物种','pairs'],details))
    if qc:
        page.append(section('FASTQ 质量（按碱基加权）',['样本','mate','碱基','Q30','GC','N','平均 Phred'],qrows))
    if (folder/'truth_metrics.tsv').is_file():
        cols='sample method rank label tp fp fn precision recall f1'.split()
        page.append(section('独立真值子集',['样本','方法','层级','标签','TP','FP','FN','precision','recall','F1'],
            ([r[k] for k in cols] for r in records(folder/'truth_metrics.tsv') if r['label'] in ('ALL','human','plant'))))
    else:
        page.append('<p>本报告没有逐对真值。归类比例和 FASTQ 质量无法估计物种归属的精确率或召回率。人源为零表示没有 read pair 达到最终人源标准。</p>')
    if (folder/'refine_status.tsv').is_file():
        import sys
        from pathlib import Path
        sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'refine'))
        from metag_refine_tables import summarize
        pooled_sources = summarize(folder)
        states = list(records(folder/'refine_status.tsv'))
        count = sum(r['state']=='complete' for r in states)
        page.append(f'<h2>来源复核与公共卫生</h2><p>已关联 {count}/{len(states)} 个样本的 refine 结果。'
                    '各样本的建库类型、BLAST 与 AMR 是否执行见下表。</p>')
        page.append(table(['样本','来源','建库','refine','BLAST','AMR'],
                          ([r[k] for k in ('sample','sample_type','library_type','state','blast','amr')] for r in states)))
        page.append('<p>食源标签表示兼容的 DNA 来源。牛源读段也可能来自乳制品、环境或加工环节。'
                    '病毒表仅核对已有种级标签。AMR 表按原始读对计数，每对最多分配给一个家族，分母是全部输入读对。'
                    'RNA 建库的 AMR 命中反映转录本证据。短读段的复核一致率受参考库范围及抽样策略影响。</p>')
        page.append('<p><a href="refine_sources_by_dataset.tsv">各队列全部来源计数和分母</a>　'
                    '<a href="refine_blast_by_bin.tsv">各 bin 的 BLAST 复核结果和分母</a></p>')
        display = [r for r in pooled_sources if r[4]=='species' and r[3] in ('plant','mammal','other_vertebrate')]
        if display:
            page.append(table(['队列','样本类型','建库','类群','原始种级标签（最多展示 20 行）','读对','有标签的样本','同组已复核非空样本'],
                              ([r[i] for i in (0,1,2,3,6,8,9,10)] for r in display[:20])))
        links = [('sources','候选来源及原有证据'),('blast','双端 BLAST 复核及抽样分母'),
                 ('amr','AMR 家族及归一化结果'),('virus','公共卫生病毒标签'),('status','执行状态与版本')]
        page.append('<ul>'+''.join(f'<li><a href="refine_{stem}.tsv">{label}</a></li>'
                                  for stem,label in links)+'</ul>')
    if (folder/'seal_status.tsv').is_file():
        states=list(records(folder/'seal_status.tsv'))
        done=sum(r['state']=='complete' for r in states)
        page.append(f'<h2>威德尔海豹序列复核</h2><p>本次选择中完成 {done}/{len(states)} 个样本。原分箱与序列复核分别报告。</p>')
        if done:
            fields=['sample','bin_taxid','tested_pairs','target_nuclear','unique_nuclear_fragments','nuclear_windows_10kb','separated_nuclear_loci_10kb','interpretation']
            page.append(table(fields,([r[k] for k in fields] for r in records(folder/'seal_bins.tsv'))))
        page.append('<p><a href="seal_bins.tsv">全量分箱复核</a> · <a href="seal_loci.tsv">跨样本位点</a> · <a href="seal_recurrent_loci.tsv">跨样本重复的参考区段</a> · <a href="seal_status.tsv">执行状态</a></p>')
    page.insert(1,'<p><a href="explore.html">交互报告与 SVG 图</a> · <a href="species_recurrence.tsv">各队列物种检出与最高贡献样本</a></p>')
    page.append('<h2>明细表</h2><ul>'+''.join(f'<li><a href="{n}">{n}</a></li>' for n in NAMES if (folder/n).is_file())+'</ul>')
    (folder/'report.html').write_text('\n'.join(page),encoding='utf-8')
    return list(samples)
