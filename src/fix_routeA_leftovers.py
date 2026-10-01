# -*- coding: utf-8 -*-
"""
fix_routeA_leftovers.py — 修复路线 A 改写后的遗留问题
========================================================
1. 摘要 Availability：改为指向已发布代码
2. 删除 "Plain language summary"（BMC 不用）
3. Table 1 脚注自指（"see Table 1 footnote"）→ 删除
4. 旧 Results 正文中的失效数字（BLOSUM / 0.392 / 0.376 / 0.364 / 0.343 / 0.164）→ 用新数字重写
5. 旧 Figure 1 图注（六族含 BLOSUM）→ 改写为与新数据一致
6. Data & Code Availability：更新
"""
from docx import Document

DOC = r'D:\GED_mutation\manuscript_routeA.docx'

NEW_AVAIL_ABS = (
    'Availability. Code, trained-model configurations and the analysis scripts that reproduce '
    'every number and figure in this paper are archived at Zenodo (DOI to be assigned on '
    'acceptance) and maintained at GitHub. All datasets are publicly available (PDB, MegaScale, '
    'ThermoMutDB, S669 and ssym).')

NEW_AVAIL_DECL = (
    'All datasets used in this study are publicly available: PDB structures (RCSB Protein Data '
    'Bank), MegaScale ([14]), ThermoMutDB ([17]), and the S669 and ssym benchmarks ([9], [3]). '
    'Code, model configurations and the scripts that regenerate every table and figure are '
    'released under an MIT licence and archived at Zenodo; the repository also contains '
    'exploratory graph-edit-distance code that is not differentiable and contributes to no result '
    'reported here.')

# 旧 Results 4.1 正文（含失效数字）→ 新版本
NEW_RES_41 = (
    'On S669, the five-feature physicochemical ridge baseline reached r = 0.390 (95% CI 0.264 to '
    '0.494). Graph encoders spanned a wide range depending on both definition and capacity: the '
    '5.8 k-parameter global-pooling GCN reached 0.026-0.190, the 5.9 k-parameter mutation-site-'
    'pooled GCN 0.310-0.381, the 50 k-parameter edge-aware GINE 0.328-0.418, the 460 k-parameter '
    'attention-pooled deep GINE 0.300-0.394, and the 509 k-parameter equivariant network -0.081 to '
    '0.393. The highest single value among graph encoders was r = 0.402, attained by the '
    '50 k-parameter GINE on the C\u03b2 definition, which exceeds the physicochemical baseline. '
    'Model performance was therefore not monotone in capacity: a 50 k-parameter encoder with '
    'explicit edge features outperformed a 460 k-parameter attention-pooled encoder on the same '
    'definition.')

NEW_RES_41B = (
    'On the secondary ssym benchmark the physicochemical baseline reached r = 0.305 (95% CI 0.141 '
    'to 0.470). The sequence baselines and the fusion exceeded it: ESM-2 zero-shot r = 0.266, the '
    'supervised ESM-2 head 0.458-0.516, and the ESM-2 plus contact-graph fusion 0.540. The relative '
    'ordering of information sources is thus benchmark-dependent, whereas the contribution of the '
    'contact-graph representation is small on both benchmarks.')

NEW_FIG1_CAP = (
    'Figure 1 | Study design and evaluation protocol. Wild-type structures were paired with '
    'FoldX-modelled mutant structures under four residue contact-graph definitions (C\u03b1, '
    'C\u03b2, side-chain centroid, all-atom) and processed identically. Training combined '
    'MegaScale with ThermoMutDB after BLAST leakage filtering, sampled 1:1 to 7,905 examples. '
    'Evaluation used the common 511-mutation S669 intersection as the exploratory benchmark and '
    'the 342-mutation ssym set as the development-independent secondary benchmark. Five encoder '
    'capacities, four graph definitions, a physicochemical ridge baseline, an ESM-2 650M sequence '
    'baseline and a sequence-structure fusion were compared under one protocol; all paired '
    'comparisons used protein-cluster bootstrap intervals.')


def main():
    d = Document(DOC)
    n = {'abs': 0, 'plain': 0, 'selfref': 0, 'res41': 0, 'res41b': 0,
         'fig1': 0, 'decl': 0}

    for p in d.paragraphs:
        t = p.text.strip()

        # 1) 摘要 Availability
        if t.startswith('Availability.'):
            p.runs[0].text = NEW_AVAIL_ABS
            for r in p.runs[1:]:
                r.text = ''
            n['abs'] += 1

        # 5) 旧 Figure 1 图注
        elif t.startswith('Figure 1 | Leakage-aware comparison'):
            p.runs[0].text = NEW_FIG1_CAP
            for r in p.runs[1:]:
                r.text = ''
            n['fig1'] += 1

        # 4) 旧 Results 4.1 正文（含失效数字）
        elif t.startswith('For S669, the five-feature physicochemical ridge baseline'):
            p.runs[0].text = NEW_RES_41
            for r in p.runs[1:]:
                r.text = ''
            n['res41'] += 1

        elif t.startswith('On the secondary ssym benchmark, local centroid pooling'):
            p.runs[0].text = NEW_RES_41B
            for r in p.runs[1:]:
                r.text = ''
            n['res41b'] += 1

        # 6) Declarations Availability
        elif t.startswith('All datasets used in this study are publicly available'):
            p.runs[0].text = NEW_AVAIL_DECL
            for r in p.runs[1:]:
                r.text = ''
            n['decl'] += 1

        # 3) Table 1 脚注自指
        elif t.startswith('Table 1 |'):
            newt = t.replace(' (see Table 1 footnote and Figure 3)', ' (Figure 3)')
            if newt != t:
                p.runs[0].text = newt
                for r in p.runs[1:]:
                    r.text = ''
                n['selfref'] += 1

    # 2) 删除 "Plain language summary" 及其后正文段
    def drop(par):
        par._p.getparent().remove(par._p)

    to_drop = []
    for i, p in enumerate(d.paragraphs):
        if p.text.strip() == 'Plain language summary':
            to_drop.append(i)
    # 该标题后通常是若干段正文，删除到下一个明显标题
    for i in to_drop:
        ps = d.paragraphs
        j = i + 1
        while j < len(ps) and ps[j].text.strip() and not ps[j].text.strip().startswith(
                ('Keywords', '1. Introduction', 'Abstract')):
            drop(ps[j])
            ps = d.paragraphs
        drop(d.paragraphs[i])
        n['plain'] += 1

    d.save(DOC)
    print('修复统计:', n)


if __name__ == '__main__':
    main()
