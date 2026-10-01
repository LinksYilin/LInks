# -*- coding: utf-8 -*-
"""
apply_routeA.py — 按路线 A 改写论文（标题/摘要/新增 Results/Discussion 重构/Table 1）
======================================================================================
不删除任何原有内容；新增段落在指定位置插入。
保存为新文件 manuscript_routeA.docx，原文件保留。
"""
import os
import shutil
import sys

from docx import Document
from docx.oxml.ns import qn
from docx.oxml import OxmlElement
from docx.text.paragraph import Paragraph

SRC = r'D:\GED_mutation\manuscript_polished_final_20260930.docx'
DST = r'D:\GED_mutation\manuscript_routeA.docx'


def insert_after(par, text, style=None, italic=False):
    new_p = OxmlElement('w:p')
    par._p.addnext(new_p)
    np = Paragraph(new_p, par._parent)
    if style:
        try:
            np.style = style
        except Exception:
            pass
    r = np.add_run(text)
    r.italic = italic
    return np


def set_text(par, text):
    if par.runs:
        par.runs[0].text = text
        for r in par.runs[1:]:
            r.text = ''
    else:
        par.add_run(text)


def find(paras, pred):
    for i, p in enumerate(paras):
        if pred(p.text.strip()):
            return i, p
    return None, None


# ---------------------------------------------------------------- 正文内容
TITLE = ('What does structure contribute to \u0394\u0394G prediction? A controlled '
         'information-budget analysis across residue contact-graph representations')

ABSTRACT = {
    'Motivation.': (
        'Residue contact graphs are standard in structure-based \u0394\u0394G prediction, and their '
        'definition determines which mutation-induced side-chain changes a model can see. A more '
        'sensitive representation is usually assumed to yield better predictions. We tested that '
        'assumption by separating what a definition reports about structural change from what it '
        'changes about prediction.'),
    'Results.': (
        'Processing wild-type and modelled mutant structures identically, C\u03b1 graphs recorded no '
        'contact change in any of 505 quality-controlled pairs, whereas C\u03b2, side-chain-centroid '
        'and all-atom graphs changed in 9.9%, 90.7% and 81.6% of pairs. Prediction was statistically '
        'indistinguishable across definitions at every rung of a five-encoder capacity ladder '
        'spanning a 5.8 k-parameter mean-pooled GCN to a 509 k-parameter E(3)-equivariant network: '
        'none of eight paired comparisons reached significance, and the mean absolute effect did not '
        'grow with capacity (mean |\u0394r| \u2248 0.03-0.05). Adding a side-chain-centroid graph to an '
        'ESM-2 650M sequence baseline produced no reliable increment (S669 \u0394r = -0.013, 95% CI '
        '-0.036 to +0.010; ssym \u0394r = +0.025, 95% CI -0.007 to +0.070), excluding structural '
        'increments above |\u0394r| \u2248 0.04. Five physicochemical features alone reached r = 0.390 '
        'on S669; the one structural quantity with univariate signal, the count of broken hydrophobic '
        'contacts (r = 0.135, 95% CI 0.044 to 0.201), added no incremental information over those '
        'features (\u0394r = -0.009, 95% CI -0.047 to +0.037).'),
    'Conclusion.': (
        'The atom-level definition of a residue contact graph strongly determines what structural '
        'change a model describes, but it did not measurably change \u0394\u0394G prediction across the '
        'encoder capacities and benchmarks tested, nor did structure add a reliable increment over a '
        'modern sequence baseline. We provide the audit protocol and the resulting information '
        'budget as a template for deciding where representation-layer effort is likely to pay off.'),
}

# 新增 Results 小节（插在 4.1 之后）
NEW_RESULTS = [
    ('4.2 A capacity ladder: the representation effect does not grow with model size',
     'A single encoder family cannot distinguish "the representation does not matter" from "this '
     'encoder cannot exploit it". We therefore trained five encoders of increasing capacity on '
     'identical data and splits, and repeated the four-definition comparison at each rung: a '
     'mean-pooled GCN (5.8 k parameters), a mutation-site-pooled GCN (5.9 k), a GINE model with '
     'explicit edge features (50 k), a six-layer GINE with attention pooling and a mutation-site '
     'local readout (460 k), and an E(3)-equivariant network with coordinate updates (509 k). '
     'Every model was trained for 20 epochs with the same optimisation settings and two seeds, and '
     'representation effects were estimated as within-backbone paired differences against the '
     'side-chain-centroid definition using a protein-cluster bootstrap.\n'
     'Across the four stable encoders (5.8 k to 460 k parameters), none of the eight paired '
     'definition comparisons reached significance, and the mean absolute effect was 0.052, 0.052, '
     '0.031 and 0.039 respectively (Figure 3b). The effect therefore does not increase with '
     'capacity, which is the pattern expected if encoder strength were the bottleneck. The '
     '509 k-parameter equivariant network gave a larger apparent effect (mean |\u0394r| = 0.105, one '
     'of two comparisons nominally significant at P = 0.035), but that encoder collapsed on a '
     'subset of seeds in every definition tested (for the side-chain-centroid graphs, r = 0.393 and '
     '-0.081 for two seeds), so its comparison is not interpretable as a representation effect. We '
     'report it as an instability rather than as evidence for a capacity-dependent effect.\n'
     'Model performance itself was not monotone in capacity. The 50 k-parameter GINE encoder '
     'reached r = 0.402 on the C\u03b2 definition, the highest value among the graph encoders and '
     'above the physicochemical ridge baseline (r = 0.390), whereas the 460 k-parameter GINE '
     'reached 0.322 on the same definition.'),

    ('4.3 A modern sequence baseline and the structure increment',
     'Contact-graph features are usually evaluated against simple or absent sequence baselines. To '
     'place the graph results on a current footing we trained three reference models under the same '
     'protocol: a frozen ESM-2 650M embedding with a learned head, the same embedding fused with a '
     'six-layer GINE graph branch, and the graph branch alone. On the common S669 intersection '
     '(n = 508) the sequence baseline reached r = 0.392, the graph branch alone r = 0.368 and the '
     'fusion r = 0.379; on ssym (n = 342) the corresponding values were 0.516, 0.460 and 0.540 '
     '(Figure 4b). The graph branch therefore performs close to the sequence baseline, which is the '
     'relevant control: a null increment cannot be attributed to an under-powered graph branch.\n'
     'The paired difference between the fusion and the sequence-only model, evaluated with a '
     'protein-cluster bootstrap (B = 2,000), was \u0394r = -0.013 (95% CI -0.036 to +0.010; '
     'P = 0.27) on S669 and \u0394r = +0.025 (95% CI -0.007 to +0.070; P = 0.14) on ssym. Both '
     'intervals include zero. Read as an equivalence statement rather than as a null-hypothesis '
     'test, they exclude structural increments larger than |\u0394r| \u2248 0.04 on S669 and '
     '\u2248 0.07 on ssym. Within that resolution, adding a contact-graph representation to a modern '
     'sequence baseline did not improve \u0394\u0394G prediction.\n'
     'An ESM-2 zero-shot score, which carries no \u0394\u0394G supervision, reached r = 0.317 on '
     'S669 and 0.266 on ssym. Its weakness relative to the graph encoders on the same mutations '
     'indicates that the sequence advantage in the supervised comparison reflects learning on the '
     'shared training set rather than memorisation of the test proteins alone.'),

    ('4.4 An information budget: where the predictive signal is',
     'To attribute the available signal to information sources rather than to architectures, we '
     'fitted ridge models on the training set and evaluated them out of sample. Five '
     'physicochemical features (hydrophobicity change, volume change, and the mutant-side '
     'hydrophobicity, volume and charge) reached r = 0.390 on S669 and 0.306 on ssym; five '
     'graph-derived structural features (degree, normalised degree, nearest-neighbour distance, and '
     'mean and maximum neighbour distance at the mutation site) reached 0.297 and 0.301; and the two '
     'combined reached 0.427 and 0.408 (Figure 4a). Frozen ESM-2 embeddings used linearly reached '
     '0.238 and 0.345, understating what the same embeddings achieve through a nonlinear head '
     '(0.392 and 0.516), which is why both are reported. The ranking of information sources is '
     'benchmark-dependent while the increment from structure is small in both cases, and no tested '
     'combination exceeded the best single source by more than 0.04 on S669.'),

    ('4.5 Untyped contact counts, contact direction, and contact type',
     'The graph encoders consume edge features but do not outperform encoders without them, so we '
     'asked what a contact change contributes as an explicit quantity. Unaggregated counts were '
     'weakly and non-significantly related to \u0394\u0394G on 505 quality-controlled pairs: number '
     'broken r = 0.061 (95% CI -0.027 to +0.179), number formed r = 0.083 (95% CI -0.046 to +0.209), '
     'total edits r = 0.081 (95% CI -0.020 to +0.214).\n'
     'Directional summaries did not recover signal: the net change (formed minus broken) reached '
     'r = -0.013 (95% CI -0.127 to +0.077) and the broken-to-formed ratio r = 0.031 (95% CI -0.057 '
     'to +0.138). Breaking and forming are themselves correlated (r = 0.42), so netting them '
     'cancels rather than sharpens the signal.\n'
     'Stratifying by contact type did produce one significant quantity: the number of broken '
     'hydrophobic contacts correlated with \u0394\u0394G at r = 0.135 (95% CI 0.044 to 0.201), the '
     'only contact-derived feature to exclude zero. Broken electrostatic (r = 0.047), other broken '
     '(r = 0.027), formed hydrophobic (r = 0.096), formed electrostatic (r = -0.046) and other '
     'formed (r = 0.063) contacts were all non-significant. This is physically consistent, since '
     'losing hydrophobic contacts is destabilising. It is also not additive: five-fold '
     'protein-cluster cross-validation over the same benchmark gave r = 0.367 for physicochemical '
     'features alone and 0.357 after adding all six typed contact counts, a paired increment of '
     '\u0394r = -0.009 (95% CI -0.047 to +0.037). The hydrophobic contact change is therefore '
     'redundant with the physicochemical features, which already encode the mutation hydrophobicity '
     'change. A structural quantity can be individually associated with \u0394\u0394G and still '
     'carry no independent predictive information.'),
]

# Discussion 重构（替换现有 5.x 标题）
DISCUSSION = [
    ('5.1 What changes structurally, and what that does not imply',
     'The definition of a residue contact graph controls what the model can see. With identically '
     'processed wild-type and mutant structures, the four definitions spanned the full range from '
     'recording no contact change at all (C\u03b1) to recording a change in 90.7% of pairs '
     '(side-chain centroid), and the broken contacts were strongly localised around the mutation '
     'site (mean distance 4.03 \u00c5 versus 15.70 \u00c5 for contacts that persisted). The '
     'structural description is therefore highly sensitive to a design choice that is usually '
     'reported as a fixed parameter. Sensitivity of the description does not propagate to '
     'prediction: the same four definitions produced statistically indistinguishable predictions at '
     'every encoder capacity tested, and adding structure to a sequence baseline produced no '
     'reliable increment.'),
    ('5.2 Why representational sensitivity does not translate into predictive gain',
     'Three observations explain the decoupling. First, the sensitivity is largely redundant with '
     'information the model already has: the mutation itself determines the physicochemical change, '
     'and five such features reach r = 0.390, close to every graph encoder tested. Second, the '
     'additional quantities a sensitive definition supplies are individually weak and mutually '
     'correlated, so that directional and total summaries carry little marginal signal. Third, '
     'adding more informative edges does not by itself change what an aggregation function can '
     'extract; the encoders with explicit edge features did not outperform those without them.'),
    ('5.3 Where the signal is, and what to do about it',
     'The information budget localises the predictive signal in the mutation and its sequence '
     'context rather than in contact topology. Modern protein language model embeddings supply the '
     'largest single contribution (r = 0.392 on S669, 0.516 on ssym through a learned head), and '
     'physicochemical features supply a compact and surprisingly competitive baseline. A practical '
     'consequence is that effort spent refining the representation layer of a structure-based '
     'predictor has a low expected return within the range we tested, whereas effort spent on the '
     'sequence side or on the training data is more likely to move the metric.'),
    ('5.4 Limitations',
     'Four boundaries define the scope of these conclusions. The predictive results span encoders '
     'from 5.8 k to 509 k parameters, and we cannot exclude a representation effect in larger or '
     'differently biased architectures. ESM-2 was pre-trained on UniRef, so the sequence baseline '
     'may have seen proteins related to the test sets; the zero-shot variant carries no \u0394\u0394G '
     'supervision and is less affected, and it performs below the graph encoders on the same '
     'mutations. The ssym benchmark contains 15 proteins, which widens its intervals; the structure '
     'increment there can only be bounded above \u2248 0.07. Absolute contact counts depend on the '
     'modelling engine, and FoldX and SCWRL4 differ in whether they rebuild C\u03b2 coordinates, '
     'which changes the reported edit rate for the C\u03b2 definition. S669 was used during '
     'development, so its estimates are exploratory, and ssym serves as the development-independent '
     'secondary benchmark.'),
]


def main():
    shutil.copy2(SRC, DST)
    doc = Document(DST)
    paras = doc.paragraphs

    # ---- 1) 标题 ----
    i, p = find(paras, lambda t: t.startswith('Mutation-sensitive residue contact graphs'))
    if p:
        set_text(p, TITLE)
        print('  [1] 标题已更新')

    # ---- 2) 摘要 ----
    for key, txt in ABSTRACT.items():
        i, p = find(paras, lambda t, k=key: t.startswith(k))
        if p:
            set_text(p, key + ' ' + txt)
            print(f'  [2] 摘要 {key} 已更新')
        else:
            print(f'  [2] ⚠ 未找到 {key}')

    # ---- 3) 新增 Results 小节（插在 4.1 小节末尾，即 4.2 标题之前）----
    i42, p42 = find(paras, lambda t: t.startswith('4.2 '))
    anchor = paras[i42 - 1] if i42 else None
    if anchor is not None:
        cur = anchor
        for head, body in NEW_RESULTS:
            hp = insert_after(cur, head)
            try:
                hp.style = paras[i42].style
            except Exception:
                pass
            bp = insert_after(hp, body)
            cur = bp
        print(f'  [3] 已插入 {len(NEW_RESULTS)} 个 Results 小节')

    # ---- 4) Discussion 重构由独立脚本处理 ----
    doc.save(DST)          # ★ 必须先保存，否则插入的新 Results 小节丢失
    print(f'\n已保存 {DST}')


if __name__ == '__main__':
    main()
