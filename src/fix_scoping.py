# -*- coding: utf-8 -*-
"""
fix_scoping.py — R2-M1：把"被操纵的因子"写进标题与摘要
==========================================================
审稿人 2 的意见：标题与摘要说的是 "structure"，实际操纵的是
FoldX 建模突变体上的接触图定义。必须在每一处主张中点名。
"""
from docx import Document

DOC = r'D:\GED_mutation\manuscript_routeA.docx'

NEW_TITLE = ('Contact-graph definitions shape what a model describes but not what it predicts: '
             'a controlled comparison across five encoder capacities')

# 摘要 Motivation：点名被操纵因子与结构流
NEW_MOTIV = (
    'Motivation. Residue contact graphs are standard in structure-based \u0394\u0394G prediction, '
    'and their definition determines which mutation-induced side-chain changes a model can see. A '
    'more sensitive definition is usually assumed to yield better predictions. We tested that '
    'assumption while holding the encoder family, training data, leakage control and evaluation '
    'protocol fixed, using contact graphs built from FoldX-modelled mutants on fixed backbones. We '
    'separate what a definition reports about structural change from what it changes about '
    'prediction.')

# 摘要 Conclusion：点名范围
NEW_CONCL = (
    'Conclusion. The atom-level definition of a residue contact graph strongly determines what '
    'structural change a model describes, but within the tested encoder capacities, training data '
    'and evaluation protocol it did not measurably change \u0394\u0394G prediction, and structure '
    'did not add a reliable increment over a modern sequence baseline. No published structure-based '
    'predictor and no structural stream other than the FoldX-modelled one were evaluated, so these '
    'bounds apply to the representation layer as tested. We provide the audit protocol and the '
    'resulting feature-block accounting as a reusable procedure for weighing representation-layer '
    'effort against sequence and physicochemical alternatives in a new setting.')


def main():
    d = Document(DOC)
    n = 0
    for p in d.paragraphs:
        if not p.runs:
            continue
        t = p.text.strip()
        if t.startswith('What does structure contribute to'):
            p.runs[0].text = NEW_TITLE
            for r in p.runs[1:]:
                r.text = ''
            n += 1
        elif t.startswith('Motivation. Residue contact graphs are standard'):
            p.runs[0].text = NEW_MOTIV
            for r in p.runs[1:]:
                r.text = ''
            n += 1
        elif t.startswith('Conclusion. The atom-level definition of a residue contact graph'):
            p.runs[0].text = NEW_CONCL
            for r in p.runs[1:]:
                r.text = ''
            n += 1
    d.save(DOC)
    print(f'已更新 {n} 处（标题 + 摘要 Motivation + 摘要 Conclusion）')
    print(f'  新标题: {NEW_TITLE}')


if __name__ == '__main__':
    main()
